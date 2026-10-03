#include "irred/gaussian_design.hpp"
#include "payload_accounting.hpp"
#include "retained_gaussian.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <numeric>
#include <utility>
namespace irred::statistics {
namespace {
void charge(detail::PayloadAccounting &b, const DesignMetadata &m) noexcept {
  b.strings(m.ordered_parameter_ids);
  b.strings(m.parameter_units);
  b.strings(m.shared_nuisance_ids);
  b.string(m.residual_unit);
  b.string(m.design_identity);
  b.string(m.dependence_identity);
}
bool unique(std::span<const std::string> ids) {
  for (std::size_t i = 0; i < ids.size(); ++i)
    if (ids[i].empty() ||
        std::find(ids.begin(), ids.begin() + i, ids[i]) != ids.begin() + i)
      return false;
  return true;
}
bool same(std::span<const std::string> a, std::span<const std::string> b) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}
bool policy_valid(DesignPolicy p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
bool supported_arithmetic_environment() noexcept {
  return std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384 &&
         std::fegetround() == FE_TONEAREST;
}
} // namespace
DesignProfile::DesignProfile(DesignProfile &&o) noexcept
    : gaussian_(std::move(o.gaussian_)),
      design_metadata_(std::move(o.design_metadata_)), x_(std::move(o.x_)),
      whitened_design_(std::move(o.whitened_design_)), qr_(std::move(o.qr_)),
      scales_(std::move(o.scales_)), tau_(std::move(o.tau_)),
      pivot_(std::move(o.pivot_)),
      status_(std::exchange(o.status_, DensityStatus::invalid_input)),
      numerical_status_(
          std::exchange(o.numerical_status_, numerics::Status::invalid_input)),
      rank_(std::exchange(o.rank_, DesignRank::unassessed)),
      preparation_whitening_attempts_(std::exchange(o.preparation_whitening_attempts_,0)),
      preparation_whitenings_completed_(std::exchange(o.preparation_whitenings_completed_,0)),
      triangular_condition_(std::exchange(o.triangular_condition_, 0)),
      transpose_triangular_condition_(
          std::exchange(o.transpose_triangular_condition_, 0)),
      preparation_sensitivity_(std::exchange(o.preparation_sensitivity_, 0)) {}
DesignProfile &DesignProfile::operator=(DesignProfile &&o) noexcept {
  if (this != &o) {
    this->~DesignProfile();
    new (this) DesignProfile(std::move(o));
  }
  return *this;
}
std::optional<std::size_t>
DesignProfile::preparation_payload_bound(const Gaussian &g, std::size_t p,
                                         const DesignMetadata &m) noexcept {
  const auto n = g.metadata().ordered_ids.size();
  if (!n || !p || n > SIZE_MAX / p || p > SIZE_MAX / p || n > SIZE_MAX / n)
    return {};
  detail::PayloadAccounting b(sizeof(DesignProfile));
  b.embedded(g.retained_payload_bound(), sizeof(Gaussian));
  charge(b, m);
  // Source factor is borrowed during preparation and transferred only on
  // success. A is retained for actual-output checks; QR is its decomposition.
  b.add(n * p, sizeof(double) + 2 * sizeof(long double));
  b.add(n, 4 * sizeof(long double));
  b.add(p, 12 * sizeof(long double) + sizeof(std::size_t));
  b.add(1, sizeof(numerics::WhiteningResult));
  return b.result();
}
std::optional<std::size_t>
DesignProfile::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  charge(b, design_metadata_);
  b.embedded(gaussian_.retained_payload_bound(), sizeof(Gaussian));
  b.vector(x_);
  b.vector(whitened_design_);
  b.vector(qr_);
  b.vector(scales_);
  b.vector(tau_);
  b.vector(pivot_);
  return b.result();
}
std::optional<std::size_t>
DesignProfile::evaluation_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(DesignResult));
  b.add(gaussian_.metadata().ordered_ids.size(), 16 * sizeof(long double));
  b.add(scales_.size(), 16 * sizeof(long double));
  return b.result();
}
DesignProfile DesignProfile::prepare(Gaussian &&g, std::span<const double> x,
                                     std::span<const std::string> ids,
                                     DesignMetadata m, DesignPolicy policy) {
  DesignProfile out;
  if (!supported_arithmetic_environment()) {
    out.status_ = DensityStatus::unsupported_domain;
    out.numerical_status_ = numerics::Status::outside_domain;
    return out;
  }
  const auto n = g.metadata().ordered_ids.size(),
             p = m.ordered_parameter_ids.size();
  if (g.status() != DensityStatus::finite) {
    out.status_ = g.status();
    out.numerical_status_ = g.numerical_status();
    return out;
  }
  if (!same(ids, g.metadata().ordered_ids) || !g.priors().empty()) {
    out.status_ = DensityStatus::incompatible_metadata;
    return out;
  }
  if (!policy_valid(policy) || p < 2 || n < p)
    return out;
  auto bound = preparation_payload_bound(g, p, m);
  if (!bound || n > policy.maximum_elements / p ||
      p > policy.maximum_elements / p ||
      *bound > policy.maximum_payload_bytes) {
    out.status_ = DensityStatus::numerical_failure;
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  if (!unique(m.ordered_parameter_ids) || m.parameter_units.size() != p ||
      m.residual_unit.empty() || m.design_identity.empty() ||
      m.dependence_identity.empty() || !unique(m.shared_nuisance_ids) ||
      std::any_of(m.parameter_units.begin(), m.parameter_units.end(),
                  [](const auto &s) { return s.empty(); }))
    return out;
  for (const auto &id : m.shared_nuisance_ids)
    if (std::find(m.ordered_parameter_ids.begin(),
                  m.ordered_parameter_ids.end(),
                  id) == m.ordered_parameter_ids.end())
      return out;
  if (x.size() != n * p || std::any_of(x.begin(), x.end(), [](double v) {
        return !std::isfinite(v);
      }))
    return out;
  auto fail = [&](numerics::Status s) {
    out.status_ = DensityStatus::numerical_failure;
    out.numerical_status_ = s;
    out.rank_ = DesignRank::unresolved;
  };
  // Return the actual failed whitening prefix even if a later retained-QR
  // allocation refuses. The source Gaussian transfers only at success.
  try {
  out.x_.assign(x.begin(), x.end());
  out.whitened_design_.resize(n * p);
  out.scales_.resize(p);
  std::vector<long double> column(n);
  for (std::size_t j = 0; j < p; ++j) {
    double maximum = 0;
    for (std::size_t i = 0; i < n; ++i)
      maximum = std::max(maximum, std::abs(x[i * p + j]));
    if (maximum == 0) {
      fail(numerics::Status::singular);
      out.rank_ = DesignRank::deficient;
      return out;
    }
    // Wide normalization preserves small entries of a binary64 column even
    // when its range is wider than binary64 normalized coordinates allow.
    for (std::size_t i = 0; i < n; ++i)
      column[i] = static_cast<long double>(x[i * p + j]) / maximum;
    ++out.preparation_whitening_attempts_;
    auto whitened =
        numerics::whiten(g.factor_, std::span<const long double>(column),
                         policy.maximum_elements, policy.maximum_payload_bytes,
                         policy.maximum_forward_sensitivity);
    if (whitened.status != numerics::Status::ok) {
      fail(whitened.status);
      return out;
    }
    ++out.preparation_whitenings_completed_;
    long double norm = 0;
    for (const auto v : whitened.value)
      norm = std::hypot(norm, v);
    out.scales_[j] = static_cast<long double>(maximum) * norm;
    if (!(norm > 0) || !std::isfinite(out.scales_[j]) ||
        std::fpclassify(out.scales_[j]) == FP_SUBNORMAL) {
      fail(numerics::Status::outside_domain);
      return out;
    }
    for (std::size_t i = 0; i < n; ++i) {
      const auto value = whitened.value[i] / norm;
      if (!std::isfinite(value) || std::fpclassify(value) == FP_SUBNORMAL ||
          (whitened.value[i] != 0 && value == 0)) {
        fail(numerics::Status::outside_domain);
        return out;
      }
      out.whitened_design_[i * p + j] = value;
    }
    out.preparation_sensitivity_ = std::max(
        out.preparation_sensitivity_, whitened.arithmetic_rounding_estimate);
  }
  out.qr_ = out.whitened_design_;
  out.tau_.resize(p);
  out.pivot_.resize(p);
  std::iota(out.pivot_.begin(), out.pivot_.end(), 0);
  const auto rank_floor = (static_cast<long double>(n) + p) *
                          std::numeric_limits<long double>::epsilon();
  for (std::size_t k = 0; k < p; ++k) {
    std::size_t best = k;
    long double largest = -1;
    // Recompute trailing norms: downdates can lose the decisive small tail.
    for (std::size_t j = k; j < p; ++j) {
      long double norm = 0;
      for (std::size_t i = k; i < n; ++i)
        norm = std::hypot(norm, out.qr_[i * p + j]);
      if (norm > largest) {
        largest = norm;
        best = j;
      }
    }
    if (!std::isfinite(largest) || largest <= rank_floor) {
      fail(numerics::Status::singular);
      return out;
    }
    if (best != k) {
      for (std::size_t i = 0; i < n; ++i)
        std::swap(out.qr_[i * p + k], out.qr_[i * p + best]);
      std::swap(out.pivot_[k], out.pivot_[best]);
    }
    const auto first = out.qr_[k * p + k];
    const auto diagonal = -std::copysign(largest, first);
    const auto leading = first - diagonal;
    out.tau_[k] = (diagonal - first) / diagonal;
    for (std::size_t i = k + 1; i < n; ++i)
      out.qr_[i * p + k] /= leading;
    out.qr_[k * p + k] = diagonal;
    for (std::size_t j = k + 1; j < p; ++j) {
      long double dot = out.qr_[k * p + j];
      for (std::size_t i = k + 1; i < n; ++i)
        dot += out.qr_[i * p + k] * out.qr_[i * p + j];
      dot *= out.tau_[k];
      out.qr_[k * p + j] -= dot;
      for (std::size_t i = k + 1; i < n; ++i)
        out.qr_[i * p + j] -= out.qr_[i * p + k] * dot;
    }
  }
  for (const auto value : out.qr_)
    if (!std::isfinite(value) || std::fpclassify(value) == FP_SUBNORMAL) {
      fail(numerics::Status::outside_domain);
      return out;
    }
  // Norm-infinity conditioning of triangular R in pivoted, equilibrated
  // coordinates. Q preserves 2-norm, not infinity norm; this is not kappa2(A).
  const auto conditions=detail::qr_conditions(out.qr_,p);
  const auto condition=conditions.triangular;
  const auto sensitivity =
      condition * (std::numeric_limits<double>::epsilon() + rank_floor);
  if (!std::isfinite(condition) || !std::isfinite(sensitivity) ||
      sensitivity > std::min(1e-8, policy.maximum_forward_sensitivity)) {
    fail(numerics::Status::conditioning_budget_exceeded);
    return out;
  }
  const auto transpose_condition=conditions.transpose;
  if (!std::isfinite(transpose_condition) ||
      transpose_condition > std::numeric_limits<double>::max()) {
    fail(numerics::Status::overflow);
    return out;
  }
  out.triangular_condition_ = static_cast<double>(condition);
  out.transpose_triangular_condition_ =
      static_cast<double>(transpose_condition);
  out.preparation_sensitivity_ =
      std::max(out.preparation_sensitivity_, static_cast<double>(sensitivity));
  out.design_metadata_ = std::move(m);
  out.gaussian_ = std::move(g);
  std::vector<double>().swap(out.gaussian_.covariance_);
  out.status_ = DensityStatus::finite;
  out.numerical_status_ = numerics::Status::ok;
  out.rank_ = DesignRank::full_within_conditioning_contract;
  return out;
  } catch (const std::bad_alloc &) {
    fail(numerics::Status::work_limit);
    return out;
  }
}
DesignResult DesignProfile::evaluate(std::span<const double> r,
                                     std::span<const std::string> ids,
                                     DesignPolicy policy) const {
  DesignResult out;
  if (!supported_arithmetic_environment()) {
    out.status = DensityStatus::unsupported_domain;
    out.numerical_status = numerics::Status::outside_domain;
    return out;
  }
  if (status_ != DensityStatus::finite) {
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!same(ids, metadata().ordered_ids)) {
    out.status = DensityStatus::incompatible_metadata;
    return out;
  }
  const auto n = metadata().ordered_ids.size(), p = scales_.size();
  if (!policy_valid(policy) || r.size() != n)
    return out;
  auto fail = [](numerics::Status s) {
    DesignResult f;
    f.status = DensityStatus::numerical_failure;
    f.numerical_status = s;
    return f;
  };
  auto bound = evaluation_payload_bound();
  if (!bound || *bound > policy.maximum_payload_bytes ||
      n > policy.maximum_elements || p > policy.maximum_elements)
    return fail(numerics::Status::work_limit);
  if (std::any_of(r.begin(), r.end(),
                  [](double v) { return !std::isfinite(v); }))
    return out;
  const auto sensitivity =
      triangular_condition_ * (std::numeric_limits<double>::epsilon() +
                               (static_cast<long double>(n) + p) *
                                   std::numeric_limits<long double>::epsilon());
  if (preparation_sensitivity_ > policy.maximum_forward_sensitivity ||
      sensitivity > std::min(1e-8, policy.maximum_forward_sensitivity))
    return fail(numerics::Status::conditioning_budget_exceeded);
  auto wr = numerics::whiten(gaussian_.factor_, r, policy.maximum_elements,
                             policy.maximum_payload_bytes,
                             policy.maximum_forward_sensitivity);
  if (wr.status != numerics::Status::ok)
    return fail(wr.status);
  auto transformed = wr.value;
  const auto view=detail::RetainedQrAccess::view(*this);
  detail::qr_transform(view,tau_,transformed);
  std::vector<long double> column_scales(p,0);
  for(std::size_t j=0;j<p;++j)
    for(std::size_t i=0;i<n;++i)
      column_scales[j]+=std::abs(whitened_design_[i*p+j]*wr.value[i]);
  return detail::qr_finish_fit(view,r,std::span<const long double>(transformed),
      std::span<const long double>(column_scales),
      [&](std::size_t i,std::size_t j){return x_[i*p+j];},
      [&](std::size_t i,std::size_t j){return whitened_design_[i*p+j];},
      [&](std::span<const double> adjusted){return numerics::whiten(gaussian_.factor_,adjusted,
          policy.maximum_elements,policy.maximum_payload_bytes,policy.maximum_forward_sensitivity);},
      policy.maximum_forward_sensitivity);
}

std::optional<std::size_t> DesignProfile::estimator_variance_payload_bound(
    const LinearFunctionalMetadata &m) const noexcept {
  detail::PayloadAccounting b(sizeof(EstimatorVarianceResult));
  b.strings(m.weight_units);
  b.string(m.functional_identity);
  b.string(m.output_unit);
  b.add(scales_.size(), 2 * sizeof(long double));
  return b.result();
}
EstimatorVarianceResult DesignProfile::estimator_variance(
    std::span<const double> w, std::span<const std::string> ids,
    LinearFunctionalMetadata m, DesignPolicy policy) const {
  EstimatorVarianceResult out;
  auto fail = [](numerics::Status ns) {
    EstimatorVarianceResult result;
    result.status = DensityStatus::numerical_failure;
    result.numerical_status = ns;
    return result;
  };
  if (!supported_arithmetic_environment()) {
    out.status = DensityStatus::unsupported_domain;
    out.numerical_status = numerics::Status::outside_domain;
    return out;
  }
  if (status_ != DensityStatus::finite) {
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!same(ids, design_metadata_.ordered_parameter_ids)) {
    out.status = DensityStatus::incompatible_metadata;
    return out;
  }
  const auto p = scales_.size();
  if (!policy_valid(policy) || w.size() != p || m.weight_units.size() != p ||
      m.functional_identity.empty() || m.output_unit.empty() ||
      std::any_of(m.weight_units.begin(), m.weight_units.end(),
                  [](const auto &unit) { return unit.empty(); }) ||
      std::any_of(w.begin(), w.end(),
                  [](double value) { return !std::isfinite(value); }))
    return out;
  const auto bound = estimator_variance_payload_bound(m);
  if (!bound || *bound > policy.maximum_payload_bytes ||
      p > policy.maximum_elements)
    return fail(numerics::Status::work_limit);
  if (preparation_sensitivity_ > policy.maximum_forward_sensitivity)
    return fail(numerics::Status::conditioning_budget_exceeded);
  out=detail::qr_variance(detail::RetainedQrAccess::view(*this),w,policy.maximum_forward_sensitivity);
  if(out.status==DensityStatus::finite) out.metadata=std::move(m);
  return out;
}

} // namespace irred::statistics
