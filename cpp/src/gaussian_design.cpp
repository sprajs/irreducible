#include "irred/gaussian_design.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
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
bool output(long double v, double &d) {
  if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<double>::max())
    return false;
  d = static_cast<double>(v);
  return v == 0 ? d == 0 : std::fpclassify(d) == FP_NORMAL;
}
bool policy_valid(DesignPolicy p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
} // namespace
DesignProfile::DesignProfile(DesignProfile &&o) noexcept
    : gaussian_(std::move(o.gaussian_)),
      design_metadata_(std::move(o.design_metadata_)),
      gram_factor_(std::move(o.gram_factor_)), x_(std::move(o.x_)),
      wx_(std::move(o.wx_)), scales_(std::move(o.scales_)),
      status_(std::exchange(o.status_, DensityStatus::invalid_input)),
      numerical_status_(
          std::exchange(o.numerical_status_, numerics::Status::invalid_input)),
      rank_(std::exchange(o.rank_, DesignRank::unassessed)),
      gram_condition_(std::exchange(o.gram_condition_, 0)),
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
  // Includes transferred source + candidate owners simultaneously; Gram
  // Cholesky inverse-condition scratch, solves, normalized columns and result
  // buffers.
  b.add(n * p, 2 * sizeof(double));
  b.add(p * p, 8 * sizeof(long double));
  b.add(n, 16 * sizeof(long double));
  b.add(p, 16 * sizeof(long double));
  return b.result();
}
std::optional<std::size_t>
DesignProfile::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  charge(b, design_metadata_);
  b.embedded(gaussian_.retained_payload_bound(), sizeof(Gaussian));
  b.embedded(gram_factor_.retained_payload_bound(), sizeof(gram_factor_));
  b.vector(x_);
  b.vector(wx_);
  b.vector(scales_);
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
  if (!policy_valid(policy) || p < 2 || n < p ||
      !unique(m.ordered_parameter_ids) || m.parameter_units.size() != p ||
      m.residual_unit.empty() || m.design_identity.empty() ||
      m.dependence_identity.empty() || !unique(m.shared_nuisance_ids) ||
      std::any_of(
          m.parameter_units.begin(), m.parameter_units.end(),
          [](const auto &s) { return s.empty(); }))
    return out;
  for (const auto &id : m.shared_nuisance_ids)
    if (std::find(m.ordered_parameter_ids.begin(),
                  m.ordered_parameter_ids.end(),
                  id) == m.ordered_parameter_ids.end())
      return out;
  auto bound = preparation_payload_bound(g, p, m);
  if (!bound || n > policy.maximum_elements / p ||
      p > policy.maximum_elements / p ||
      *bound > policy.maximum_payload_bytes) {
    out.status_ = DensityStatus::numerical_failure;
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  if (x.size() != n * p || std::any_of(x.begin(), x.end(), [](double v) {
        return !std::isfinite(v);
      }))
    return out;
  auto fail = [&](numerics::Status s) {
    out.status_ = DensityStatus::numerical_failure;
    out.numerical_status_ = s;
    out.rank_ = DesignRank::unresolved;
  };
  out.x_.assign(x.begin(), x.end());
  out.wx_.resize(n * p);
  out.scales_.resize(p);
  std::vector<double> column(n), gram(p * p);
  // First normalize by max absolute entry, avoiding overflow in the first
  // solve.
  for (std::size_t j = 0; j < p; ++j) {
    double maximum = 0;
    for (std::size_t i = 0; i < n; ++i)
      maximum = std::max(maximum, std::abs(x[i * p + j]));
    if (maximum == 0) {
      fail(numerics::Status::singular);
      out.rank_ = DesignRank::deficient;
      return out;
    }
    for (std::size_t i = 0; i < n; ++i)
      column[i] = x[i * p + j] / maximum;
    auto solved =
        numerics::solve(g.factor_, column, policy.maximum_forward_sensitivity);
    if (solved.status != numerics::Status::ok) {
      fail(solved.status);
      return out;
    }
    long double q = 0;
    for (std::size_t i = 0; i < n; ++i)
      q += static_cast<long double>(column[i]) * solved.value[i];
    if (!(q > 0) || !output(static_cast<long double>(maximum) * std::sqrt(q),
                            out.scales_[j])) {
      fail(numerics::Status::overflow);
      return out;
    }
    const long double root = std::sqrt(q);
    for (std::size_t i = 0; i < n; ++i)
      if (!output(static_cast<long double>(solved.value[i]) / root,
                  out.wx_[i * p + j])) {
        fail(numerics::Status::outside_domain);
        return out;
      }
    out.preparation_sensitivity_ = std::max(
        out.preparation_sensitivity_, solved.estimated_forward_sensitivity);
  }
  for (std::size_t j = 0; j < p; ++j)
    for (std::size_t k = 0; k <= j; ++k) {
      long double v = 0;
      for (std::size_t i = 0; i < n; ++i)
        v += (static_cast<long double>(x[i * p + j]) / out.scales_[j]) *
             out.wx_[i * p + k];
      if (!output(v, gram[j * p + k])) {
        fail(numerics::Status::overflow);
        return out;
      }
      gram[k * p + j] = gram[j * p + k];
    }
  out.gram_factor_ =
      numerics::cholesky(gram, p, policy.maximum_elements, g.arithmetic());
  if (out.gram_factor_.status() != numerics::Status::ok) {
    fail(out.gram_factor_.status());
    return out;
  }
  out.gram_condition_ = out.gram_factor_.condition_estimate_inf();
  if (out.gram_condition_ * std::numeric_limits<double>::epsilon() >
      std::min(1e-8, policy.maximum_forward_sensitivity)) {
    fail(numerics::Status::conditioning_budget_exceeded);
    return out;
  }
  out.design_metadata_ = std::move(m);
  out.gaussian_ = std::move(g);
  std::vector<double>().swap(out.gaussian_.covariance_);
  out.status_ = DensityStatus::finite;
  out.numerical_status_ = numerics::Status::ok;
  out.rank_ = DesignRank::full_within_conditioning_contract;
  return out;
}
DesignResult DesignProfile::evaluate(std::span<const double> r,
                                     std::span<const std::string> ids,
                                     DesignPolicy policy) const {
  DesignResult out;
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
  if (!policy_valid(policy) || r.size() != n ||
      std::any_of(r.begin(), r.end(),
                  [](double v) { return !std::isfinite(v); }))
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
  if (preparation_sensitivity_ > policy.maximum_forward_sensitivity ||
      gram_condition_ * std::numeric_limits<double>::epsilon() >
          std::min(1e-8, policy.maximum_forward_sensitivity))
    return fail(numerics::Status::conditioning_budget_exceeded);
  auto wr =
      numerics::solve(gaussian_.factor_, r, policy.maximum_forward_sensitivity);
  if (wr.status != numerics::Status::ok)
    return fail(wr.status);
  std::vector<double> rhs(p);
  for (std::size_t j = 0; j < p; ++j) {
    long double v = 0;
    for (std::size_t i = 0; i < n; ++i)
      v += (static_cast<long double>(x_[i * p + j]) / scales_[j]) * wr.value[i];
    if (!output(v, rhs[j]))
      return fail(numerics::Status::overflow);
  }
  auto beta =
      numerics::solve(gram_factor_, rhs, policy.maximum_forward_sensitivity);
  if (beta.status != numerics::Status::ok)
    return fail(beta.status);
  out.coefficients.resize(p);
  out.adjusted_residuals.resize(n);
  for (std::size_t j = 0; j < p; ++j)
    if (!output(static_cast<long double>(beta.value[j]) / scales_[j],
                out.coefficients[j]))
      return fail(numerics::Status::outside_domain);
  for (std::size_t i = 0; i < n; ++i) {
    long double v = r[i];
    for (std::size_t j = 0; j < p; ++j)
      v -= static_cast<long double>(x_[i * p + j]) * out.coefficients[j];
    if (!output(v, out.adjusted_residuals[i]))
      return fail(numerics::Status::outside_domain);
  }
  auto we = numerics::solve(gaussian_.factor_, out.adjusted_residuals,
                            policy.maximum_forward_sensitivity);
  if (we.status != numerics::Status::ok)
    return fail(we.status);
  long double q = 0, defect = 0, denom = 0;
  for (std::size_t i = 0; i < n; ++i)
    q += static_cast<long double>(out.adjusted_residuals[i]) * we.value[i];
  for (std::size_t j = 0; j < p; ++j) {
    long double d = 0, scale = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const long double z =
          static_cast<long double>(x_[i * p + j]) / scales_[j];
      d += z * we.value[i];
      scale += std::abs(z * wr.value[i]);
    }
    defect = std::max(defect, std::abs(d));
    denom = std::max(denom, scale);
  }
  // Absolute fallback at zero RHS: normalized coordinates are dimensionless.
  const long double normal = denom == 0 ? defect : defect / denom;
  if (!std::isfinite(normal) || normal > policy.maximum_forward_sensitivity)
    return fail(numerics::Status::conditioning_budget_exceeded);
  if (q < 0 || !output(q, out.quadratic) ||
      !output(-q / 2, out.relative_log_score))
    return fail(numerics::Status::outside_domain);
  out.normalized_normal_equation_residual = static_cast<double>(normal);
  out.covariance_solve_backward_residual = we.backward_residual;
  out.covariance_solve_forward_sensitivity = we.estimated_forward_sensitivity;
  out.status = DensityStatus::finite;
  out.numerical_status = numerics::Status::ok;
  return out;
}
} // namespace irred::statistics
