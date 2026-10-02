#include "irred/gaussian_posterior.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <utility>
namespace irred::statistics {
namespace {
using S = numerics::Status;
bool normal(double x) { return x == 0 || std::fpclassify(x) == FP_NORMAL; }
bool cast(long double x, double &v) {
  if (!std::isfinite(x) || std::abs(x) > std::numeric_limits<double>::max())
    return false;
  v = static_cast<double>(x);
  return x == 0 ? v == 0 : std::fpclassify(v) == FP_NORMAL;
}
bool same(std::span<const std::string> a, std::span<const std::string> b) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}
bool valid(PosteriorPolicy p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
bool environment() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384;
}
void charge(detail::PayloadAccounting &b, const ParameterPrior &p) {
  b.strings(p.ordered_parameter_ids);
  b.strings(p.parameter_units);
  b.strings(p.shared_nuisance_ids);
  b.vector(p.mean);
  b.vector(p.covariance);
  for (const auto *s : {&p.prior_identity, &p.design_identity, &p.residual_unit,
                        &p.parameter_measure, &p.dependence_identity})
    b.string(*s);
}
std::optional<std::size_t> work(std::size_t n, std::size_t p) {
  if (!n || p < 2 || n > SIZE_MAX / n || p > SIZE_MAX / p || n > SIZE_MAX / p)
    return {};
  if (p > SIZE_MAX / 8)
    return {};
  detail::PayloadAccounting b(0);
  b.add(n * n, p);
  b.add(n * p, p);
  b.add(p * p, 8 * p);
  const auto units = b.result();
  if (!units || *units > SIZE_MAX / 32)
    return {};
  return *units * 32;
}
} // namespace
GaussianPosterior::GaussianPosterior(GaussianPosterior &&o) noexcept
    : source_(std::move(o.source_)), posterior_(std::move(o.posterior_)),
      prior_(std::move(o.prior_)), design_(std::move(o.design_)),
      precision_(std::move(o.precision_)),
      whitened_design_(std::move(o.whitened_design_)),
      covariance_error_(o.covariance_error_),
      source_whitening_error_(o.source_whitening_error_),
      source_inverse_norm_estimate_(o.source_inverse_norm_estimate_),
      status_(std::exchange(o.status_, DensityStatus::invalid_input)),
      numerical_status_(std::exchange(o.numerical_status_, S::invalid_input)) {}
GaussianPosterior &
GaussianPosterior::operator=(GaussianPosterior &&o) noexcept {
  if (this != &o) {
    this->~GaussianPosterior();
    new (this) GaussianPosterior(std::move(o));
  }
  return *this;
}
std::optional<std::size_t>
GaussianPosterior::preparation_payload_bound(const Gaussian &g,
                                             const ParameterPrior &p) noexcept {
  const auto n = g.metadata().ordered_ids.size(),
             k = p.ordered_parameter_ids.size();
  if (!work(n, k))
    return {};
  detail::PayloadAccounting b(sizeof(GaussianPosterior));
  b.embedded(g.retained_payload_bound(), sizeof(Gaussian));
  charge(b, p);
  charge(b, p); // candidate metadata and transferred prior capacities
  b.add(n * k, 4 * sizeof(long double));
  b.add(k * k, 32 * sizeof(long double));
  b.add(n + k, 16 * sizeof(long double));
  return b.result();
}
std::optional<std::size_t>
GaussianPosterior::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  b.embedded(source_.retained_payload_bound(), sizeof(Gaussian));
  b.embedded(posterior_.retained_payload_bound(), sizeof(Gaussian));
  charge(b, prior_);
  b.vector(design_);
  b.vector(precision_);
  b.vector(whitened_design_);
  return b.result();
}
std::optional<std::size_t>
GaussianPosterior::evaluation_payload_bound() const noexcept {
  const auto n = source_.metadata().ordered_ids.size(), p = prior_.mean.size();
  detail::PayloadAccounting b(sizeof(PosteriorMean) + sizeof(GaussianResult));
  b.add(n + p, 24 * sizeof(long double));
  const auto e = posterior_.evaluation_payload_bound(p, false);
  if (!e)
    return {};
  b.add(*e, 1);
  return b.result();
}
GaussianPosterior GaussianPosterior::prepare(Gaussian &&g,
                                             std::span<const double> x,
                                             std::span<const std::string> ids,
                                             ParameterPrior prior,
                                             PosteriorPolicy policy) {
  GaussianPosterior out;
  if (!environment()) {
    out.status_ = DensityStatus::unsupported_domain;
    return out;
  }
  const auto n = g.metadata().ordered_ids.size(),
             p = prior.ordered_parameter_ids.size();
  if (g.status() != DensityStatus::finite || !valid(policy) || !work(n, p))
    return out;
  const auto peak = preparation_payload_bound(g, prior), units = work(n, p);
  if (n > policy.maximum_elements / n || p > policy.maximum_elements / p ||
      n > policy.maximum_elements / p || !peak ||
      *peak > policy.maximum_payload_bytes ||
      *units > policy.maximum_work_units) {
    out.numerical_status_ = S::work_limit;
    return out;
  }
  if (!same(ids, g.metadata().ordered_ids)) {
    out.status_ = DensityStatus::incompatible_metadata;
    return out;
  }
  if (x.size() != n * p || prior.mean.size() != p ||
      prior.covariance.size() != p * p || prior.parameter_units.size() != p ||
      !g.priors().empty() || !g.mean_shift().empty() ||
      !prior.noise_independence_declared || prior.prior_identity.empty() ||
      prior.design_identity.empty() || prior.residual_unit.empty() ||
      prior.parameter_measure.empty() || prior.dependence_identity.empty())
    return out;
  for (std::size_t j = 0; j < p; ++j)
    if (prior.ordered_parameter_ids[j].empty() ||
        prior.parameter_units[j].empty() ||
        std::find(prior.ordered_parameter_ids.begin(),
                  prior.ordered_parameter_ids.begin() + j,
                  prior.ordered_parameter_ids[j]) !=
            prior.ordered_parameter_ids.begin() + j)
      return out;
  for (std::size_t j = 0; j < prior.shared_nuisance_ids.size(); ++j)
    if (std::find(prior.ordered_parameter_ids.begin(),
                  prior.ordered_parameter_ids.end(),
                  prior.shared_nuisance_ids[j]) ==
            prior.ordered_parameter_ids.end() ||
        std::find(prior.shared_nuisance_ids.begin(),
                  prior.shared_nuisance_ids.begin() + j,
                  prior.shared_nuisance_ids[j]) !=
            prior.shared_nuisance_ids.begin() + j)
      return out;
  for (auto values : {g.covariance(), x, std::span<const double>(prior.mean),
                      std::span<const double>(prior.covariance)})
    for (auto v : values)
      if (!normal(v))
        return out;
  const auto sf = numerics::cholesky(prior.covariance, p,
                                     policy.maximum_elements, g.arithmetic());
  if (sf.status() != S::ok) {
    out.numerical_status_ = sf.status();
    return out;
  }
  std::vector<long double> w(n * p), u(p * p);
  std::vector<double> rhs(std::max(n, p), 0);
  long double c_error = 0, s_error = 0;
  const long double eps = std::numeric_limits<long double>::epsilon();
  // These retained-factor estimates are empirical forward diagnostics, not
  // certified coordinates. Bound each Gram term before the binary64 cast.
  const auto coordinate_error = [&](const numerics::WhiteningResult &v,
                                    const numerics::Factorization &f,
                                    std::size_t size) {
    return std::sqrt(static_cast<long double>(size) *
                     f.condition_estimate_inf()) *
           (v.backward_residual + v.arithmetic_rounding_estimate +
            64 * eps * size);
  };
  for (std::size_t j = 0; j < p; ++j) {
    for (std::size_t i = 0; i < n; ++i)
      rhs[i] = x[i * p + j];
    auto v =
        numerics::whiten(g.factor_, std::span<const double>(rhs.data(), n),
                         policy.maximum_elements, policy.maximum_payload_bytes,
                         policy.maximum_forward_sensitivity);
    if (v.status != S::ok) {
      out.numerical_status_ = v.status;
      return out;
    }
    c_error = std::max(c_error, coordinate_error(v, g.factor_, n));
    for (std::size_t i = 0; i < n; ++i)
      w[i * p + j] = v.value[i];
    std::fill(rhs.begin(), rhs.begin() + p, 0);
    rhs[j] = 1;
    v = numerics::whiten(sf, std::span<const double>(rhs.data(), p),
                         policy.maximum_elements, policy.maximum_payload_bytes,
                         policy.maximum_forward_sensitivity);
    if (v.status != S::ok) {
      out.numerical_status_ = v.status;
      return out;
    }
    s_error = std::max(s_error, coordinate_error(v, sf, p));
    for (std::size_t i = 0; i < p; ++i)
      u[i * p + j] = v.value[i];
  }
  std::vector<double> precision(p * p);
  long double error_norm = 0;
  for (std::size_t i = 0; i < p; ++i) {
    long double row_error = 0;
    for (std::size_t j = 0; j < p; ++j) {
      long double v = 0, ca = 0, sa = 0;
      for (std::size_t k = 0; k < n; ++k) {
        const auto t = w[k * p + i] * w[k * p + j];
        v += t;
        ca += std::abs(t);
      }
      for (std::size_t k = 0; k < p; ++k) {
        const auto t = u[k * p + i] * u[k * p + j];
        v += t;
        sa += std::abs(t);
      }
      if (!cast(v, precision[i * p + j])) {
        out.numerical_status_ = S::overflow;
        return out;
      }
      row_error += (2 * c_error + c_error * c_error) * ca +
                   (2 * s_error + s_error * s_error) * sa +
                   (4 * (n + p) + 4) * eps * (ca + sa) +
                   std::abs(v - precision[i * p + j]);
    }
    error_norm = std::max(error_norm, row_error);
  }
  Metadata md;
  md.ordered_ids = prior.ordered_parameter_ids;
  md.measure = prior.parameter_measure;
  md.table_identity =
      "conditional parameter coordinates; " + prior.prior_identity;
  md.uncertainty_identity = "proper Gaussian posterior covariance";
  md.ordering_provenance = "explicit prior parameter order";
  md.calibration_provenance = prior.prior_identity;
  md.dependence_provenance = prior.dependence_identity;
  md.source_semantics =
      "derived conditional Gaussian, not observed measurements";
  auto posterior = prepare_gaussian(
      precision, MatrixKind::precision, std::move(md), policy.maximum_elements,
      policy.maximum_forward_sensitivity, g.arithmetic());
  if (posterior.status() != DensityStatus::finite) {
    out.status_ = posterior.status();
    out.numerical_status_ = posterior.numerical_status();
    return out;
  }
  long double inverse_norm = 0;
  for (std::size_t i = 0; i < p; ++i) {
    long double row = 0;
    for (std::size_t j = 0; j < p; ++j)
      row += std::abs(posterior.covariance()[i * p + j]);
    inverse_norm = std::max(inverse_norm, row);
  }
  const auto eta = inverse_norm * error_norm;
  if (!std::isfinite(eta) || eta >= .01L ||
      eta / (1 - eta) > policy.maximum_forward_sensitivity) {
    out.status_ = DensityStatus::numerical_failure;
    out.numerical_status_ = S::conditioning_budget_exceeded;
    return out;
  }
  out.design_.assign(x.begin(), x.end());
  out.precision_ = std::move(precision);
  out.whitened_design_ = std::move(w);
  out.prior_ = std::move(prior);
  out.posterior_ = std::move(posterior);
  out.covariance_error_ = static_cast<double>(eta / (1 - eta));
  out.source_whitening_error_ = static_cast<double>(c_error);
  long double source_norm = 0;
  for (std::size_t i = 0; i < n; ++i) {
    long double row = 0;
    for (std::size_t j = 0; j < n; ++j)
      row += std::abs(g.covariance()[i * n + j]);
    source_norm = std::max(source_norm, row);
  }
  out.source_inverse_norm_estimate_ =
      g.factor_.condition_estimate_inf() / source_norm;
  if (!std::isfinite(out.source_inverse_norm_estimate_)) {
    out.status_ = DensityStatus::numerical_failure;
    out.numerical_status_ = S::conditioning_budget_exceeded;
    return out;
  }
  out.source_ = std::move(g);
  out.status_ = DensityStatus::finite;
  out.numerical_status_ = S::ok;
  return out;
}
PosteriorMean GaussianPosterior::condition(std::span<const double> r,
                                           std::span<const std::string> ids,
                                           PosteriorPolicy policy) const {
  PosteriorMean out;
  if (!environment() || !valid(policy))
    return out;
  if (status_ != DensityStatus::finite) {
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!same(ids, source_.metadata().ordered_ids)) {
    out.status = DensityStatus::incompatible_metadata;
    return out;
  }
  const auto n = ids.size(), p = prior_.mean.size();
  const auto peak = evaluation_payload_bound();
  if (!peak || *peak > policy.maximum_payload_bytes ||
      n * p > policy.maximum_elements || p * p > policy.maximum_elements) {
    out.numerical_status = S::work_limit;
    return out;
  }
  if (r.size() != n)
    return out;
  std::vector<long double> centered(n), rhs(p), delta(p), rhs_error(p);
  long double centering_error = 0;
  const auto eps = std::numeric_limits<long double>::epsilon();
  for (std::size_t i = 0; i < n; ++i) {
    if (!normal(r[i]))
      return out;
    long double v = r[i], absolute = std::abs(v);
    for (std::size_t j = 0; j < p; ++j) {
      const auto t =
          static_cast<long double>(design_[i * p + j]) * prior_.mean[j];
      v -= t;
      absolute += std::abs(t);
    }
    centering_error = std::max(centering_error, (4 * p + 4) * eps * absolute);
    if (!std::isfinite(v)) {
      out.numerical_status = S::overflow;
      return out;
    }
    centered[i] = v;
  }
  const auto z = numerics::whiten(
      source_.factor_, centered, policy.maximum_elements,
      policy.maximum_payload_bytes, policy.maximum_forward_sensitivity);
  if (z.status != S::ok) {
    out.numerical_status = z.status;
    return out;
  }
  const long double z_relative =
      std::sqrt(static_cast<long double>(n) *
                source_.factor_.condition_estimate_inf()) *
      (z.backward_residual + z.arithmetic_rounding_estimate + 64 * eps * n);
  const long double z_absolute =
      std::sqrt(static_cast<long double>(n) * source_inverse_norm_estimate_) *
      centering_error;
  for (std::size_t j = 0; j < p; ++j) {
    long double absolute = 0, w_sum = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const auto t = whitened_design_[i * p + j] * z.value[i];
      rhs[j] += t;
      absolute += std::abs(t);
      w_sum += std::abs(whitened_design_[i * p + j]);
    }
    rhs_error[j] = (source_whitening_error_ + z_relative +
                    source_whitening_error_ * z_relative + (4 * n + 4) * eps) *
                       absolute +
                   w_sum * z_absolute;
  }
  out.value.resize(p);
  out.absolute_error_estimates.resize(p);
  long double magnitude = 0, inverse_norm = 0, rhs_error_norm = 0;
  for (auto v : rhs_error)
    rhs_error_norm = std::max(rhs_error_norm, v);
  for (std::size_t j = 0; j < p; ++j) {
    long double row = 0;
    for (std::size_t k = 0; k < p; ++k) {
      delta[j] +=
          static_cast<long double>(posterior_.covariance()[j * p + k]) * rhs[k];
      row += std::abs(posterior_.covariance()[j * p + k]);
    }
    inverse_norm = std::max(inverse_norm, row);
    if (!cast(prior_.mean[j] + delta[j], out.value[j])) {
      out.value.clear();
      out.numerical_status = S::overflow;
      return out;
    }
    magnitude = std::max(magnitude, std::abs(delta[j]));
  }
  long double residual = 0, denominator = 0;
  for (std::size_t i = 0; i < p; ++i) {
    long double v = -rhs[i], norm = std::abs(rhs[i]);
    for (std::size_t j = 0; j < p; ++j) {
      const auto t = static_cast<long double>(precision_[i * p + j]) *
                     (static_cast<long double>(out.value[j]) - prior_.mean[j]);
      v += t;
      norm += std::abs(t);
    }
    residual = std::max(residual, std::abs(v));
    denominator = std::max(denominator, norm);
  }
  const long double backward = denominator ? residual / denominator : 0;
  long double sensitivity = 0;
  for (std::size_t j = 0; j < p; ++j) {
    const long double absolute =
        inverse_norm * rhs_error_norm + covariance_error_ * magnitude +
        posterior_.factor_.condition_estimate_inf() * backward * magnitude +
        std::abs(prior_.mean[j] + delta[j] - out.value[j]) +
        (4 * p + 4) * eps * (std::abs(prior_.mean[j]) + magnitude);
    if (!cast(absolute, out.absolute_error_estimates[j])) {
      out.value.clear();
      out.absolute_error_estimates.clear();
      out.numerical_status = S::conditioning_budget_exceeded;
      return out;
    }
    sensitivity = std::max(
        sensitivity,
        absolute / (1 + std::abs(static_cast<long double>(out.value[j]))));
  }
  if (!std::isfinite(sensitivity) ||
      sensitivity > policy.maximum_forward_sensitivity) {
    out.value.clear();
    out.absolute_error_estimates.clear();
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = S::conditioning_budget_exceeded;
    return out;
  }
  out.backward_residual = static_cast<double>(backward);
  out.estimated_forward_sensitivity = static_cast<double>(sensitivity);
  out.status = DensityStatus::finite;
  out.numerical_status = S::ok;
  return out;
}
GaussianResult GaussianPosterior::log_density(
    std::span<const double> r, std::span<const std::string> rows,
    std::span<const double> beta, std::span<const std::string> parameters,
    PosteriorPolicy policy) const {
  GaussianResult out;
  if (!same(parameters, prior_.ordered_parameter_ids)) {
    out.density.status = DensityStatus::incompatible_metadata;
    return out;
  }
  auto mean = condition(r, rows, policy);
  if (mean.status != DensityStatus::finite) {
    out.density.status = mean.status;
    out.density.numerical_status = mean.numerical_status;
    return out;
  }
  if (beta.size() != mean.value.size())
    return out;
  std::vector<double> centered(beta.size());
  std::vector<long double> errors(beta.size());
  for (std::size_t i = 0; i < beta.size(); ++i) {
    if (!normal(beta[i]))
      return out;
    const auto wide = static_cast<long double>(beta[i]) - mean.value[i];
    if (!cast(wide, centered[i])) {
      out.density.numerical_status = S::overflow;
      return out;
    }
    errors[i] = mean.absolute_error_estimates[i] + std::abs(wide - centered[i]);
  }
  out = posterior_.evaluate(centered, parameters,
                            policy.maximum_forward_sensitivity);
  if (out.density.status != DensityStatus::finite)
    return out;
  long double projection = 0;
  for (std::size_t i = 0; i < beta.size(); ++i)
    for (std::size_t j = 0; j < beta.size(); ++j)
      projection +=
          std::abs(static_cast<long double>(precision_[i * beta.size() + j])) *
          (std::abs(static_cast<long double>(centered[i])) * errors[j] +
           .5L * errors[i] * errors[j]);
  const long double estimate =
      projection +
      .5L * (out.quadratic * covariance_error_ / (1 - covariance_error_) +
             beta.size() *
                 -std::log1p(-static_cast<long double>(covariance_error_)));
  if (!std::isfinite(estimate) ||
      estimate > policy.maximum_forward_sensitivity *
                     (1 + std::abs(out.density.log_value))) {
    out.density.status = DensityStatus::numerical_failure;
    out.density.numerical_status = S::conditioning_budget_exceeded;
    return out;
  }
  out.estimated_forward_sensitivity = std::max(
      out.estimated_forward_sensitivity,
      static_cast<double>(estimate / (1 + std::abs(out.density.log_value))));
  return out;
}
} // namespace irred::statistics
