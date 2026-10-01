#include "irred/correlated_calibration.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <utility>
namespace irred::statistics {
namespace {
void charge(detail::PayloadAccounting &b, const CalibrationPrior &p) noexcept {
  b.strings(p.ordered_parameter_ids);
  b.strings(p.parameter_units);
  b.vector(p.mean);
  b.vector(p.covariance);
  b.string(p.prior_identity);
  b.string(p.response_identity);
  b.string(p.residual_unit);
  b.string(p.calibration_identity);
  b.string(p.dependence_identity);
  b.string(p.measure_identity);
}
bool environment() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384;
}
bool normal(double v) { return v == 0 || std::fpclassify(v) == FP_NORMAL; }
bool cast(long double v, double &d) {
  if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<double>::max())
    return false;
  d = static_cast<double>(v);
  return v == 0 ? d == 0 : std::fpclassify(d) == FP_NORMAL;
}
bool same(std::span<const std::string> a, std::span<const std::string> b) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}
bool valid(CalibrationPolicy p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
} // namespace
CorrelatedCalibration::CorrelatedCalibration(CorrelatedCalibration &&o) noexcept
    : inverse_norm_estimate_(o.inverse_norm_estimate_),
      covariance_rounding_eta_(o.covariance_rounding_eta_),
      mean_rounding_inf_(o.mean_rounding_inf_),
      mean_rounding_l1_(o.mean_rounding_l1_), source_(std::move(o.source_)),
      effective_(std::move(o.effective_)), prior_(std::move(o.prior_)),
      response_(std::move(o.response_)), shift_(std::move(o.shift_)),
      status_(std::exchange(o.status_, DensityStatus::invalid_input)),
      numerical_status_(std::exchange(o.numerical_status_,
                                      numerics::Status::invalid_input)) {}
CorrelatedCalibration &
CorrelatedCalibration::operator=(CorrelatedCalibration &&o) noexcept {
  if (this != &o) {
    this->~CorrelatedCalibration();
    new (this) CorrelatedCalibration(std::move(o));
  }
  return *this;
}
std::optional<std::size_t> CorrelatedCalibration::preparation_payload_bound(
    const Gaussian &g, const CalibrationPrior &p) noexcept {
  auto n = g.metadata().ordered_ids.size(), k = p.ordered_parameter_ids.size();
  if (!n || k < 2 || n > SIZE_MAX / n || k > SIZE_MAX / k || n > SIZE_MAX / k)
    return {};
  detail::PayloadAccounting b(sizeof(CorrelatedCalibration));
  b.embedded(g.retained_payload_bound(), sizeof(Gaussian));
  charge(b, p);
  b.embedded(gaussian_preparation_payload_bound(n, MatrixKind::covariance,
                                                g.arithmetic(), g.metadata()),
             sizeof(Gaussian));
  // Candidate covariance, response, shift/centered residual; prior factor and
  // its condition-estimation scratch. Source remains simultaneously owned.
  b.add(n * n, 12 * sizeof(long double));
  b.add(n * k, sizeof(double));
  b.add(n, 4 * sizeof(long double));
  b.add(k * k, 12 * sizeof(long double));
  return b.result();
}
std::optional<std::size_t>
CorrelatedCalibration::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  b.embedded(source_.retained_payload_bound(), sizeof(Gaussian));
  b.embedded(effective_.retained_payload_bound(), sizeof(Gaussian));
  charge(b, prior_);
  b.vector(response_);
  b.vector(shift_);
  return b.result();
}
std::optional<std::size_t>
CorrelatedCalibration::evaluation_payload_bound() const noexcept {
  auto n = shift_.size();
  detail::PayloadAccounting b(sizeof(GaussianResult));
  b.add(n, sizeof(double));
  auto e = effective_.evaluation_payload_bound(n, false);
  if (!e)
    return {};
  b.add(*e, 1);
  return b.result();
}
CorrelatedCalibration
CorrelatedCalibration::prepare(Gaussian &&g, std::span<const double> x,
                               std::span<const std::string> ids,
                               CalibrationPrior p, CalibrationPolicy policy) {
  CorrelatedCalibration out;
  if (!environment()) {
    out.status_ = DensityStatus::unsupported_domain;
    return out;
  }
  auto n = g.metadata().ordered_ids.size(), k = p.ordered_parameter_ids.size();
  if (g.status() != DensityStatus::finite || !valid(policy) || k < 2 || !n ||
      n > policy.maximum_elements / n || k > policy.maximum_elements / k ||
      n > policy.maximum_elements / k)
    return out;
  auto peak = preparation_payload_bound(g, p);
  if (!peak || *peak > policy.maximum_payload_bytes) {
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  if (!same(ids, g.metadata().ordered_ids) ||
      p.measure_identity != g.metadata().measure) {
    out.status_ = DensityStatus::incompatible_metadata;
    return out;
  }
  if (x.size() != n * k || p.mean.size() != k || p.covariance.size() != k * k ||
      p.parameter_units.size() != k || !g.priors().empty() ||
      !g.mean_shift().empty() || !p.noise_independence_declared ||
      p.prior_identity.empty() || p.response_identity.empty() ||
      p.residual_unit.empty() || p.calibration_identity.empty() ||
      p.dependence_identity.empty() || p.measure_identity.empty())
    return out;
  for (std::size_t j = 0; j < k; ++j)
    if (p.ordered_parameter_ids[j].empty() || p.parameter_units[j].empty() ||
        std::find(p.ordered_parameter_ids.begin(),
                  p.ordered_parameter_ids.begin() + j,
                  p.ordered_parameter_ids[j]) !=
            p.ordered_parameter_ids.begin() + j)
      return out;
  for (auto v : x)
    if (!normal(v))
      return out;
  for (auto v : p.mean)
    if (!normal(v))
      return out;
  for (auto v : p.covariance)
    if (!normal(v))
      return out;
  auto prior_factor = numerics::cholesky(
      p.covariance, k, policy.maximum_elements, g.arithmetic());
  if (prior_factor.status() != numerics::Status::ok) {
    out.numerical_status_ = prior_factor.status();
    return out;
  }
  if (prior_factor.condition_estimate_inf() *
          std::numeric_limits<double>::epsilon() >
      policy.maximum_forward_sensitivity) {
    out.numerical_status_ = numerics::Status::conditioning_budget_exceeded;
    return out;
  }
  std::vector<double> c(n * n), shift(n);
  std::vector<long double> row_error(n, 0);
  const auto eps = std::numeric_limits<long double>::epsilon();
  const long double gamma = (4.L * k * k + 4) * eps;
  if (gamma >= 0.01L) {
    out.status_ = DensityStatus::unsupported_domain;
    return out;
  }
  for (std::size_t i = 0; i < n; ++i) {
    long double m = 0, ma = 0;
    for (std::size_t a = 0; a < k; ++a) {
      auto t = static_cast<long double>(x[i * k + a]) * p.mean[a];
      m += t;
      ma += std::abs(t);
    }
    if (!cast(m, shift[i])) {
      out.status_ = DensityStatus::numerical_failure;
      return out;
    }
    const auto me =
        gamma * ma + std::abs(m - static_cast<long double>(shift[i]));
    out.mean_rounding_inf_ = std::max(out.mean_rounding_inf_, me);
    out.mean_rounding_l1_ += me;
    for (std::size_t j = 0; j <= i; ++j) {
      long double v = g.covariance()[i * n + j], va = std::abs(v);
      for (std::size_t a = 0; a < k; ++a)
        for (std::size_t b = 0; b < k; ++b) {
          auto t = static_cast<long double>(x[i * k + a]) *
                   p.covariance[a * k + b] * x[j * k + b];
          v += t;
          va += std::abs(t);
        }
      if (!cast(v, c[i * n + j])) {
        out.status_ = DensityStatus::numerical_failure;
        return out;
      }
      const auto ce =
          gamma * va + std::abs(v - static_cast<long double>(c[i * n + j]));
      row_error[i] += ce;
      if (i != j)
        row_error[j] += ce;
      c[j * n + i] = c[i * n + j];
    }
  }
  auto effective = prepare_gaussian(
      c, MatrixKind::covariance, g.metadata(), policy.maximum_elements,
      policy.maximum_forward_sensitivity, g.arithmetic());
  if (effective.status() != DensityStatus::finite) {
    out.status_ = effective.status();
    out.numerical_status_ = effective.numerical_status();
    return out;
  }
  long double norm = 0;
  for (std::size_t i = 0; i < n; ++i) {
    long double row = 0;
    for (std::size_t j = 0; j < n; ++j)
      row += std::abs(c[i * n + j]);
    norm = std::max(norm, row);
  }
  out.inverse_norm_estimate_ =
      static_cast<double>(effective.factor_.condition_estimate_inf() / norm);
  out.covariance_rounding_eta_ = static_cast<double>(
      out.inverse_norm_estimate_ *
      *std::max_element(row_error.begin(), row_error.end()));
  if (!std::isfinite(out.inverse_norm_estimate_) ||
      !std::isfinite(out.covariance_rounding_eta_) ||
      out.covariance_rounding_eta_ >= 0.01) {
    out.status_ = DensityStatus::numerical_failure;
    return out;
  }
  // Finish every allocation before consuming the caller's source.
  out.response_.assign(x.begin(), x.end());
  out.shift_ = std::move(shift);
  out.prior_ = std::move(p);
  out.effective_ = std::move(effective);
  out.source_ = std::move(g);
  out.status_ = DensityStatus::finite;
  out.numerical_status_ = numerics::Status::ok;
  return out;
}
CalibrationResult
CorrelatedCalibration::evaluate(std::span<const double> r,
                                std::span<const std::string> ids,
                                CalibrationPolicy policy) const {
  CalibrationResult out;
  if (!environment()) {
    out.density.status = DensityStatus::unsupported_domain;
    return out;
  }
  if (status_ != DensityStatus::finite || !valid(policy) ||
      r.size() != shift_.size())
    return out;
  if (!same(ids, source_.metadata().ordered_ids)) {
    out.density.status = DensityStatus::incompatible_metadata;
    return out;
  }
  auto peak = evaluation_payload_bound();
  if (!peak || *peak > policy.maximum_payload_bytes ||
      r.size() > policy.maximum_elements) {
    out.density.numerical_status = numerics::Status::work_limit;
    return out;
  }
  std::vector<double> centered(r.size());
  long double mean_inf = mean_rounding_inf_, mean_l1 = mean_rounding_l1_;
  for (std::size_t i = 0; i < r.size(); ++i) {
    if (!normal(r[i]) ||
        !cast(static_cast<long double>(r[i]) - shift_[i], centered[i])) {
      out.density.status = DensityStatus::numerical_failure;
      return out;
    }
    const long double ce =
        std::abs((static_cast<long double>(r[i]) - shift_[i]) - centered[i]) +
        std::numeric_limits<long double>::epsilon() *
            (std::abs(static_cast<long double>(r[i])) + std::abs(shift_[i]));
    mean_inf = std::max(mean_inf, mean_rounding_inf_ + ce);
    mean_l1 += ce;
  }
  static_cast<GaussianResult &>(out) =
      effective_.evaluate(centered, ids, policy.maximum_forward_sensitivity);
  if (out.density.status != DensityStatus::finite)
    return out;
  long double einf = 0;
  for (auto v : centered)
    einf = std::max(einf, std::abs(static_cast<long double>(v)));
  const long double eta = covariance_rounding_eta_;
  const long double qe =
      (eta * out.quadratic +
       inverse_norm_estimate_ * (2 * einf * mean_l1 + mean_inf * mean_l1)) /
      (1 - eta);
  const long double de = centered.size() * (-std::log1p(-eta));
  const long double le = (qe + de) / 2;
  if (qe > policy.maximum_forward_sensitivity * (1 + std::abs(out.quadratic)) ||
      de > policy.maximum_forward_sensitivity *
               (1 + std::abs(out.log_determinant)) ||
      le > policy.maximum_forward_sensitivity *
               (1 + std::abs(out.density.log_value))) {
    out = {};
    out.density.status = DensityStatus::numerical_failure;
    out.density.numerical_status =
        numerics::Status::conditioning_budget_exceeded;
  }
  if (out.density.status == DensityStatus::finite) {
    out.quadratic_rounding_estimate = static_cast<double>(qe);
    out.log_determinant_rounding_estimate = static_cast<double>(de);
    out.log_density_rounding_estimate = static_cast<double>(le);
  }
  return out;
}
} // namespace irred::statistics
