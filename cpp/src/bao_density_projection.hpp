#pragma once
#include "irred/bao.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace irred::bao::detail {
inline long double
inverse_covariance_norm_estimate(const numerics::Factorization &factor,
                                 std::span<const double> covariance, size_t n) {
  long double norm = 0;
  for (size_t j = 0; j < n; ++j) {
    long double row = 0;
    for (size_t k = 0; k < n; ++k)
      row += std::abs((long double)covariance[j * n + k]);
    norm = std::max(norm, row);
  }
  return factor.condition_estimate_inf() / norm;
}
inline numerics::Status subtract_observations(std::span<const double> observed,
                                              std::span<const double> mu,
                                              std::vector<double> &r,
                                              std::vector<long double> &eps) {
  for (size_t j = 0; j < observed.size(); ++j) {
    const long double wide = (long double)observed[j] - mu[j];
    const double cast = static_cast<double>(wide);
    if (!std::isfinite(cast) ||
        (cast != 0 && std::fpclassify(cast) != FP_NORMAL) ||
        (cast == 0 && wide != 0))
      return numerics::Status::overflow;
    r.push_back(cast);
    if (!eps.empty())
      eps[j] +=
          std::abs(wide - cast) + std::numeric_limits<long double>::epsilon() *
                                      (std::abs((long double)observed[j]) +
                                       std::abs((long double)mu[j]));
  }
  return numerics::Status::ok;
}
inline numerics::Status project_density(
    const numerics::Factorization &factor, const statistics::Gaussian &gaussian,
    std::span<const std::string> ordered_ids, std::span<const double> r,
    std::span<const long double> eps, long double inverse_norm_estimate,
    double maximum_forward_sensitivity,
    double maximum_projection_log_density_error,
    std::optional<double> &projection_estimate,
    std::optional<statistics::GaussianResult> &result) {
  const size_t n = r.size();
  auto cause = numerics::Status::ok;
  const auto solved = numerics::solve(factor, r, maximum_forward_sensitivity);
  cause = solved.status;
  if (cause == numerics::Status::ok) {
    long double ainf = 0, einf = 0, e1 = 0;
    for (size_t j = 0; j < n; ++j) {
      ainf = std::max(ainf, std::abs((long double)solved.value[j]));
      einf = std::max(einf, eps[j]);
      e1 += eps[j];
    }
    const long double estimate =
        ainf * e1 + inverse_norm_estimate * einf * e1 / 2;
    const double reported = static_cast<double>(estimate);
    if (!std::isfinite(estimate) || estimate < 0 || !std::isfinite(reported))
      cause = numerics::Status::overflow;
    else if (estimate != 0 &&
             (reported == 0 || std::fpclassify(reported) != FP_NORMAL))
      cause = numerics::Status::outside_domain;
    else {
      projection_estimate = reported;
      if (estimate > maximum_projection_log_density_error)
        cause = numerics::Status::conditioning_budget_exceeded;
    }
  }
  if (cause == numerics::Status::ok) {
    result = gaussian.evaluate(r, ordered_ids, maximum_forward_sensitivity);
    cause = result->density.numerical_status;
  }
  return cause;
}
} // namespace irred::bao::detail
