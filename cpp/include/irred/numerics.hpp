#pragma once
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
namespace irred::numerics {
enum class Status {
  ok,
  invalid_input,
  nonfinite_input,
  overflow,
  work_limit,
  outside_domain,
  singular,
  not_positive_definite,
  conditioning_budget_exceeded
};
// error_estimate is an empirical estimator or arithmetic diagnostic, never a
// proven bound.
struct ScalarResult {
  Status status = Status::invalid_input;
  double value = 0;
  double error_estimate = 0;
  std::size_t evaluations = 0;
};
ScalarResult compensated_sum(std::span<const double>) noexcept;
ScalarResult
    log_sum_exp(std::span<const double>) noexcept; // nonempty finite input
ScalarResult log1p_checked(double) noexcept;       // x > -1
ScalarResult expm1_checked(double) noexcept;       // finite result required
ScalarResult
log_gamma_positive(double) noexcept; // bounded first slice [0.125,100]
ScalarResult interpolate_linear(
    std::span<const double> x, std::span<const double> y,
    double query) noexcept; // strictly increasing, no extrapolation
struct IntegrationPolicy {
  double absolute_tolerance;
  double relative_tolerance;
  std::size_t max_evaluations;
  unsigned max_depth;
};
using Integrand = double (*)(double, const void *);
// Caller supplies smooth resolved interval; split known discontinuities before
// calls. Callback must not throw. Narrow unobserved features can defeat
// empirical estimate.
ScalarResult integrate(Integrand, const void *, double lower, double upper,
                       IntegrationPolicy);
enum class Arithmetic { binary64_legacy_v1 = 0, longdouble_cpu_v1 = 1 };
class Factorization {
public:
  Factorization() = default;
  Factorization(const Factorization &) = default;
  Factorization &operator=(const Factorization &) = default;
  // Moving preserves the destination and leaves the source explicitly invalid.
  // Self move assignment is a no-op and retains the original valid state.
  Factorization(Factorization &&) noexcept;
  Factorization &operator=(Factorization &&) noexcept;

  Status status() const noexcept { return status_; }
  Arithmetic arithmetic() const noexcept { return arithmetic_; }
  const char *arithmetic_id() const noexcept {
    if (arithmetic_ == Arithmetic::longdouble_cpu_v1)
      return "F02/longdouble-cpu/v1";
    if (arithmetic_ == Arithmetic::binary64_legacy_v1)
      return "F02/binary64-legacy/v1";
    return "F02/unsupported-arithmetic";
  }
  std::size_t size() const noexcept { return n_; }
  // Owner header plus retained vector capacities; excludes allocator/RSS.
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  double log_determinant() const noexcept { return log_determinant_; }
  double condition_estimate_inf() const noexcept {
    return condition_estimate_inf_;
  }

private:
  void swap(Factorization &) noexcept;
  Status status_ = Status::invalid_input;
  std::size_t n_ = 0;
  Arithmetic arithmetic_ = Arithmetic::binary64_legacy_v1;
  std::vector<double> lower_, original_;
  std::vector<long double> wide_lower_;
  double log_determinant_ = 0, norm_inf_ = 0, condition_estimate_inf_ = 0;
  long double wide_norm_inf_ = 0, wide_condition_inf_ = 0;
  friend Factorization cholesky(std::span<const double>, std::size_t,
                                std::size_t, Arithmetic);
  friend struct FactorAccess;
};
// Finite exactly symmetric row-major matrix; checked shape/resource limit, no
// repair/jitter. Initial qualification covers tested n<=16 only; size support
// is not qualification.
Factorization cholesky(std::span<const double>, std::size_t n,
                       std::size_t maximum_elements,
                       Arithmetic arithmetic = Arithmetic::binary64_legacy_v1);
struct SolveResult {
  Status status = Status::invalid_input;
  std::vector<double> value;
  double backward_residual = 0;
  double estimated_forward_sensitivity = 0;
  double pre_cast_backward_residual = 0, pre_cast_sensitivity_estimate = 0,
         output_rounding_error_inf = 0, output_rounding_error_relative = 0,
         post_cast_backward_residual = 0;
};
// Requested sensitivity budget is caller/consumer policy, not a global epsilon.
SolveResult solve(const Factorization &, std::span<const double> rhs,
                  double maximum_forward_sensitivity);
// One retained-factor forward triangular solve, L z = rhs. Wide outputs are
// intentional internal numerical coordinates, not rounded binary64 observables.
// The diagnostic concerns stored L, not a covariance forward-error certificate.
struct WhiteningResult {
  Status status = Status::invalid_input;
  std::vector<long double> value;
  double backward_residual = 0, arithmetic_rounding_estimate = 0;
};
std::optional<std::size_t> whitening_payload_bound(std::size_t n) noexcept;
// Borrow rhs/factor; no factor copy or matrix readback. No vector payload on
// failure. maximum_elements bounds n; payload bounds scratch/result only.
WhiteningResult whiten(const Factorization &, std::span<const double> rhs,
                       std::size_t maximum_elements,
                       std::size_t maximum_payload_bytes,
                       double maximum_backward_error);
WhiteningResult whiten(const Factorization &, std::span<const long double> rhs,
                       std::size_t maximum_elements,
                       std::size_t maximum_payload_bytes,
                       double maximum_backward_error);
// Multiply the retained triangle L*z, without copying L or reading original C.
// Errors concern accumulation against stored L, not covariance reconstruction.
struct ColouringResult {
  Status status = Status::invalid_input;
  std::vector<long double> value, absolute_error_estimates;
};
std::optional<std::size_t> colouring_payload_bound(std::size_t) noexcept;
ColouringResult colour(const Factorization &, std::span<const long double>,
                       std::size_t maximum_elements,
                       std::size_t maximum_payload_bytes,
                       double maximum_scaled_arithmetic_error);
} // namespace irred::numerics
