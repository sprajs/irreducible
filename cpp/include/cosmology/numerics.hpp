#pragma once
#include <cstddef>
#include <span>
#include <vector>
namespace cosmology::numerics {
enum class Status { ok, invalid_input, nonfinite_input, overflow, work_limit, outside_domain, singular, not_positive_definite, conditioning_budget_exceeded };
// error_estimate is an empirical estimator or arithmetic diagnostic, never a proven bound.
struct ScalarResult { Status status=Status::invalid_input; double value=0; double error_estimate=0; std::size_t evaluations=0; };
ScalarResult compensated_sum(std::span<const double>) noexcept;
ScalarResult log_sum_exp(std::span<const double>) noexcept; // nonempty finite input
ScalarResult log1p_checked(double) noexcept; // x > -1
ScalarResult expm1_checked(double) noexcept; // finite result required
ScalarResult log_gamma_positive(double) noexcept; // bounded first slice [0.125,100]
ScalarResult interpolate_linear(std::span<const double> x,std::span<const double> y,double query) noexcept; // strictly increasing, no extrapolation
struct IntegrationPolicy { double absolute_tolerance; double relative_tolerance; std::size_t max_evaluations; unsigned max_depth; };
using Integrand = double (*)(double, const void*);
// Caller supplies smooth resolved interval; split known discontinuities before calls.
// Callback must not throw. Narrow unobserved features can defeat empirical estimate.
ScalarResult integrate(Integrand,const void*,double lower,double upper,IntegrationPolicy);
class Factorization {
public:
 Status status() const noexcept { return status_; }
 std::size_t size() const noexcept { return n_; }
 double log_determinant() const noexcept { return log_determinant_; }
 double condition_estimate_inf() const noexcept { return condition_estimate_inf_; }
private:
 Status status_=Status::invalid_input;
 std::size_t n_=0;
 std::vector<double> lower_,original_;
 double log_determinant_=0,norm_inf_=0,condition_estimate_inf_=0;
 friend Factorization cholesky(std::span<const double>,std::size_t,std::size_t);
 friend struct FactorAccess;
};
// Finite exactly symmetric row-major matrix; checked shape/resource limit, no repair/jitter.
// Initial qualification covers tested n<=16 only; size support is not qualification.
Factorization cholesky(std::span<const double>,std::size_t n,std::size_t maximum_elements);
struct SolveResult { Status status=Status::invalid_input; std::vector<double> value; double backward_residual=0; double estimated_forward_sensitivity=0; };
// Requested sensitivity budget is caller/consumer policy, not a global epsilon.
SolveResult solve(const Factorization&,std::span<const double> rhs,double maximum_forward_sensitivity);
} // namespace cosmology::numerics
