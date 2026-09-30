#include "irred/numerics.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numeric>
#include <type_traits>
namespace irred::numerics {
namespace {
ScalarResult finite_result(long double x, double error = 0,
                           std::size_t evaluations = 0) noexcept {
  if (!std::isfinite(x) || std::abs(x) > std::numeric_limits<double>::max())
    return {Status::overflow, 0, error, evaluations};
  return {Status::ok, static_cast<double>(x), error, evaluations};
}
} // namespace
ScalarResult compensated_sum(std::span<const double> xs) noexcept {
  long double total = 0, correction = 0, absolute = 0;
  for (double x : xs) {
    if (!std::isfinite(x))
      return {Status::nonfinite_input};
    const long double next = total + x;
    correction += std::abs(total) >= std::abs(x) ? (total - next) + x
                                                 : (x - next) + total;
    total = next;
    absolute += std::abs(static_cast<long double>(x));
  }
  return finite_result(
      total + correction,
      static_cast<double>(absolute * std::numeric_limits<double>::epsilon()),
      xs.size());
}
ScalarResult log_sum_exp(std::span<const double> xs) noexcept {
  if (xs.empty())
    return {Status::invalid_input};
  double maximum = -std::numeric_limits<double>::infinity();
  for (double x : xs) {
    if (!std::isfinite(x))
      return {Status::nonfinite_input};
    maximum = std::max(maximum, x);
  }
  long double sum = 0;
  for (double x : xs)
    sum += std::exp(static_cast<long double>(x) - maximum);
  return finite_result(static_cast<long double>(maximum) + std::log(sum), 0,
                       xs.size());
}
ScalarResult log1p_checked(double x) noexcept {
  if (!std::isfinite(x))
    return {Status::nonfinite_input};
  if (x <= -1)
    return {Status::outside_domain};
  return finite_result(std::log1p(x));
}
ScalarResult expm1_checked(double x) noexcept {
  if (!std::isfinite(x))
    return {Status::nonfinite_input};
  return finite_result(std::expm1(x));
}
ScalarResult log_gamma_positive(double x) noexcept {
  if (!std::isfinite(x))
    return {Status::nonfinite_input};
  if (x < .125 || x > 100)
    return {Status::outside_domain};
  return finite_result(std::lgamma(x));
}
ScalarResult interpolate_linear(std::span<const double> x,
                                std::span<const double> y,
                                double query) noexcept {
  if (x.size() != y.size() || x.size() < 2)
    return {Status::invalid_input};
  if (!std::isfinite(query))
    return {Status::nonfinite_input};
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
      return {Status::nonfinite_input};
    if (i && x[i] <= x[i - 1])
      return {Status::invalid_input};
  }
  if (query < x.front() || query > x.back())
    return {Status::outside_domain};
  auto upper = std::lower_bound(x.begin(), x.end(), query);
  const auto i = static_cast<std::size_t>(upper - x.begin());
  if (*upper == query)
    return {Status::ok, y[i]};
  const auto j = i - 1;
  const long double fraction = (static_cast<long double>(query) - x[j]) /
                               (static_cast<long double>(x[i]) - x[j]);
  return finite_result((1 - fraction) * y[j] + fraction * y[i]);
}
namespace {
struct Adaptive {
  Integrand function;
  const void *context;
  IntegrationPolicy policy;
  std::size_t evaluations = 0;
  Status status = Status::ok;
  double evaluate(double x) {
    if (status != Status::ok)
      return 0;
    if (evaluations >= policy.max_evaluations) {
      status = Status::work_limit;
      return 0;
    }
    ++evaluations;
    double y = function(x, context);
    if (!std::isfinite(y))
      status = Status::nonfinite_input;
    return y;
  }
  struct Part {
    long double value = 0, error = 0;
  };
  Part recurse(double a, double b, double fa, double fm, double fb,
               long double coarse, long double allowance, unsigned depth) {
    if (status != Status::ok)
      return {};
    const double m = std::midpoint(a, b), left = std::midpoint(a, m),
                 right = std::midpoint(m, b);
    if (left == a || left == m || right == m || right == b) {
      status = Status::work_limit;
      return {};
    }
    const auto fl = evaluate(left), fr = evaluate(right);
    if (status != Status::ok)
      return {};
    const long double first =
        (static_cast<long double>(m) - a) * (fa + 4.L * fl + fm) / 6;
    const long double second =
        (static_cast<long double>(b) - m) * (fm + 4.L * fr + fb) / 6;
    const long double refined = first + second, delta = refined - coarse,
                      error = std::abs(delta) / 15;
    if (!std::isfinite(refined) || !std::isfinite(error)) {
      status = Status::overflow;
      return {};
    }
    if (error <= allowance)
      return {refined + delta / 15, error};
    if (depth >= policy.max_depth) {
      status = Status::work_limit;
      return {};
    }
    const auto l = recurse(a, m, fa, fl, fm, first, allowance / 2, depth + 1);
    const auto r = recurse(m, b, fm, fr, fb, second, allowance / 2, depth + 1);
    return {l.value + r.value, l.error + r.error};
  }
};
} // namespace
ScalarResult integrate(Integrand f, const void *context, double a, double b,
                       IntegrationPolicy p) {
  if (!f || !std::isfinite(a) || !std::isfinite(b) || a > b ||
      !std::isfinite(p.absolute_tolerance) ||
      !std::isfinite(p.relative_tolerance) || p.absolute_tolerance < 0 ||
      p.relative_tolerance < 0 ||
      (p.absolute_tolerance == 0 && p.relative_tolerance == 0) ||
      p.max_evaluations < 3 || p.max_depth > 60)
    return {Status::invalid_input};
  if (a == b)
    return {Status::ok, 0, 0, 0};
  Adaptive runner{f, context, p};
  double m = std::midpoint(a, b);
  const auto fa = runner.evaluate(a), fm = runner.evaluate(m),
             fb = runner.evaluate(b);
  if (runner.status != Status::ok)
    return {runner.status, 0, 0, runner.evaluations};
  const long double coarse =
      (static_cast<long double>(b) - a) * (fa + 4.L * fm + fb) / 6;
  if (!std::isfinite(coarse))
    return {Status::overflow, 0, 0, runner.evaluations};
  long double allowance =
      p.absolute_tolerance + p.relative_tolerance * std::abs(coarse);
  for (;;) {
    const auto r = runner.recurse(a, b, fa, fm, fb, coarse, allowance, 0);
    if (runner.status != Status::ok)
      return {runner.status, 0, 0, runner.evaluations};
    // Relative accuracy belongs to the final integral, including cancellation.
    // The initial coarse value supplies a work proposal, never an acceptance
    // scale.
    const long double final_budget =
        p.absolute_tolerance + p.relative_tolerance * std::abs(r.value);
    if (r.error <= final_budget)
      return finite_result(r.value, static_cast<double>(r.error),
                           runner.evaluations);
    allowance = std::min(allowance / 2, final_budget);
  }
}
struct WideRangeGuard {
  std::fenv_t saved{};
  bool supported = false;
  WideRangeGuard() { supported = std::feholdexcept(&saved) == 0; }
  ~WideRangeGuard() {
    if (supported)
      std::fesetenv(&saved);
  }
  Status status() const {
    auto flags = std::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID |
                                   FE_DIVBYZERO);
    if (flags & FE_UNDERFLOW)
      return Status::outside_domain;
    if (flags)
      return Status::overflow;
    return Status::ok;
  }
};
struct FactorAccess {
  static bool wide_environment() {
    return std::numeric_limits<long double>::radix == 2 &&
           std::numeric_limits<long double>::digits >= 64 &&
           std::numeric_limits<long double>::min_exponent <
               std::numeric_limits<double>::min_exponent &&
           std::numeric_limits<long double>::max_exponent >
               std::numeric_limits<double>::max_exponent &&
           std::fegetround() == FE_TONEAREST;
  }
  template <class Real> static auto &lower(Factorization &f) {
    if constexpr (std::is_same_v<Real, double>)
      return f.lower_;
    else
      return f.wide_lower_;
  }
  template <class Real> static const auto &lower(const Factorization &f) {
    if constexpr (std::is_same_v<Real, double>)
      return f.lower_;
    else
      return f.wide_lower_;
  }
  template <class Real> static bool valid(const Factorization &f) {
    const auto &l = lower<Real>(f);
    if (f.status_ != Status::ok || !f.n_ ||
        f.n_ > std::numeric_limits<std::size_t>::max() / f.n_ ||
        l.size() != f.n_ * f.n_ || f.original_.size() != l.size())
      return false;
    if (!std::isfinite(f.norm_inf_) ||
        !std::isfinite(f.condition_estimate_inf_) ||
        f.condition_estimate_inf_ < 1)
      return false;
    for (auto x : l)
      if (!std::isfinite(x))
        return false;
    for (double x : f.original_)
      if (!std::isfinite(x))
        return false;
    for (std::size_t i = 0; i < f.n_; ++i)
      if (l[i * f.n_ + i] <= 0)
        return false;
    return true;
  }
  template <class Real>
  static std::vector<Real> triangular(const Factorization &f,
                                      std::span<const double> rhs) {
    const auto &l = lower<Real>(f);
    std::vector<Real> x(rhs.begin(), rhs.end());
    const auto n = f.n_;
    for (std::size_t i = 0; i < n; ++i) {
      long double sum = x[i];
      for (std::size_t j = 0; j < i; ++j)
        sum -= static_cast<long double>(l[i * n + j]) * x[j];
      x[i] = static_cast<Real>(sum / l[i * n + i]);
      if constexpr (!std::is_same_v<Real, double>)
        if (std::fpclassify(x[i]) == FP_SUBNORMAL)
          return {};
    }
    for (std::size_t k = n; k > 0; --k) {
      const auto i = k - 1;
      long double sum = x[i];
      for (std::size_t j = i + 1; j < n; ++j)
        sum -= static_cast<long double>(l[i * n + j]) * x[j];
      x[i] = static_cast<Real>(sum / l[i * n + i]);
      if constexpr (!std::is_same_v<Real, double>)
        if (std::fpclassify(x[i]) == FP_SUBNORMAL)
          return {};
    }
    return x;
  }
  template <class Real>
  static SolveResult solve_impl(const Factorization &f,
                                std::span<const double> rhs, double budget) {
    if constexpr (!std::is_same_v<Real, double>)
      if (!wide_environment())
        return {Status::outside_domain, {}};
    if (!valid<Real>(f) || rhs.size() != f.n_ || !std::isfinite(budget) ||
        budget <= 0)
      return {Status::invalid_input, {}};
    for (double x : rhs)
      if (!std::isfinite(x))
        return {Status::nonfinite_input, {}};
    auto work = triangular<Real>(f, rhs);
    if (work.size() != f.n_)
      return {Status::outside_domain, {}};
    SolveResult out;
    out.value.resize(f.n_);
    double normx = 0, normb = 0;
    long double residual = 0, pre_residual = 0, pre_norm = 0, cast = 0;
    for (std::size_t i = 0; i < f.n_; ++i) {
      if constexpr (!std::is_same_v<Real, double>)
        if (std::fpclassify(work[i]) == FP_SUBNORMAL)
          return {Status::outside_domain, {}};
      if (!std::isfinite(work[i]) ||
          std::abs(static_cast<long double>(work[i])) >
              std::numeric_limits<double>::max())
        return {Status::overflow, {}};
      out.value[i] = static_cast<double>(work[i]);
      if constexpr (!std::is_same_v<Real, double>)
        if (work[i] != 0 && (out.value[i] == 0 ||
                             std::fpclassify(out.value[i]) == FP_SUBNORMAL))
          return {Status::outside_domain, {}};
      normx = std::max(normx, std::abs(out.value[i]));
      normb = std::max(normb, std::abs(rhs[i]));
      pre_norm =
          std::max(pre_norm, std::abs(static_cast<long double>(work[i])));
      cast = std::max(
          cast, std::abs(static_cast<long double>(work[i]) - out.value[i]));
      long double sum = -rhs[i], pre_sum = -rhs[i];
      for (std::size_t j = 0; j < f.n_; ++j) {
        sum += static_cast<long double>(f.original_[i * f.n_ + j]) *
               static_cast<double>(work[j]);
        pre_sum +=
            static_cast<long double>(f.original_[i * f.n_ + j]) * work[j];
      }
      residual = std::max(residual, std::abs(sum));
      pre_residual = std::max(pre_residual, std::abs(pre_sum));
    }
    const long double scale =
        static_cast<long double>(f.norm_inf_) * normx + normb;
    out.backward_residual =
        scale == 0 ? 0 : static_cast<double>(residual / scale);
    out.post_cast_backward_residual = out.backward_residual;
    if constexpr (std::is_same_v<Real, double>) {
      out.estimated_forward_sensitivity =
          f.condition_estimate_inf_ *
          (out.backward_residual +
           f.n_ * std::numeric_limits<double>::epsilon());
      out.pre_cast_backward_residual = out.backward_residual;
      out.pre_cast_sensitivity_estimate = out.estimated_forward_sensitivity;
    } else {
      const auto pre_scale = f.wide_norm_inf_ * pre_norm + normb;
      const auto eta = pre_scale == 0 ? 0 : pre_residual / pre_scale;
      const auto sensitivity =
          f.wide_condition_inf_ *
          (eta + f.n_ * std::numeric_limits<long double>::epsilon());
      out.pre_cast_backward_residual = static_cast<double>(eta);
      out.pre_cast_sensitivity_estimate = static_cast<double>(sensitivity);
      out.output_rounding_error_inf = static_cast<double>(cast);
      out.output_rounding_error_relative =
          normx == 0 ? 0 : static_cast<double>(cast / normx);
      out.estimated_forward_sensitivity = static_cast<double>(
          normx == 0 ? sensitivity
                     : sensitivity * (pre_norm / normx) + cast / normx);
    }
    out.status = out.estimated_forward_sensitivity <= budget
                     ? Status::ok
                     : Status::conditioning_budget_exceeded;
    return out;
  }
  template <class Real>
  static Factorization factor(std::span<const double> matrix, std::size_t n,
                              std::size_t cap, Arithmetic arithmetic) {
    Factorization f;
    f.arithmetic_ = arithmetic;
    if constexpr (!std::is_same_v<Real, double>)
      if (!wide_environment()) {
        f.status_ = Status::outside_domain;
        return f;
      }
    if (!n || n > std::numeric_limits<std::size_t>::max() / n ||
        matrix.size() != n * n)
      return f;
    if (n * n > cap) {
      f.status_ = Status::work_limit;
      return f;
    }
    for (double x : matrix)
      if (!std::isfinite(x)) {
        f.status_ = Status::nonfinite_input;
        return f;
      }
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < i; ++j)
        if (matrix[i * n + j] != matrix[j * n + i])
          return f;
    f.n_ = n;
    f.original_.assign(matrix.begin(), matrix.end());
    auto &l = lower<Real>(f);
    l.assign(n * n, 0);
    long double logdet = 0;
    double norm = 0;
    long double wide_norm = 0;
    for (std::size_t i = 0; i < n; ++i) {
      long double row = 0;
      for (std::size_t j = 0; j < n; ++j)
        row += std::abs(static_cast<long double>(matrix[i * n + j]));
      if (row > std::numeric_limits<double>::max()) {
        f.status_ = Status::overflow;
        return f;
      }
      norm = std::max(norm, static_cast<double>(row));
      wide_norm = std::max(wide_norm, row);
      for (std::size_t j = 0; j <= i; ++j) {
        long double sum = matrix[i * n + j];
        for (std::size_t k = 0; k < j; ++k)
          sum -= static_cast<long double>(l[i * n + k]) * l[j * n + k];
        if (i == j) {
          if (sum <= 0) {
            f.status_ = Status::not_positive_definite;
            return f;
          }
          l[i * n + j] = static_cast<Real>(std::sqrt(sum));
          logdet += std::log(sum);
        } else
          l[i * n + j] = static_cast<Real>(sum / l[j * n + j]);
        if constexpr (!std::is_same_v<Real, double>)
          if (std::fpclassify(l[i * n + j]) == FP_SUBNORMAL) {
            f.status_ = Status::outside_domain;
            return f;
          }
        if (!std::isfinite(l[i * n + j]) || (i == j && l[i * n + j] == 0)) {
          f.status_ = Status::overflow;
          return f;
        }
      }
    }
    // Reuse the already allocated upper triangle as the transpose of L.
    // Back substitution keeps identical operands and ascending j order while
    // reading contiguous rows; both arithmetic policies share this layout.
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j)
        l[i * n + j] = l[j * n + i];
    f.norm_inf_ = norm;
    f.wide_norm_inf_ = wide_norm;
    f.log_determinant_ = static_cast<double>(logdet);
    std::vector<long double> row_sums(n, 0);
    std::vector<double> basis(n, 0);
    for (std::size_t j = 0; j < n; ++j) {
      basis[j] = 1;
      const auto col = triangular<Real>(f, basis);
      if (col.size() != n) {
        f.status_ = Status::outside_domain;
        return f;
      }
      basis[j] = 0;
      for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(col[i])) {
          f.status_ = Status::overflow;
          return f;
        }
        row_sums[i] += std::abs(static_cast<long double>(col[i]));
      }
    }
    const long double condition =
        (std::is_same_v<Real, double> ? static_cast<long double>(norm)
                                      : wide_norm) *
        (*std::max_element(row_sums.begin(), row_sums.end()));
    if (!std::isfinite(condition) ||
        condition > std::numeric_limits<double>::max()) {
      f.status_ = Status::overflow;
      return f;
    }
    f.condition_estimate_inf_ = std::max(1., static_cast<double>(condition));
    f.wide_condition_inf_ = std::max(1.L, condition);
    f.status_ = Status::ok;
    return f;
  }
};
Factorization cholesky(std::span<const double> matrix, std::size_t n,
                       std::size_t cap, Arithmetic arithmetic) {
  if (arithmetic == Arithmetic::binary64_legacy_v1)
    return FactorAccess::factor<double>(matrix, n, cap, arithmetic);
  if (arithmetic == Arithmetic::longdouble_cpu_v1) {
    WideRangeGuard guard;
    if (!guard.supported) {
      Factorization out;
      out.arithmetic_ = arithmetic;
      out.status_ = Status::outside_domain;
      return out;
    }
    auto out = FactorAccess::factor<long double>(matrix, n, cap, arithmetic);
    auto status = guard.status();
    if (status != Status::ok)
      out.status_ = status;
    return out;
  }
  Factorization invalid;
  invalid.arithmetic_ = arithmetic;
  return invalid;
}
SolveResult solve(const Factorization &f, std::span<const double> rhs,
                  double budget) {
  if (f.arithmetic() == Arithmetic::longdouble_cpu_v1) {
    WideRangeGuard guard;
    if (!guard.supported)
      return {Status::outside_domain, {}};
    auto out = FactorAccess::solve_impl<long double>(f, rhs, budget);
    auto status = guard.status();
    if (status != Status::ok)
      return {status, {}};
    return out;
  }
  if (f.arithmetic() != Arithmetic::binary64_legacy_v1)
    return {Status::invalid_input, {}};
  return FactorAccess::solve_impl<double>(f, rhs, budget);
}
} // namespace irred::numerics
