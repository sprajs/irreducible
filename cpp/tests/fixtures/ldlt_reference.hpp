#pragma once
// Test-only independent unit-lower LDL^T recurrence, long-double arithmetic.
// Inputs are CANONICAL binary64 promoted to long double, not reinterpreted
// exact decimal source. D_i=A_ii-sum_{k<i} L_ik^2 D_k;
// L_ij=(A_ij-sum_{k<j} L_ik L_jk D_k)/D_j. No repair/pivot/jitter.
// Production uses a separately implemented binary64 Cholesky kernel.
// This is empirical reference arithmetic, not interval certification; logs
// share the system libm. n<=1701 bounded developer comparison only.
#include <cmath>
#include <span>
#include <vector>
namespace irred::test_reference {
struct LDLT {
  bool valid = false;
  std::size_t n = 0;
  std::vector<long double> lower, diagonal;
  long double logdet = 0;
};
inline LDLT factor(std::span<const double> a, std::size_t n) {
  LDLT out;
  if (!n || n > 1701 || a.size() != n * n)
    return out;
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      if (!std::isfinite(a[i * n + j]) || a[i * n + j] != a[j * n + i])
        return out;
  out.n = n;
  out.lower.assign(n * n, 0);
  out.diagonal.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      long double sum = a[i * n + j];
      for (std::size_t k = 0; k < j; ++k)
        sum -= out.lower[i * n + k] * out.diagonal[k] * out.lower[j * n + k];
      out.lower[i * n + j] = sum / out.diagonal[j];
    }
    long double d = a[i * n + i];
    for (std::size_t k = 0; k < i; ++k)
      d -= out.lower[i * n + k] * out.lower[i * n + k] * out.diagonal[k];
    if (!(d > 0) || !std::isfinite(d))
      return out;
    out.diagonal[i] = d;
    out.lower[i * n + i] = 1;
    out.logdet += std::log(d);
  }
  out.valid = true;
  return out;
}
template <class Real>
inline std::vector<long double> solve_work(const LDLT &f,
                                           std::span<const Real> b) {
  if (!f.valid || b.size() != f.n || f.lower.size() != f.n * f.n ||
      f.diagonal.size() != f.n)
    return {};
  auto n = f.n;
  std::vector<long double> y(n), x(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(b[i]))
      return {};
    y[i] = b[i];
    for (std::size_t j = 0; j < i; ++j)
      y[i] -= f.lower[i * n + j] * y[j];
  }
  for (std::size_t i = 0; i < n; ++i)
    x[i] = y[i] / f.diagonal[i];
  for (std::size_t i = n; i-- > 0;)
    for (std::size_t j = i + 1; j < n; ++j)
      x[i] -= f.lower[j * n + i] * x[j];
  return x;
}
inline std::vector<long double> solve(const LDLT &f,
                                      std::span<const double> b) {
  return solve_work<double>(f, b);
}
inline std::vector<long double>
solve_longdouble(const LDLT &f, std::span<const long double> b) {
  return solve_work<long double>(f, b);
}
inline long double quadratic(std::span<const double> b,
                             std::span<const long double> x) {
  long double q = 0;
  for (std::size_t i = 0; i < b.size(); ++i)
    q += static_cast<long double>(b[i]) * x[i];
  return q;
}
} // namespace irred::test_reference
