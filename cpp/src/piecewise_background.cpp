#include "irred/piecewise_background.hpp"
#include <algorithm>
#include <cmath>
namespace irred::cosmology::detail {
long double exprel(long double t) {
  if (t == 0)
    return 1;
  if (std::abs(t) > .0001L)
    return std::expm1(t) / t;
  long double sum = 1, term = 1;
  for (int k = 1; k <= 12; ++k) {
    term *= t / (k + 1);
    sum += term;
  }
  return sum;
}
std::size_t piecewise_visits(double z) noexcept {
  std::size_t n = 0;
  for (std::size_t i = 0; i < 5; ++i)
    if (z > piecewise_q_edges[i])
      ++n;
  return n;
}
PiecewiseIntegral piecewise_integral(const FixedFiveBinQ &p,
                                     const std::array<long double, 5> &start,
                                     double z, bool radial, bool clock) {
  PiecewiseIntegral out;
  out.bin = std::min(std::size_t(4),
                     std::size_t(std::upper_bound(piecewise_q_edges.begin() + 1,
                                                  piecewise_q_edges.end(), z) -
                                 piecewise_q_edges.begin() - 1));
  for (std::size_t b = 0; b < piecewise_visits(z); ++b) {
    const auto lo = (long double)piecewise_q_edges[b];
    const auto hi =
        std::min((long double)z, (long double)piecewise_q_edges[b + 1]);
    const auto L = std::log1p((hi - lo) / (1 + lo));
    const auto q = (long double)p.q[b];
    if (radial)
      out.radial += (1 + lo) / start[b] * L * exprel(-q * L);
    if (clock)
      out.clock += 1 / start[b] * L * exprel(-(1 + q) * L);
  }
  const auto lo = (long double)piecewise_q_edges[out.bin];
  const auto L = std::log1p(((long double)z - lo) / (1 + lo));
  out.E = start[out.bin] * std::exp((1 + (long double)p.q[out.bin]) * L);
  return out;
}
} // namespace irred::cosmology::detail
