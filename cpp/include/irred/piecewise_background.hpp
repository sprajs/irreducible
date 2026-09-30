#pragma once
// Internal analytic kernel used by the current compiled Expansion. The fixed
// edges are model identity, not runtime input or a separate provider interface.
#include "irred/background.hpp"
namespace irred::cosmology::detail {
struct PiecewiseIntegral {
  long double E = 0, radial = 0, clock = 0;
  std::size_t bin = 0;
};
long double exprel(long double);
PiecewiseIntegral piecewise_integral(const FixedFiveBinQ &,
                                     const std::array<long double, 5> &, double,
                                     bool, bool);
std::size_t piecewise_visits(double) noexcept;
} // namespace irred::cosmology::detail
