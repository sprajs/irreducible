#include "thermal_cc_sample.hpp"
#include <cmath>
namespace irred::cosmology::detail {
ThermalCCSample thermal_cc_sample(double q, long double y, long double scale,
                                 bool need_pressure) noexcept {
  if (q == 0) return {0, 0};
  const long double x = q, energy = std::hypot(x / scale, y / scale);
  const long double e = std::exp(-x), fd = e / (1 + e);
  return {static_cast<double>(x*x*energy*fd),
          need_pressure ? static_cast<double>(x*x*x*x*fd/(3*energy)) : 0};
}
} // namespace irred::cosmology::detail
