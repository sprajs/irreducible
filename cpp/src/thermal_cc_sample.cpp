#include "thermal_cc_sample.hpp"
#include "thermal_fd_sample.hpp"
namespace irred::cosmology::detail {
ThermalCCSample thermal_cc_sample(double q, long double y, long double scale,
                                 bool need_pressure) noexcept {
  if (q == 0) return {0, 0};
  const auto state=thermal_fd_state(q,y,scale);
  return {thermal_fd_value(state,false),
          need_pressure ? thermal_fd_value(state,true) : 0};
}
} // namespace irred::cosmology::detail
