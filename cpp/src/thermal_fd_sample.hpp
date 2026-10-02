#pragma once
#include <cmath>
namespace irred::cosmology::detail {
// One physical FD equation and operation order for both numerical algorithms.
struct ThermalFDState { long double x, energy, occupation; };
inline ThermalFDState thermal_fd_state(double q, long double y, long double scale) {
  const long double x=q;
  const long double energy=std::hypot(x/scale,y/scale);
  const long double e=std::exp(-x), fd=e/(1+e);
  return {x,energy,fd};
}
inline double thermal_fd_value(const ThermalFDState &s, bool pressure) {
  if(s.x==0) return 0;
  const long double x=s.x, energy=s.energy, fd=s.occupation;
  const long double result=pressure ? x*x*x*x*fd/(3*energy) : x*x*energy*fd;
  return static_cast<double>(result);
}
} // namespace irred::cosmology::detail
