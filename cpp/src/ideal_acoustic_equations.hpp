#pragma once
#include <span>
namespace irred::cosmology::detail {
// Model-specific variables, not a replacement background or photon hierarchy.
// L and Rbg are the actual supplied shared-epoch coefficients. Exact compatible
// source algebra has Rbg=0; the numerical defect is never silently deleted.
struct IdealAcousticCoefficients {
  long double x2, L, F, B, loading_over_one_plus_loading, sound_speed_squared;
  long double acceleration_defect;
};
inline void ideal_acoustic_derivative(std::span<const long double,5> y,
                                     const IdealAcousticCoefficients &e,
                                     std::span<long double,5> dy) noexcept {
  dy[0]=-e.x2*y[3]+4.5L*e.B*y[2]+3*e.acceleration_defect*y[3];
  dy[1]=e.x2*y[2];
  dy[2]=(e.L-e.loading_over_one_plus_loading)*y[2]+
        e.sound_speed_squared*(y[0]-y[1]);
  dy[3]=(e.L-1)*y[3]+y[4];
  dy[4]=-y[4]+1.5L*e.F*y[3]+1.5L*e.B*y[2];
}
} // namespace irred::cosmology::detail
