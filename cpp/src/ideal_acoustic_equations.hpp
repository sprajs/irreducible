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
// Each component returns the original rounded Wide expression. This lets a
// signed-response consumer combine it directly with its source component
// without materializing another physical-state vector or changing grouping.
template<unsigned I>
inline long double ideal_acoustic_derivative_coordinate(
    std::span<const long double,5> y,
    const IdealAcousticCoefficients &e) noexcept {
  static_assert(I<5);
  if constexpr (I==0)
    return -e.x2*y[3]+4.5L*e.B*y[2]+3*e.acceleration_defect*y[3];
  else if constexpr (I==1)
    return e.x2*y[2];
  else if constexpr (I==2)
    return (e.L-e.loading_over_one_plus_loading)*y[2]+
           e.sound_speed_squared*(y[0]-y[1]);
  else if constexpr (I==3)
    return (e.L-1)*y[3]+y[4];
  else
    return -y[4]+1.5L*e.F*y[3]+1.5L*e.B*y[2];
}
inline void ideal_acoustic_derivative(std::span<const long double,5> y,
                                     const IdealAcousticCoefficients &e,
                                     std::span<long double,5> dy) noexcept {
  dy[0]=ideal_acoustic_derivative_coordinate<0>(y,e);
  dy[1]=ideal_acoustic_derivative_coordinate<1>(y,e);
  dy[2]=ideal_acoustic_derivative_coordinate<2>(y,e);
  dy[3]=ideal_acoustic_derivative_coordinate<3>(y,e);
  dy[4]=ideal_acoustic_derivative_coordinate<4>(y,e);
}
} // namespace irred::cosmology::detail
