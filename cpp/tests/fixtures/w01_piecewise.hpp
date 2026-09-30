#pragma once
#include <array>
namespace irred::test_reference {
// Same original canonical binary64 table/covariance and strict zHD>.01 order
// as w01_historical.hpp. Analytic five-bin provider + retained wide Cholesky
// compared to split scale-factor GL8 32/64/128 and long-double LDLT. Shared
// libm and owner ancestry: empirical reference, not an interval certificate.
// Six declared points only; no priors/posterior/campaign qualification.
// Absolute score allocation1e-6: reference2e-7, solve3e-7, background5e-7;
// magnitude1e-8, explicit longdouble_cpu_v1 guard1e-10, segment cap50000.
inline constexpr std::array<std::array<double, 5>, 6> w01_piecewise_q{
    {{0, 0, 0, 0, 0},
     {-1, -1, -1, -1, -1},
     {.5, .5, .5, .5, .5},
     {-.4, -.4, -.2, .1, .3},
     {-1, 0, -1, 0, -1},
     {-3, 2, -3, 2, -3}}};
inline constexpr std::array<double, 6> w01_piecewise_scores{
    -754.19381072963506, -1029.9347255770519, -1050.5727707147407,
    -702.29430438551958, -746.32743999130741, -982.69593471471467};
// Read-only original lib/cosmology.py qbins/Pantheon fixedpoint route:
// Python3.12.14/NumPy2.2.6/SciPy1.15.3/pandas2.2.3, BLAS1, no priors.
// Source lib/cosmology.py SHA256
// 6c1d9c8cf13e7126aa31b6dc7976b39ee6c34b030dc4735cc60f3ffbf2209e1b Original
// table/covariance SHA256 pinned in w01_historical.hpp. Historical compressed
// spline is separately approximate: five points fail the fixed2e-7 reference
// allocation; retain discrepancies, do not promote.
inline constexpr std::array<double, 6> w01_piecewise_historical_stable{
    -754.193810729634, -1029.93472557705,  -1050.5727707147403,
    -702.294304385518, -746.3274399913068, -982.6959347147131};
inline constexpr std::array<double, 6> w01_piecewise_historical_full_inverse{
    -754.193810732776,  -1029.9347255796993, -1050.5727707173892,
    -702.2943043883897, -746.3274399938423,  -982.6959347171035};
inline constexpr std::array<double, 6> w01_piecewise_historical_compressed{
    -754.1938104647413, -1029.9347256074907, -1050.5727701741653,
    -702.2943055977524, -746.3274375139254,  -982.69592412781};
} // namespace irred::test_reference
