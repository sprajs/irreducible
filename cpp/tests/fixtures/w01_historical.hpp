// Frozen read-only historical comparison; not a production dependency.
// Python3.12.14 / NumPy2.2.6 / SciPy1.15.3 / pandas2.2.3.
// lib/cosmology.py SHA256
// 6c1d9c8cf13e7126aa31b6dc7976b39ee6c34b030dc4735cc60f3ffbf2209e1b. zHD>0.01,
// 1590 ordered rows; zHD expansion integral, (1+zHEL) prefactor.
// mu=5log10((1+zHEL)*I(zHD)); free additive offset; H0 absorbed.
// Unchanged GL48 integration and inverse-projection full score; stable score
// separately uses SciPy Cholesky solves. Spline40 compressed is approximate,
// and fails the 2e-7 precise-reference allocation at some frozen points.
#pragma once
#include <array>
#include <string_view>
namespace irred::test_reference {
inline constexpr std::string_view w01_table_sha256 =
    "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8";
inline constexpr std::string_view w01_cov_sha256 =
    "abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc";
inline constexpr std::string_view w01_selected_f64le_sha256 =
    "64339b81cfd28998ca04f589394c9209b44c93dd63112704eb21bdeb6f300ee4";
struct HistoricalPoint {
  bool lcdm;
  double parameter, full_score, stable_score, compressed_score;
};
inline constexpr std::array<HistoricalPoint, 7> w01_historical{{
    {true, 0.0, -1029.9347255796138, -1029.9347255770522, -1029.9347256074911},
    {true, 0.3, -703.0339824813782, -703.0339824794266, -703.0339824816658},
    {true, 1.0, -1050.5727707172978, -1050.5727707147396, -1050.5727701741625},
    {false, -1.0, -1029.9347255796138, -1029.9347255770522,
     -1029.9347256074911},
    {false, -0.5, -732.9247611716803, -732.9247611699597, -732.9247610890417},
    {false, 0.0, -754.193810732856, -754.1938107296344, -754.1938104647393},
    {false, 0.5, -1050.5727707173141, -1050.5727707147396, -1050.5727701741625},
}};
} // namespace irred::test_reference
