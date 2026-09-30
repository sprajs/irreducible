#pragma once
#include <array>
#include <string_view>
namespace irred::test_reference {
// Four explicitly frozen radiation-free flat CPL SN points, selected zHD>0.01,
// zHD radial integral and zHEL luminosity prefactor; H0 absorbed in free
// offset. Fine native wide-Cholesky/profile scores compared to independent
// fixed GL8 z32/64 + longdouble LDLT on the SAME canonical selected binary64
// covariance. Shared libm/equations and owner authorship disclosed; separate
// peer tests. Deterministic score budget1e-6: reference2e-7, factor/solve3e-7,
// background5e-7; magnitude1e-8; explicit sensitivity guard1e-10.
// Optional original-data tests MUST SHA-check both sources before/after.
// No fitting/posterior, broad domain, or historical campaign reproduction
// claim.
inline constexpr std::string_view w01_cpl_table_sha256 =
    "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8";
inline constexpr std::string_view w01_cpl_covariance_sha256 =
    "abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc";
struct W01CplPoint {
  double omega_m, constant_q, w0, wa, relative_score;
};
inline constexpr std::array<W01CplPoint, 4> w01_cpl_points{
    {{.3, 0, -1, 0, -703.03398247940595},
     {0, 0, -0.66666666666666663, 0, -732.92476116995874},
     {.3, 0, -.9, .4, -703.21948165929484},
     {.3, 0, -1.1, -.4, -715.39258058819257}}};
} // namespace irred::test_reference
