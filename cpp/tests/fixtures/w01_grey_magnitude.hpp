#pragma once
#include <array>
#include <string_view>
// Extracted unchanged historical IntegratedLuminosity fixed-amplitude facts.
// Same current Pantheon1590, NOT original Dovekie1820 campaign replay.
// Geometry supplied by native C++; original class whitened zHD basis also
// exercised independently of native B. Two score columns retain both routes.
// NumPy2.2.6,
// SciPy1.15.3, Python3.12.14. Relative -q/2 removes constant half class
// lognorm; prior-integrated ratios unused, no absolute evidence normalization
// claim. Source:
// studies/unified_cosmology/code/inference/luminosity_sensitivity.py original
// nodes lines28-38,41-93, no source edits. Native optional test requires
// external original SHA guard. Frozen reference allocation abs2e-7 per score.
namespace irred::test_reference {
inline constexpr std::string_view grey_class_source_sha256 =
    "3f831be63b3776b02c0b547dbb6d9273ddbad68b272c6cb372eda390abafaab7";
inline constexpr std::string_view grey_class_node_sha256 =
    "bd049c686b11ddca3fb3d896b48026b60a309fa935da73f88fd9abf039c032ed";
inline constexpr std::string_view grey_class_selected_covariance_sha256 =
    "64339b81cfd28998ca04f589394c9209b44c93dd63112704eb21bdeb6f300ee4";
struct GreyMagnitudeReference {
  unsigned model;
  double omega_m, constant_q, w0, wa, epsilon_mag, relative_score,
      original_basis_score;
};
inline constexpr std::array<GreyMagnitudeReference, 12>
    grey_magnitude_reference{{
        {0, 0.29999999999999999, 0, -1, 0, -0.20000000000000001,
         -727.19945172038319, -727.19945172038308},
        {0, 0.29999999999999999, 0, -1, 0, 0, -703.03398247940606,
         -703.03398247940606},
        {0, 0.29999999999999999, 0, -1, 0, 0.20000000000000001,
         -758.50345207021394, -758.50345207021303},
        {1, 0, 0, -1, 0, -0.20000000000000001, -882.14774011293116,
         -882.14774011293161},
        {1, 0, 0, -1, 0, 0, -754.19381072963699, -754.19381072963699},
        {1, 0, 0, -1, 0, 0.20000000000000001, -705.87482017812704,
         -705.87482017812681},
        {1, 0, 0.5, -1, 0, -0.20000000000000001, -1324.634705976197,
         -1324.634705976197},
        {1, 0, 0.5, -1, 0, 0, -1050.572770714746, -1050.572770714746},
        {1, 0, 0.5, -1, 0, 0.20000000000000001, -856.14577428507937,
         -856.14577428507937},
        {2, 0.34999999999999998, 0, -0.90000000000000002, 0.20000000000000001,
         -0.20000000000000001, -780.00355638111887, -780.00355638111932},
        {2, 0.34999999999999998, 0, -0.90000000000000002, 0.20000000000000001,
         0, -707.79383429753193, -707.79383429753193},
        {2, 0.34999999999999998, 0, -0.90000000000000002, 0.20000000000000001,
         0.20000000000000001, -715.21905104572841, -715.21905104572886},
    }};
} // namespace irred::test_reference
