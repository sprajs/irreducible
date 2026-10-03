#pragma once
// Original synthetic source frozen before runtime. This is a proof-caller
// fixture, not a physical abundance/mass/temperature probability law.
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <string_view>
namespace abundance_history_cohort {
inline constexpr std::string_view model_id =
    "synthetic-two-equal-state-H1-He4-neutral-mass-photon-temperature-"
    "retained-HHe-working-history-law/v1";
inline constexpr std::string_view method_id =
    "two-public-history-owners-pooled-Tm-opacity-exact-half-population-moments/v1";
inline constexpr std::string_view role = "synthetic_control";
inline constexpr std::string_view temperature_role = "supplied-photon-temperature-today";
inline constexpr std::string_view distribution_origin = "original synthetic two-equal-state history support/v1";
inline constexpr std::string_view dependence_origin = "one complete support state shared by all ordered rows and outputs/v1";
inline constexpr std::string_view source_origin_a = "original synthetic complete thermal-abundance state A/v1";
inline constexpr std::string_view source_origin_b = "original synthetic complete thermal-abundance state B/v1";
inline constexpr std::string_view mass_origin = "supplied synthetic neutral H1-He4 effective masses; existing installed fixture ancestry/v1";
inline constexpr std::string_view photon_origin_a = "supplied synthetic photon T_CMB,0=2.7 K/v1";
inline constexpr std::string_view photon_origin_b = "supplied synthetic photon T_CMB,0=2.75 K/v1";
inline constexpr std::array<std::string_view, 2> state_ids{"synthetic-A", "synthetic-B"};
inline constexpr std::array<std::string_view, 3> row_ids{"initial-2700", "interior-1300", "late-300"};
inline constexpr double h0_a = 63, h0_b = 77, omega_b = .02237,
    omega_cdm = .12, tcmb_a = 2.7, tcmb_b = 2.75, omega_other = 1.7e-5,
    y_a = .24, y_b = .26,
    neutral_hydrogen_mass_kg = 1.6735328383153192e-27,
    neutral_helium_mass_kg = 6.646479071583153e-27,
    relative_mass = 1, probability = .5,
    initial_redshift = 2700, interior_redshift = 1300, late_redshift = 300,
    covariance_resolution = .01;
inline constexpr unsigned outputs = 24;
inline constexpr std::size_t preparation_work_cap = 8000006,
    evaluation_work_cap = 1024, two_evaluation_work_cap = 8002054,
    whole_live_payload_cap = 128 * 1024 * 1024,
    child_work_cap = 4000000, child_payload_cap = 32 * 1024 * 1024,
    default_base_intervals = 8192, maximum_fine_intervals = 65536,
    states = 2, rows = 3, axes = 6, covariance_cells = 36,
    maximum_origin_bytes = 4096;
constexpr std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
static_assert(bits(h0_a) == 0x404f800000000000ULL);
static_assert(bits(h0_b) == 0x4053400000000000ULL);
static_assert(bits(omega_b) == 0x3f96e82949a56580ULL);
static_assert(bits(omega_cdm) == 0x3fbeb851eb851eb8ULL);
static_assert(bits(tcmb_a) == 0x400599999999999aULL);
static_assert(bits(tcmb_b) == 0x4006000000000000ULL);
static_assert(bits(omega_other) == 0x3ef1d3671ac14c66ULL);
static_assert(bits(y_a) == 0x3fceb851eb851eb8ULL);
static_assert(bits(y_b) == 0x3fd0a3d70a3d70a4ULL);
static_assert(bits(neutral_hydrogen_mass_kg) == 0x3a6092e8e991bb4cULL);
static_assert(bits(neutral_helium_mass_kg) == 0x3a8074b4e33ac660ULL);
static_assert(bits(relative_mass) == 0x3ff0000000000000ULL);
static_assert(bits(probability) == 0x3fe0000000000000ULL);
static_assert(bits(initial_redshift) == 0x40a5180000000000ULL);
static_assert(bits(interior_redshift) == 0x4094500000000000ULL);
static_assert(bits(late_redshift) == 0x4072c00000000000ULL);
static_assert(bits(covariance_resolution) == 0x3f847ae147ae147bULL);
static_assert(two_evaluation_work_cap == preparation_work_cap + 2 * evaluation_work_cap);
} // namespace abundance_history_cohort
