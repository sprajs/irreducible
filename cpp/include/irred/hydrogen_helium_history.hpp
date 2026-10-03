#pragma once
#include "irred/thermal_neutrino.hpp"
#include <limits>
#include <string>
namespace irred::cosmology {
namespace detail {
struct HydrogenHeliumCellAccess;
struct HydrogenHeliumPreparationAccess;
}
// A distinct bounded singlet-only approximation: source-code-convention He
// detailed balance/escape, shared NIST central ionization energies, F_H=1,
// T_photo=T_m, and Compton/adiabatic temperature. Not literal SSS1999/RECFAST.
inline constexpr std::string_view hydrogen_helium_history_model_id =
    "HII-HeII-singlet-RecfastCLASS-convention-NIST-central-FH1-TphotoTm-"
    "Compton-adiabatic-bounded/v1";
inline constexpr std::string_view hydrogen_helium_history_method_id =
    "shared-charge-coupled-BE-temperature-elimination-quadratic-z-"
    "Richardson2-refinement3/v1";
inline constexpr std::string_view hydrogen_helium_supplied_history_model_id =
    "HII-HeII-singlet-RecfastCLASS-convention-NIST-central-FH1-TphotoTm-"
    "Compton-adiabatic-supplied-binary64-boundary-bounded/v1";
enum class HydrogenHeliumHistoryBoundary {
  restricted_two_stage_saha,
  supplied_binary64
};
struct HydrogenHeliumHistoryRequest {
  ThermalPhysicalModel model;
  // Independently supplied total physical nuclei densities today, m^-3.
  // No abundance, mass, or consistency with omega_b is inferred here.
  double hydrogen_nuclei_today_per_cubic_metre = 0;
  double helium_nuclei_today_per_cubic_metre = 0;
  std::string nuclei_origin;
  double initial_redshift = 2700, late_redshift = 300;
};
struct HydrogenHeliumHistoryPolicy {
  std::size_t base_intervals = 8192, maximum_fine_intervals = 65536;
  std::size_t maximum_total_work = 4000000;
  std::size_t maximum_native_bytes = 32 * 1024 * 1024;
  double absolute_fraction_tolerance = 1e-8, relative_fraction_tolerance = 2e-6;
  double absolute_temperature_tolerance_kelvin = 2e-4,
         relative_temperature_tolerance = 2e-6;
  double absolute_opacity_tolerance = 2e-9, relative_opacity_tolerance = 3e-6;
  ThermalPolicy thermal;
};
// These emitted binary64 values are the exact boundary of a conditional IVP.
// Earlier wide-to-double casts and upstream numerical/physical uncertainty
// are not propagated. Fractions use their own nuclei denominator.
struct HydrogenHeliumSuppliedInitialState {
  double hydrogen_ionized_fraction = 0;
  double helium_singly_ionized_fraction = 0;
  double matter_temperature_kelvin = 0;
  std::string origin;
  std::string source_identity;
};
struct HydrogenHeliumSuppliedHistoryRequest {
  HydrogenHeliumHistoryRequest history;
  HydrogenHeliumSuppliedInitialState initial;
};
struct HydrogenHeliumInitialImportWitness {
  // Order: HII/H, HeII/He, Tm in K. Promotion is exact on this wide profile.
  std::array<long double, 3> promoted_values{};
  std::array<long double, 3> measured_absolute_promotion_loss{};
};
struct HydrogenHeliumRateDomainWitness {
  // Includes every attempted residual/rate temperature, including failed
  // line-search/mesh work. Absence is not a zero-temperature range.
  std::optional<std::array<long double, 2>> attempted_kelvin_range;
  bool invalid_temperature_attempted = false;
  // Present only after complete central Richardson node/cell construction.
  // Encloses its linear cell T, not the exact ODE or its uncertainty.
  std::optional<std::array<long double, 2>> retained_kelvin_range;
};
struct HydrogenHeliumHistoryValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0;
};
inline constexpr unsigned history_hydrogen_ionized = 1,
                          history_helium_ionized = 2,
                          history_electron_density = 4,
                          history_matter_temperature = 8,
                          history_thomson_opacity = 16;
struct HydrogenHeliumHistoryRow {
  double redshift = 0;
  HydrogenHeliumHistoryValue hydrogen_ionized_fraction,
      helium_singly_ionized_fraction, electron_number_density_per_cubic_metre,
      matter_temperature_kelvin, thomson_opacity_per_redshift;
};
struct HydrogenHeliumHistoryBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested_outputs = 0;
  std::vector<HydrogenHeliumHistoryRow> rows;
};
struct HydrogenHeliumHistoryWork {
  std::size_t background_evaluations = 0, momentum_callbacks = 0,
              initial_charge_evaluations = 0, rhs_evaluations = 0,
              initial_boundary_evaluations = 0;
  std::size_t total() const noexcept {
    std::size_t sum = 0;
    for (const auto n : {background_evaluations, momentum_callbacks,
                        initial_charge_evaluations, rhs_evaluations,
                        initial_boundary_evaluations}) {
      if (n > std::numeric_limits<std::size_t>::max() - sum)
        return std::numeric_limits<std::size_t>::max();
      sum += n;
    }
    return sum;
  }
};
class HydrogenHeliumHistory {
public:
  HydrogenHeliumHistory() = default;
  HydrogenHeliumHistory(const HydrogenHeliumHistory &) = default;
  HydrogenHeliumHistory &operator=(const HydrogenHeliumHistory &);
  HydrogenHeliumHistory(HydrogenHeliumHistory &&) noexcept;
  HydrogenHeliumHistory &operator=(HydrogenHeliumHistory &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const HydrogenHeliumHistoryRequest *source() const noexcept {
    return source_ ? &*source_ : nullptr;
  }
  const ThermalBackground *background() const noexcept {
    return source_ ? &background_ : nullptr;
  }
  HydrogenHeliumHistoryWork work() const noexcept { return work_; }
  // A tag on a refused owner describes the attempted profile, not admission.
  std::optional<HydrogenHeliumHistoryBoundary> boundary_kind() const noexcept {
    return boundary_;
  }
  std::string_view model_identity() const noexcept;
  std::string_view method_identity() const noexcept {
    return boundary_ ? hydrogen_helium_history_method_id : std::string_view{};
  }
  const HydrogenHeliumSuppliedInitialState *supplied_initial_state() const noexcept {
    return supplied_initial_ ? &*supplied_initial_ : nullptr;
  }
  const HydrogenHeliumInitialImportWitness *initial_import_witness() const noexcept {
    return initial_import_ ? &*initial_import_ : nullptr;
  }
  const HydrogenHeliumRateDomainWitness *rate_domain_witness() const noexcept {
    return rate_domain_ ? &*rate_domain_ : nullptr;
  }
  // Original map metadata is retained once, including on later preparation
  // refusal. No remapping of the supplied physical source is performed.
  const std::array<ThermalScalarMapWitness, 4> *thermal_mapping_witnesses() const noexcept {
    return mapping_witnesses_ ? &*mapping_witnesses_ : nullptr;
  }
  // Positive omitted-stage activity Qe exp(-chiHeII/kT)/ne at initialization.
  // A physical scope witness, excluded from numerical diagnostics. <=1e-12.
  std::optional<double> excluded_initial_heiii_activity() const noexcept {
    return excluded_heiii_activity_;
  }
  // Maximum log(HeIII/HeII Saha activity) over the retained nodes. This is a
  // model-scope witness, not a bound on a full kinetic HeIII population.
  std::optional<double> maximum_log_excluded_heiii_activity() const noexcept {
    return maximum_log_heiii_activity_;
  }
  HydrogenHeliumHistoryBatch evaluate(
      std::span<const double>, unsigned outputs,
      std::size_t maximum_output_points = 4096,
      std::size_t maximum_native_bytes = 32 * 1024 * 1024) const;
private:
  struct Node {
    long double z = 0, hydrogen = 0, helium = 0, temperature = 0;
    long double hydrogen_error = 0, helium_error = 0, temperature_error = 0;
    long double opacity_coefficient = 0, opacity_coefficient_error = 0;
  };
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<HydrogenHeliumHistoryRequest> source_;
  std::optional<HydrogenHeliumHistoryBoundary> boundary_;
  std::optional<HydrogenHeliumSuppliedInitialState> supplied_initial_;
  std::optional<HydrogenHeliumInitialImportWitness> initial_import_;
  std::optional<HydrogenHeliumRateDomainWitness> rate_domain_;
  ThermalBackground background_;
  HydrogenHeliumHistoryPolicy policy_;
  HydrogenHeliumHistoryWork work_;
  std::optional<std::array<ThermalScalarMapWitness, 4>> mapping_witnesses_;
  std::optional<double> excluded_heiii_activity_;
  std::optional<double> maximum_log_heiii_activity_;
  std::vector<Node> nodes_;
  friend struct detail::HydrogenHeliumCellAccess;
  friend struct detail::HydrogenHeliumPreparationAccess;
  friend std::optional<std::size_t> hydrogen_helium_history_payload_bound(
      std::size_t, std::size_t, std::size_t, std::size_t, std::size_t,
      std::size_t) noexcept;
};
HydrogenHeliumHistory prepare_hydrogen_helium_history(
    const HydrogenHeliumHistoryRequest &, HydrogenHeliumHistoryPolicy = {});
HydrogenHeliumHistory prepare_hydrogen_helium_supplied_history(
    const HydrogenHeliumSuppliedHistoryRequest &, HydrogenHeliumHistoryPolicy = {});
std::optional<std::size_t> hydrogen_helium_history_payload_bound(
    std::size_t fine_intervals, std::size_t output_points,
    std::size_t species, std::size_t nuclei_origin_bytes,
    std::size_t supplied_origin_bytes = 0,
    std::size_t supplied_source_identity_bytes = 0) noexcept;
} // namespace irred::cosmology
