#pragma once
#include "irred/thermal_neutrino.hpp"
#include <string>
namespace irred::cosmology {
// A distinct bounded singlet-only approximation: source-code-convention He
// detailed balance/escape, shared NIST central ionization energies, F_H=1,
// T_photo=T_m, and Compton/adiabatic temperature. Not literal SSS1999/RECFAST.
inline constexpr std::string_view hydrogen_helium_history_model_id =
    "HII-HeII-singlet-RecfastCLASS-convention-NIST-central-FH1-TphotoTm-"
    "Compton-adiabatic-bounded/v1";
inline constexpr std::string_view hydrogen_helium_history_method_id =
    "shared-charge-coupled-BE-temperature-elimination-quadratic-z-"
    "Richardson2-refinement3/v1";
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
              initial_charge_evaluations = 0, rhs_evaluations = 0;
  std::size_t total() const noexcept {
    return background_evaluations + momentum_callbacks +
           initial_charge_evaluations + rhs_evaluations;
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
  ThermalBackground background_;
  HydrogenHeliumHistoryPolicy policy_;
  HydrogenHeliumHistoryWork work_;
  std::optional<double> excluded_heiii_activity_;
  std::optional<double> maximum_log_heiii_activity_;
  std::vector<Node> nodes_;
  friend HydrogenHeliumHistory prepare_hydrogen_helium_history(
      const HydrogenHeliumHistoryRequest &, HydrogenHeliumHistoryPolicy);
  friend std::optional<std::size_t> hydrogen_helium_history_payload_bound(
      std::size_t, std::size_t, std::size_t, std::size_t) noexcept;
};
HydrogenHeliumHistory prepare_hydrogen_helium_history(
    const HydrogenHeliumHistoryRequest &, HydrogenHeliumHistoryPolicy = {});
std::optional<std::size_t> hydrogen_helium_history_payload_bound(
    std::size_t fine_intervals, std::size_t output_points,
    std::size_t species, std::size_t nuclei_origin_bytes) noexcept;
} // namespace irred::cosmology
