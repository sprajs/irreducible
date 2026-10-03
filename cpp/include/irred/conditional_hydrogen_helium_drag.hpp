// Conditional supplied-depth sensitivity under exact emitted working inputs.
// Physical present-day drag requires a separately qualified late-tail law.
#pragma once
#include <irred/baryon_abundance.hpp>
#include <irred/hydrogen_helium_history.hpp>
#include <array>
#include <cstdint>

namespace irred::cosmology {
namespace detail { struct ConditionalDragAccess; }
inline constexpr std::string_view conditional_hydrogen_helium_drag_model_id =
    "exact-emitted-H1-He4-nuclei-and-retained-thermal-state-supplied-D-"
    "conditional-drag-ruler/v1";
struct ConditionalHydrogenHeliumDragRequest {
  ThermalPhysicalModel model; // sole omega_b source
  double helium4_neutral_mass_fraction = 0;
  double hydrogen1_neutral_effective_mass_kg = 0;
  double helium4_neutral_effective_mass_kg = 0;
  std::string source_origin, mass_origin;
  double initial_redshift = 2700, late_redshift = 300;
};
struct ConditionalDragInterval {
  std::string id;
  double lower_depth = 0, upper_depth = 0;
  std::string origin;
};
struct ConditionalHydrogenHeliumDragPolicy {
  HydrogenHeliumHistoryPolicy history;
  std::size_t maximum_total_work = 4000000;
  std::size_t maximum_native_bytes = 32 * 1024 * 1024;
  std::size_t maximum_intervals = 4, maximum_root_endpoints = 8;
  unsigned maximum_endpoint_trials = 128, maximum_bracket_proposals = 8;
  std::size_t maximum_outer_callbacks_per_integral = 200000;
  unsigned maximum_depth = 30;
  double total_ruler_absolute_tolerance_mpc = 1e-3;
  double total_ruler_relative_tolerance = 0;
  double capacity_absolute_tolerance = 1e-6, capacity_relative_tolerance = 1e-8;
};
struct ConditionalDragEndpoint {
  numerics::Status root_status = numerics::Status::invalid_input;
  numerics::Status ruler_status = numerics::Status::invalid_input;
  std::optional<double> redshift, comoving_ruler_mpc;
  double redshift_numerical_estimate = 0, total_ruler_numerical_estimate_mpc = 0;
};
struct ConditionalDragIntervalRow {
  ConditionalDragInterval source;
  // first=min z (D_upper), second=max z (D_lower); ruler order reverses.
  std::array<ConditionalDragEndpoint, 2> endpoints;
};
struct ConditionalDragWork {
  // OUTPUT ONLY; never an input ledger. checked_total() cannot wrap.
  std::size_t source_map = 0, background = 0, momentum = 0, charge = 0, rhs = 0;
  std::size_t cell = 0, prefix = 0, primitive = 0, root = 0, bound = 0;
  std::size_t capacity_outer = 0, ruler_outer = 0;
  std::optional<std::size_t> checked_total() const noexcept;
};
struct ConditionalDragSourceSnapshot {
  BaryonAbundanceRow nuclei_today;
  std::array<std::uint64_t, 2> emitted_nuclei_binary64_bits{};
  ThermalFlatModel fixed_mapped_thermal_source;
  // Order: photon, baryon, CDM, other massless fractional density.
  std::array<ThermalScalarMapWitness, 4> thermal_source_map;
  double relative_mass_reconstruction_residual = 0;
  double relative_photon_reconstruction_residual = 0;
  double reconstruction_arithmetic_estimate = 0;
  double baryon_loading = 0, loading_arithmetic_estimate = 0;
  // All source-to-emitted witnesses above are provenance, unpropagated.
  // Loading arithmetic is propagated through drag and ruler diagnostics.
};
// ADDITIVE diagnostic needed by Repro N16-R3. Fixed-size owned audit metadata;
// no public ledger mutation, no physical error/probability interpretation.
struct ConditionalDragBudgetDiagnostics {
  std::size_t requested_total_work = 0, requested_native_bytes = 0;
  std::size_t requested_history_work = 0, requested_history_native_bytes = 0;
  // Engaged when the history admission phase chooses effective ceilings.
  // They remain available on refusal; never fabricated if this phase is absent.
  std::optional<std::size_t> served_history_work, served_history_native_bytes;
  std::optional<std::size_t> earned_work_before_history, parent_live_bytes_before_history;
  bool history_preparation_attempted = false;
};
struct ConditionalDragCapacity {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> upper_depth_capacity;
  double numerical_estimate = 0;
};
struct ConditionalDragBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::vector<ConditionalDragIntervalRow> rows;
  ConditionalDragWork work;
  ConditionalDragCapacity late_charge_capacity;
  std::size_t retained_payload_bytes = 0, peak_payload_bytes = 0;
};
class ConditionalHydrogenHeliumDrag {
public:
  ConditionalHydrogenHeliumDrag() = default;
  ConditionalHydrogenHeliumDrag(const ConditionalHydrogenHeliumDrag &) = default;
  ConditionalHydrogenHeliumDrag &operator=(const ConditionalHydrogenHeliumDrag &);
  ConditionalHydrogenHeliumDrag(ConditionalHydrogenHeliumDrag &&) noexcept;
  ConditionalHydrogenHeliumDrag &operator=(ConditionalHydrogenHeliumDrag &&) noexcept;
  numerics::Status status() const noexcept;
  const ConditionalHydrogenHeliumDragRequest *source() const noexcept;
  // Exact emitted values/bits, mapper diagnostics and reconstruction witnesses
  // belong in a read-only source snapshot; source-map errors are unpropagated.
  const ConditionalDragSourceSnapshot *source_snapshot() const noexcept;
  const HydrogenHeliumHistory *history() const noexcept;
  ConditionalDragWork preparation_work() const noexcept;
  // Supplied/effective cap record survives preparation refusal/copy/move.
  // A moved-from/never-admitted owner returns nullptr.
  const ConditionalDragBudgetDiagnostics *budget_diagnostics() const noexcept;
  ConditionalDragBatch evaluate(std::span<const ConditionalDragInterval>) const;
private:
  struct Prefix { long double value = 0, error = 0; };
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<ConditionalHydrogenHeliumDragRequest> source_;
  BaryonAbundance abundance_;
  HydrogenHeliumHistory history_;
  std::optional<ConditionalDragSourceSnapshot> snapshot_;
  std::optional<ConditionalDragBudgetDiagnostics> budget_;
  ConditionalHydrogenHeliumDragPolicy policy_;
  ConditionalDragWork work_;
  ConditionalDragCapacity capacity_;
  long double loading_ = 0, loading_error_ = 0, loading_numerator_ = 0, loading_denominator_ = 0;
  std::vector<Prefix> prefix_;
  std::size_t spent_ = 0, peak_payload_ = 0;
  std::optional<std::size_t> payload() const noexcept;
  friend struct detail::ConditionalDragAccess;
  friend ConditionalHydrogenHeliumDrag prepare_conditional_hydrogen_helium_drag(
      const ConditionalHydrogenHeliumDragRequest &,
      ConditionalHydrogenHeliumDragPolicy);
};
ConditionalHydrogenHeliumDrag prepare_conditional_hydrogen_helium_drag(
    const ConditionalHydrogenHeliumDragRequest &,
    ConditionalHydrogenHeliumDragPolicy = {});
} // namespace irred::cosmology
