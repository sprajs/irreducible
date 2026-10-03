// Fixed synthetic public-API consumer. No observational or reference admission.
#include "irred/finite_opacity_source.hpp"
#include "irred/quantities.hpp"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <optional>
#include <span>
#include <utility>

namespace {
namespace c = irred::cosmology;
namespace p = irred::projection;
using W = long double;
using S = irred::numerics::Status;
const char *status(S value) {
  switch (value) {
  case S::ok: return "ok";
  case S::invalid_input: return "invalid_input";
  case S::nonfinite_input: return "nonfinite_input";
  case S::overflow: return "overflow";
  case S::work_limit: return "work_limit";
  case S::outside_domain: return "outside_domain";
  case S::singular: return "singular";
  case S::not_positive_definite: return "not_positive_definite";
  case S::conditioning_budget_exceeded: return "conditioning_budget_exceeded";
  }
  return "unsupported_status";
}
template<class T> void optional(const std::optional<T> &value) {
  if (value) std::cout << *value;
  else std::cout << "absent";
}
void values(std::span<const W> row) {
  for (W value : row) std::cout << '\t' << value;
}
// SAME first-trial thermal ancestry and Simpson16384 input-axis construction.
// This fixture supplies emitted eta knots; it is no independent clock reference.
struct Age {
  S code = S::ok;
  std::optional<W> value;
  std::size_t calls = 0;
  W partial_sum = 0, last_a = 0, last_p = 0;
};
Age age(const c::ThermalBackground &background, W upper) {
  Age out;
  constexpr unsigned panels = 16384;
  for (unsigned i = 0; i <= panels; ++i) {
    out.last_a = upper * i / panels;
    ++out.calls;
    const auto query = background.scaled_expansion(out.last_a);
    out.last_p = query.a4_e2;
    if (query.status != S::ok) { out.code = query.status; return out; }
    if (!(out.last_p > 0) || !std::isfinite(out.last_p)) {
      out.code = S::outside_domain; return out;
    }
    const W weight = i == 0 || i == panels ? 1 : i % 2 ? 4 : 2;
    out.partial_sum += weight * (W(irred::speed_of_light_m_per_s) / 1000) /
                       background.source().h0_km_s_mpc / std::sqrt(out.last_p);
    if (!std::isfinite(out.partial_sum)) { out.code = S::overflow; return out; }
  }
  const W computed = upper * out.partial_sum / (3 * panels);
  if (!(computed > 0) || !std::isfinite(computed)) { out.code = S::overflow; return out; }
  out.value = computed;
  return out;
}
void print_age(const Age &value, W scale_factor) {
  std::cout << "fixture_age\ta=" << scale_factor << "\tstatus=" << status(value.code)
            << "\tcalls=" << value.calls << "\teta_wide=";
  optional(value.value);
  std::cout << "\tpartial_sum=" << value.partial_sum << "\tlast_a=" << value.last_a
            << "\tlast_p=" << value.last_p << "\tage_error=absent\n";
}
void print_boundary(const c::FiniteOpacityIdentity &id) {
  const auto &b = id.boundary;
  // Preparation retains its identity BEFORE completing clock/seeds/boundary.
  // eta_i is assigned with tau_i/survival_i only after all seeds are saved.
  const bool scalar_boundary_computed = b.eta_i_mpc > 0;
  for (std::size_t k = 0; k < id.original.k_mpc_inverse.size(); ++k) {
    std::cout << "boundary\tk_index=" << k << "\tk=" << id.original.k_mpc_inverse[k]
              << "\tscalars_computed=" << scalar_boundary_computed << "\teta_i=";
    if (scalar_boundary_computed) std::cout << b.eta_i_mpc; else std::cout << "absent";
    std::cout << "\ttau_i=";
    if (scalar_boundary_computed) std::cout << b.tau_i; else std::cout << "absent";
    std::cout << "\tw_i=";
    if (scalar_boundary_computed) std::cout << b.survival_i; else std::cout << "absent";
    std::cout << "\tmonopole=";
    if (k < b.photon_monopole.size()) std::cout << b.photon_monopole[k]; else std::cout << "absent";
    std::cout << "\tdipole=";
    if (k < b.photon_dipole_theta_over_k.size()) std::cout << b.photon_dipole_theta_over_k[k]; else std::cout << "absent";
    std::cout << "\tomitted_temperature_bound=";
    if (k < b.omitted_temperature_absolute_bound.size()) std::cout << b.omitted_temperature_absolute_bound[k]; else std::cout << "absent";
    std::cout << "\tinitial_e=";
    if (k < id.seeds.size()) std::cout << 0; else std::cout << "absent";
    std::cout << "\tboundary_temperature=absent\tfull_temperature=absent"
                 "\tboundary_radial_error=absent\n";
  }
}
void print_source_receipt(const c::FiniteOpacitySourceResult &result) {
  const auto &w = result.work;
  std::cout << "producer\tstatus=" << status(result.status)
            << "\tsource_admission=" << result.source_numerically_admitted
            << "\tcomplete_grid=" << result.source.has_value() << "\tpeak_native_bytes=";
  optional(result.peak_owned_payload_bound);
  std::cout << "\nbegun_steps\t" << w.attempted_steps << "\tdenied_steps\t" << w.denied_step_requests
            << "\tbackground_calls\t" << w.background_clock_calls
            << "\tcoupled_solves\t" << w.coupled_stage_solves
            << "\ttail_inversions\t" << w.tail_block_inversions
            << "\tcore_factors\t" << w.core_factorizations
            << "\tlogical_writes\t" << w.destination_writes
            << "\treserved_prefix_writes\t" << w.reserved_refusal_writes
            << "\trefused_write_request\t" << w.refused_write_request << '\n';
  for (std::size_t ai = 0; ai < result.attempts.size(); ++ai) {
    const auto &a = result.attempts[ai];
    std::cout << "attempt\t" << ai << "\tL=" << a.hierarchy << "\ttime_refinement=" << a.time_refinement
              << "\tstatus=" << status(a.status) << "\tbegun=" << a.attempted_steps
              << "\tcompleted=" << a.completed_steps << "\tdenied=" << a.denied_step_requests
              << "\tnodes=" << a.nodes.size() << "\treached_k=" << a.reached_wavenumbers
              << "\tstage_residual=" << a.maximum_stage_residual << "\tminimum_pivot=" << a.minimum_scaled_pivot
              << "\tclock_residual=" << a.maximum_clock_stage_residual
              << "\tobserver_a_mismatch=" << a.maximum_observer_scale_factor_mismatch
              << "\tobserver_a_consistent=" << a.observer_scale_factor_consistent << '\n';
    if (result.source) continue;
    // Retain all actual source-refusal nodes and matching reached core/tails.
    // Unreached allocated slots are not emitted as zero-valued physical states.
    for (const auto &node : a.nodes) {
      std::cout << "refused_node\t" << ai << '\t' << node.k_index << '\t' << node.eta_index
                << '\t' << node.eta_mpc << '\t' << node.scale_factor;
      values(node.core); values(node.raw_channels); values(node.channel_cast_loss);
      std::cout << '\t' << node.psi << '\t' << node.phi_prime_mpc_inverse
                << '\t' << node.hamiltonian_residual << '\t' << node.hamiltonian_term_scale
                << '\t' << node.normalized_hamiltonian_residual << '\t' << node.closure_defect
                << '\t' << node.actual_p_minus_shadow << '\t' << node.conditional_hcal_estimate
                << '\t' << node.survival << '\t' << node.visibility_mpc_inverse << '\n';
    }
    const std::size_t stride = a.hierarchy >= 3 ? a.hierarchy - 2 : 0;
    for (std::size_t k = 0; k < a.reached_wavenumbers; ++k) {
      std::cout << "refused_state\t" << ai << '\t' << k;
      if (k < a.final_core.size() && k < a.final_eta_mpc.size() && k < a.final_scale_factor.size()) {
        std::cout << '\t' << a.final_eta_mpc[k] << '\t' << a.final_scale_factor[k];
        values(a.final_core[k]);
      } else std::cout << "\tstate=absent";
      std::cout << '\n';
      for (bool temperature : {true, false}) {
        const auto &tail = temperature ? a.final_temperature_tail : a.final_polarization_tail;
        std::cout << (temperature ? "refused_F_tail" : "refused_G_tail") << '\t' << ai << '\t' << k;
        if (stride && k < tail.size() / stride) values(std::span<const W>(tail).subspan(k * stride, stride));
        else std::cout << "\tabsent";
        std::cout << '\n';
      }
    }
  }
  for (const auto &d : result.diagnostics) {
    std::cout << "source_witness\t" << d.eta_index << '\t' << d.k_index << '\t' << d.channel
              << '\t' << d.allocation << '\t' << d.time_difference << '\t' << d.hierarchy_difference
              << '\t' << d.previous_time_difference << '\t' << d.previous_hierarchy_difference
              << '\t' << d.measured_cast_loss << "\tcommon_clock_error=";
    optional(d.common_background_clock_error); std::cout << "\tarithmetic_linear_error=";
    optional(d.arithmetic_linear_error); std::cout << "\tsource_grid_error=";
    optional(d.source_grid_error); std::cout << '\n';
  }
}
void print_grid(const p::ContinuousCmbSource &grid) {
  std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "grid_identity\t" << grid.producer_id << '\t' << grid.signed_mode_id << '\t' << grid.normalization_id
            << "\nsupport\teta_first=" << grid.eta_mpc.front() << "\teta_last=" << grid.eta_mpc.back()
            << "\tobserver_eta=" << grid.observer_eta_mpc
            << "\tinterpolation=within_cell_linear\tk_interpolation=absent\toutside_extension=absent\n";
  for (std::size_t ti = 0; ti < grid.eta_mpc.size(); ++ti)
    for (std::size_t ki = 0; ki < grid.k_mpc_inverse.size(); ++ki) {
      const auto i = ti * grid.k_mpc_inverse.size() + ki;
      std::cout << "raw_grid\t" << ti << '\t' << ki << '\t' << grid.eta_mpc[ti] << '\t' << grid.k_mpc_inverse[ki]
                << '\t' << grid.t0[i] << '\t' << grid.t1[i] << '\t' << grid.t2[i] << '\t' << grid.polarization[i]
                << "\tsource_admission=false\tchannels_unit=Mpc_inverse\n";
    }
  std::cout << std::setprecision(std::numeric_limits<W>::max_digits10);
}
} // namespace

int main() {
  std::cout.imbue(std::locale::classic());
  std::cout << std::boolalpha << std::scientific << std::setprecision(std::numeric_limits<W>::max_digits10)
            << "contract\tsynthetic=true\tphysical_qualification=false\tfull_temperature=absent"
               "\treference_admission=false\tmissing_integrated_errors=absent\n";
  c::FiniteOpacityRequest input;
  input.model = {67.4, .0224, .12, 2.7255, 0, {}};
  input.initial_scale_factor = .001; input.observer_scale_factor = .00101;
  input.k_mpc_inverse = {.005, .01};
  const auto map = c::map_thermal_physical_model(input.model);
  std::cout << "fixture_map\tstatus=" << status(map.status) << '\n';
  if (map.status != S::ok || !map.model) return 1;
  const auto background = c::prepare_thermal_background(*map.model);
  std::cout << "fixture_background\tstatus=" << status(background.status()) << '\n';
  if (background.status() != S::ok) return 1;
  const auto first = age(background, input.initial_scale_factor); print_age(first, input.initial_scale_factor);
  if (first.code != S::ok || !first.value) return 1;
  const auto last = age(background, input.observer_scale_factor); print_age(last, input.observer_scale_factor);
  if (last.code != S::ok || !last.value) return 1;
  input.opacity.eta_mpc = {static_cast<double>(*first.value), static_cast<double>(*last.value)};
  input.opacity.differential_opacity_mpc_inverse = {.1, .2};
  input.opacity.origin = "synthetic-no-electron-history";
  input.opacity.law_id = "emitted-linear-conformal-opacity/v1";
  input.opacity.exact_member_digest = "synthetic-fixture-axis-construction-not-content-certification";
  input.source_origin = "synthetic-finite-opacity-source-grid-consumer/v1";
  c::FiniteOpacityPolicy policy;
  policy.maximum_attempted_steps = 2000; policy.maximum_background_clock_calls = 400000;
  policy.maximum_coupled_stage_solves = 2000; policy.maximum_destination_writes = 100000000;
  policy.maximum_native_bytes = 8 * 1024 * 1024;
  auto producer = c::prepare_finite_opacity_source(std::move(input), policy);
  std::cout << "source_preparation\tstatus=" << status(producer.status()) << '\n';
  if (producer.identity()) {
    const auto &id = *producer.identity();
    std::cout << "identity\t" << id.model_id << '\t' << id.mode_id << '\t' << id.method_id << '\t' << id.arithmetic_id << '\n';
    std::cout << "original_input\tsource_origin=" << id.original.source_origin
              << "\topacity_origin=" << id.original.opacity.origin << "\topacity_law=" << id.original.opacity.law_id
              << "\topacity_digest=" << id.original.opacity.exact_member_digest
              << "\ta_initial=" << id.original.initial_scale_factor << "\ta_observer=" << id.original.observer_scale_factor << '\n';
    for (std::size_t i = 0; i < id.original.opacity.eta_mpc.size(); ++i)
      std::cout << "original_opacity\t" << i << '\t' << id.original.opacity.eta_mpc[i]
                << '\t' << id.original.opacity.differential_opacity_mpc_inverse[i] << '\n';
    std::cout << "preparation_prefix\tbackground_calls=" << id.preparation_work.background_clock_calls
              << "\tlogical_writes=" << id.preparation_work.destination_writes
              << "\trefused_write_request=" << id.preparation_work.refused_write_request
              << "\tretained_seeds=" << id.seeds.size() << "\tclock_consistent=" << id.clock.endpoint_clock_consistent << '\n';
    std::cout << "stored_clock_witness\tinitial_age=" << id.clock.initial_age_mpc
              << "\tobserver_age=" << id.clock.observer_age_mpc
              << "\tinitial_difference=" << id.clock.initial_quadrature_difference_mpc
              << "\tobserver_difference=" << id.clock.observer_quadrature_difference_mpc
              << "\tinitial_mismatch=" << id.clock.initial_eta_mismatch_mpc
              << "\tobserver_mismatch=" << id.clock.observer_eta_mismatch_mpc
              << "\tendpoint_consistent=" << id.clock.endpoint_clock_consistent
              << "\tcomplete_clock_error=absent\n";
    for (std::size_t k = 0; k < id.seeds.size(); ++k) {
      const auto &seed = id.seeds[k];
      std::cout << "initial_seed\t" << k;
      for (double value : seed.emitted_core) std::cout << '\t' << value;
      values(seed.absolute_cast_loss);
      std::cout << "\thamiltonian_residual=" << seed.hamiltonian_residual << "\tzeta_residual=" << seed.zeta_residual << '\n';
    }
    print_boundary(id);
  }
  if (producer.status() != S::ok || !producer.identity()) return 1;
  const auto *original_identity = producer.identity();
  const auto *original_boundary = &original_identity->boundary;
  auto receipt = producer.produce(policy); print_source_receipt(receipt);
  const bool same_identity = receipt.identity.get() == original_identity && receipt.boundary() == original_boundary;
  std::cout << "producer_ownership\tsame_immutable_identity_and_boundary=" << same_identity << '\n';
  if (receipt.status != S::conditioning_budget_exceeded || receipt.source_numerically_admitted || !receipt.source || !same_identity) return 1;
  const bool expected_grid = receipt.source->k_mpc_inverse.size() == 2 && receipt.source->eta_mpc.size() == 5 &&
      receipt.source->t0.size() == 10 && receipt.source->t1.size() == 10 &&
      receipt.source->t2.size() == 10 && receipt.source->polarization.size() == 10 && receipt.attempts.size() == 9;
  std::cout << "consumer_grid\texpected_2k_5eta_10cells_9attempts=" << expected_grid << '\n';
  if (!expected_grid) return 1;
  print_grid(*receipt.source);
  std::array<const double *, 6> original_buffers{};
  {
    const auto &s = *receipt.source;
    original_buffers = {s.k_mpc_inverse.data(), s.eta_mpc.data(), s.t0.data(),
                        s.t1.data(), s.t2.data(), s.polarization.data()};
  }
  p::ContinuousCmbPreparationPolicy acquisition;
  acquisition.maximum_k = 2; acquisition.maximum_eta = 5; acquisition.maximum_cells = 10;
  acquisition.maximum_payload_bytes = 8 * 1024 * 1024;
  auto projection = p::prepare_continuous_cmb_projection(std::move(*receipt.source), acquisition);
  std::cout << "projection_preparation\tstatus=" << status(projection.status()) << '\n';
  if (projection.status() != S::ok || !projection.source()) {
    std::cout << "projection_preparation_refusal\toriginal_complete_grid_retained="
              << (receipt.source->k_mpc_inverse.size() == 2 && receipt.source->eta_mpc.size() == 5 &&
                  receipt.source->t0.size() == 10 && receipt.source->t1.size() == 10 &&
                  receipt.source->t2.size() == 10 && receipt.source->polarization.size() == 10) << '\n';
    return 2;
  }
  receipt.source.reset(); // only AFTER successful ownership transfer
  const auto &moved = *projection.source();
  const std::array<const double *, 6> moved_buffers{moved.k_mpc_inverse.data(), moved.eta_mpc.data(), moved.t0.data(),
                                                 moved.t1.data(), moved.t2.data(), moved.polarization.data()};
  const bool same_buffers = original_buffers == moved_buffers;
  std::cout << "ownership\tsame_six_grid_buffers=" << same_buffers << '\n';
  if (!same_buffers) return 2;
  producer = c::FiniteOpacitySourceProducer{};
  std::cout << "lifetime\tproducer_destroyed=true\tsame_identity_retained=" << (receipt.identity.get() == original_identity)
            << "\tsame_boundary_retained=" << (receipt.boundary() == original_boundary)
            << "\tsource_admission=" << receipt.source_numerically_admitted << '\n';
  p::ContinuousCmbPolicy geometry;
  geometry.maximum_multipoles = 1; geometry.maximum_kernel_evaluations = 2000;
  geometry.maximum_depth = 12; geometry.maximum_payload_bytes = 8 * 1024 * 1024;
  constexpr std::array<unsigned, 1> multipoles{2};
  const auto forced = p::project_continuous_cmb(projection, multipoles, p::continuous_temperature | p::continuous_e_mode, geometry);
  const bool same_owner = forced.source_owner.get() == projection.source();
  std::cout << "projection\tstatus=" << status(forced.status) << "\tkernel_evaluations=" << forced.kernel_evaluations
            << "\tpayload_bound=" << forced.payload_bound << "\tsame_grid_owner="
            << same_owner << "\trow_count=" << forced.rows.size() << '\n';
  for (const auto &row : forced.rows) {
    std::cout << "forced_only\t" << row.multipole_index << '\t' << row.k_index << '\t' << row.ell << '\t' << row.k_mpc_inverse
              << "\tstatus=" << status(row.status) << "\tT="; optional(row.temperature);
    std::cout << "\tE="; optional(row.e_mode);
    std::cout << "\tT_quadrature=" << row.temperature_quadrature_estimate << "\tT_arithmetic=" << row.temperature_arithmetic_estimate
              << "\tE_quadrature=" << row.e_quadrature_estimate << "\tE_arithmetic=" << row.e_arithmetic_estimate
              << "\tradial_error="; optional(row.radial_error_estimate);
    std::cout << "\tsource_grid_error="; optional(row.source_grid_error_estimate);
    std::cout << "\tkernel_evaluations=" << row.kernel_evaluations << "\tcompleted_source_cells=" << row.completed_source_cells
              << "\tsource_admission=false\tfull_temperature=absent\tboundary_composition=absent\n";
  }
  return forced.status == S::ok && same_owner && forced.rows.size() == 2 && std::cout.good() ? 0 : 2;
}
