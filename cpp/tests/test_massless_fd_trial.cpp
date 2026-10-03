// Deliberate private-source regression: exercise actual early refusal/rounding
// code without adding a public control route or repeating a full campaign.
// The production implementation is included once in this test translation unit;
// irred_core supplies only its existing thermal/numerical dependencies.
#include "../src/massless_fd_transfer.cpp"
#include <cstdlib>
#include <iostream>
using namespace irred::cosmology;
void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL " << message << '\n';
    std::exit(1);
  }
}
int main() {
  const auto background =
      prepare_thermal_background({70, 0, 0, 0, .3, {{0, .0002, 2}}});
  require(background.status() == S::ok, "exact retained source");
  const auto radiation = background.scaled_expansion(0);
  require(radiation.status == S::ok, "retained early radiation");
  MasslessFDTransferPolicy policy;
  // Initial138 plus complete RK4 physical assemblies8*133 plus three
  // reduced stage casts3*6; the final six-component cast is refused.
  policy.maximum_scalar_updates_per_point = 1220;
  policy.maximum_total_scalar_updates = 1220;
  MasslessFDTransferWork total, row;
  Budget budget{policy, total, &row};
  Control control;
  control.kind = 5;
  control.round_core = true;
  const W start = W(1e-14) / 4;
  const auto trial = evolve(background, .01L, start, 1e-4L, radiation.a4_e2,
                            radiation.error_estimate, 128, 4, control, budget);
  require(trial.status == S::work_limit && total.scalar_updates == 1220 &&
              total.rhs == 4 && trial.state_initialized,
          "kind5 final cast cap retained");
  require(
      trial.state_a > start && trial.metric_a == start &&
          trial.state_a != trial.metric_a && trial.y[6] > trial.initial.eta_mpc,
      "new combined state epoch survives refusal and cannot use stale lapse");
  require(trial.maximum_phase_increment > 0 &&
              trial.maximum_phase_increment <= policy.maximum_phase_step &&
              trial.maximum_stage_phase <= policy.maximum_phase_step,
          "actual attempted clock and stage phase bounds retained");
  // In radiation x grows as exp(N): h*x_start alone misses actual phase.
  const W old_h = .04L;
  require(std::exp(old_h) - 1 > old_h, "original start-only phase law fails");
  const W conservative_h = old_h * std::exp(-old_h);
  require(conservative_h * std::exp(conservative_h) < old_h &&
              clock_phase_upper(1, 1, std::exp(conservative_h)) < old_h,
          "conservative stage and whole-step phase allocations");
  require(
      !requested_cross_failure(massless_fd_comoving_cdm, 1, 0, 1, .1L) &&
          requested_cross_failure(massless_fd_spatial_potential, 1, 0, 1,
                                  .1L) &&
          requested_cross_failure(7, 1, 0, 1, .1L) &&
          !requested_cross_failure(massless_fd_lapse_potential, 0, 0, 1, .1L) &&
          requested_cross_failure(massless_fd_comoving_cdm, 0, 0, 1, .1L),
      "only requested coordinate cross discrepancy refuses its result");
  for (unsigned j = 1; j <= 32; ++j) {
    MasslessFDTransferValue value;
    std::array<W, 5> terms;
    for (unsigned i = 0; i < terms.size(); ++i)
      terms[i] = W(1e-10) * (1 + W(j + i) * eps_w);
    accept(value, .5L, terms, MasslessFDTransferPolicy{});
    const W published = W(value.time_refinement) + value.initial_refinement +
                        value.hierarchy_refinement +
                        value.background_age_sensitivity +
                        value.arithmetic_cast_sensitivity;
    require(value.status == S::ok && value.value &&
                W(value.absolute_error_estimate) >= published,
            "total bounds every outward published contributor");
  }
  std::cout << "massless_fd_trial_contract passed "
               "epoch/cap/phase/published-error regressions\n";
}
