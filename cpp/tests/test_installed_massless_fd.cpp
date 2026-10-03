#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <irred/massless_fd_transfer.hpp>
// Actual independently installed SDK caller, also compiled by GCC/Clang native
// CI. No test-only reference header, ABI/CLI route, or fitted input is
// involved.
int main() {
  using namespace irred::cosmology;
  using S = irred::numerics::Status;
  const auto background =
      prepare_thermal_background({70, 0, 0, 0, .3, {{0, .0002, 2}}});
  const auto owner = prepare_massless_fd_transfer(background, 1e-14);
  const std::array<double, 2> k{1e-7, .01};
  const auto out =
      owner.evaluate(k, 1e-4,
                     massless_fd_comoving_cdm | massless_fd_spatial_potential |
                         massless_fd_lapse_potential);
  std::cout << "model=" << massless_fd_transfer_model_id
            << " method=" << massless_fd_transfer_method_id
            << " arithmetic=" << massless_fd_transfer_arithmetic_id << '\n';
  std::cout << std::setprecision(21) << "batch_status=" << int(out.status)
            << " rhs=" << out.work.rhs << " updates=" << out.work.scalar_updates
            << " background=" << out.work.background_queries
            << " age_nodes=" << out.work.age_quadrature_evaluations << '\n';
  if (owner.status() != S::ok || out.status != S::ok ||
      out.rows.size() != k.size())
    return 1;
  for (std::size_t j = 0; j < k.size(); ++j) {
    const auto &r = out.rows[j];
    std::cout << "k=" << r.wavenumber_mpc_inverse
              << " attempts=" << r.attempts_recorded
              << " constraint=" << r.maximum_constraint_residual
              << " chi=" << r.maximum_closure_defect
              << " Rbg=" << r.maximum_background_identity_defect << '\n';
    for (unsigned n = 0; n < r.attempts_recorded; ++n) {
      const auto &t = r.attempts[n];
      std::cout << "attempt=" << n << " status=" << int(t.status)
                << " L=" << t.hierarchy_l << " time_divisor=" << t.time_divisor
                << " control=" << t.control_kind << " sign=" << t.control_sign
                << " a=" << t.endpoint_scale_factor
                << " metric_a=" << t.metric_epoch_scale_factor
                << " lapse=" << t.lapse_available << " rhs=" << t.work.rhs
                << " updates=" << t.work.scalar_updates
                << " stage_phase=" << t.maximum_stage_phase_bound
                << " clock_phase=" << t.maximum_phase_increment_upper << '\n';
      if (t.status == S::ok &&
          (!t.lapse_available ||
           t.endpoint_scale_factor != t.metric_epoch_scale_factor ||
           t.maximum_stage_phase_bound > static_cast<long double>(.04) ||
           t.maximum_phase_increment_upper > static_cast<long double>(.04)))
        return 4;
    }
    if (r.wavenumber_mpc_inverse != k[j] || r.attempts_recorded != 17 ||
        r.maximum_constraint_residual > 1e-7)
      return 2;
    for (const auto *v :
         {&r.comoving_cdm, &r.spatial_potential, &r.lapse_potential}) {
      std::cout << "field_status=" << int(v->status);
      if (v->value)
        std::cout << " value=" << *v->value;
      std::cout << " error=" << v->absolute_error_estimate
                << " time=" << v->time_refinement
                << " start=" << v->initial_refinement
                << " L=" << v->hierarchy_refinement
                << " background_age=" << v->background_age_sensitivity
                << " cast=" << v->arithmetic_cast_sensitivity << '\n';
      if (v->status != S::ok || !v->value || !std::isfinite(*v->value) ||
          !std::isnormal(v->absolute_error_estimate) ||
          v->absolute_error_estimate >
              5 * (1e-7 + 1e-4 * std::abs(*v->value)) / 6)
        return 3;
    }
  }
  std::cout << "installed pure explicit massless-FD control passed; full "
               "standard transfer/CMB remain open\n";
}
