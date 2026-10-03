#include "irred/massless_fd_transfer.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>
using namespace irred::cosmology;
using S = irred::numerics::Status;
void need(bool v, const char *label) {
  if (!v) {
    std::cerr << "FAIL " << label << '\n';
    std::exit(1);
  }
}
int main() {
  const auto background =
      prepare_thermal_background({70, 0, 0, 0, .3, {{0, .0002, 2}}});
  const auto owner = prepare_massless_fd_transfer(background, 1e-14);
  need(owner.status() == S::ok &&
           owner.background()->source().species.size() == 1,
       "explicit massless source retained");
  need(prepare_massless_fd_transfer(
           prepare_thermal_background({70, 0, 0, 0, .3, {}}), 1e-14)
               .status() == S::outside_domain,
       "zero-species EdS is a separate mathematical limit");
  need(prepare_massless_fd_transfer(
           prepare_thermal_background({70, .00001, 0, 0, .3, {{0, .0002, 2}}}),
           1e-14)
               .status() == S::outside_domain,
       "photons excluded");
  need(prepare_massless_fd_transfer(
           prepare_thermal_background({70, 0, .00001, 0, .3, {{0, .0002, 2}}}),
           1e-14)
               .status() == S::outside_domain,
       "additional massless density excluded");
  need(prepare_massless_fd_transfer(
           prepare_thermal_background({70, 0, 0, .04, .3, {{0, .0002, 2}}}),
           1e-14)
               .status() == S::outside_domain,
       "baryons excluded");
  need(prepare_massless_fd_transfer(
           prepare_thermal_background({70, 0, 0, 0, .3, {{.01, .0002, 2}}}),
           1e-14)
               .status() == S::outside_domain,
       "positive mass excluded");
  need(prepare_massless_fd_transfer(background, 3e-20).status() ==
           S::outside_domain,
       "start boundary");
  auto copied = owner;
  auto moved = std::move(copied);
  need(copied.status() == S::invalid_input && !copied.background() &&
           moved.status() == S::ok,
       "moved owner");
  auto *same = &moved;
  moved = std::move(*same);
  need(moved.status() == S::ok, "self move");
  const double k[]{.01};
  need(owner.evaluate(k, 1e-4, 0).status == S::invalid_input, "empty mask");
  need(owner.evaluate(k, 1e-4, 8).status == S::invalid_input, "unknown mask");
  need(owner.evaluate(k, 1e-5, 1).status == S::outside_domain, "target domain");
  need(owner.evaluate(k, std::numeric_limits<double>::quiet_NaN(), 1).status ==
           S::nonfinite_input,
       "target finite");
  MasslessFDTransferPolicy p;
  p.maximum_native_bytes = 1;
  need(owner.evaluate(k, 1e-4, 1, p).status == S::work_limit,
       "payload admission");
  p = {};
  p.absolute_tolerance = 2e-7;
  need(owner.evaluate(k, 1e-4, 1, p).status == S::invalid_input,
       "budget cannot widen");
  p = {};
  p.maximum_rhs_per_point = 0;
  const double ordered[]{std::numeric_limits<double>::quiet_NaN(), 0, .01, .01};
  const auto trials =
      owner.evaluate(ordered, 1e-4, massless_fd_comoving_cdm, p);
  need(trials.status == S::ok && trials.rows.size() == 4,
       "original order retained");
  need(trials.rows[0].comoving_cdm.status == S::nonfinite_input &&
           trials.rows[1].comoving_cdm.status == S::outside_domain,
       "invalid original rows retained");
  for (unsigned j = 2; j < 4; ++j) {
    const auto &row = trials.rows[j];
    need(row.comoving_cdm.status == S::work_limit && !row.comoving_cdm.value &&
             !row.spatial_potential.value && row.work.rhs == 0 &&
             row.attempts_recorded == 1,
         "failed first attempt and requested mask");
    const auto &attempt = row.attempts[0];
    const auto &w = attempt.initial;
    need(attempt.endpoint_available && attempt.lapse_available &&
             w.status == S::ok && attempt.hierarchy_l == 128 &&
             attempt.time_divisor == 4 &&
             attempt.endpoint_scale_factor == attempt.metric_epoch_scale_factor,
         "actual projected initializer survives refusal");
    const long double x = ordered[j] * w.eta_mpc;
    need(std::abs(w.projected_delta / (x * x) - 7.L / 19) < 2e-8L &&
             std::abs(w.leading_delta_c - 15.L / 19) < 1e-18L &&
             w.projected_delta > 0 && w.leading_psi < 0 && w.projected_phi < 0,
         "projected 7/19 is distinct from signed Newtonian delta and regular "
         "1/4");
    need(std::abs(w.hamiltonian_after) < 1e-18L &&
             row.work.scalar_updates == 138,
         "actual L128 initialization assignments counted");
  }
  need(trials.rows[2].attempts[0].initial.projected_delta ==
           trials.rows[3].attempts[0].initial.projected_delta,
       "duplicate projected state replay");
  p = {};
  p.maximum_total_rhs = 4;
  const auto capped = owner.evaluate(k, 1e-4, 7, p);
  need(capped.rows[0].comoving_cdm.status == S::work_limit &&
           capped.work.rhs == 4 && capped.rows[0].attempts_recorded == 1 &&
           capped.rows[0].attempts[0].endpoint_available,
       "all attempted stage work retained at global cap");
  for (unsigned mask : {massless_fd_comoving_cdm, massless_fd_spatial_potential,
                        massless_fd_lapse_potential}) {
    const auto required = owner.evaluate(k, 1e-4, mask, p);
    const auto &r = required.rows[0];
    const auto *field =
        mask == massless_fd_comoving_cdm
            ? &r.comoving_cdm
            : (mask == massless_fd_spatial_potential ? &r.spatial_potential
                                                     : &r.lapse_potential);
    need(field->status == S::work_limit && !field->value &&
             required.work.rhs == capped.work.rhs &&
             required.work.scalar_updates == capped.work.scalar_updates &&
             required.work.background_queries == capped.work.background_queries,
         "single requested coordinate retains mandatory coupled work refusal");
    auto invalid = p;
    invalid.maximum_constraint_residual = 0;
    need(owner.evaluate(k, 1e-4, mask, invalid).status == S::invalid_input,
         "output mask cannot disable shared constraint admission");
  }
  std::cout << "massless_fd_transfer_contract passed "
               "source/modes/lifetime/masks/refusal/work controls\n";
}
