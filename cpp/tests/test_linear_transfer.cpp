#include "irred/linear_transfer.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
using namespace irred::cosmology;
using S = irred::numerics::Status;
void need(bool value, const char *message) {
  if (!value) {
    std::cerr << "FAIL " << message << '\n';
    std::exit(1);
  }
}
void close(double value, long double expected, double relative,
           const char *message) {
  if (std::abs(value - expected) >
      relative * std::max(1.L, std::abs(expected))) {
    std::cerr << "value=" << value << " expected=" << (double)expected << '\n';
    need(false, message);
  }
}
int main() {
  auto dust = prepare_thermal_background({70, 0, 0, 0, 1, {}});
  auto owner = prepare_perfect_fluid_transfer(dust, 1e-9);
  need(owner.status() == S::ok, "dust owner");
  const double ks[]{.003, .001, .003};
  auto output = owner.evaluate(ks, 1, 3);
  need(output.status == S::ok && output.rows.size() == 3, "ordered rows");
  for (std::size_t j = 0; j < 3; ++j) {
    auto &r = output.rows[j];
    need(r.comoving_cdm.status == S::ok && r.metric.status == S::ok,
         "EdS accepted");
    close(*r.comoving_cdm.value,
          .4L * ks[j] * ks[j] * std::pow(299792.458L / 70, 2), 2e-7,
          "EdS comoving transfer");
    close(*r.metric.value, -.6L, 1e-12, "EdS metric");
    need(r.wavenumber_mpc_inverse == ks[j] &&
             r.maximum_constraint_residual < 1e-6,
         "order and Einstein gate");
  }
  need(output.rows[0].comoving_cdm.value == output.rows[2].comoving_cdm.value &&
           output.rows[0].callbacks == output.rows[2].callbacks,
       "duplicate ordered replay");
  auto masked = owner.evaluate(std::span(ks, 1), .5, 1);
  need(masked.rows[0].comoving_cdm.value && !masked.rows[0].metric.value,
       "requested mask");
  const double small[]{1e-15};
  auto tiny = owner.evaluate(small, .1, 1);
  need(tiny.rows[0].comoving_cdm.status == S::ok &&
           *tiny.rows[0].comoving_cdm.value > 0,
       "small k cancellation-safe positive transfer");
  const long double tiny_expected =
      .04L * small[0] * small[0] * std::pow(299792.458L / 70, 2);
  need(std::abs(*tiny.rows[0].comoving_cdm.value / tiny_expected - 1) < 2e-7,
       "small k relative analytic limit");
  auto copied = owner;
  auto moved = std::move(copied);
  need(copied.status() == S::invalid_input && !copied.background() &&
           moved.status() == S::ok,
       "move owner");
  auto *self = &moved;
  moved = std::move(*self);
  need(moved.status() == S::ok, "self move");
  auto baryon = prepare_thermal_background({70, 0, 0, .1, .9, {}});
  need(prepare_perfect_fluid_transfer(baryon, 1e-9).status() ==
           S::outside_domain,
       "baryons excluded");
  auto relic = prepare_thermal_background({70, 0, 0, 0, .3, {{0, .0001, 2}}});
  need(prepare_perfect_fluid_transfer(relic, 1e-9).status() ==
           S::outside_domain,
       "relic perturbations excluded");
  need(owner.evaluate(ks, 1, 0).status == S::invalid_input,
       "empty mask invalid");
  need(owner.evaluate(ks, 0, 1).status == S::outside_domain, "target domain");
  const double bad[]{0, std::numeric_limits<double>::quiet_NaN(), 2, .001};
  auto badrows = owner.evaluate(bad, .1, 1);
  need(badrows.rows.size() == 4 &&
           badrows.rows[0].comoving_cdm.status == S::outside_domain &&
           badrows.rows[1].comoving_cdm.status == S::nonfinite_input &&
           badrows.rows[2].comoving_cdm.status == S::outside_domain &&
           badrows.rows[3].comoving_cdm.status == S::ok,
       "refusals retain order");
  TransferPolicy p;
  p.maximum_total_callbacks = 7;
  auto capped = owner.evaluate(ks, .1, 1, p);
  need(capped.callbacks == 7 && capped.rows[0].callbacks == 7 &&
           capped.rows[0].comoving_cdm.status == S::work_limit &&
           capped.rows[1].callbacks == 0 && !capped.rows[0].comoving_cdm.value,
       "global failed work counted");
  p = {};
  p.maximum_native_bytes = 1;
  need(owner.evaluate(ks, .1, 1, p).status == S::work_limit,
       "payload admission");
  auto rad = prepare_perfect_fluid_transfer(
      prepare_thermal_background({70, 1, 0, 0, 0, {}}), 1e-9);
  need(rad.evaluate(std::span(ks, 1), .1, 1).rows[0].comoving_cdm.status ==
           S::ok,
       "radiation tracer allowed");
  PrimordialBand law{2.1e-9, .965, .05, .001, .002};
  BandVariancePolicy vp;
  vp.base_panels = 16;
  auto amplitude = owner.sigma8_band(law, 1, vp);
  need(amplitude.status == S::ok && amplitude.variance && amplitude.sigma &&
           amplitude.primordial_support_is_finite_band &&
           amplitude.primordial.maximum_mpc_inverse == .002,
       "finite band amplitude");
  close(amplitude.radius_mpc, 8 / .7, 1e-14, "8/h Mpc window units");
  need(std::isnormal(amplitude.sigma_absolute_error_estimate) &&
           amplitude.absolute_error_estimate > 0,
       "sigma diagnostic has separate units");
  law.amplitude *= 4;
  auto scaled = owner.sigma8_band(law, 1, vp);
  need(scaled.status == S::ok, "scaled band accepted");
  close(*scaled.variance, 4 * (*amplitude.variance), 1e-13,
        "primordial variance scaling");
  close(*scaled.sigma, 2 * (*amplitude.sigma), 1e-13,
        "primordial sigma scaling");
  need(rad.sigma8_band(law, .1, vp).status != S::ok,
       "no matter variance without CDM");
  TransferPolicy strict;
  strict.absolute_tolerance = 0;
  strict.relative_tolerance = 1e-18;
  auto unattainable = owner.evaluate(std::span(ks, 1), .1, 1, strict);
  need(unattainable.rows[0].comoving_cdm.status ==
               S::conditioning_budget_exceeded &&
           !unattainable.rows[0].comoving_cdm.value,
       "unattainable accuracy withheld");
  auto late_start = prepare_perfect_fluid_transfer(
      prepare_thermal_background({70, .000085, 0, 0, .3, {}}), 1e-3);
  need(late_start.evaluate(std::span(ks, 1), .1, 1)
               .rows[0]
               .comoving_cdm.status == S::outside_domain,
       "unqualified radiation start refused without changing start");
  law.minimum_mpc_inverse = .5;
  law.maximum_mpc_inverse = 1;
  need(owner.band_variance(law, 1, 1000, vp).status == S::outside_domain,
       "unresolved window phase refused");
  std::cout << "linear_transfer_contract passed "
               "EdS/order/masks/units/lifetime/invalid/work/band controls\n";
}
