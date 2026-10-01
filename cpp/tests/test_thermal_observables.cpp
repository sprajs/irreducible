#include "irred/thermal_observables.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void require(bool b, const char *message) {
  ++checks;
  if (!b)
    throw std::runtime_error(message);
}
unsigned all = (thermal_ruler_mask << 1) - 1;
const EarlyLateValue &value(const ThermalObservableBatch &b, std::size_t row,
                            EarlyLateOutput x) {
  return b.rows[row].outputs[static_cast<unsigned>(x)];
}
void close(long double actual, long double expected, long double abs = 1e-8L,
           long double rel = 2e-10L) {
  require(std::abs(actual - expected) <= abs + rel * std::abs(expected),
          "frozen numerical allocation");
}
} // namespace
int main() {
  ThermalObservablePolicy p;
  ThermalObservableRequest source{};
  source.model.h0_km_s_mpc = 67.4;
  source.model.physical_baryon_density = .0224;
  source.model.physical_cdm_density = .12;
  source.model.tcmb_kelvin = 2.7255;
  source.model.physical_massless_nonphoton_density = 1.68e-5;
  source.model.species.push_back({.06, 1.95, 2});
  source.z_drag = 1059.95;
  source.drag_origin = "supplied-test-drag";
  source.source_origin = "explicit-synthetic-test-source";
  const auto mapping = map_thermal_physical_model(source.model);
  require(mapping.status == S::ok && mapping.model.has_value(),
          "physical source maps explicitly");
  close(mapping.model->omega_b, .0224L / (.674L * .674L), 0, 2e-15L);
  close(mapping.model->omega_cdm, .12L / (.674L * .674L), 0, 2e-15L);
  close(mapping.model->species[0].temperature_today_ev,
        1.95L * 1.380649e-23L / 1.602176634e-19L, 0, 2e-15L);
  auto owner = prepare_thermal_observables(source, p);
  require(owner.status() == S::ok, "thermal observable owner prepared");
  const double zs[]{0, .1, 1, 10};
  auto batch = owner.evaluate(zs, all, p);
  require(batch.status == S::ok && batch.ruler && batch.ruler->status == S::ok,
          "shared ruler accepted");
  for (const auto &row : batch.rows)
    for (const auto &v : row.outputs)
      require(v.status == S::ok && v.value.has_value(),
              "all requested coordinates accepted");
  require(batch.callbacks == batch.outer_callbacks + batch.momentum_callbacks &&
              batch.outer_callbacks > 0 && batch.momentum_callbacks > 0 &&
              batch.callbacks <= p.maximum_total_callbacks,
          "nested callbacks accounted separately and together");
  for (auto id : {EarlyLateOutput::dm_mpc, EarlyLateOutput::dl_mpc,
                  EarlyLateOutput::dv_mpc, EarlyLateOutput::dm_over_rs,
                  EarlyLateOutput::dv_over_rs})
    require(*value(batch, 0, id).value == 0, "present geometric zero controls");
  require(*value(batch, 0, EarlyLateOutput::e).value == 1,
          "present exact E normalization");
  require(*value(batch, 0, EarlyLateOutput::h_km_s_mpc).value == 67.4,
          "present exact H normalization");
  for (std::size_t i = 1; i < 4; ++i) {
    close(*value(batch, i, EarlyLateOutput::dl_mpc).value,
          (1 + zs[i]) * *value(batch, i, EarlyLateOutput::dm_mpc).value);
    close(*value(batch, i, EarlyLateOutput::dm_over_rs).value,
          *value(batch, i, EarlyLateOutput::dm_mpc).value / *batch.ruler->value,
          1e-10L, 5e-10L);
  }
  auto massless_source = source;
  massless_source.model.species[0].mass_ev = 0;
  const auto massless_owner = prepare_thermal_observables(massless_source, p);
  const auto &mapped = massless_owner.background().source();
  const SoundHorizonRequest old_source{
      {mapped.h0_km_s_mpc, mapped.omega_b + mapped.omega_cdm,
       mapped.omega_gamma + mapped.omega_massless_nonphoton +
           *massless_owner.background().omega_species_today(),
       mapped.omega_b, mapped.omega_gamma},
      massless_source.z_drag,
      massless_source.drag_origin};
  const EarlyLatePolicy old_policy{
      p.absolute_tolerance_mpc,
      p.relative_tolerance,
      p.absolute_tolerance_ratio,
      p.relative_tolerance_ratio,
      p.maximum_callbacks_per_point,
      p.maximum_total_callbacks,
      p.maximum_depth,
      p.maximum_points,
      p.maximum_native_bytes,
      {p.absolute_tolerance_mpc, p.relative_tolerance,
       p.maximum_callbacks_per_point, p.maximum_depth, p.maximum_points,
       p.maximum_total_callbacks, p.maximum_native_bytes}};
  const auto old = evaluate_early_late(
      old_source, zs, (1u << early_late_output_count) - 1, old_policy);
  const auto massless = massless_owner.evaluate(zs, all, p);
  massless_source.model.species[0].mass_ev = 1e-12;
  const auto small_mass =
      prepare_thermal_observables(massless_source, p).evaluate(zs, all, p);
  require(old.ruler && old.ruler->sound_horizon_mpc && massless.ruler &&
              massless.ruler->value && small_mass.ruler &&
              small_mass.ruler->value,
          "shared massless and tiny-relic rulers available");
  close(*massless.ruler->value, *old.ruler->sound_horizon_mpc);
  close(*small_mass.ruler->value, *massless.ruler->value);
  for (std::size_t row = 0; row < 4; ++row)
    for (unsigned i = 0; i < early_late_output_count; ++i) {
      const auto &prior = old.rows[row].outputs[i];
      const auto &now = massless.rows[row].outputs[i];
      const auto &tiny = small_mass.rows[row].outputs[i];
      require(prior.status == S::ok && now.status == S::ok &&
                  tiny.status == S::ok,
              "massless and tiny-relic shared-state coordinates admitted");
      const long double absolute = i < 2 ? 0 : (i < 6 ? 1e-8L : 1e-10L);
      const long double relative = i < 6 ? 2e-10L : 5e-10L;
      close(*now.value, *prior.value, absolute, relative);
      close(*tiny.value, *now.value, absolute, relative);
    }
  auto no_ruler =
      owner.evaluate(zs, early_late_mask(EarlyLateOutput::dm_mpc), p);
  require(!no_ruler.ruler && !value(no_ruler, 1, EarlyLateOutput::e).value,
          "unrequested scalar and ruler absent");
  const auto dh_only =
      owner.evaluate({zs + 1, 1}, early_late_mask(EarlyLateOutput::dh_mpc), p);
  require(dh_only.outer_callbacks == 0 && dh_only.momentum_callbacks > 0,
          "DH only has no distance quadrature");
  const auto ruler_only = owner.evaluate(zs, thermal_ruler_mask, p);
  const auto ruler_empty = owner.evaluate({}, thermal_ruler_mask, p);
  require(ruler_only.callbacks == ruler_empty.callbacks &&
              ruler_only.ruler->value == ruler_empty.ruler->value &&
              ruler_only.rows[1].callbacks == 0,
          "ruler only has no extraneous redshift work");
  auto with_e = owner.evaluate({zs + 1, 1},
                               early_late_mask(EarlyLateOutput::dm_mpc) |
                                   early_late_mask(EarlyLateOutput::e),
                               p);
  auto dm_only =
      owner.evaluate({zs + 1, 1}, early_late_mask(EarlyLateOutput::dm_mpc), p);
  require(dm_only.momentum_callbacks < with_e.momentum_callbacks &&
              dm_only.outer_callbacks == with_e.outer_callbacks,
          "DM only omits point expansion work");
  auto scalar = owner.evaluate({zs, 1}, early_late_mask(EarlyLateOutput::e), p);
  require(scalar.callbacks == 0 && scalar.outer_callbacks == 0,
          "present E only needs retained normalization");
  auto q = p;
  q.thermal.maximum_total_callbacks = 4000000;
  const auto inherited = owner.evaluate(zs, all, q);
  require(inherited.momentum_callbacks <= 4000000 &&
              value(inherited, 3, EarlyLateOutput::dm_mpc).status ==
                  S::work_limit,
          "explicit smaller inherited background quota refuses full consumer "
          "batch");
  q = p;
  q.maximum_total_callbacks = 0;
  auto partial =
      owner.evaluate({zs, 1},
                     early_late_mask(EarlyLateOutput::e) |
                         early_late_mask(EarlyLateOutput::dm_over_rs),
                     q);
  require(value(partial, 0, EarlyLateOutput::e).status == S::ok &&
              partial.ruler->status == S::work_limit &&
              value(partial, 0, EarlyLateOutput::dm_over_rs).status ==
                  S::work_limit,
          "ruler failure retains independent present E");
  q = p;
  q.maximum_callbacks_per_point = 7;
  auto limited = owner.evaluate({zs + 1, 1}, all, q);
  require(limited.callbacks <= 14 && limited.ruler->status == S::work_limit,
          "per-point nested quota includes outer work");
  q = p;
  q.thermal.maximum_total_callbacks = 13;
  limited = owner.evaluate({zs + 1, 1}, all, q);
  require(limited.momentum_callbacks <= 13,
          "momentum-only global quota enforced across ruler and rows");
  const auto payload = thermal_observables_payload_bound(
      4, 1, source.source_origin.size() + source.drag_origin.size() + 2);
  require(payload.has_value(), "payload bound exists");
  q = p;
  q.maximum_native_bytes = *payload - 1;
  require(owner.evaluate(zs, all, q).status == S::work_limit,
          "payload admission precedes rows");
  q = p;
  q.maximum_native_bytes = 0;
  require(prepare_thermal_observables(source, q).status() == S::work_limit,
          "preparation admission before source copy");
  q = p;
  q.maximum_points = 0;
  require(owner.evaluate(zs, all, q).rows.empty(),
          "point quota precedes allocation");
  q = p;
  q.thermal.maximum_species = 0;
  require(owner.evaluate(zs, all, q).status == S::work_limit,
          "retained species quota still applies");
  const double invalid[]{-1, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()};
  limited =
      owner.evaluate(invalid, early_late_mask(EarlyLateOutput::dm_mpc), p);
  require(limited.rows.size() == 3 &&
              value(limited, 0, EarlyLateOutput::dm_mpc).status ==
                  S::outside_domain &&
              value(limited, 1, EarlyLateOutput::dm_mpc).status ==
                  S::nonfinite_input,
          "invalid rows retained in source order");
  require(owner.evaluate(zs, 0, p).status == S::invalid_input &&
              owner.evaluate(zs, all + 1, p).status == S::invalid_input,
          "invalid masks");
  q = p;
  q.relative_tolerance = 0;
  q.absolute_tolerance_mpc = 1e-30;
  limited =
      owner.evaluate({zs + 1, 1}, early_late_mask(EarlyLateOutput::dh_mpc), q);
  require(value(limited, 0, EarlyLateOutput::dh_mpc).status ==
              S::conditioning_budget_exceeded,
          "unattainable tolerance refused");
  for (int field = 0; field < 4; ++field) {
    auto bad = source;
    if (field == 0)
      bad.source_origin.clear();
    if (field == 1)
      bad.drag_origin.clear();
    if (field == 2)
      bad.z_drag = -1;
    if (field == 3)
      bad.model.tcmb_kelvin = 0;
    require(prepare_thermal_observables(bad, p).status() == S::outside_domain,
            "invalid source identity/physical domain refused");
  }
  auto overflow_source = source;
  overflow_source.model = {1e100, 0, 0, 1e50, 0, {}};
  auto extreme = prepare_thermal_observables(overflow_source, p);
  require(extreme.status() == S::ok,
          "extreme dimensional source is explicit and admitted");
  const double extreme_z[]{1e110};
  const auto projected =
      extreme.evaluate(extreme_z,
                       early_late_mask(EarlyLateOutput::e) |
                           early_late_mask(EarlyLateOutput::h_km_s_mpc) |
                           early_late_mask(EarlyLateOutput::dh_mpc),
                       p);
  require(
      value(projected, 0, EarlyLateOutput::e).status == S::ok &&
          value(projected, 0, EarlyLateOutput::h_km_s_mpc).status ==
              S::outside_domain &&
          value(projected, 0, EarlyLateOutput::dh_mpc).status ==
              S::outside_domain,
      "dimensionless E survives dimensional overflow and positive underflow");
  const double tiny_z[]{std::numeric_limits<double>::denorm_min()};
  const auto tiny =
      owner.evaluate(tiny_z, early_late_mask(EarlyLateOutput::dm_mpc), p);
  require(value(tiny, 0, EarlyLateOutput::dm_mpc).status == S::outside_domain,
          "positive distance underflow refused");
  auto tiny_h_source = source;
  tiny_h_source.model = {
      std::numeric_limits<double>::denorm_min(), 0, 0, 1e-162, 0, {}};
  auto tiny_h_owner = prepare_thermal_observables(tiny_h_source, p);
  require(tiny_h_owner.status() == S::ok,
          "subnormal H source has an explicit admitted photon fraction");
  const auto tiny_h =
      tiny_h_owner.evaluate({zs, 1},
                            early_late_mask(EarlyLateOutput::e) |
                                early_late_mask(EarlyLateOutput::h_km_s_mpc) |
                                early_late_mask(EarlyLateOutput::dh_mpc),
                            p);
  require(value(tiny_h, 0, EarlyLateOutput::e).status == S::ok &&
              *value(tiny_h, 0, EarlyLateOutput::e).value == 1 &&
              value(tiny_h, 0, EarlyLateOutput::h_km_s_mpc).status ==
                  S::outside_domain &&
              value(tiny_h, 0, EarlyLateOutput::dh_mpc).status ==
                  S::outside_domain,
          "exact present identities preserve dimensional representability");
  auto failed_source = source;
  failed_source.model.physical_baryon_density = 2;
  auto failed_owner = prepare_thermal_observables(failed_source, p);
  require(failed_owner.status() == S::outside_domain &&
              failed_owner.source().source_origin == source.source_origin,
          "failed physical closure preserves admitted source provenance");
  auto radiation = prepare_thermal_background({70, 1, 0, 0, 0, {}});
  const auto origin = radiation.scaled_expansion(0);
  require(origin.status == S::ok && origin.a4_e2 == 1 && origin.callbacks == 0,
          "finite exact radiation endpoint");
  require(radiation.evaluate(std::array<double, 1>{0}, thermal_e)
                  .rows[0]
                  .e.status == S::outside_domain,
          "origin is not fabricated finite E");
  for (double a : {0., 1e-300, .01, 1.}) {
    const auto v = radiation.scaled_expansion(a);
    require(v.status == S::ok && v.a4_e2 == 1,
            "scaled massless radiation limit");
  }
  for (double matter : {0., 1.}) {
    auto zero = prepare_thermal_background({70, 0, 0, matter, 0, {}});
    const auto endpoint = zero.scaled_expansion(0);
    require(endpoint.status == S::ok && endpoint.a4_e2 == 0 &&
                endpoint.error_estimate == 0,
            "radiation-free scaled coordinate has finite exact zero endpoint");
    const double a[]{.5};
    const auto finite = zero.evaluate(a, thermal_e);
    require(finite.rows[0].e.status == S::ok,
            "radiation-free E remains finite at positive a");
    close(*finite.rows[0].e.value, matter == 0 ? 1.L : std::sqrt(8.L), 0,
          2e-10L);
  }
  const auto pure_lambda = prepare_thermal_background({70, 0, 0, 0, 0, {}});
  for (long double a : {1e-1233L, 1e-2000L})
    require(pure_lambda.scaled_expansion(a).status == S::outside_domain,
            "positive wide-coordinate/diagnostic underflow refused");
  static_assert(std::is_nothrow_move_constructible_v<ThermalObservables>);
  static_assert(std::is_nothrow_move_assignable_v<ThermalObservables>);
  auto copy = owner;
  auto moved = std::move(owner);
  require(owner.status() == S::invalid_input &&
              owner.source().model.species.empty() &&
              owner.background().status() == S::invalid_input,
          "moved source fully invalidated");
  require(owner.evaluate(zs, all, p).rows.empty(),
          "moved source emits no rows");
  auto compare_owner = [&](const ThermalObservables &o) {
    const auto result =
        o.evaluate({zs + 1, 1}, early_late_mask(EarlyLateOutput::dh_mpc), p);
    require(value(result, 0, EarlyLateOutput::dh_mpc).value ==
                    value(batch, 1, EarlyLateOutput::dh_mpc).value &&
                o.source().source_origin == source.source_origin,
            "owner lifetime retains physical identity and projections");
  };
  compare_owner(copy);
  compare_owner(moved);
  ThermalObservables assigned;
  assigned = std::move(moved);
  require(moved.status() == S::invalid_input,
          "move assignment invalidates old source");
  auto *alias = &assigned;
  assigned = std::move(*alias);
  compare_owner(assigned);
  const int rounding = std::fegetround();
  std::fesetround(FE_UPWARD);
  const auto unsupported = assigned.evaluate(zs, all, p);
  std::fesetround(rounding);
  require(unsupported.rows.empty() && unsupported.status == S::invalid_input,
          "unsupported arithmetic refused");
  std::cout << "thermal_observables owner checks=" << checks
            << " callbacks=" << batch.callbacks << '\n';
}
