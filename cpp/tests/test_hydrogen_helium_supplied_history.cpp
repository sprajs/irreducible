#include "hydrogen_helium_supplied_cases.hpp"
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace cases = hydrogen_helium_supplied_cases;
namespace {
unsigned checks = 0;
void need(bool condition, const char *message) {
  ++checks; if (!condition) throw std::runtime_error(message);
}
bool same_bits(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}
double value(const HydrogenHeliumHistoryValue &group) {
  need(group.status == S::ok && group.value.has_value(), "accepted group");
  return *group.value;
}
const HydrogenHeliumHistoryValue *group(const HydrogenHeliumHistoryRow &row, unsigned k) {
  const HydrogenHeliumHistoryValue *items[]{&row.hydrogen_ionized_fraction,
      &row.helium_singly_ionized_fraction, &row.electron_number_density_per_cubic_metre,
      &row.matter_temperature_kelvin, &row.thomson_opacity_per_redshift};
  return items[k];
}
void retained(const HydrogenHeliumHistory &owner, const HydrogenHeliumSuppliedHistoryRequest &source) {
  need(owner.source() && owner.background() && owner.supplied_initial_state() &&
       owner.initial_import_witness() && owner.thermal_mapping_witnesses(), "earned source fields");
  need(owner.source()->nuclei_origin == source.history.nuclei_origin &&
       owner.source()->model.species.size() == source.history.model.species.size() &&
       owner.supplied_initial_state()->origin == source.initial.origin &&
       owner.supplied_initial_state()->source_identity == source.initial.source_identity,
       "owned immutable source identity");
  const auto &w = *owner.initial_import_witness();
  const double expected[]{source.initial.hydrogen_ionized_fraction,
      source.initial.helium_singly_ionized_fraction, source.initial.matter_temperature_kelvin};
  for (unsigned k = 0; k < 3; ++k)
    need(w.promoted_values[k] == static_cast<long double>(expected[k]) &&
         w.measured_absolute_promotion_loss[k] == 0, "exact binary64 import");
}
void refusal_controls() {
  auto source = cases::source();
  auto policy = HydrogenHeliumHistoryPolicy{};
  policy.maximum_total_work = 0;
  auto none = prepare_hydrogen_helium_supplied_history(source, policy);
  need(none.status() == S::work_limit && none.work().total() == 0 &&
       !none.supplied_initial_state() && !none.initial_import_witness(), "zero work before import");
  policy.maximum_total_work = 1;
  auto imported = prepare_hydrogen_helium_supplied_history(source, policy);
  need(imported.status() == S::work_limit && imported.work().total() == 1 &&
       imported.work().background_evaluations == 0 && imported.work().momentum_callbacks == 0 &&
       imported.supplied_initial_state() && imported.initial_import_witness() &&
       !imported.source() && !imported.thermal_mapping_witnesses() && !imported.rate_domain_witness(),
       "earned import survives exhausted remaining work before thermal");
  policy.maximum_total_work = 10;
  auto coefficient_prefix = prepare_hydrogen_helium_supplied_history(source, policy);
  need(coefficient_prefix.status() == S::work_limit && coefficient_prefix.work().total() <= 10 &&
       coefficient_prefix.supplied_initial_state() && coefficient_prefix.initial_import_witness() &&
       coefficient_prefix.source() && coefficient_prefix.thermal_mapping_witnesses() &&
       !coefficient_prefix.rate_domain_witness(), "coefficient failure retains only earned witnesses");
  auto bad_thermal = policy; bad_thermal.maximum_total_work = 4000000;
  bad_thermal.thermal.absolute_tolerance = -1;
  auto map_refusal = prepare_hydrogen_helium_supplied_history(source, bad_thermal);
  need(map_refusal.status() == S::invalid_input && map_refusal.supplied_initial_state() &&
       map_refusal.initial_import_witness() && !map_refusal.source() &&
       !map_refusal.rate_domain_witness(), "supplied source survives later thermal-policy refusal");
  policy = {}; policy.base_intervals = 2; policy.maximum_total_work = 12;
  auto rate_prefix = prepare_hydrogen_helium_supplied_history(source, policy);
  need(rate_prefix.status() != S::ok && rate_prefix.work().total() <= 12 &&
       rate_prefix.rate_domain_witness() && rate_prefix.work().rhs_evaluations > 0 &&
       rate_prefix.rate_domain_witness()->attempted_kelvin_range &&
       !rate_prefix.rate_domain_witness()->retained_kelvin_range,
       "attempted rate range survives partial solver refusal without complete range");
  auto moved_prefix = std::move(rate_prefix);
  need(moved_prefix.rate_domain_witness() && !rate_prefix.rate_domain_witness() &&
       !rate_prefix.boundary_kind() && rate_prefix.model_identity().empty(), "partial move lifetime");
  auto invalid = source;
  for (unsigned coordinate = 0; coordinate < 2; ++coordinate) {
    for (double bad : {0., 1., -std::numeric_limits<double>::denorm_min(), std::nextafter(1., 2.)}) {
      invalid = source;
      (coordinate ? invalid.initial.helium_singly_ionized_fraction :
                    invalid.initial.hydrogen_ionized_fraction) = bad;
      auto rejected = prepare_hydrogen_helium_supplied_history(invalid);
      need(rejected.status() == S::outside_domain && rejected.work().initial_boundary_evaluations == 1 &&
           !rejected.initial_import_witness(), "open fraction boundary refusal");
    }
  }
  for (unsigned coordinate = 0; coordinate < 3; ++coordinate) {
    for (double bad : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
      invalid = source;
      double *items[]{&invalid.initial.hydrogen_ionized_fraction,
          &invalid.initial.helium_singly_ionized_fraction, &invalid.initial.matter_temperature_kelvin};
      *items[coordinate] = bad;
      need(prepare_hydrogen_helium_supplied_history(invalid).status() == S::nonfinite_input,
           "nonfinite supplied coordinate");
    }
  }
  for (double bad : {std::nextafter(4000., 0.), 10000.}) {
    invalid = source; invalid.initial.matter_temperature_kelvin = bad;
    need(prepare_hydrogen_helium_supplied_history(invalid).status() == S::outside_domain,
         "initial temperature scope");
  }
  invalid = source; invalid.initial.origin.clear();
  need(prepare_hydrogen_helium_supplied_history(invalid).status() == S::invalid_input, "required origin");
  invalid = source; invalid.initial.source_identity.clear();
  need(prepare_hydrogen_helium_supplied_history(invalid).status() == S::invalid_input, "required boundary ID");
  invalid = source; invalid.initial.origin.assign(65536, 'o');
  need(prepare_hydrogen_helium_supplied_history(invalid).status() == S::work_limit, "aggregate metadata bound");
  invalid = source; invalid.history.model.species.push_back({.06, 1.95, 2});
  need(prepare_hydrogen_helium_supplied_history(invalid).status() == S::outside_domain, "same massless domain");
  policy = {}; policy.maximum_native_bytes = 1;
  need(prepare_hydrogen_helium_supplied_history(source, policy).status() == S::work_limit, "byte preflight");
  policy = {}; policy.maximum_fine_intervals = 4;
  need(prepare_hydrogen_helium_supplied_history(source, policy).status() == S::work_limit, "unchanged mesh quota");
  policy = {}; policy.absolute_fraction_tolerance = -1;
  need(prepare_hydrogen_helium_supplied_history(source, policy).status() == S::invalid_input, "unchanged allocation policy");
  HydrogenHeliumHistoryWork overflow;
  overflow.rhs_evaluations = std::numeric_limits<std::size_t>::max();
  overflow.initial_boundary_evaluations = 1;
  need(overflow.total() == std::numeric_limits<std::size_t>::max(), "checked public work total");
  need(!hydrogen_helium_history_payload_bound(SIZE_MAX, 0, 0, 1, 1, 1) &&
       !hydrogen_helium_history_payload_bound(1, SIZE_MAX, 0, 1, 1, 1) &&
       !hydrogen_helium_history_payload_bound(1, 0, SIZE_MAX, 1, 1, 1) &&
       !hydrogen_helium_history_payload_bound(1, 0, 0, 1, SIZE_MAX, 1), "checked payload overflow");
  const auto bytes = hydrogen_helium_history_payload_bound(8, 0, 0,
      source.history.nuclei_origin.size(), source.initial.origin.size(), source.initial.source_identity.size());
  need(bool(bytes), "finite payload boundary");
  policy = {}; policy.base_intervals = 2; policy.maximum_total_work = 1;
  policy.maximum_native_bytes = *bytes - 1;
  need(!prepare_hydrogen_helium_supplied_history(source, policy).initial_import_witness(), "one byte below preflight");
  policy.maximum_native_bytes = *bytes;
  need(prepare_hydrogen_helium_supplied_history(source, policy).initial_import_witness(), "exact preflight includes import");
  const int old_round = std::fegetround();
  need(std::fesetround(FE_DOWNWARD) == 0, "set adversarial rounding");
  const auto rounded = prepare_hydrogen_helium_supplied_history(source);
  std::fesetround(old_round);
  need(rounded.status() == S::invalid_input && rounded.work().total() == 0, "strict arithmetic refusal");
  // The exact wide upper bound is tested, never repaired with clipping.
  invalid = source;
  const long double upper = static_cast<long double>(source.history.model.tcmb_kelvin) *
                            (1 + static_cast<long double>(source.history.initial_redshift));
  invalid.initial.matter_temperature_kelvin = static_cast<double>(upper);
  policy = {}; policy.maximum_total_work = 1;
  auto upper_owner = prepare_hydrogen_helium_supplied_history(invalid, policy);
  need(static_cast<long double>(invalid.initial.matter_temperature_kelvin) <= upper
       ? upper_owner.initial_import_witness() != nullptr
       : upper_owner.status() == S::outside_domain, "literal emitted upper-bound meaning");
  invalid.initial.matter_temperature_kelvin = 4000;
  need(prepare_hydrogen_helium_supplied_history(invalid, policy).initial_import_witness(), "4000K initial endpoint admitted");
}
} // namespace
int main() {
  try {
    refusal_controls();
    for (unsigned id = 0; id < 3; ++id) {
      auto source = cases::source(id), saved = source;
      HydrogenHeliumHistoryPolicy policy; policy.base_intervals = 16384;
      auto owner = prepare_hydrogen_helium_supplied_history(source, policy);
      std::cerr << "case=" << id << " status=" << int(owner.status()) << " work=" << owner.work().total() << '\n';
      need(owner.status() == S::ok, "named supplied case prepared");
      retained(owner, saved);
      need(owner.boundary_kind() == HydrogenHeliumHistoryBoundary::supplied_binary64 &&
           owner.model_identity() == hydrogen_helium_supplied_history_model_id &&
           owner.method_identity() == hydrogen_helium_history_method_id &&
           owner.work().initial_boundary_evaluations == 1 && owner.work().initial_charge_evaluations == 0 &&
           owner.work().total() <= 4000000, "distinct boundary identity and same charged kernel");
      need(owner.rate_domain_witness() && owner.rate_domain_witness()->attempted_kelvin_range &&
           owner.rate_domain_witness()->retained_kelvin_range &&
           !owner.rate_domain_witness()->invalid_temperature_attempted, "complete rate witnesses");
      const auto z = cases::queries(id, false);
      auto result = owner.evaluate(z, 31);
      need(result.status == S::ok && result.rows.size() == z.size(), "coarse ordered consumer");
      const auto &first = result.rows.front();
      need(same_bits(value(first.hydrogen_ionized_fraction), saved.initial.hydrogen_ionized_fraction) &&
           same_bits(value(first.helium_singly_ionized_fraction), saved.initial.helium_singly_ionized_fraction) &&
           same_bits(value(first.matter_temperature_kelvin), saved.initial.matter_temperature_kelvin), "exact initial source recovery");
      for (std::size_t j = 0; j < z.size(); ++j) {
        const auto &r = result.rows[j];
        for (unsigned k = 0; k < 5; ++k) value(*group(r, k));
        const long double u = 1 + static_cast<long double>(z[j]),
            ne = u * u * u * (saved.history.hydrogen_nuclei_today_per_cubic_metre * value(r.hydrogen_ionized_fraction) +
                            saved.history.helium_nuclei_today_per_cubic_metre * value(r.helium_singly_ionized_fraction));
        need(r.redshift == z[j] && std::abs(value(r.electron_number_density_per_cubic_metre) - ne) <= 8e-15L * ne,
             "ordered shared-charge row");
        const auto range = *owner.rate_domain_witness()->retained_kelvin_range;
        need(value(r.matter_temperature_kelvin) >= range[0] - 1e-11L &&
             value(r.matter_temperature_kelvin) <= range[1] + 1e-11L, "retained linear-cell temperature enclosure");
      }
      for (unsigned mask = 1; mask < 32; ++mask) {
        const double query[]{saved.history.initial_redshift};
        const auto single = owner.evaluate(query, mask);
        need(single.status == S::ok && single.rows.size() == 1, "all output masks");
        for (unsigned k = 0; k < 5; ++k)
          need(bool(group(single.rows[0], k)->value) == bool(mask & (1u << k)), "masked output identity");
      }
      const auto before = owner.work().total();
      const double arbitrary[]{1300, saved.history.initial_redshift, 1300,
          std::numeric_limits<double>::quiet_NaN(), 299, 2801};
      const auto ordered = owner.evaluate(arbitrary, 31);
      need(ordered.status == S::ok && ordered.rows.size() == 6 &&
           ordered.rows[0].redshift == 1300 && ordered.rows[2].redshift == 1300 &&
           ordered.rows[3].hydrogen_ionized_fraction.status == S::nonfinite_input &&
           ordered.rows[4].hydrogen_ionized_fraction.status == S::outside_domain &&
           ordered.rows[5].hydrogen_ionized_fraction.status == S::outside_domain &&
           owner.work().total() == before, "arbitrary source order and invalid rows without ODE work");
      need(owner.evaluate({}, 31).status == S::ok && owner.evaluate(z, 0).status == S::invalid_input &&
           owner.evaluate(z, 32).status == S::invalid_input && owner.evaluate(z, 31, 1).status == S::work_limit &&
           owner.evaluate(z, 31, 4096, 1).status == S::work_limit, "batch refusals and empty request");
      source.initial = {}; source.history.nuclei_origin.clear(); source.history.model.species.clear();
      auto copied = owner;
      owner = HydrogenHeliumHistory{};
      auto moved = std::move(copied);
      retained(moved, saved);
      need(!copied.boundary_kind() && copied.model_identity().empty() && !copied.supplied_initial_state() &&
           !copied.initial_import_witness() && !copied.rate_domain_witness() && !copied.source() &&
           !copied.thermal_mapping_witnesses(), "complete moved-from invalidation");
      auto *self = &moved; moved = std::move(*self); retained(moved, saved);
      auto replacement = prepare_hydrogen_helium_supplied_history(cases::source((id + 1) % 3), policy);
      replacement = moved; retained(replacement, saved);
      const auto held = replacement.evaluate(z, 31);
      need(held.rows.size() == result.rows.size(), "copy replacement independent storage");
      for (std::size_t j = 0; j < held.rows.size(); ++j) for (unsigned k = 0; k < 5; ++k)
        need(group(held.rows[j], k)->value == group(result.rows[j], k)->value &&
             group(held.rows[j], k)->absolute_error_estimate == group(result.rows[j], k)->absolute_error_estimate,
             "copy/move exact retained values and diagnostics");
    }
    std::cout << "hydrogen_helium_supplied_history checks=" << checks << " failures=0\n";
  } catch (const std::exception &error) { std::cerr << "FAILED " << error.what() << '\n'; return 1; }
}
