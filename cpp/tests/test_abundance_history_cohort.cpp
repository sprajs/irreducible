#include "abundance_history_cohort_wire.hpp"
#include "hydrogen_helium_history_reference.hpp"
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string_view>
namespace allocation_observation {
struct alignas(std::max_align_t) Header { void *base; std::size_t bytes; };
inline bool armed = false;
inline long fail_at = -1;
inline std::size_t calls = 0, live = 0, peak = 0;
inline bool fired = false;
void *allocate(std::size_t n, std::size_t alignment) {
  if (armed) {
    const auto index = calls++;
    if (fail_at >= 0 && index == std::size_t(fail_at)) { fired = true; throw std::bad_alloc(); }
  }
  alignment = std::max(alignment, alignof(Header));
  if (n > SIZE_MAX - sizeof(Header) - alignment) throw std::bad_alloc();
  void *base = std::malloc((n ? n : 1) + sizeof(Header) + alignment);
  if (!base) throw std::bad_alloc();
  const auto start = reinterpret_cast<std::uintptr_t>(base) + sizeof(Header);
  const auto address = (start + alignment - 1) & ~(std::uintptr_t(alignment) - 1);
  auto *header = reinterpret_cast<Header *>(address) - 1;
  header->base = base; header->bytes = n;
  live += n; if (armed) peak = std::max(peak, live);
  return reinterpret_cast<void *>(address);
}
void release(void *p) noexcept {
  if (!p) return;
  auto *header = static_cast<Header *>(p) - 1;
  live -= header->bytes; std::free(header->base);
}
void arm(long site = -1) { calls = 0; peak = live; fired = false; fail_at = site; armed = true; }
void disarm() { armed = false; fail_at = -1; }
}
#if defined(__GNUC__) || defined(__clang__)
#define COHORT_NOINLINE __attribute__((noinline))
#else
#define COHORT_NOINLINE
#endif
COHORT_NOINLINE void *operator new(std::size_t n) { return allocation_observation::allocate(n, alignof(std::max_align_t)); }
COHORT_NOINLINE void *operator new[](std::size_t n) { return ::operator new(n); }
COHORT_NOINLINE void *operator new(std::size_t n, std::align_val_t a) { return allocation_observation::allocate(n, std::size_t(a)); }
COHORT_NOINLINE void *operator new[](std::size_t n, std::align_val_t a) { return ::operator new(n, a); }
COHORT_NOINLINE void operator delete(void *p) noexcept { allocation_observation::release(p); }
COHORT_NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
COHORT_NOINLINE void operator delete(void *p, std::size_t) noexcept { ::operator delete(p); }
COHORT_NOINLINE void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }
COHORT_NOINLINE void operator delete(void *p, std::align_val_t) noexcept { ::operator delete(p); }
COHORT_NOINLINE void operator delete[](void *p, std::align_val_t) noexcept { ::operator delete(p); }
COHORT_NOINLINE void operator delete(void *p, std::size_t, std::align_val_t) noexcept { ::operator delete(p); }
COHORT_NOINLINE void operator delete[](void *p, std::size_t, std::align_val_t) noexcept { ::operator delete(p); }
namespace {
namespace h = abundance_history_cohort;
namespace ref = hydrogen_helium_history_reference;
using S = h::S; using W = h::W;
unsigned checks = 0;
void need(bool value, const char *message) { ++checks; if (!value) throw std::runtime_error(message); }
void algebra() {
  std::array<std::array<W, h::axes>, 2> f{{{2, 8, 4, 10, 6, 12}, {6, 4, 8, 6, 10, 8}}}, e{};
  std::array<bool, h::axes> constant{};
  h::EvaluationWork work; std::optional<h::Moments> m;
  need(h::reduce(f, e, constant, .01, work, 1024, m) == S::ok && m.has_value(), "exact two-half toy law");
  // These exact-rational expected integers are independently specified; no
  // engine history/map output or native reducer defines the expectation.
  const std::array<double, h::axes> expected_mean{4, 6, 6, 8, 8, 10};
  for (std::size_t i = 0; i < h::axes; ++i) {
    need(m->mean[i] == expected_mean[i], "exact rational half mean");
    for (std::size_t j = 0; j < h::axes; ++j) {
      const double expected = i % 2 == j % 2 ? 4 : -4;
      need(m->covariance[i * h::axes + j] == expected, "signed full cross-row population covariance");
      need(h::bits(m->covariance[i * h::axes + j]) == h::bits(m->covariance[j * h::axes + i]), "one stored symmetric triangle");
      need(m->covariance_error[i * h::axes + j] > 0 && m->covariance_error[i * h::axes + j] < .04, "operation/cast error separate from spread");
    }
  }
  need(work.means == 6 && work.differences == 6 && work.products == 21 && work.projections == 27, "actual bounded reduction categories");
  auto reverse = f; std::swap(reverse[0], reverse[1]);
  h::EvaluationWork reverse_work; std::optional<h::Moments> reversed;
  need(h::reduce(reverse, e, constant, .01, reverse_work, 1024, reversed) == S::ok, "support reversal");
  need(reversed->mean == m->mean && reversed->covariance == m->covariance, "support reversal preserves the same law");
  const W dyadic_error = std::ldexp(W(1), -20);
  auto noisy = e; for (auto &state : noisy) state.fill(dyadic_error);
  h::EvaluationWork noisy_work; std::optional<h::Moments> noisy_moments;
  need(h::reduce(f, noisy, constant, .01, noisy_work, 1024, noisy_moments) == S::ok,
       "nonconstant exact-dyadic inherited diagnostics");
  const W inherited_covariance = 4 * dyadic_error + dyadic_error * dyadic_error;
  for (double x : noisy_moments->mean_error) need(W(x) >= dyadic_error, "half-law inherited mean error retained");
  for (double x : noisy_moments->covariance_error) need(W(x) >= inherited_covariance,
       "signed delta-four law independently requires four-epsilon plus epsilon-squared");
  auto mean_cast = f;
  mean_cast[0].fill(1); mean_cast[1].fill(1 + std::ldexp(W(1), -52));
  h::EvaluationWork mean_cast_work; std::optional<h::Moments> cast_moments;
  need(h::reduce(mean_cast, e, constant, .01, mean_cast_work, 1024, cast_moments) == S::ok &&
       cast_moments->mean[0] == 1 && W(cast_moments->mean_cast_and_round[0]) >= std::ldexp(W(1), -53),
       "independent midpoint mean output cast loss two-to-minus53");
  auto covariance_cast = f;
  covariance_cast[0].fill(2); covariance_cast[1].fill(1 + std::ldexp(W(1), -30));
  h::EvaluationWork covariance_cast_work;
  need(h::reduce(covariance_cast, e, constant, .01, covariance_cast_work, 1024, cast_moments) == S::ok &&
       W(cast_moments->covariance_cast_and_round[0]) >= std::ldexp(W(1), -62),
       "independent squared-dyadic covariance output cast loss two-to-minus62");
  auto unresolved = f; unresolved[0].fill(2); unresolved[1].fill(2 + std::ldexp(W(1), -40));
  h::EvaluationWork unresolved_work;
  need(h::reduce(unresolved, noisy, constant, .01, unresolved_work, 1024, cast_moments) == S::conditioning_budget_exceeded &&
       !cast_moments, "nonconstant positive spread unresolved by inherited errors refuses");
  h::EvaluationWork exact_work; std::optional<h::Moments> exact;
  const auto measured = *work.checked_total();
  need(h::reduce(f, e, constant, .01, exact_work, measured, exact) == S::ok, "exact reduction work threshold");
  h::EvaluationWork below_work;
  need(h::reduce(f, e, constant, .01, below_work, measured - 1, exact) == S::work_limit && !exact, "one-below reduction work refusal");
  for (auto &state : f) state.fill(10);
  for (auto &state : e) state.fill(1e-6L);
  constant.fill(true); h::EvaluationWork zero_work;
  need(h::reduce(f, e, constant, .01, zero_work, 1024, m) == S::ok, "mathematical constant witness");
  for (double x : m->covariance) need(x == 0, "constant physical covariance exactly zero");
  for (double x : m->covariance_error) need(x == 0, "shared deterministic error creates no spread");
  need(m->mean_error[0] > 0, "constant mean still has forward error");
  constant.fill(false); h::EvaluationWork equal_work;
  need(h::reduce(f, e, constant, .01, equal_work, 1024, m) == S::conditioning_budget_exceeded && !m, "equal rounded varying-source values lack witness");
  for (auto &state : e) state.fill(0);
  f[0].fill(std::ldexp(W(1), -529)); f[1].fill(std::ldexp(W(1), -530));
  h::EvaluationWork subnormal_work;
  need(h::reduce(f, e, constant, .01, subnormal_work, 1024, m) == S::ok &&
       m->covariance[0] > 0 && std::fpclassify(m->covariance[0]) == FP_SUBNORMAL, "resolved positive subnormal covariance");
  f[0].fill(std::ldexp(W(1), -600)); f[1].fill(std::ldexp(W(1), -601));
  h::EvaluationWork underflow_work;
  need(h::reduce(f, e, constant, .01, underflow_work, 1024, m) == S::conditioning_budget_exceeded && !m, "positive covariance never projects to zero");
  need(h::outward(std::ldexp(W(1), -1200)).value_or(0) == std::numeric_limits<double>::denorm_min(), "known positive diagnostic projects upward");
  h::Bytes bytes; bytes.count(SIZE_MAX, 2); need(!bytes.result(), "checked product overflow");
  std::size_t maximum = SIZE_MAX; need(!h::add(maximum, 1) && maximum == SIZE_MAX, "checked work overflow leaves destination");
  h::PreparationWork bad_ledger; bad_ledger.overflowed = true;
  need(!bad_ledger.checked_total(), "counter overflow cannot produce an accepted total");
  f[0].fill(std::numeric_limits<W>::infinity()); h::EvaluationWork infinity_work;
  need(h::reduce(f, e, constant, .01, infinity_work, 1024, m) == S::conditioning_budget_exceeded && !m,
       "nonfinite required parent withholds all moments");
}
void source_controls() {
  auto in = h::original_input();
  need(h::structural(in), "complete original source");
  need(!h::constant_axis(in, 0) && !h::constant_axis(in, 2), "original supplied temperature varies");
  auto common = in; common.support[1].model.tcmb_kelvin = common.support[0].model.tcmb_kelvin;
  need(h::constant_axis(common, 0) && !h::constant_axis(common, 2), "initial matter-T has narrow raw photon witness only");
  auto duplicate = in; duplicate.support[1] = duplicate.support[0]; duplicate.support[1].id = "distinct duplicated-support ID";
  for (std::size_t i = 0; i < h::axes; ++i) need(h::constant_axis(duplicate, i), "full raw-input constant witness ignores labels");
  auto malformed = in; malformed.support[1].mass = 2;
  auto rejected = h::prepare(malformed);
  need(rejected.status() == S::invalid_input && !rejected.source(), "unequal masses are outside exact-half intake");
  malformed = in; malformed.requested_outputs = 31;
  need(h::prepare(malformed).status() == S::invalid_input, "fixed mask cannot change after failure");
  malformed = in; malformed.queries[1].id = malformed.queries[0].id;
  need(h::prepare(malformed).status() == S::invalid_input, "duplicate row ID refusal");
  auto role = in; role.temperature_role = "supplied-matter-temperature-vector";
  auto no_adapter = h::prepare(role);
  need(no_adapter.status() == S::outside_domain && no_adapter.source() &&
       no_adapter.preparation_work().background == 0, "matter-temperature role is not silently relabeled");
  auto p = h::Policy{}; p.history.absolute_temperature_tolerance_kelvin *= 2;
  need(h::prepare(in, p).status() == S::invalid_input, "inherited numerical allocation cannot relax");
  p = {}; p.maximum_preparation_work = h::preparation_work_cap + 1;
  need(h::prepare(in, p).status() == S::invalid_input, "scope cap cannot silently rise");
  const auto bound = h::Cohort::whole_live_bound(in, {});
  need(bound && *bound <= h::whole_live_payload_cap, "explicit full two-envelope layout bound");
  p = {}; p.maximum_whole_live_bytes = *bound - 1;
  auto below = h::prepare(in, p);
  need(below.status() == S::work_limit && !below.source() && below.preparation_work().checked_total() == 0, "one-below bytes refuse before acquisition");
  auto endpoint = in; endpoint.support[0].helium_fraction = 0; endpoint.support[1].helium_fraction = 0;
  auto density_only = h::prepare(endpoint);
  need(density_only.status() == S::outside_domain && density_only.source(), "both zero-He states retained without substitute");
  for (std::size_t s = 0; s < 2; ++s) {
    need(density_only.slot(s).snapshot.nuclei && density_only.slot(s).snapshot.nuclei->helium_nuclei.value == 0 &&
         !density_only.slot(s).budget.history_attempted, "endpoint map survives history refusal");
  }
  auto refused = density_only.evaluate();
  need(refused.status == S::outside_domain && !refused.moments && refused.source->probabilities[0] == .5 &&
       refused.source->probabilities[1] == .5 && refused.attempts[0].rows.empty() && refused.attempts[1].rows.empty(), "whole-law refusal with no fabricated native rows");
  p = {}; p.maximum_preparation_work = 2;
  auto limited = h::prepare(in, p);
  need(limited.status() == S::work_limit && limited.source() && limited.preparation_work().checked_total() == 2,
       "actual preparation quota and later retained refusal");
  need(limited.slot(1).status == S::work_limit && !limited.slot(1).budget.history_attempted, "later state is retained unattempted");
  const int rounding = std::fegetround(); std::fesetround(FE_DOWNWARD);
  auto wrong_rounding = h::prepare(in); std::fesetround(rounding);
  need(wrong_rounding.status() == S::invalid_input && !wrong_rounding.source(), "strict nearest profile gate");
}
void emitted_constructor_controls() {
  ref::Model model; model.H0 = 63; model.T0 = 2.7; model.massless_species.clear();
  const ref::EmittedThermalSource source{.0001, .05, .25, .00002};
  ref::Physics working(model, source);
  need(working.emitted_source.has_value() && working.fixed_photon_energy_today.has_value(), "explicit emitted-source role");
  auto changed = model; changed.baryon = .8; changed.cdm = .7; changed.other = .5;
  ref::Physics same_working(changed, source);
  need(working.hubble(2) == same_working.hubble(2) &&
       working.fixed_photon_energy_today == same_working.fixed_photon_energy_today,
       "physical density provenance never remaps the emitted law");
  auto colder = model; colder.T0 = 2.75;
  ref::Physics same_photons(colder, source);
  need(working.fixed_photon_energy_today == same_photons.fixed_photon_energy_today,
       "fixed mapped photons drive Compton rather than raw T0^4 remap");
  ref::Physics physical(model);
  need(!physical.emitted_source && !physical.fixed_photon_energy_today, "existing physical-reference route retained");
  ref::Result sorted;
  sorted.rows = {{2700, 1, 2, 3, 4, 5}, {300, 6, 7, 8, 9, 10}};
  sorted.stats.steps = 17;
  auto restored = ref::restore_original_rows(std::move(sorted), {300, 2700, 300});
  need(restored.rows.size() == 3 && restored.rows[0].z == 300 && restored.rows[1].z == 2700 &&
       restored.rows[2].z == 300 && restored.rows[0].opacity == 10 && restored.rows[2].opacity == 10 &&
       restored.stats.steps == 17, "reference restores original duplicate order and unchanged counters");
}
void wire_controls() {
  h::Wire w;
  const auto baseline = allocation_observation::live;
  allocation_observation::arm();
  w.scalar("", "bits", h::tcmb_a); w.wide("", "wide", W(1.5)); w.finish();
  const auto calls = allocation_observation::calls;
  allocation_observation::disarm();
  const std::string_view serialized{w.bytes.data(), w.used};
  need(!w.failed && calls == 0 && allocation_observation::live == baseline &&
       serialized.find("bits=400599999999999a\n") != std::string_view::npos &&
       serialized.find("wide=0:f:00000000000000001:c000000000000000\n") != std::string_view::npos &&
       w.fields == 7 && w.wide_nibbles == 16 && w.wide_decompositions == 1,
       "fixed zero-allocation exact double bits and canonical wide significand");
  h::Wire limited; limited.limit_work(5); limited.scalar("", "unadmitted", h::tcmb_a); limited.finish();
  need(limited.failed && limited.fields == 5 && limited.used < h::maximum_wire_bytes &&
       std::string_view(limited.bytes.data(), limited.used).find("unadmitted=") == std::string_view::npos,
       "serialization quota preserves a complete capped failure trailer without partial field");
  auto in = h::original_input();
  h::Wire original; h::wire_input(original, in, {}); original.finish();
  need(!original.failed && original.fields < h::maximum_wire_fields && original.work() <= h::evaluation_work_cap,
       "original offered support fits fixed wire field/work bounds");
}
void allocation_sweep(bool histories) {
  auto in = h::original_input();
  if (!histories) for (auto &s : in.support) s.helium_fraction = 0;
  const auto baseline = allocation_observation::live;
  allocation_observation::arm();
  std::size_t sites = 0, peak = 0, bound = 0;
  {
    auto owner = h::prepare(in);
    sites = allocation_observation::calls; peak = allocation_observation::peak - baseline;
    bound = owner.conservative_whole_live_bytes();
  }
  allocation_observation::disarm();
  need(allocation_observation::live == baseline && sites > 0 && peak <= bound, "measured preparation sites/layout clean baseline");
  for (std::size_t site = 0; site < sites; ++site) {
    bool valid = false, fired = false;
    allocation_observation::arm(long(site));
    try {
      {
        auto owner = h::prepare(in);
        valid = owner.status() != S::ok && (!owner.source() || owner.source()->probabilities == std::array<double, 2>{.5, .5});
      }
    } catch (const std::bad_alloc &) { valid = true; }
    fired = allocation_observation::fired;
    allocation_observation::disarm();
    need(valid && fired && allocation_observation::live == baseline, "each measured allocation failure retains law/refusal and RAII");
  }
  std::cout << "allocation_scope=" << (histories ? "complete_history_preparation" : "source_and_abundance_refusal")
            << " measured_sites=" << sites << " observed_peak_dynamic=" << peak << " conservative_bound=" << bound << '\n';
}
void evaluation_allocation_sweep() {
  auto in = h::original_input();
  const auto baseline = allocation_observation::live;
  std::size_t sites = 0, peak = 0, bound = 0;
  {
    auto owner = h::prepare(in);
    need(owner.status() == S::ok, "evaluation fault baseline original histories");
    allocation_observation::arm();
    {
      auto first = owner.evaluate(); auto second = owner.evaluate();
      sites = allocation_observation::calls / 2;
      peak = allocation_observation::peak - baseline;
      bound = owner.conservative_whole_live_bytes();
      need(first.status == S::ok && second.status == S::ok && first.moments && second.moments &&
           allocation_observation::calls == 2 * sites, "two-receipt allocation baseline");
    }
    allocation_observation::disarm();
  }
  need(allocation_observation::live == baseline && sites > 0 && peak <= bound,
       "both complete query receipts fit observed dynamic envelope");
  for (std::size_t site = 0; site < sites; ++site) {
    bool valid = false, fired = false;
    {
      auto owner = h::prepare(in);
      need(owner.status() == S::ok, "evaluation fault original preparation");
      const auto work = owner.preparation_work().checked_total();
      allocation_observation::arm(long(site));
      {
        auto refused = owner.evaluate();
        valid = refused.status == S::work_limit && !refused.moments && refused.source == owner.source() &&
            refused.source->probabilities == std::array<double, 2>{.5, .5} &&
            refused.work.dispatches == 2 && refused.work.queried_rows == 6 &&
            owner.preparation_work().checked_total() == work;
        fired = allocation_observation::fired;
      }
      allocation_observation::disarm();
    }
    need(valid && fired && allocation_observation::live == baseline,
         "each query allocation refusal preserves complete law/work and RAII");
  }
  std::cout << "allocation_scope=two_receipts_and_each_query_site measured_sites=" << sites
            << " observed_peak_dynamic=" << peak << " conservative_bound=" << bound << '\n';
}
void history_controls() {
  auto in = h::original_input(); auto owner = h::prepare(in);
  need(owner.status() == S::ok && owner.source(), "fixed original support prepared");
  const auto work = owner.preparation_work().checked_total();
  auto first = owner.evaluate(); auto second = owner.evaluate();
  need(first.status == S::ok && second.status == S::ok && first.moments && second.moments,
       "two pooled full-law evaluations");
  need(first.moments->mean == second.moments->mean && first.moments->covariance == second.moments->covariance &&
       owner.preparation_work().checked_total() == work, "retention gives identical results with no repeated trajectory work");
  need(first.work.dispatches == 2 && first.work.queried_rows == 6 && first.work.collections == 12 &&
       first.work.witnesses == 6 && first.work.means == 6 && first.work.differences == 6 &&
       first.work.products == 21 && first.work.projections == 27, "complete admitted logical evaluation graph");
  h::Wire wire; h::wire_input(wire, in, {}); h::wire_owner(wire, owner);
  wire.limit_work(owner.serialization_budget(second)); h::wire_evaluations(wire, first, second); wire.finish();
  need(!wire.failed && wire.fields == h::successful_wire_fields && wire.work() == h::successful_wire_work &&
       owner.record_serialization(second, wire.fields, wire.wide_nibbles, wire.wide_decompositions) == S::ok &&
       second.work.checked_total() == 964 &&
       owner.record_serialization(second, wire.fields, wire.wide_nibbles, wire.wide_decompositions) == S::work_limit,
       "complete frozen wire graph charges once to original final-evaluation cap");
  auto third = owner.evaluate();
  need(third.status == S::work_limit && !third.moments && third.attempts[0].rows.empty() &&
       third.state_status[0] == S::ok && third.state_status[1] == S::ok &&
       third.snapshots[0].nuclei && third.snapshots[1].thermal_witnesses &&
       third.snapshots[0].native_history_source_present,
       "third-call quota retains earned preparation metadata without query rows");
  const auto held = second.source;
  auto moved = std::move(owner);
  need(!owner.source() && owner.status() == S::invalid_input, "move clears original owner");
  moved = std::move(moved); need(moved.source() == held, "self move preserves source");
  moved = {}; need(second.source == held && second.moments.has_value(), "receipt survives native owner release");
  auto bad_row = in; bad_row.queries[0].redshift = 4000;
  auto row_owner = h::prepare(bad_row); auto rows = row_owner.evaluate();
  need(row_owner.status() == S::ok && rows.status == S::outside_domain && !rows.moments,
       "required invalid row withholds whole moments");
  for (const auto &attempt : rows.attempts) need(attempt.rows.size() == 3 &&
      attempt.rows[0].matter_temperature_kelvin.status == S::outside_domain &&
      attempt.rows[1].matter_temperature_kelvin.value && attempt.rows[2].thomson_opacity_per_redshift.value,
      "successful partner rows retained after required refusal");
  auto duplicate_rows = in;
  duplicate_rows.queries = {{{"late-first", 300}, {"initial-next", 2700}, {"late-duplicate", 300}}};
  auto ordered_owner = h::prepare(duplicate_rows); auto ordered = ordered_owner.evaluate();
  need(ordered.status == S::ok && ordered.moments, "unsorted duplicate coordinates with distinct IDs");
  for (const auto &attempt : ordered.attempts) need(attempt.rows[0].redshift == 300 &&
      attempt.rows[1].redshift == 2700 && attempt.rows[2].redshift == 300, "native exact original order/duplicates");
  auto causal = in;
  causal.support[1].hydrogen_mass_kg = 0;
  causal.queries[0].redshift = std::numeric_limits<double>::quiet_NaN();
  auto causal_owner = h::prepare(causal); auto causal_rows = causal_owner.evaluate();
  need(causal_owner.status() == S::outside_domain && causal_rows.status == S::nonfinite_input &&
       !causal_rows.moments && causal_rows.state_status[1] == S::outside_domain &&
       causal_rows.attempts[0].rows.size() == 3 && causal_rows.attempts[1].rows.empty(),
       "state-A query cause precedes state-B preparation cause without erasing either");
  auto rounding_owner = h::prepare(in);
  const int rounding = std::fegetround(); std::fesetround(FE_DOWNWARD);
  auto rejected_profile = rounding_owner.evaluate(); std::fesetround(rounding);
  need(rejected_profile.status == S::invalid_input && !rejected_profile.moments &&
       rejected_profile.attempts[0].rows.empty() && rejected_profile.state_status[0] == S::ok &&
       rejected_profile.state_status[1] == S::ok && rejected_profile.snapshots[0].nuclei &&
       rejected_profile.snapshots[1].emitted_thermal && rejected_profile.snapshots[0].thermal_witnesses,
       "unsupported evaluation profile preserves all earned preparation metadata");
}
}
int main(int argc, char **argv) {
  try {
#ifdef IRRED_COHORT_EXPECT_REJECTED_FP
    static_cast<void>(argc); static_cast<void>(argv);
    need(!h::profile(), "compiled unsafe FP profile is refused");
    auto in = h::original_input();
    auto owner = h::prepare(in);
    need(owner.status() == S::invalid_input && !owner.source() && owner.preparation_work().checked_total() == 0,
         "unsafe caller reaches no native physical work");
#else
    const std::string_view mode = argc == 1 ? "source" : argc == 2 ? std::string_view(argv[1]) : "invalid";
    if (mode == "source") { algebra(); source_controls(); emitted_constructor_controls(); wire_controls(); allocation_sweep(false); }
    else if (mode == "--history") history_controls();
    else if (mode == "--fault-history") { allocation_sweep(true); evaluation_allocation_sweep(); }
    else throw std::runtime_error("unrecognized proof-control mode");
#endif
    std::cout << "abundance_history_cohort checks=" << checks
              << " complete_reference=withheld physical_law=unqualified\n";
  } catch (const std::exception &e) {
    allocation_observation::disarm(); std::cerr << "FAIL after " << checks << ": " << e.what() << '\n'; return 1;
  }
}
