// Native owner/resource controls; separate peer has independent time/frequency
// references.
#include "irred/temporal_photometry.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <numbers>
#include <stdexcept>
#include <utility>
static long fail_at = -1, allocations = 0, live = 0;
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(std::size_t n) {
  if (fail_at >= 0 && allocations++ == fail_at)
    throw std::bad_alloc();
  if (void *p = std::malloc(n ? n : 1)) {
    ++live;
    return p;
  }
  throw std::bad_alloc();
}
NOINLINE void *operator new[](std::size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
NOINLINE void operator delete(void *p, std::size_t) noexcept {
  ::operator delete(p);
}
NOINLINE void operator delete[](void *p, std::size_t) noexcept {
  ::operator delete(p);
}
namespace {
namespace p = irred::photometry;
using S = irred::numerics::Status;
using W = long double;
unsigned checks = 0;
void check(bool b, const char *what) {
  ++checks;
  if (!b)
    throw std::runtime_error(what);
}
void near(const p::Outcome &v, W expected, const char *what) {
  check(v.availability == p::Availability::available && v.value.has_value(),
        what);
  if (expected > 0)
    check(*v.value > 0, what);
  else
    check(*v.value == 0, what);
  check(std::abs(W(*v.value) - expected) <=
            2e-12L * std::abs(expected) + 1e-300L,
        what);
}
p::TemporalExposure exposure() { return {0, 0, 2, 1, 3, 10, 10, 30}; }
void same(const p::Outcome &a, const p::Outcome &b) {
  check(a.availability == b.availability &&
            a.numerical_status == b.numerical_status,
        "shared scalar state");
  if (a.value && b.value)
    check(std::abs(W(*a.value) - *b.value) <=
              2e-12L * std::abs(W(*b.value)) + 1e-300L,
          "shared scalar value");
}
} // namespace
int main() {
  try {
    constexpr double m = 0x1p-20;
    std::array times{0., 10.};
    std::array wave{m, 2 * m};
    std::array lum{4., 4., 4., 4.};
    std::array bw{2 * m, 4 * m};
    std::array trans{.5, .5};
    std::array grids{p::TemporalGrid{
        "grid-a", "declared synthetic constant grid", times, wave, lum}};
    std::array bands{p::TemporalBand{
        "band-a", "fixed synthetic optical calibration", {bw, trans}}};
    auto prepared = p::prepare_temporal(grids, bands);
    check(prepared.status() == S::ok, "prepare constant grid");
    check(prepared.grid_count() == 1 && prepared.band_count() == 1,
          "pool counts");
    auto gv = *prepared.grid(0);
    auto bv = *prepared.band(0);
    check(gv.id == grids[0].id && bv.id == bands[0].id,
          "source identities retained");
    for (auto pair :
         {std::pair{gv.rest_time_second, std::span<const double>(times)},
          std::pair{gv.rest_wavelength_metre, std::span<const double>(wave)},
          std::pair{gv.luminosity_watt_per_metre, std::span<const double>(lum)},
          std::pair{bv.passband.wavelength_metre, std::span<const double>(bw)},
          std::pair{bv.passband.optical_transmission,
                    std::span<const double>(trans)}}) {
      check(pair.first.data() != pair.second.data(),
            "immutable separate source owner");
      for (std::size_t j = 0; j < pair.first.size(); ++j)
        check(std::bit_cast<std::uint64_t>(pair.first[j]) ==
                  std::bit_cast<std::uint64_t>(pair.second[j]),
              "owned input exact bits/order");
    }
    std::array rows{exposure()};
    auto out = p::evaluate_temporal(prepared, rows);
    check(out.status == S::ok && out.rows.size() == 1 &&
              out.rows[0].admission_status == S::ok,
          "constant exposure admission");
    const auto &r = out.rows[0];
    const W flux = W(m) / (4 * std::numbers::pi_v<W>);
    near(r.mean_flux_watt_per_square_metre, flux, "constant incident mean");
    near(r.energy_joule, 30 * flux, "constant energy");
    near(r.expected_photons,
         30 * flux * 3 * m / (6.62607015e-34L * 299792458.L),
         "constant photons");
    check(r.coverage == p::TemporalCoverage::full && r.covered_fraction == 1 &&
              r.observer_duration_second == 20,
          "full support observer convention");
    const std::array steady{4., 4.};
    auto native = p::evaluate_sampled(
        {{wave, steady}, {bw, trans}, 2, 1, 3, 20}, {7, 64, 64});
    same(r.mean_flux_watt_per_square_metre, native.flux_watt_per_square_metre);
    same(r.energy_joule, native.energy_joule);
    same(r.expected_photons, native.expected_photons);
    lum.fill(99);
    times[0] = -99;
    trans.fill(0);
    check(prepared.grid(0)->luminosity_watt_per_metre[0] == 4 &&
              prepared.grid(0)->rest_time_second[0] == 0 &&
              prepared.band(0)->passband.optical_transmission[0] == .5,
          "caller mutation does not alter source");
    auto copy = prepared;
    prepared = p::PreparedTemporal{};
    auto moved = std::move(copy);
    check(copy.status() == S::invalid_input && !copy.grid(0) &&
              moved.status() == S::ok,
          "move leaves failed source");
    auto *self = &moved;
    moved = std::move(*self);
    check(moved.status() == S::ok, "self move retains source");
    check(p::evaluate_temporal(copy, rows).status == S::invalid_input,
          "failed moved owner cannot evaluate");
    prepared = moved;
    rows[0].observer_lower_second = 0;
    out = p::evaluate_temporal(prepared, rows);
    const auto &partial = out.rows[0];
    check(partial.coverage == p::TemporalCoverage::partial &&
              std::abs(partial.covered_fraction - 2.L / 3) < 1e-18L,
          "incomplete support explicit");
    near(partial.mean_flux_watt_per_square_metre, 2.L / 3 * flux,
         "full exposure mean includes zero outside");
    near(partial.energy_joule, 30 * flux, "empty observer time adds no energy");
    rows[0].observer_lower_second = 30;
    rows[0].observer_upper_second = 40;
    out = p::evaluate_temporal(prepared, rows);
    check(out.rows[0].coverage == p::TemporalCoverage::no_overlap,
          "exact temporal touching zero");
    near(out.rows[0].mean_flux_watt_per_square_metre, 0, "touching mean zero");
    near(out.rows[0].energy_joule, 0, "touching energy zero");
    near(out.rows[0].expected_photons, 0, "touching photons zero");
    rows[0] = exposure();
    rows[0].source_epoch_observer_second = -0.;
    rows[0].observer_lower_second = 0;
    rows[0].observer_upper_second = 20;
    out = p::evaluate_temporal(prepared, rows);
    check(std::bit_cast<std::uint64_t>(
              out.rows[0].source.source_epoch_observer_second) ==
              std::bit_cast<std::uint64_t>(-0.),
          "epoch exact signed-zero source bits");
    // Independent linear-time mean fact; shared wavelength reference checks
    // composition only.
    times = {0, 10};
    lum = {2, 2, 6, 6};
    trans = {.5, .5};
    prepared = p::prepare_temporal(grids, bands);
    rows[0] = exposure();
    rows[0].observer_lower_second = 20;
    out = p::evaluate_temporal(prepared, rows);
    near(out.rows[0].mean_flux_watt_per_square_metre, 5.L / 4 * flux,
         "linear time clipped mean");
    near(out.rows[0].energy_joule, 15 * 5.L / 4 * flux,
         "linear time exposure integral");
    std::array badrows{exposure(), exposure(), exposure(), exposure()};
    badrows[1].grid_index = 9;
    badrows[2].redshift = -1;
    badrows[3].observer_upper_second = badrows[3].observer_lower_second;
    out = p::evaluate_temporal(prepared, badrows);
    check(out.status == S::ok && out.rows.size() == 4,
          "mixed coarse batch keeps order");
    for (std::size_t i = 0; i < 4; ++i)
      check(out.rows[i].source.grid_index == badrows[i].grid_index &&
                out.rows[i].source.redshift == badrows[i].redshift,
            "row source order");
    check(out.rows[0].admission_status == S::ok &&
              out.rows[1].admission_status == S::invalid_input &&
              out.rows[2].admission_status == S::outside_domain &&
              out.rows[3].admission_status == S::outside_domain,
          "typed mixed failures");
    p::TemporalPolicy policy;
    policy.requested_outputs = p::collected_energy;
    out = p::evaluate_temporal(prepared, badrows, policy);
    check(out.rows[0].mean_flux_watt_per_square_metre.availability ==
                  p::Availability::omitted &&
              out.rows[0].expected_photons.availability ==
                  p::Availability::omitted,
          "omitted outputs");
    check(out.rows[1].energy_joule.availability == p::Availability::failed &&
              out.rows[1].mean_flux_watt_per_square_metre.availability ==
                  p::Availability::omitted,
          "failed requested group only");
    policy = {};
    policy.maximum_segment_work = 3;
    std::array two{exposure(), exposure()};
    out = p::evaluate_temporal(prepared, two, policy);
    check(out.status == S::ok && out.segment_work == 3 &&
              out.rows[0].admission_status == S::ok &&
              out.rows[1].admission_status == S::work_limit,
          "global work budget and prior row preserved");
    policy.maximum_native_bytes = 0;
    check(p::evaluate_temporal(prepared, two, policy).status == S::work_limit,
          "evaluation byte admission");
    policy = {};
    policy.maximum_rows = 0;
    check(p::evaluate_temporal(prepared, two, policy).status == S::work_limit,
          "row quota admission");
    policy = {};
    policy.maximum_segment_work = 4000001;
    check(p::evaluate_temporal(prepared, two, policy).status ==
              S::invalid_input,
          "hard global segment-work policy refuses");
    policy = {};
    policy.maximum_rows = 65537;
    check(p::evaluate_temporal(prepared, two, policy).status ==
              S::invalid_input,
          "hard row policy refuses");
    policy = {};
    policy.requested_outputs = 0;
    check(p::evaluate_temporal(prepared, two, policy).status ==
              S::invalid_input,
          "empty output mask refuses");
    {
      std::vector<double> too_many_times(65537, NAN);
      auto hostile = grids[0];
      hostile.rest_time_second = too_many_times;
      const std::array hostile_grids{hostile};
      allocations = 0;
      fail_at = 0;
      auto no_scan = p::prepare_temporal(hostile_grids, bands);
      fail_at = -1;
      check(no_scan.status() == S::work_limit && allocations == 0,
            "hard axis cap before nonfinite scans/allocation");
    }
    {
      std::vector<double> large_times(4096), large_luminosity(8192, 4);
      for (std::size_t i = 0; i < large_times.size(); ++i)
        large_times[i] = i;
      const p::TemporalGrid unused{"unused-grid",
                                   "synthetic unrequested source", large_times,
                                   wave, large_luminosity};
      const std::array pooled{grids[0], unused};
      auto pooled_owner = p::prepare_temporal(pooled, bands);
      policy = {};
      policy.maximum_native_bytes = 1024;
      const auto subset = p::evaluate_temporal(pooled_owner, rows, policy);
      check(subset.status == S::ok && subset.rows[0].admission_status == S::ok,
            "unrequested grid imposes no batch scratch or byte refusal");
    }
    const auto payload = *prepared.retained_payload_bytes();
    p::TemporalPreparationPolicy prep;
    prep.maximum_native_bytes = payload;
    check(p::prepare_temporal(grids, bands, prep).status() == S::ok,
          "exact retained payload admitted");
    prep.maximum_native_bytes = payload - 1;
    auto refused = p::prepare_temporal(grids, bands, prep);
    check(refused.status() == S::work_limit &&
              !refused.retained_payload_bytes() && !refused.grid(0),
          "one under payload no owner");
    prep = {};
    prep.maximum_grid_values = 3;
    check(p::prepare_temporal(grids, bands, prep).status() == S::work_limit,
          "grid-cell quota");
    prep = {};
    prep.maximum_total_knots = 5;
    check(p::prepare_temporal(grids, bands, prep).status() == S::work_limit,
          "pooled axes quota");
    times[1] = times[0];
    check(p::prepare_temporal(grids, bands).status() == S::outside_domain,
          "duplicate time refused");
    times[1] = 10;
    lum[0] = -1;
    check(p::prepare_temporal(grids, bands).status() == S::outside_domain,
          "negative luminosity refused");
    lum[0] = 2;
    lum[0] = NAN;
    check(p::prepare_temporal(grids, bands).status() == S::nonfinite_input,
          "nonfinite luminosity refused");
    lum[0] = 2;
    trans[0] = 2;
    check(p::prepare_temporal(grids, bands).status() == S::outside_domain,
          "invalid optical calibration refused");
    trans[0] = .5;
    auto original = grids[0];
    grids[0].luminosity_watt_per_metre = std::span<const double>(lum).first(3);
    check(p::prepare_temporal(grids, bands).status() == S::invalid_input,
          "rectangular shape refused");
    grids[0] = original;
    std::array duplicates{grids[0], grids[0]};
    check(p::prepare_temporal(duplicates, bands).status() == S::invalid_input,
          "duplicate source IDs refused");
    // Binary-exact adjacent time knots must retain strictly positive exposure
    // energy.
    times = {1, std::nextafter(1., 2.)};
    bw = {m, 2 * m};
    lum.fill(1);
    prepared = p::prepare_temporal(grids, bands);
    rows[0] = {0, 0, 2, 0, 3, 0, times[0], times[1]};
    out = p::evaluate_temporal(prepared, rows);
    check(out.rows[0].admission_status == S::ok,
          "adjacent time support admitted");
    check(out.rows[0].energy_joule.value && *out.rows[0].energy_joule.value > 0,
          "adjacent positive time energy");
    // True positive temporal overlap cannot become a fabricated
    // no-overlap/zero.
    times = {.5, 1};
    prepared = p::prepare_temporal(grids, bands);
    rows[0] = {0, 0, 2, std::ldexp(1., -65), 3, 0, 1, 2};
    out = p::evaluate_temporal(prepared, rows);
    check(out.rows[0].admission_status == S::conditioning_budget_exceeded &&
              out.rows[0].coverage == p::TemporalCoverage::unassessed &&
              !out.rows[0].energy_joule.value,
          "unresolved positive temporal edge refused");
    times = {0, 10};
    bw = {2 * m, 4 * m};
    prepared = p::prepare_temporal(grids, bands);
    rows[0] = exposure();
    rows[0].source_epoch_observer_second = 1e308;
    rows[0].observer_lower_second = 1e308;
    rows[0].observer_upper_second = std::nextafter(1e308, INFINITY);
    out = p::evaluate_temporal(prepared, rows);
    check(out.rows[0].admission_status == S::conditioning_budget_exceeded &&
              !out.rows[0].energy_joule.value,
          "large epoch lost positive interval refused");
    rows[0] = exposure();
    rows[0].observer_lower_second = -1e308;
    rows[0].observer_upper_second = 1e308;
    out = p::evaluate_temporal(prepared, rows);
    check(out.rows[0].admission_status == S::overflow,
          "unrepresentable observer duration refused");
    rows[0] = exposure();
    rows[0].luminosity_distance_metre = 1e308;
    rows[0].collecting_area_square_metre = 0;
    out = p::evaluate_temporal(prepared, rows);
    check(out.rows[0].mean_flux_watt_per_square_metre.availability ==
                  p::Availability::failed &&
              out.rows[0].energy_joule.availability ==
                  p::Availability::available &&
              out.rows[0].expected_photons.availability ==
                  p::Availability::available,
          "shared spectral partial failures independent");
    // Scoped failure sweeps prove native ownership cleanup; no C ABI status
    // exists here.
    times = {0, 10};
    lum.fill(4);
    rows[0] = exposure();
    long prep_sites = 0, eval_sites = 0;
    {
      const auto baseline = live;
      allocations = 0;
      fail_at = 1000000;
      {
        auto p0 = p::prepare_temporal(grids, bands);
        prep_sites = allocations;
        fail_at = -1;
        check(p0.status() == S::ok, "allocation prepare baseline");
      }
      check(live == baseline, "prepare baseline cleanup");
      for (long i = 0; i < prep_sites; ++i) {
        allocations = 0;
        fail_at = i;
        bool threw = false;
        try {
          auto v = p::prepare_temporal(grids, bands);
          (void)v;
        } catch (const std::bad_alloc &) {
          threw = true;
        }
        fail_at = -1;
        check(threw && live == baseline,
              "every preparation allocation cleanup");
      }
    }
    prepared = p::prepare_temporal(grids, bands);
    {
      const auto baseline = live;
      allocations = 0;
      fail_at = 1000000;
      {
        auto v = p::evaluate_temporal(prepared, rows);
        eval_sites = allocations;
        fail_at = -1;
        check(v.status == S::ok, "evaluation allocation baseline");
      }
      check(live == baseline, "evaluation baseline cleanup");
      for (long i = 0; i < eval_sites; ++i) {
        allocations = 0;
        fail_at = i;
        bool threw = false;
        try {
          auto v = p::evaluate_temporal(prepared, rows);
          (void)v;
        } catch (const std::bad_alloc &) {
          threw = true;
        }
        fail_at = -1;
        check(threw && live == baseline, "every evaluation allocation cleanup");
      }
    }
    std::cout
        << "PASS " << checks
        << " temporal photometry owner controls preparation_allocation_sites="
        << prep_sites << " evaluation_allocation_sites=" << eval_sites
        << " prepared_payload=" << payload << '\n';
  } catch (const std::exception &e) {
    fail_at = -1;
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
