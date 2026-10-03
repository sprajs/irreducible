#include "hydrogen_helium_supplied_cases.hpp"
#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
// Single-thread allocation campaign. Counters cover requested allocation bytes
// including original/replacement owners simultaneously. Allocation headers,
// allocator metadata and RSS are deliberately outside this byte witness.
namespace {
std::atomic<bool> armed = false;
std::atomic<std::size_t> calls = 0, fail_on = 0, live_bytes = 0, peak_bytes = 0;
std::atomic<long> live_count = 0;
struct alignas(std::max_align_t) Header { std::size_t bytes; };
void start(std::size_t failure = 0) {
  calls = 0; fail_on = failure; peak_bytes = live_bytes.load(); armed = true;
}
} // namespace
#if defined(__GNUC__) || defined(__clang__)
#define HHE_NOINLINE __attribute__((noinline))
#else
#define HHE_NOINLINE
#endif
HHE_NOINLINE void *operator new(std::size_t bytes) {
  if (armed && ++calls == fail_on) throw std::bad_alloc();
  if (bytes > std::numeric_limits<std::size_t>::max() - sizeof(Header)) throw std::bad_alloc();
  auto *header = static_cast<Header *>(std::malloc(sizeof(Header) + (bytes ? bytes : 1)));
  if (!header) throw std::bad_alloc();
  header->bytes = bytes;
  ++live_count;
  const auto current = live_bytes.fetch_add(bytes) + bytes;
  if (armed && current > peak_bytes) peak_bytes = current;
  return header + 1;
}
HHE_NOINLINE void *operator new[](std::size_t bytes) { return ::operator new(bytes); }
HHE_NOINLINE void operator delete(void *pointer) noexcept {
  if (!pointer) return;
  auto *header = static_cast<Header *>(pointer) - 1;
  live_bytes -= header->bytes; --live_count; std::free(header);
}
HHE_NOINLINE void operator delete[](void *pointer) noexcept { ::operator delete(pointer); }
HHE_NOINLINE void operator delete(void *pointer, std::size_t) noexcept { ::operator delete(pointer); }
HHE_NOINLINE void operator delete[](void *pointer, std::size_t) noexcept { ::operator delete(pointer); }
int main() {
  using namespace irred::cosmology;
  using S = irred::numerics::Status;
  auto source = hydrogen_helium_supplied_cases::source(2);
  source.initial.origin.assign(128, 'o'); source.initial.source_identity.assign(128, 's');
  source.history.nuclei_origin.assign(128, 'n');
  HydrogenHeliumHistoryPolicy policy; policy.base_intervals = 16384;
  const auto before_preparation = live_bytes.load();
  start();
  auto owner = prepare_hydrogen_helium_supplied_history(source, policy);
  armed = false;
  const auto preparation_calls = calls.load(), preparation_peak = peak_bytes.load() - before_preparation;
  if (owner.status() != S::ok || preparation_calls == 0 || preparation_calls > 256 ||
      preparation_peak > policy.maximum_native_bytes) return 1;
  const double queries[]{2800, 1300, 300};
  for (std::size_t site = 1; site <= preparation_calls; ++site) {
    std::cout << "preparation_fault_site=" << site << std::endl;
    const auto old_bytes = live_bytes.load(); const auto old_count = live_count.load();
    {
      start(site);
      auto rejected = prepare_hydrogen_helium_supplied_history(source, policy);
      armed = false;
      std::cout << "preparation_fault_result=" << site << " status=" << int(rejected.status())
                << " work=" << rejected.work().total() << " import=" << bool(rejected.initial_import_witness())
                << " common_source=" << bool(rejected.source()) << " rate=" << bool(rejected.rate_domain_witness()) << '\n';
      if (rejected.status() != S::work_limit || rejected.work().total() > policy.maximum_total_work ||
          !rejected.evaluate(queries, 31).rows.empty()) return 2;
      if (rejected.initial_import_witness() && (!rejected.supplied_initial_state() ||
          rejected.supplied_initial_state()->origin != source.initial.origin ||
          rejected.supplied_initial_state()->source_identity != source.initial.source_identity)) return 3;
      if (rejected.rate_domain_witness() && rejected.rate_domain_witness()->retained_kelvin_range) return 3;
    }
    if (live_bytes != old_bytes || live_count != old_count) return 4;
  }
  const auto old_query_bytes = live_bytes.load();
  start(); auto observed = owner.evaluate(queries, 31); armed = false;
  const auto query_calls = calls.load(), query_peak = peak_bytes.load() - old_query_bytes;
  if (observed.rows.size() != 3 || query_calls == 0 || query_calls > 256) return 5;
  observed = {};
  for (std::size_t site = 1; site <= query_calls; ++site) {
    std::cout << "query_fault_site=" << site << std::endl;
    const auto old_bytes = live_bytes.load(); const auto old_count = live_count.load();
    {
      start(site); auto refused = owner.evaluate(queries, 31); armed = false;
      if (refused.status != S::work_limit || !refused.rows.empty()) return 5;
    }
    if (live_bytes != old_bytes || live_count != old_count) return 6;
  }
  std::size_t copy_calls = 0, copy_peak = 0, copy_total_peak = 0;
  {
    const auto before = live_bytes.load();
    start(); auto copy = owner; armed = false;
    copy_calls = calls; copy_peak = peak_bytes.load() - before;
    copy_total_peak = peak_bytes;
    if (!copy.supplied_initial_state() || copy.evaluate(queries, 31).rows.size() != 3 ||
        copy_calls == 0 || copy_calls > 256) return 7;
  }
  for (std::size_t site = 1; site <= copy_calls; ++site) {
    std::cout << "copy_fault_site=" << site << std::endl;
    const auto old_bytes = live_bytes.load(); const auto old_count = live_count.load();
    bool threw = false;
    start(site);
    try { auto copy = owner; armed = false; if (copy.evaluate(queries, 31).rows.size() != 3) return 7; }
    catch (const std::bad_alloc &) { threw = true; }
    armed = false;
    if (!threw || live_bytes != old_bytes || live_count != old_count) return 8;
  }
  auto destination_source = hydrogen_helium_supplied_cases::source(0);
  auto destination = prepare_hydrogen_helium_supplied_history(destination_source, policy);
  if (destination.status() != S::ok) return 9;
  const double destination_queries[]{2700, 1300, 300};
  const auto original = destination.evaluate(destination_queries, 31);
  const auto old_work = destination.work();
  const auto old_import = *destination.initial_import_witness();
  const auto old_range = *destination.rate_domain_witness();
  auto same = [](const HydrogenHeliumHistoryValue &a, const HydrogenHeliumHistoryValue &b) {
    return a.status == b.status && a.value == b.value && a.absolute_error_estimate == b.absolute_error_estimate;
  };
  for (std::size_t site = 1; site <= copy_calls; ++site) {
    std::cout << "replacement_fault_site=" << site << std::endl;
    const auto old_bytes = live_bytes.load(); const auto old_count = live_count.load();
    bool threw = false;
    start(site);
    try { destination = owner; } catch (const std::bad_alloc &) { threw = true; }
    armed = false;
    if (!threw || live_bytes != old_bytes || live_count != old_count || !destination.source() ||
        !destination.supplied_initial_state() || destination.source()->initial_redshift != 2700 ||
        destination.supplied_initial_state()->source_identity != destination_source.initial.source_identity ||
        destination.model_identity() != hydrogen_helium_supplied_history_model_id ||
        destination.work().total() != old_work.total() ||
        destination.initial_import_witness()->promoted_values != old_import.promoted_values ||
        destination.rate_domain_witness()->attempted_kelvin_range != old_range.attempted_kelvin_range ||
        destination.rate_domain_witness()->retained_kelvin_range != old_range.retained_kelvin_range) return 10;
    const auto current = destination.evaluate(destination_queries, 31);
    if (current.status != original.status || current.rows.size() != original.rows.size()) return 10;
    for (std::size_t j = 0; j < current.rows.size(); ++j) {
      const auto &a = current.rows[j], &b = original.rows[j];
      if (a.redshift != b.redshift || !same(a.hydrogen_ionized_fraction, b.hydrogen_ionized_fraction) ||
          !same(a.helium_singly_ionized_fraction, b.helium_singly_ionized_fraction) ||
          !same(a.electron_number_density_per_cubic_metre, b.electron_number_density_per_cubic_metre) ||
          !same(a.matter_temperature_kelvin, b.matter_temperature_kelvin) ||
          !same(a.thomson_opacity_per_redshift, b.thomson_opacity_per_redshift)) return 10;
    }
  }
  const auto before_assignment = live_bytes.load();
  start(); destination = owner; armed = false;
  const auto assignment_calls = calls.load(), assignment_peak = peak_bytes.load() - before_assignment;
  const auto assignment_total_peak = peak_bytes.load();
  if (assignment_calls != copy_calls || destination.supplied_initial_state()->source_identity != source.initial.source_identity ||
      destination.evaluate(queries, 31).rows.size() != 3) return 11;
  std::cout << "supplied_allocation_sites preparation=" << preparation_calls << " query=" << query_calls
            << " copy=" << copy_calls << " replacement=" << assignment_calls
            << " requested_delta_peak_bytes=" << preparation_peak << ',' << query_peak << ',' << copy_peak << ',' << assignment_peak
            << " retained_all_owner_requested_bytes=" << live_bytes.load()
            << " simultaneous_copy_replacement_total_peaks=" << copy_total_peak << ',' << assignment_total_peak
            << " prep_cap_only=" << policy.maximum_native_bytes << " allocator_metadata_RSS_included=false\n";
  // The copy/replacement peaks are measured, not declared a general cap on
  // arbitrary caller-owned multiplicity or on total process memory.
  return 0;
}
