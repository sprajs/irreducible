#include "conditional_drag_reference_ownership.hpp"
#include <cstdio>
#include <type_traits>
#include <utility>

namespace ref = conditional_drag_reference_certificate;
namespace {
unsigned checks = 0, failures = 0;
void need(bool answer, const char *label) {
  ++checks;
  if (!answer) { ++failures; std::fprintf(stderr, "FAILED %s\n", label); }
}
// These controls exercise the conditional policy ledger WITHOUT allocating.
// Declaring this profile does not earn actual allocator/frame/library coverage.
constexpr ref::Profile conditional{true, true, true, true, false, false};
struct ProposedNodeLayout { ref::W entries[13]; };
void early_refusal() {
  ref::Owner unsupported;
  const auto initial = unsupported.live_bytes();
  ref::PayloadReservation frame;
  need(!ref::PayloadAccess::frame(unsupported, frame, ref::Payload::stage, 4096) &&
       !frame && unsupported.failure().cause == ref::Cause::unsupported_profile &&
       unsupported.live_bytes() == initial && unsupported.pending_bytes() == 0,
       "unknown frame profile refuses before payload acquisition");
  ref::Owner no_allocator(ref::Profile{true, true, true, false, false, false});
  ref::PayloadReservation block;
  need(!ref::PayloadAccess::allocation(no_allocator, block,
       {1, 4096, 16, 64, ref::Payload::source, ref::Storage::ancillary}) &&
       !block && no_allocator.failure().cause == ref::Cause::unsupported_profile &&
       no_allocator.pending_bytes() == 0,
       "unknown allocator metadata profile refuses before pending reservation");
}
void old_pending_peak_and_unwind() {
  ref::Owner owner(conditional);
  const auto initial = owner.live_bytes();
  constexpr ref::U n = ref::maximum_node_capacity;
  const ref::U each = n * sizeof(ProposedNodeLayout) + 64 + 15 + ref::metadata_per_block;
  ref::Counts first_prefix{};
  {
    ref::PayloadReservation old, pending;
    const ref::AllocationRequest request{n, sizeof(ProposedNodeLayout), 16, 64,
                                        ref::Payload::stage, ref::Storage::reference_nodes};
    need(ref::PayloadAccess::allocation(owner, old, request) &&
         owner.live_bytes() == initial && owner.pending_bytes() == each,
         "requested old node block is pending before external allocator entry");
    need(old.commit() && owner.live_bytes() == initial + each &&
         owner.pending_bytes() == 0, "successful external result commits pending to live");
    need(ref::PayloadAccess::allocation(owner, pending, request) &&
         owner.live_bytes() == initial + each && owner.pending_bytes() == each &&
         owner.peak_bytes() == initial + 2 * each,
         "old and actual requested pending block contribute simultaneously once");
    owner.fail(ref::Cause::allocation_failure); // synthetic external outcome
    first_prefix = owner.failure().served_prefix;
    need(pending.reset() && owner.pending_bytes() == 0 &&
         owner.peak_bytes() == initial + 2 * each &&
         owner.failure().served_prefix == first_prefix,
         "failed pending result unwinds without reducing prior peak/causal prefix");
  }
  need(owner.live_bytes() == initial && owner.pending_bytes() == 0 &&
       owner.cleanup_reserved() == ref::Counts{} &&
       owner.failure().served_prefix == first_prefix &&
       owner.served()[unsigned(ref::Counter::transitions)] ==
       first_prefix[unsigned(ref::Counter::transitions)] + 2,
       "both real cleanup exits consume pre-reserved credits after first failure");
}
void ownership_move_boundaries() {
  ref::Owner owner(conditional);
  const auto initial = owner.live_bytes();
  {
    ref::PayloadReservation a;
    need(ref::PayloadAccess::frame(owner, a, ref::Payload::stage, 4096),
         "first fixed owned frame acquired");
    ref::PayloadReservation b(std::move(a));
    need(!a && bool(b) && owner.live_bytes() == initial + 4096,
         "move construction transfers reservation once and clears the source");
    b = std::move(b);
    need(bool(b) && owner.live_bytes() == initial + 4096,
         "self move leaves one owned reservation");
    ref::PayloadReservation c;
    need(ref::PayloadAccess::frame(owner, c, ref::Payload::stage, 2048),
         "assignment source acquired while destination stays live");
    b = std::move(c);
    need(!c && bool(b) && owner.live_bytes() == initial + 2048 &&
         owner.peak_bytes() == initial + 6144,
         "move assignment releases destination after preserving coexistence peak");
  }
  need(owner.live_bytes() == initial && owner.cleanup_reserved() == ref::Counts{},
       "source-cleared moves leave no duplicate resource ownership");

  ref::Owner first_owner(conditional), second_owner(conditional);
  const auto first_initial = first_owner.live_bytes(), second_initial = second_owner.live_bytes();
  ref::PayloadReservation first, second;
  need(ref::PayloadAccess::frame(first_owner, first, ref::Payload::stage, 4096) &&
       ref::PayloadAccess::frame(second_owner, second, ref::Payload::stage, 2048),
       "two independent request ownership controls acquired");
  second = std::move(first);
  need(bool(first) && bool(second) && first_owner.failure().cause == ref::Cause::invalid_ownership &&
       first_owner.live_bytes() == first_initial + 4096 &&
       second_owner.live_bytes() == second_initial + 2048,
       "cross-request assignment refuses before damaging either owned frame");
  ref::PayloadReservation absent(std::move(first));
  need(!absent && bool(first),
       "post-refusal move cannot acquire new ownership or silently empty source");
  need(first.reset() && second.reset() && first_owner.live_bytes() == first_initial &&
       second_owner.live_bytes() == second_initial,
       "both acquired lifetimes unwind after refused transfer");
}
void request_and_phase_refusals() {
  ref::Owner excessive(conditional);
  ref::PayloadReservation rejected;
  const auto initial = excessive.live_bytes();
  need(!ref::PayloadAccess::allocation(excessive, rejected,
       {ref::maximum_node_capacity + 1, sizeof(ProposedNodeLayout), 16, 64,
        ref::Payload::stage, ref::Storage::reference_nodes}) &&
       excessive.failure().cause == ref::Cause::payload_limit && !rejected &&
       excessive.failure().requested_count == ref::maximum_node_capacity + 1 &&
       excessive.live_bytes() == initial && excessive.pending_bytes() == 0,
       "allocator over-request refuses BEFORE allocation and retains exact count");
  ref::Owner arithmetic_overflow(conditional);
  ref::PayloadReservation overflow;
  need(!ref::PayloadAccess::allocation(arithmetic_overflow, overflow,
       {std::numeric_limits<ref::U>::max(), 2, 16, 64,
        ref::Payload::source, ref::Storage::ancillary}) &&
       arithmetic_overflow.failure().cause == ref::Cause::counter_overflow &&
       arithmetic_overflow.failure().requested_element_bytes == 2 && !overflow,
       "unrepresentable actual count-times-element request cannot wrap into a fit");
  ref::Owner category(conditional);
  ref::PayloadReservation double_category;
  need(!ref::PayloadAccess::allocation(category, double_category,
       {1, ref::payload_caps[unsigned(ref::Payload::allocator_metadata)], 16, 64,
        ref::Payload::allocator_metadata, ref::Storage::ancillary}) &&
       category.failure().cause == ref::Cause::payload_limit && !double_category,
       "data plus metadata in one category share one checked remaining cap");

  ref::Owner native(conditional);
  ref::PayloadReservation child, reference;
  need(ref::PayloadAccess::allocation(native, child,
       {1, ref::native_child_byte_cap, 16, 64,
        ref::Payload::native_batches, ref::Storage::native_internal}) && child.commit(),
       "active native raw storage keeps its distinct original32MiB child fence");
  need(!ref::PayloadAccess::allocation(native, reference,
       {1, sizeof(ProposedNodeLayout), 16, 64,
        ref::Payload::stage, ref::Storage::reference_nodes}) && !reference &&
       native.failure().cause == ref::Cause::payload_limit,
       "reference history cannot start while native internals still live");
  ref::Owner too_large_native(conditional);
  ref::PayloadReservation larger;
  need(!ref::PayloadAccess::allocation(too_large_native, larger,
       {1, ref::native_child_byte_cap + 1, 16, 64,
        ref::Payload::native_batches, ref::Storage::native_internal}) && !larger &&
       too_large_native.failure().cause == ref::Cause::payload_limit,
       "whole128MiB policy does not enlarge original native32MiB fence");
}
void refused_entry_context_and_local_cleanup() {
  auto no_guard = ref::policy_caps;
  no_guard[unsigned(ref::Counter::guards)] = 0;
  ref::Owner early(conditional, no_guard);
  ref::PayloadReservation absent;
  const ref::AllocationRequest original{ref::maximum_node_capacity, sizeof(ProposedNodeLayout),
                                        16, 64, ref::Payload::stage, ref::Storage::reference_nodes};
  need(!ref::PayloadAccess::allocation(early, absent, original) && !absent &&
       early.failure().cause == ref::Cause::work_limit &&
       early.failure().requested_count == original.count &&
       early.failure().requested_element_bytes == original.element_bytes &&
       early.failure().requested_alignment == original.alignment &&
       early.failure().requested_header_bytes == original.header_bytes &&
       early.failure().requested_category == original.category &&
       early.failure().requested_storage == unsigned(original.storage) &&
       early.failure().served_prefix == ref::Counts{},
       "zero-guard unserved entry retains the exact fixed integer payload request");
  need(!ref::PayloadAccess::frame(early, absent, ref::Payload::source, 128) &&
       early.failure().requested_count == original.count &&
       early.failure().requested_header_bytes == original.header_bytes,
       "another request cannot overwrite the first unserved payload context");

  auto late_caps = ref::policy_caps;
  late_caps[unsigned(ref::Counter::transitions)] = 33;
  ref::Owner late(conditional, late_caps);
  ref::PayloadReservation unacquired;
  // This path has23 semantic guards and10 actual iterations before the fixed
  // cleanup reservation; it owns33 served transitions. The reservation is
  // refused BEFORE acquiring a block, and has its own literalU1024 increment.
  need(!ref::PayloadAccess::frame(late, unacquired, ref::Payload::stage, 4096) &&
       !unacquired && late.failure().cause == ref::Cause::work_limit &&
       late.failure().requested_count == 1 && late.failure().requested_element_bytes == 4096 &&
       late.failure().requested_bytes == 4096 &&
       late.failure().served_prefix[unsigned(ref::Counter::transitions)] == 33 &&
       late.failure().requested[unsigned(ref::Counter::integer_upper)] == 1024 &&
       late.cleanup_reserved() == ref::Counts{},
       "late cleanup-credit exhaustion retains request and exact unserved cost");

  ref::Owner local(conditional);
  {
    ref::MacroScope relocation(local, ref::Macro::relocation);
    need(bool(relocation), "original relocation macro admitted");
    const ref::U guard_cap = 8 * ref::macro_destinations[unsigned(ref::Macro::relocation)];
    const ref::U prefix = guard_cap - 23;
    // Synthetic served-prefix control, not executed/scientific guard evidence.
    need(local.action(ref::Action{{0, 0, 0, prefix, prefix, 0, 0, 256 * prefix}}),
         "synthetic local prefix leaves exactly23 entry guards");
    ref::PayloadReservation refused;
    need(!ref::PayloadAccess::frame(local, refused, ref::Payload::stage, 4096) && !refused &&
         local.failure().cause == ref::Cause::macro_limit &&
         local.failure().requested_bytes == 4096 &&
         local.relocation_served()[unsigned(ref::Counter::guards)] == guard_cap &&
         local.failure().requested[unsigned(ref::Counter::integer_upper)] == 1024,
         "payload cleanup cannot escape the original local relocation guard fence");
  }
  ref::Owner completed(conditional);
  {
    ref::MacroScope source(completed, ref::Macro::source);
    ref::PayloadReservation frame;
    need(bool(source) && ref::PayloadAccess::frame(completed, frame, ref::Payload::stage, 4096),
         "frame acquired inside another named macro");
    const auto before = completed.relocation_served();
    need(frame.reset() &&
         completed.relocation_served()[unsigned(ref::Counter::guards)] -
         before[unsigned(ref::Counter::guards)] == 1 &&
         completed.relocation_served()[unsigned(ref::Counter::calls)] -
         before[unsigned(ref::Counter::calls)] == 1 &&
         completed.relocation_served()[unsigned(ref::Counter::integer_upper)] -
         before[unsigned(ref::Counter::integer_upper)] == 1024,
         "actual nested cleanup is served once in the shared relocation owner");
  }
}
}
static_assert(sizeof(ProposedNodeLayout) == 208);
static_assert(!std::is_copy_constructible_v<ref::PayloadReservation>);
static_assert(std::is_nothrow_move_constructible_v<ref::PayloadReservation>);
int main() {
  early_refusal(); old_pending_peak_and_unwind(); ownership_move_boundaries();
  request_and_phase_refusals();
  refused_entry_context_and_local_cleanup();
  std::printf("conditional_drag_reference_ownership checks=%u failures=%u "
              "actual_allocator_called=false whole_payload_qualified=false\n",
              checks, failures);
  return failures ? 1 : 0;
}
