#pragma once
// Conditional request-owned payload reservations. These counters do not prove
// libm/ABI/allocator/exception frame sizes. Default profiles refuse; actual
// allocation/callback wiring and immutable measured/profile admission remain
// prerequisites for whole_payload_qualified (still false).
#include "conditional_drag_reference_certificate.hpp"

namespace conditional_drag_reference_certificate {
enum class Storage : unsigned { ancillary, reference_nodes, native_internal };
// This is the original CHILD native byte ceiling, confined to that phase.
// It is not a reference byte/scalar ceiling or a category-work conversion.
inline constexpr U native_child_byte_cap = U{32} << 20;
struct AllocationRequest {
  U count = 0, element_bytes = 0, alignment = 0, header_bytes = 0;
  Payload category = Payload::other;
  Storage storage = Storage::ancillary;
};
class PayloadReservation {
  Owner *owner_ = nullptr;
  U slot_ = 0, generation_ = 0;
  friend struct PayloadAccess;
public:
  PayloadReservation() noexcept = default;
  PayloadReservation(const PayloadReservation &) = delete;
  PayloadReservation &operator=(const PayloadReservation &) = delete;
  PayloadReservation(PayloadReservation &&other) noexcept;
  PayloadReservation &operator=(PayloadReservation &&other) noexcept;
  ~PayloadReservation();
  explicit operator bool() const noexcept { return owner_ != nullptr; }
  bool commit() noexcept;
  bool reset() noexcept;
};

struct PayloadAccess {
private:
  static bool refuse(Owner &o, Cause cause, const AllocationRequest &r,
                     U requested_bytes = 0) noexcept {
    if (!o.ok()) return false;
    o.fail(cause);
    o.failure_.requested_count = r.count;
    o.failure_.requested_element_bytes = r.element_bytes;
    o.failure_.requested_alignment = r.alignment;
    o.failure_.requested_header_bytes = r.header_bytes;
    o.failure_.requested_category = r.category;
    o.failure_.requested_storage = unsigned(r.storage);
    o.failure_.requested_bytes = requested_bytes;
    return false;
  }
  static bool acquire(Owner &o, PayloadReservation &out,
                      const AllocationRequest &r, bool heap) noexcept {
    // Fixed integer context precedes the first counter guard, just as the
    // first-cause snapshot exists before any floating object. This bounded
    // refusal path annotates even an unserved entry, without changing an older
    // causal receipt. It performs no floating arithmetic or formatting.
    if (!o.ok()) return false;
    o.failure_.requested_count = r.count;
    o.failure_.requested_element_bytes = r.element_bytes;
    o.failure_.requested_alignment = r.alignment;
    o.failure_.requested_header_bytes = r.header_bytes;
    o.failure_.requested_category = r.category;
    o.failure_.requested_storage = unsigned(r.storage);
    o.failure_.requested_bytes = 0;
    if (!o.guard()) return false;
    if (out.owner_ || unsigned(r.category) >= unsigned(Payload::count))
      return refuse(o, Cause::invalid_ownership, r);
    if (!o.guard()) return false;
    if (!(o.profile_.frames && (!heap || o.profile_.allocator)))
      return refuse(o, Cause::unsupported_profile, r);
    if (!o.guard()) return false;
    if (!r.count || !r.element_bytes || !r.alignment ||
        r.alignment > 16 || (r.alignment & (r.alignment - 1)))
      return refuse(o, Cause::invalid_input, r);
    if (!o.guard()) return false;
    if (unsigned(r.storage) > unsigned(Storage::native_internal) ||
        (!heap && r.storage != Storage::ancillary))
      return refuse(o, Cause::invalid_input, r);
    const bool node = r.storage == Storage::reference_nodes;
    const bool native = r.storage == Storage::native_internal;
    if (!o.guard()) return false;
    if (node && (r.count > maximum_node_capacity || r.element_bytes > 208 ||
                 o.node_blocks_ >= 2 || o.native_data_ != 0))
      return refuse(o, Cause::payload_limit, r);
    if (!o.guard()) return false;
    if (r.count > std::numeric_limits<U>::max() / r.element_bytes)
      return refuse(o, Cause::counter_overflow, r);
    const U data = r.count * r.element_bytes;
    if (!o.guard()) return false;
    if (native && (o.node_blocks_ != 0 ||
                   !Owner::fits(o.native_data_, data, native_child_byte_cap)))
      return refuse(o, Cause::payload_limit, r, data);
    if (!o.guard()) return false;
    if (heap && r.header_bytes > std::numeric_limits<U>::max() -
                               (r.alignment - 1) - metadata_per_block)
      return refuse(o, Cause::counter_overflow, r, data);
    const U overhead = heap ? r.header_bytes + r.alignment - 1 + metadata_per_block : 0;
    if (!o.guard()) return false;
    if (data > std::numeric_limits<U>::max() - overhead)
      return refuse(o, Cause::counter_overflow, r, data);
    const U total = data + overhead;
    o.failure_.requested_bytes = total;
    // Treat two increments into allocator_metadata as one checked category
    // increment, rather than testing each independently against the same room.
    std::array<U, unsigned(Payload::count)> delta{};
    if (!node && !native) delta[unsigned(r.category)] = data;
    if (!o.guard()) return false;
    if (delta[unsigned(Payload::allocator_metadata)] >
        std::numeric_limits<U>::max() - overhead)
      return refuse(o, Cause::counter_overflow, r, total);
    delta[unsigned(Payload::allocator_metadata)] += overhead;
    const U ancillary = node || native ? overhead : total;
    for (unsigned i = 0; i < delta.size(); ++i) {
      if (!o.guard() || !o.iteration()) return false;
      if (!Owner::fits(o.ancillary_[i], delta[i], payload_caps[i]))
        return refuse(o, Cause::payload_limit, r, total);
    }
    if (!o.guard()) return false;
    if (!Owner::fits(o.ancillary_live_, ancillary, ancillary_byte_cap) ||
        !Owner::fits(o.live_, o.pending_, whole_byte_cap) ||
        !Owner::fits(o.live_ + o.pending_, total, whole_byte_cap))
      return refuse(o, Cause::payload_limit, r, total);
    U slot = 0;
    for (;;) {
      if (!o.guard()) return false;
      if (slot == maximum_blocks) return refuse(o, Cause::payload_limit, r, total);
      if (!o.iteration() || !o.guard()) return false;
      if (!o.blocks_[slot].used) break;
      ++slot;
    }
    if (!o.guard()) return false;
    if (o.next_generation_ == std::numeric_limits<U>::max())
      return refuse(o, Cause::counter_overflow, r, total);
    if (!o.reserve_payload_cleanup()) return false;
    const U generation = o.next_generation_++;
    o.blocks_[slot] = {data, overhead, generation, r.category, true, heap, node, heap, native};
    for (unsigned i = 0; i < delta.size(); ++i) o.ancillary_[i] += delta[i];
    o.ancillary_live_ += ancillary;
    if (node) ++o.node_blocks_;
    if (native) o.native_data_ += data;
    if (heap) o.pending_ += total; else o.live_ += total;
    const U simultaneous = o.live_ + o.pending_;
    if (simultaneous > o.peak_) o.peak_ = simultaneous;
    out.owner_ = &o; out.slot_ = slot; out.generation_ = generation;
    // A later unrelated arithmetic refusal must not report this successful
    // payload request as though it were the failed attempt.
    o.failure_.requested_count = o.failure_.requested_element_bytes = 0;
    o.failure_.requested_alignment = o.failure_.requested_header_bytes = 0;
    o.failure_.requested_category = Payload::other;
    o.failure_.requested_storage = 0; o.failure_.requested_bytes = 0;
    return true;
  }
  static bool valid(const PayloadReservation &r) noexcept {
    return r.owner_ && r.slot_ < maximum_blocks &&
      r.owner_->blocks_[r.slot_].used &&
      r.owner_->blocks_[r.slot_].generation == r.generation_;
  }
public:
  static bool frame(Owner &o, PayloadReservation &out, Payload category,
                    U actual_bytes, U alignment = 16) noexcept {
    return acquire(o, out, {1, actual_bytes, alignment, 0, category, Storage::ancillary}, false);
  }
  static bool allocation(Owner &o, PayloadReservation &out,
                         const AllocationRequest &request) noexcept {
    return acquire(o, out, request, true);
  }
  static bool commit(PayloadReservation &r) noexcept {
    if (!r.owner_) return false;
    Owner &o = *r.owner_;
    if (!o.guard()) return false;
    if (!valid(r)) return o.fail(Cause::invalid_ownership);
    auto &b = o.blocks_[r.slot_];
    if (!o.guard()) return false;
    if (!b.heap || !b.pending) return o.fail(Cause::invalid_ownership);
    const U total = b.bytes + b.overhead;
    o.pending_ -= total; o.live_ += total; b.pending = false;
    return true;
  }
  static bool release(PayloadReservation &r) noexcept {
    if (!r.owner_) return true;
    Owner &o = *r.owner_;
    // Fixed integer cleanup is served from a previously reserved exit. This
    // function calls no allocator/destructor and performs no floating work.
    if (!valid(r)) return o.fail(Cause::invalid_ownership);
    o.consume_payload_cleanup();
    auto &b = o.blocks_[r.slot_];
    const U total = b.bytes + b.overhead;
    if (b.pending) o.pending_ -= total; else o.live_ -= total;
    if (!b.node && !b.native) o.ancillary_[unsigned(b.category)] -= b.bytes;
    o.ancillary_[unsigned(Payload::allocator_metadata)] -= b.overhead;
    o.ancillary_live_ -= b.node || b.native ? b.overhead : total;
    if (b.node) --o.node_blocks_;
    if (b.native) o.native_data_ -= b.bytes;
    b.used = false;
    r.owner_ = nullptr; r.slot_ = 0; r.generation_ = 0;
    return true;
  }
  static bool transfer(PayloadReservation &target, PayloadReservation &source) noexcept {
    if (&target == &source) return !source.owner_ || source.owner_->guard();
    if (!source.owner_) return release(target);
    Owner &o = *source.owner_;
    if (!o.guard()) return false;
    if (!valid(source) || (target.owner_ && target.owner_ != &o))
      return o.fail(Cause::invalid_ownership);
    if (!release(target)) return false;
    target.owner_ = source.owner_; target.slot_ = source.slot_;
    target.generation_ = source.generation_;
    source.owner_ = nullptr; source.slot_ = 0; source.generation_ = 0;
    return true;
  }
};
inline PayloadReservation::PayloadReservation(PayloadReservation &&other) noexcept {
  PayloadAccess::transfer(*this, other);
}
inline PayloadReservation &PayloadReservation::operator=(PayloadReservation &&other) noexcept {
  PayloadAccess::transfer(*this, other); return *this;
}
inline PayloadReservation::~PayloadReservation() { PayloadAccess::release(*this); }
inline bool PayloadReservation::commit() noexcept { return PayloadAccess::commit(*this); }
inline bool PayloadReservation::reset() noexcept { return PayloadAccess::release(*this); }
} // namespace conditional_drag_reference_certificate
