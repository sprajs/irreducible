#pragma once
// Test-only arithmetic/resource owner for the ORIGINAL conditional reference.
// This is neither a physical tail law nor an earned complete error certificate.
// Its strict binary80/basic/adjacency/frame/allocator assumptions need their own
// immutable built admission. The native consumer's 4M/32MiB limits are separate.
#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#if defined(__x86_64__)
#include <xmmintrin.h>
#endif

namespace conditional_drag_reference_certificate {
using W = long double;
using U = std::uint64_t;
inline constexpr bool complete_arithmetic_qualified = false;
inline constexpr bool whole_payload_qualified = false;
enum class Cause : unsigned {
  none, unsupported_profile, invalid_input, nonfinite, unresolved_positive,
  counter_overflow, work_limit, payload_limit, invalid_ownership, macro_limit,
  allocation_failure, unavailable
};
enum class Counter : unsigned {
  destinations, primitives, copies, guards, transitions, iterations, calls,
  integer_upper, count
};
inline constexpr unsigned counter_count = unsigned(Counter::count);
using Counts = std::array<U, counter_count>;
inline constexpr Counts policy_caps{
    U{1} << 44, U{1} << 44, U{1} << 44, U{1} << 47,
    U{1} << 49, U{1} << 40, U{1} << 34, U{1} << 58};
enum class Macro : unsigned {
  rate, quadrature, newton, damping, solve, certificate, cell, partial,
  endpoint, source, saha, clock, knot, peer, relocation, count
};
inline constexpr std::array<U, unsigned(Macro::count)> macro_destinations{
    U{1} << 18, U{1} << 16, U{1} << 12, U{1} << 8, U{1} << 10,
    U{1} << 16, U{1} << 14, U{1} << 14, U{1} << 17, U{1} << 17,
    U{1} << 10, U{1} << 18, U{1} << 12, U{1} << 24, U{1} << 28};
enum class Payload : unsigned {
  native_batches, endpoints, source, stage, elementary, receipts, encoder,
  allocator_metadata, other, count
};
inline constexpr std::array<U, unsigned(Payload::count)> payload_caps{
    U{2} << 20, U{1} << 20, U{512} << 10, U{1} << 20, U{16} << 10,
    U{512} << 10, U{128} << 10, U{1} << 20, U{512} << 10};
inline constexpr U whole_byte_cap = U{128} << 20;
inline constexpr U ancillary_byte_cap = U{8} << 20;
inline constexpr U metadata_per_block = 512;
inline constexpr U maximum_blocks = 2048;
inline constexpr U maximum_node_capacity = 200001;

// These are explicit conditional assumptions, not runtime tests proving them.
// Default construction never fabricates an arithmetic or opaque-frame admission.
struct Profile {
  bool basic_rounding = false, adjacent = false, frames = false;
  bool allocator = false, exceptions = false, formatter = false;
};
struct Location {
  U source = 0, level = 0, endpoint = 0, step = 0, newton = 0, trial = 0;
};
struct Failure {
  Cause cause = Cause::none;
  Location location;
  Counts requested{}, served_prefix{};
  U live_bytes = 0, pending_bytes = 0, peak_bytes = 0;
  U requested_bytes = 0;
};
struct Action {
  Counts count{};
  static constexpr Action floating(U destinations, U primitives,
                                   U copies = 0) noexcept {
    return {{{destinations, primitives, copies, 0, 1, 0, 0, 128}}};
  }
  static constexpr Action guard() noexcept {
    return {{{0, 0, 0, 1, 1, 0, 0, 256}}};
  }
  static constexpr Action call() noexcept {
    return {{{0, 0, 0, 0, 1, 0, 1, 640}}};
  }
};

class Owner {
  struct MacroFrame { Macro id; Counts spent; };
  struct Block {
    U bytes = 0, generation = 0;
    Payload category = Payload::other;
    bool used = false, pending = false, node = false;
  };
  Profile profile_;
  Counts caps_, served_{}, cleanup_reserved_{};
  std::array<MacroFrame, 16> macros_{};
  unsigned depth_ = 0;
  std::array<Block, maximum_blocks> blocks_{};
  std::array<U, unsigned(Payload::count)> ancillary_{};
  U ancillary_live_ = 0, live_ = 0, pending_ = 0, peak_ = 0;
  U next_generation_ = 1;
  Location location_;
  Failure failure_;
  bool arithmetic_admitted_ = false;
  static constexpr U at(Counter c) noexcept { return unsigned(c); }
  static constexpr Counts cleanup_cost{0, 0, 0, 1, 1, 0, 0, 512};
  static bool fits(U spent, U increment, U cap) noexcept {
    return spent <= cap && increment <= cap - spent;
  }
  Counts macro_caps(Macro id) const noexcept {
    const U c = macro_destinations[unsigned(id)];
    return {c, c, c, 8 * c, 32 * c, caps_[at(Counter::iterations)],
            caps_[at(Counter::calls)], caps_[at(Counter::integer_upper)]};
  }
  bool reserve_cleanup() noexcept {
    for (unsigned i = 0; i < counter_count; ++i) {
      if (!fits(served_[i], cleanup_reserved_[i], caps_[i]) ||
          !fits(served_[i] + cleanup_reserved_[i], cleanup_cost[i], caps_[i]))
        return fail(Cause::work_limit, cleanup_cost);
    }
    for (unsigned i = 0; i < counter_count; ++i)
      cleanup_reserved_[i] += cleanup_cost[i];
    return true;
  }
  void consume_cleanup() noexcept {
    // Admission reserved this fixed integer-only exit before acquiring the
    // scope. It remains executable after failure without a floating epilogue.
    // The first failure's served_prefix stays the pre-failure snapshot; the
    // final served ledger also includes actual cleanup on the way out.
    for (unsigned i = 0; i < counter_count; ++i) {
      cleanup_reserved_[i] -= cleanup_cost[i];
      served_[i] += cleanup_cost[i];
    }
  }
public:
  // No Wide member is initialized by this owner. All failure metadata is fixed.
  explicit Owner(Profile p = {}, Counts caps = policy_caps) noexcept
      : profile_(p), caps_(caps) {
    for (unsigned i = 0; i < counter_count; ++i)
      if (caps_[i] > policy_caps[i]) {
        fail(Cause::invalid_input); return;
      }
    // The selected profile must later earn actual caller/frame coexistence.
    // Account the complete fixed owner, including all 2048 bookkeeping slots.
    const U bytes = sizeof(*this);
    if (bytes > payload_caps[unsigned(Payload::receipts)]) {
      fail(Cause::payload_limit); return;
    }
    ancillary_[unsigned(Payload::receipts)] = bytes;
    ancillary_live_ = live_ = peak_ = bytes;
  }
  Owner(const Owner &) = delete;
  Owner &operator=(const Owner &) = delete;
  Owner(Owner &&) = delete;
  Owner &operator=(Owner &&) = delete;
  bool ok() const noexcept { return failure_.cause == Cause::none; }
  const Failure &failure() const noexcept { return failure_; }
  const Counts &served() const noexcept { return served_; }
  const Counts &caps() const noexcept { return caps_; }
  U live_bytes() const noexcept { return live_; }
  U pending_bytes() const noexcept { return pending_; }
  U peak_bytes() const noexcept { return peak_; }
  const Counts &cleanup_reserved() const noexcept { return cleanup_reserved_; }
  void location(Location v) noexcept { location_ = v; }
  bool fail(Cause why, Counts requested = {}) noexcept {
    if (ok()) {
      failure_.cause = why; failure_.location = location_;
      failure_.requested = requested; failure_.served_prefix = served_;
      failure_.live_bytes = live_; failure_.pending_bytes = pending_;
      failure_.peak_bytes = peak_;
    }
    return false;
  }
  bool action(const Action &a) noexcept {
    if (!ok()) return false;
    if (a.count[at(Counter::copies)] > a.count[at(Counter::destinations)])
      return fail(Cause::invalid_input, a.count);
    // Transactional: no field is committed before ALL local/global checks.
    const Counts local = depth_ ? macro_caps(macros_[depth_ - 1].id) : caps_;
    for (unsigned i = 0; i < counter_count; ++i) {
      if (served_[i] > std::numeric_limits<U>::max() - a.count[i])
        return fail(Cause::counter_overflow, a.count);
      if (!fits(served_[i], a.count[i], caps_[i]))
        return fail(Cause::work_limit, a.count);
      if (!fits(served_[i] + a.count[i], cleanup_reserved_[i], caps_[i]))
        return fail(Cause::work_limit, a.count);
      if (depth_ && !fits(macros_[depth_ - 1].spent[i], a.count[i], local[i]))
        return fail(Cause::macro_limit, a.count);
    }
    for (unsigned i = 0; i < counter_count; ++i) {
      served_[i] += a.count[i];
      if (depth_) macros_[depth_ - 1].spent[i] += a.count[i];
    }
    return true;
  }
  bool guard() noexcept { return action(Action::guard()); }
  bool iteration() noexcept {
    return action(Action{{0, 0, 0, 0, 1, 1, 0, 32}});
  }
  bool enter(Macro id) noexcept {
    if (!action(Action::call())) return false;
    if (unsigned(id) >= unsigned(Macro::count) || depth_ == macros_.size())
      return fail(Cause::invalid_input);
    const Counts reservation = macro_caps(id);
    for (unsigned i = 0; i < counter_count; ++i)
      if ((i <= at(Counter::transitions)) &&
          (!fits(served_[i], cleanup_reserved_[i], caps_[i]) ||
           !fits(served_[i] + cleanup_reserved_[i], reservation[i], caps_[i]) ||
           !fits(served_[i] + cleanup_reserved_[i] + reservation[i],
                 cleanup_cost[i], caps_[i])))
        return fail(Cause::work_limit, reservation);
    if (!reserve_cleanup()) return false;
    macros_[depth_++] = {id, {}};
    return true;
  }
  void leave() noexcept {
    // Integer-only fixed cleanup, also available after a causal refusal.
    if (depth_) { consume_cleanup(); --depth_; }
  }
  bool admit_arithmetic() noexcept {
    if (!guard()) return false;
    if (!(profile_.basic_rounding && profile_.adjacent && profile_.frames))
      return fail(Cause::unsupported_profile);
    if (!action(Action::call())) return false;
    if (!(sizeof(W) == 16 && alignof(W) <= 16 &&
          std::numeric_limits<W>::radix == 2 &&
          std::numeric_limits<W>::digits == 64 &&
          std::numeric_limits<W>::min_exponent == -16381 &&
          std::numeric_limits<W>::max_exponent == 16384 &&
          std::fegetround() == FE_TONEAREST))
      return fail(Cause::unsupported_profile);
#if defined(__x86_64__)
    unsigned short control = 0;
    __asm__ __volatile__("fnstcw %0" : "=m"(control));
    const unsigned mxcsr = _mm_getcsr();
    if (!guard()) return false;
    if ((control & 0x0f3f) != 0x033f || (mxcsr & 0xffc0) != 0x1f80)
      return fail(Cause::unsupported_profile);
#else
    return fail(Cause::unsupported_profile);
#endif
    arithmetic_admitted_ = true;
    return true;
  }
  bool arithmetic_ready() noexcept {
    if (!guard()) return false;
    return arithmetic_admitted_ || fail(Cause::unsupported_profile);
  }
  bool opaque_science_ready() noexcept {
    if (!arithmetic_ready() || !guard()) return false;
    return (profile_.allocator && profile_.exceptions && profile_.formatter) ||
           fail(Cause::unsupported_profile);
  }
};

class MacroScope {
  Owner *owner_ = nullptr;
public:
  MacroScope(Owner &o, Macro id) noexcept { if (o.enter(id)) owner_ = &o; }
  ~MacroScope() { if (owner_) owner_->leave(); }
  MacroScope(const MacroScope &) = delete;
  MacroScope &operator=(const MacroScope &) = delete;
  MacroScope(MacroScope &&) = delete;
  MacroScope &operator=(MacroScope &&) = delete;
  explicit operator bool() const noexcept { return owner_ != nullptr; }
};

// Unassigned Wide storage: brace construction performs no hidden FP zeroing.
// All outputs are in place. No Interval can cross a layer by value or copy.
struct Interval {
  W lower, upper;
  bool valid;
  Interval() noexcept : valid(false) {}
  Interval(const Interval &) = delete;
  Interval &operator=(const Interval &) = delete;
  Interval(Interval &&) = delete;
  Interval &operator=(Interval &&) = delete;
};
enum class Direction : unsigned { lower, upper };
enum class Operation : unsigned { add, subtract, multiply, divide };

class Arithmetic {
  Owner &owner_;
  W positive_direction_, negative_direction_;
  bool directions_prepared_ = false;
  bool normal(const W &v) noexcept {
    if (!owner_.guard()) return false;
    return (v == 0 || std::isnormal(v)) || owner_.fail(Cause::nonfinite);
  }
  bool inputs(const Interval &a, const Interval &b) noexcept {
    if (!owner_.arithmetic_ready() || !owner_.guard()) return false;
    if (!(a.valid && b.valid)) return owner_.fail(Cause::invalid_input);
    if (!normal(a.lower) || !normal(a.upper) || !normal(b.lower) ||
        !normal(b.upper) || !owner_.guard()) return false;
    return (a.lower <= a.upper && b.lower <= b.upper) ||
           owner_.fail(Cause::invalid_input);
  }
  bool output(Interval &out, const Interval &a, const Interval &b) noexcept {
    if (!owner_.guard()) return false;
    if (&out == &a || &out == &b) return owner_.fail(Cause::invalid_ownership);
    out.valid = false;
    return inputs(a, b);
  }
  bool directed(W &out, const W &a, const W &b, Operation op,
                Direction direction) noexcept {
    if (!owner_.guard()) return false;
    if (!directions_prepared_) return owner_.fail(Cause::unavailable);
    if (!owner_.action(Action::floating(1, 1))) return false;
    W rounded; // assigned only AFTER its action was admitted
    switch (op) {
    case Operation::add: rounded = a + b; break;
    case Operation::subtract: rounded = a - b; break;
    case Operation::multiply: rounded = a * b; break;
    case Operation::divide: rounded = a / b; break;
    }
    if (!normal(rounded) || !owner_.guard()) return false;
    if (rounded == 0) {
      // Witness EXACT arithmetic of the stored endpoints. This does not label
      // an entire uncertain physical function zero (e.g. cold exp lower0).
      bool exact = false;
      switch (op) {
      case Operation::add: {
        if (!owner_.action(Action::floating(1, 1))) return false;
        W opposite; opposite = -b;
        if (!owner_.guard()) return false;
        exact = a == opposite;
        break;
      }
      case Operation::subtract: exact = a == b; break;
      case Operation::multiply: exact = a == 0 || b == 0; break;
      case Operation::divide: exact = a == 0 && b != 0; break;
      }
      if (!exact) return owner_.fail(Cause::unresolved_positive);
      if (!owner_.action(Action::floating(1, 0, 1))) return false;
      out = rounded;
      return true;
    }
    const W &toward = direction == Direction::lower ? negative_direction_ :
                                                      positive_direction_;
    if (!owner_.guard()) return false;
    if (rounded == toward) return owner_.fail(Cause::nonfinite);
    if (!owner_.action(Action::floating(1, 1))) return false;
    out = std::nextafter(rounded, toward);
    return normal(out);
  }
public:
  explicit Arithmetic(Owner &o) noexcept : owner_(o) {}
  Arithmetic(const Arithmetic &) = delete;
  Arithmetic &operator=(const Arithmetic &) = delete;
  Arithmetic(Arithmetic &&) = delete;
  Arithmetic &operator=(Arithmetic &&) = delete;
  // Two admitted source destinations ONCE, outside individual interval DAGs.
  // Subsequent nextafter direction arguments borrow these owned exact values.
  bool prepare_directions() noexcept {
    if (!owner_.arithmetic_ready() || !owner_.guard()) return false;
    if (directions_prepared_) return owner_.fail(Cause::invalid_input);
    if (!owner_.action(Action::floating(1, 0, 1))) return false;
    positive_direction_ = std::numeric_limits<W>::max();
    if (!owner_.action(Action::floating(1, 1))) return false;
    negative_direction_ = -positive_direction_;
    directions_prepared_ = true;
    return true;
  }
  bool enclose(Interval &out, const W &stored_lower, const W &stored_upper) noexcept {
    if (!owner_.guard()) return false;
    if (&stored_lower == &out.lower || &stored_lower == &out.upper ||
        &stored_upper == &out.lower || &stored_upper == &out.upper)
      return owner_.fail(Cause::invalid_ownership);
    out.valid = false;
    if (!owner_.arithmetic_ready() || !normal(stored_lower) || !normal(stored_upper) ||
        !owner_.guard()) return false;
    if (stored_lower > stored_upper) return owner_.fail(Cause::invalid_input);
    if (!owner_.action(Action::floating(2, 0, 2))) return false;
    out.lower = stored_lower; out.upper = stored_upper; out.valid = true;
    return true;
  }
  bool point(Interval &out, const W &exact_stored_value) noexcept {
    return enclose(out, exact_stored_value, exact_stored_value);
  }
  bool copy(Interval &out, const Interval &in) noexcept {
    if (!owner_.guard()) return false;
    if (&out == &in) return owner_.fail(Cause::invalid_ownership);
    out.valid = false;
    if (!inputs(in, in)) return false;
    if (!owner_.action(Action::floating(2, 0, 2))) return false;
    out.lower = in.lower; out.upper = in.upper; out.valid = true;
    return true;
  }
  bool add(Interval &out, const Interval &a, const Interval &b) noexcept {
    if (!output(out, a, b) ||
        !directed(out.lower, a.lower, b.lower, Operation::add, Direction::lower) ||
        !directed(out.upper, a.upper, b.upper, Operation::add, Direction::upper))
      return false;
    out.valid = true;
    return true;
  }
  bool subtract(Interval &out, const Interval &a, const Interval &b) noexcept {
    if (!output(out, a, b) ||
        !directed(out.lower, a.lower, b.upper, Operation::subtract, Direction::lower) ||
        !directed(out.upper, a.upper, b.lower, Operation::subtract, Direction::upper))
      return false;
    out.valid = true;
    return true;
  }
  bool multiply(Interval &out, const Interval &a, const Interval &b) noexcept {
    if (!output(out, a, b)) return false;
    std::array<W, 4> lo, hi; // uninitialized until each admitted action
    unsigned i = 0;
    for (;;) {
      if (!owner_.guard()) return false;
      if (i == 4) break;
      if (!owner_.iteration()) return false;
      const W &x = i < 2 ? a.lower : a.upper;
      const W &y = i % 2 ? b.upper : b.lower;
      if (!directed(lo[i], x, y, Operation::multiply, Direction::lower) ||
          !directed(hi[i], x, y, Operation::multiply, Direction::upper)) return false;
      ++i;
    }
    unsigned il = 0, ih = 0;
    i = 1;
    for (;;) {
      if (!owner_.guard()) return false;
      if (i == 4) break;
      if (!owner_.iteration()) return false;
      if (!owner_.guard()) return false;
      if (lo[i] < lo[il]) il = i;
      if (!owner_.guard()) return false;
      if (hi[i] > hi[ih]) ih = i;
      ++i;
    }
    if (!owner_.action(Action::floating(2, 0, 2))) return false;
    out.lower = lo[il]; out.upper = hi[ih]; out.valid = true;
    return true;
  }
  bool divide(Interval &out, const Interval &a, const Interval &b) noexcept {
    if (!output(out, a, b) || !owner_.guard()) return false;
    if (!(b.lower > 0 || b.upper < 0)) return owner_.fail(Cause::invalid_input);
    if (!owner_.action(Action::floating(1, 0, 1))) return false;
    W one; one = 1;
    Interval reciprocal;
    if (!directed(reciprocal.lower, one, b.upper, Operation::divide, Direction::lower) ||
        !directed(reciprocal.upper, one, b.lower, Operation::divide, Direction::upper))
      return false;
    reciprocal.valid = true;
    return multiply(out, a, reciprocal);
  }
};
} // namespace conditional_drag_reference_certificate
