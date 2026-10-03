#include "conditional_drag_reference_certificate.hpp"
#include <cstdio>
#include <type_traits>

namespace ref = conditional_drag_reference_certificate;
namespace {
unsigned checks = 0, failures = 0;
void need(bool value, const char *label) {
  ++checks;
  if (!value) { ++failures; std::fprintf(stderr, "FAILED %s\n", label); }
}
// Finite controls declare the selected arithmetic assumptions conditionally.
// They do NOT earn the compiler/adjacency/frame/whole-payload profile by fiat.
constexpr ref::Profile conditional_arithmetic{true, true, true, false, false, false};
void unsupported_profile() {
  ref::Owner owner;
  need(!owner.admit_arithmetic() &&
       owner.failure().cause == ref::Cause::unsupported_profile,
       "missing arithmetic assumptions refuse before arithmetic");
  need(owner.served()[unsigned(ref::Counter::destinations)] == 0 &&
       owner.served()[unsigned(ref::Counter::primitives)] == 0,
       "unsupported early profile has no invented floating prefix");
  const auto original = owner.failure();
  need(!owner.action(ref::Action::floating(2, 1)) &&
       owner.failure().cause == original.cause &&
       owner.failure().served_prefix == original.served_prefix,
       "first causal failure and served prefix survive later attempted work");
}
void transactional_limits() {
  auto caps = ref::policy_caps;
  caps[unsigned(ref::Counter::destinations)] = 1;
  ref::Owner owner(conditional_arithmetic, caps);
  const auto before = owner.served();
  need(!owner.action(ref::Action::floating(2, 1, 2)) &&
       owner.failure().cause == ref::Cause::work_limit &&
       owner.served() == before,
       "unserved multi-field action never partially commits");
  need(owner.failure().requested[unsigned(ref::Counter::destinations)] == 2 &&
       owner.failure().requested[unsigned(ref::Counter::copies)] == 2,
       "exact refused increment is retained");
  ref::Owner overflow(conditional_arithmetic);
  need(overflow.action(ref::Action::floating(1, 0, 1)), "nonzero prefix admitted");
  const auto prefix = overflow.served();
  auto increment = ref::Action::floating(std::numeric_limits<ref::U>::max(), 0);
  need(!overflow.action(increment) &&
       overflow.failure().cause == ref::Cause::counter_overflow &&
       overflow.served() == prefix,
       "nonzero plus UINT64_MAX cannot wrap into freshness or commit");
  ref::Owner local(conditional_arithmetic);
  ref::Counts refusal_prefix{};
  {
    ref::MacroScope scope(local, ref::Macro::damping);
    need(bool(scope), "fixed damping macro reservation admitted");
    need(local.action(ref::Action::floating(256, 0, 256)),
         "local fence serves its exact boundary");
    const auto local_prefix = local.served();
    need(!local.action(ref::Action::floating(1, 0, 1)) &&
         local.failure().cause == ref::Cause::macro_limit &&
         local.served() == local_prefix,
         "local exhaustion refuses despite ample aggregate remaining");
    refusal_prefix = local.failure().served_prefix;
    need(local.cleanup_reserved()[unsigned(ref::Counter::transitions)] == 1,
         "scope owns its integer cleanup before the terminal refusal");
  }
  need(local.failure().served_prefix == refusal_prefix &&
       local.served()[unsigned(ref::Counter::transitions)] ==
       refusal_prefix[unsigned(ref::Counter::transitions)] + 1 &&
       local.cleanup_reserved() == ref::Counts{},
       "actual post-refusal cleanup counts without erasing first failure");
}
void signed_enclosures() {
  ref::Owner owner(conditional_arithmetic);
  need(owner.admit_arithmetic(), "selected conditional arithmetic environment observed");
  ref::Arithmetic arithmetic(owner);
  const auto setup_before = owner.served();
  need(arithmetic.prepare_directions(), "once-owned direction source admitted");
  need(owner.served()[unsigned(ref::Counter::destinations)] -
       setup_before[unsigned(ref::Counter::destinations)] == 2 &&
       owner.served()[unsigned(ref::Counter::primitives)] -
       setup_before[unsigned(ref::Counter::primitives)] == 1 &&
       owner.served()[unsigned(ref::Counter::copies)] -
       setup_before[unsigned(ref::Counter::copies)] == 1,
       "once-owned direction setup reports S2/primitive1/C1 exactly");
  ref::Interval a, b, product, quotient, sum, difference;
  const ref::W amin = -3, amax = 2, bmin = -2, bmax = 4;
  need(arithmetic.enclose(a, amin, amax) && arithmetic.enclose(b, bmin, bmax),
       "mixed-sign source intervals retained in place");
  const auto product_before = owner.served();
  need(arithmetic.multiply(product, a, b) && product.lower <= -12 && product.upper >= 8,
       "all four corners enclose exact mixed-sign product extrema");
  need(owner.served()[unsigned(ref::Counter::destinations)] -
       product_before[unsigned(ref::Counter::destinations)] <= 18,
       "actual product respects admitted18-destination graph fence");
  need(owner.served()[unsigned(ref::Counter::iterations)] -
       product_before[unsigned(ref::Counter::iterations)] == 7,
       "four corner and three selection iterations are actually charged");
  need(arithmetic.add(sum, a, b) && sum.lower <= -5 && sum.upper >= 6,
       "signed sum encloses independent integer extrema");
  need(arithmetic.subtract(difference, a, b) &&
       difference.lower <= -7 && difference.upper >= 4,
       "signed subtraction uses cross endpoints");
  ref::Interval positive;
  const ref::W two = 2, four = 4;
  need(arithmetic.enclose(positive, two, four), "denominator retained");
  const auto division_before = owner.served();
  need(arithmetic.divide(quotient, a, positive) &&
       quotient.lower <= -1.5L && quotient.upper >= 1,
       "division encloses independent rational extrema");
  need(owner.served()[unsigned(ref::Counter::destinations)] -
       division_before[unsigned(ref::Counter::destinations)] <= 32,
       "actual division respects admitted32-destination graph fence");
  need(owner.served()[unsigned(ref::Counter::copies)] <=
       owner.served()[unsigned(ref::Counter::destinations)],
       "copies remain a destination subcount");
  need(!owner.opaque_science_ready() &&
       owner.failure().cause == ref::Cause::unsupported_profile,
       "unearned allocator/exception/formatter blocks scientific execution");
}
void zero_and_positive_underflow() {
  ref::Owner zero_owner(conditional_arithmetic);
  need(zero_owner.admit_arithmetic(), "zero control environment");
  ref::Arithmetic zero_arithmetic(zero_owner);
  need(zero_arithmetic.prepare_directions(), "once-owned direction source admitted");
  ref::Interval a, b, zero;
  const ref::W plus = 0.5L, minus = -0.5L;
  need(zero_arithmetic.point(a, plus) && zero_arithmetic.point(b, minus) &&
       zero_arithmetic.add(zero, a, b) && zero.lower == 0 && zero.upper == 0,
       "exact stored dyadic cancellation retains genuine zero witness");
  ref::Owner small_owner(conditional_arithmetic);
  need(small_owner.admit_arithmetic(), "positive underflow control environment");
  ref::Arithmetic small_arithmetic(small_owner);
  need(small_arithmetic.prepare_directions(), "once-owned direction source admitted");
  ref::Interval minimum, positive_product;
  const ref::W least_normal = std::numeric_limits<ref::W>::min();
  need(small_arithmetic.point(minimum, least_normal), "positive normal input present");
  need(!small_arithmetic.multiply(positive_product, minimum, minimum) &&
       small_owner.failure().cause == ref::Cause::unresolved_positive &&
       !positive_product.valid,
       "mathematically positive product underflow cannot become accepted zero");
}
void domain_and_alias_refusals() {
  ref::Owner owner(conditional_arithmetic);
  need(owner.admit_arithmetic(), "alias control environment");
  ref::Arithmetic arithmetic(owner);
  need(arithmetic.prepare_directions(), "once-owned direction source admitted");
  ref::Interval a, b;
  const ref::W one = 1, two = 2;
  need(arithmetic.point(a, one) && arithmetic.point(b, two), "alias inputs");
  need(!arithmetic.add(a, a, b) && owner.failure().cause == ref::Cause::invalid_ownership &&
       a.valid && a.lower == one && a.upper == one,
       "output alias refuses before damaging the acquired input");
  ref::Owner denominator_owner(conditional_arithmetic);
  need(denominator_owner.admit_arithmetic(), "denominator control environment");
  ref::Arithmetic denominator_arithmetic(denominator_owner);
  need(denominator_arithmetic.prepare_directions(), "once-owned direction source admitted");
  ref::Interval numerator, denominator, output;
  const ref::W minus = -1;
  need(denominator_arithmetic.point(numerator, one) &&
       denominator_arithmetic.enclose(denominator, minus, one), "zero-crossing denominator");
  need(!denominator_arithmetic.divide(output, numerator, denominator) &&
       denominator_owner.failure().cause == ref::Cause::invalid_input && !output.valid,
       "zero-crossing denominator has typed refusal");
}
void isolated_copy_refusals() {
  // Corrupted inputs are intentionally manufactured by this adversarial
  // fixture. They are not admitted science objects or an earned profile.
  const ref::W one = 1, two = 2;
  ref::Interval input, copied;
  input.lower = one; input.upper = two; input.valid = true;
  ref::Owner unsupported;
  ref::Arithmetic unprepared(unsupported);
  need(!unprepared.copy(copied, input) &&
       unsupported.failure().cause == ref::Cause::unsupported_profile && !copied.valid &&
       unsupported.served()[unsigned(ref::Counter::destinations)] == 0,
       "isolated copy cannot bypass missing arithmetic admission");
  for (unsigned i = 0; i != 3; ++i) {
    ref::Owner owner(conditional_arithmetic);
    need(owner.admit_arithmetic(), "corrupt copy conditional environment");
    ref::Arithmetic arithmetic(owner);
    ref::Interval corrupt, out;
    corrupt.lower = i == 0 ? std::numeric_limits<ref::W>::quiet_NaN() :
                    i == 1 ? std::numeric_limits<ref::W>::denorm_min() : two;
    corrupt.upper = one; corrupt.valid = true;
    const auto before = owner.served()[unsigned(ref::Counter::destinations)];
    need(!arithmetic.copy(out, corrupt) && !out.valid &&
         owner.served()[unsigned(ref::Counter::destinations)] == before &&
         owner.failure().cause == (i == 2 ? ref::Cause::invalid_input : ref::Cause::nonfinite),
         "isolated copy rejects NaN/subnormal/reversed endpoints before writes");
  }
}
}
static_assert(!std::is_copy_constructible_v<ref::Owner>);
static_assert(!std::is_move_constructible_v<ref::Owner>);
static_assert(!std::is_copy_constructible_v<ref::Interval>);
static_assert(!std::is_move_constructible_v<ref::Interval>);
int main() {
  unsupported_profile(); transactional_limits(); signed_enclosures();
  zero_and_positive_underflow(); domain_and_alias_refusals();
  isolated_copy_refusals();
  std::printf("conditional_drag_reference_certificate checks=%u failures=%u "
              "complete_arithmetic_qualified=false whole_payload_qualified=false\n",
              checks, failures);
  return failures ? 1 : 0;
}
