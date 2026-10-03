#include "conditional_drag_reference_constants.hpp"
#include "hydrogen_helium_history_reference.hpp"
#include <type_traits>
#include <cstdio>

namespace ref = conditional_drag_reference_certificate;
namespace {
unsigned checks = 0, failures = 0;
void need(bool answer, const char *label) {
  ++checks;
  if (!answer) { ++failures; std::fprintf(stderr, "FAILED %s\n", label); }
}
constexpr ref::Profile assumptions{true, true, true, false, false, false};
void rational_decimal_controls() {
  ref::Owner owner(assumptions);
  need(owner.admit_arithmetic(), "conditional arithmetic observed");
  ref::Arithmetic arithmetic(owner);
  need(arithmetic.prepare_directions(), "owned direction sources prepared once");
  ref::SourceNumbers numbers(arithmetic);
  ref::Interval half, plus, minus, zero, power, base, mpc;
  need(numbers.rational(half, 1, 2) && half.lower <= 0.5L && half.upper >= 0.5L,
       "rational1/2 contains independent exact dyadic value");
  need(numbers.decimal(plus, 125, -2) && plus.lower <= 1.25L && plus.upper >= 1.25L &&
       numbers.decimal(minus, -125, -2) && minus.lower <= -1.25L && minus.upper >= -1.25L,
       "both signed decimal125/100 contain independent dyadic value");
  need(numbers.decimal(zero, 0, -64) && zero.lower == 0 && zero.upper == 0,
       "exact integer zero remains a structural witness");
  need(numbers.integer(base, -2) && numbers.power(power, base, 5) &&
       power.lower <= -32 && power.upper >= -32,
       "finite signed product DAG encloses exact odd integer power");
  const auto original_mpc = hydrogen_helium_history_reference::mpc();
  need(numbers.mpc(mpc) && mpc.lower <= original_mpc && mpc.upper >= original_mpc,
       "Mpc encloses actual unchanged ancestor AU/pi grouping");
}
void alternating_machin_controls() {
  ref::Owner owner(assumptions);
  need(owner.admit_arithmetic(), "Machin conditional environment");
  ref::Arithmetic arithmetic(owner);
  need(arithmetic.prepare_directions(), "Machin direction source");
  ref::SourceNumbers numbers(arithmetic);
  ref::MacroScope source(owner, ref::Macro::source);
  need(bool(source), "fixed source macro reserved");
  ref::Interval atan5, pi;
  // Independently: q-q^3/3=74/375 >19/100; atan(q)<q=1/5.
  need(numbers.arctan_inverse(atan5, 5) && atan5.lower > 0.19L && atan5.upper < 0.2L,
       "32-term alternating enclosure respects elementary rational limits");
  need(numbers.pi(pi) && pi.lower > 3 && pi.upper < 3.25L,
       "Machin16atan1/5-4atan1/239 respects independent dyadic pi limits");
  need(owner.served()[unsigned(ref::Counter::iterations)] > 64,
       "actual series and signed-product loops charge iterations");
}
void domain_and_first_failure_controls() {
  auto one_destination = ref::policy_caps;
  one_destination[unsigned(ref::Counter::destinations)] = 1;
  ref::Owner one_owner(assumptions, one_destination);
  need(one_owner.admit_arithmetic(), "single-owner limited environment");
  ref::Arithmetic same_owner(one_owner);
  ref::SourceNumbers single(same_owner);
  ref::Interval pending_integer;
  need(!single.integer(pending_integer, 1) && !pending_integer.valid &&
       one_owner.failure().cause == ref::Cause::work_limit &&
       one_owner.failure().served_prefix[unsigned(ref::Counter::destinations)] == 1 &&
       one_owner.failure().served_prefix[unsigned(ref::Counter::copies)] == 0,
       "constant import and retained copies share the same causal one-destination cap");

  ref::Owner zero_denominator(assumptions);
  need(zero_denominator.admit_arithmetic(), "denominator conditional environment");
  ref::Arithmetic a(zero_denominator);
  need(a.prepare_directions(), "denominator direction source");
  ref::SourceNumbers numbers(a);
  ref::Interval out;
  need(!numbers.rational(out, 1, 0) && !out.valid &&
       zero_denominator.failure().cause == ref::Cause::invalid_input,
       "exact zero rational denominator refuses before candidate arithmetic");
  const auto first = zero_denominator.failure();
  need(!numbers.pi(out) && zero_denominator.failure().cause == first.cause &&
       zero_denominator.failure().served_prefix == first.served_prefix,
       "later constant request cannot reset first failure");

  ref::Owner invalid_exponent(assumptions);
  need(invalid_exponent.admit_arithmetic(), "exponent conditional environment");
  ref::Arithmetic b(invalid_exponent);
  need(b.prepare_directions(), "exponent direction source");
  ref::SourceNumbers bounded(b);
  ref::Interval unsupported;
  const auto before = invalid_exponent.served()[unsigned(ref::Counter::destinations)];
  need(!bounded.decimal(unsupported, 1, -65) && !unsupported.valid &&
       invalid_exponent.failure().cause == ref::Cause::invalid_input &&
       invalid_exponent.served()[unsigned(ref::Counter::destinations)] == before,
       "unsupported fixed decimal domain is refused before floating writes");

  auto caps = ref::policy_caps;
  caps[unsigned(ref::Counter::primitives)] = 8;
  ref::Owner work_limited(assumptions, caps);
  need(work_limited.admit_arithmetic(), "limited conditional environment");
  ref::Arithmetic c(work_limited);
  need(c.prepare_directions(), "limited direction source");
  ref::SourceNumbers limited(c);
  ref::Interval pending;
  need(!limited.decimal(pending, 1, -23) && !pending.valid &&
       work_limited.failure().cause == ref::Cause::work_limit &&
       work_limited.failure().served_prefix[unsigned(ref::Counter::destinations)] > 0,
       "mid-constant exhaustion retains served work and withholds output");
}
}
static_assert(!std::is_constructible_v<ref::SourceNumbers, ref::Owner &, ref::Arithmetic &>);
static_assert(std::is_constructible_v<ref::SourceNumbers, ref::Arithmetic &>);
int main() {
  rational_decimal_controls(); alternating_machin_controls();
  domain_and_first_failure_controls();
  std::printf("conditional_drag_reference_constants checks=%u failures=%u "
              "complete_arithmetic_qualified=false whole_payload_qualified=false\n",
              checks, failures);
  return failures ? 1 : 0;
}
