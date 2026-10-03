#include "conditional_drag_reference_saha.hpp"
#include <cstdio>

namespace ref = conditional_drag_reference_certificate;
namespace {
unsigned checks = 0, failures = 0;
void need(bool answer, const char *label) {
  ++checks;
  if (!answer) { ++failures; std::fprintf(stderr, "FAILED %s\n", label); }
}
constexpr ref::Profile assumptions{true, true, true, false, false, false};
void analytic_neutrality_control() {
  ref::Owner owner(assumptions);
  need(owner.admit_arithmetic(), "Saha conditional arithmetic environment");
  ref::Arithmetic arithmetic(owner);
  need(arithmetic.prepare_directions(), "Saha once-owned direction source");
  ref::SourceNumbers numbers(arithmetic);
  ref::Interval nH, nHe, A, B;
  need(numbers.integer(nH, 1) && numbers.integer(nHe, 1) &&
       numbers.integer(A, 1) && numbers.integer(B, 1), "exact synthetic shared-charge source");
  const ref::W midpoint = 0.625L, returned = 1 / (1 + midpoint), theta = 1;
  ref::InitialSahaCertificate result;
  const bool acquired = ref::certify_initial_saha(arithmetic, result, nH, nHe, A, B,
                                                 midpoint, returned, returned, theta) && result.valid;
  need(acquired,
       "actual-midpoint residual certificate acquired");
  if (!acquired) return; // Wide outputs remain intentionally unassigned on refusal.
  // Exact independent control: for A=B=nH=nHe=1, neutrality is t^2+t-1=0.
  // At632/1024 and633/1024 its numerators are -1984 and305, respectively.
  // The unique root (and both ion fractions) lie between those exact dyadics.
  const ref::W lower = 0.6171875L, upper = 0.6181640625L;
  need(result.electron_ratio.lower <= lower && result.electron_ratio.upper >= upper &&
       result.hydrogen_fraction.lower <= lower && result.hydrogen_fraction.upper >= upper &&
       result.helium_fraction.lower <= lower && result.helium_fraction.upper >= upper,
       "residual-derived intervals cover independent integer-polynomial root limits");
  need(result.electron_error >= midpoint - lower &&
       result.hydrogen_error >= upper - returned && result.helium_error >= upper - returned &&
       result.temperature_ratio_error == 0,
       "returned-field errors cover exact-source root and retain theta1 exact-zero witness");
}
void uncertainty_is_retained() {
  ref::Owner owner(assumptions);
  need(owner.admit_arithmetic(), "uncertain Saha conditional environment");
  ref::Arithmetic arithmetic(owner);
  need(arithmetic.prepare_directions(), "uncertain Saha direction source");
  ref::SourceNumbers numbers(arithmetic);
  ref::Interval nH, nHe, A, B;
  const ref::W lo = 0.999L, hi = 1.001L;
  need(numbers.integer(nH, 1) && numbers.integer(nHe, 1) &&
       arithmetic.enclose(A, lo, hi) && arithmetic.enclose(B, lo, hi),
       "positive whole-source A/B uncertainty supplied to auxiliary certificate");
  const ref::W midpoint = 0.625L, returned = 1 / (1 + midpoint), theta = 1;
  ref::InitialSahaCertificate result;
  need(ref::certify_initial_saha(arithmetic, result, nH, nHe, A, B,
                                midpoint, returned, returned, theta) &&
       result.valid && result.electron_error > 0.01L,
       "source interval uncertainty contributes without a Newton/epsilon proxy");
}
void bad_input_and_alias_controls() {
  ref::Owner owner(assumptions);
  need(owner.admit_arithmetic(), "bad-source Saha conditional environment");
  ref::Arithmetic arithmetic(owner);
  need(arithmetic.prepare_directions(), "bad-source Saha direction source");
  ref::SourceNumbers numbers(arithmetic);
  ref::Interval nH, nHe, A, B;
  need(numbers.integer(nH, 0) && numbers.integer(nHe, 1) &&
       numbers.integer(A, 1) && numbers.integer(B, 1), "zero independent nuclei control");
  const ref::W t = 0.625L, p = 0.625L, theta = 1;
  ref::InitialSahaCertificate refused;
  need(!ref::certify_initial_saha(arithmetic, refused, nH, nHe, A, B, t, p, p, theta) &&
       !refused.valid && owner.failure().cause == ref::Cause::invalid_input,
       "zero independently supplied nuclei refuses without dropping its state");

  ref::Owner alias_owner(assumptions);
  need(alias_owner.admit_arithmetic(), "alias Saha conditional environment");
  ref::Arithmetic alias_arithmetic(alias_owner);
  need(alias_arithmetic.prepare_directions(), "alias Saha direction source");
  ref::SourceNumbers alias_numbers(alias_arithmetic);
  ref::Interval h, he, positive;
  ref::InitialSahaCertificate acquired;
  need(alias_numbers.integer(h, 1) && alias_numbers.integer(he, 1) &&
       alias_numbers.integer(positive, 1) &&
       alias_numbers.integer(acquired.electron_ratio, 1), "acquired alias input retained");
  need(!ref::certify_initial_saha(alias_arithmetic, acquired, h, he,
       acquired.electron_ratio, positive, t, p, p, theta) &&
       alias_owner.failure().cause == ref::Cause::invalid_ownership &&
       acquired.electron_ratio.valid && acquired.electron_ratio.lower == 1,
       "output ownership overlap refuses before altering acquired input");
}
}
int main() {
  analytic_neutrality_control(); uncertainty_is_retained(); bad_input_and_alias_controls();
  std::printf("conditional_drag_reference_saha checks=%u failures=%u "
              "complete_arithmetic_qualified=false whole_payload_qualified=false\n",
              checks, failures);
  return failures ? 1 : 0;
}
