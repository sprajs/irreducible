#include "conditional_drag_reference_affine.hpp"
#include <cstdio>
namespace ref = conditional_drag_reference_certificate;
namespace {
unsigned checks = 0, failures = 0;
bool need(bool answer, const char *label) {
  ++checks;
  if (!answer) { ++failures; std::fprintf(stderr, "FAILED %s\n", label); }
  return answer;
}
constexpr ref::Profile assumptions{true, true, true, false, false, false};
struct Inputs {
  ref::Interval hi, lo, H, He, query;
  std::array<ref::Interval, 2> shi, slo;
};
bool source(ref::Owner &o, ref::Arithmetic &a, Inputs &in, bool wide_query) {
  if (!need(o.admit_arithmetic(), "affine conditional arithmetic environment") ||
      !need(a.prepare_directions(), "affine once-owned direction source")) return false;
  ref::SourceNumbers numbers(a);
  const ref::W qlo = 0.5L, qhi = 1.5L;
  return need(numbers.integer(in.hi, 2) && numbers.integer(in.lo, 0) &&
    numbers.integer(in.H, 2) && numbers.integer(in.He, 4) &&
    numbers.rational(in.shi[0], 3, 4) && numbers.rational(in.slo[0], 1, 4) &&
    numbers.rational(in.shi[1], 1, 2) && numbers.rational(in.slo[1], 1, 2) &&
    (wide_query ? a.enclose(in.query, qlo, qhi) : numbers.integer(in.query, 1)),
    "same independent nuclei and exact represented affine source acquired");
}
void polynomial_and_uniform_controls() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  Inputs inputs;
  if (!source(owner, arithmetic, inputs, false)) return;
  ref::AffineCellCertificate result;
  if (!need(ref::certify_affine_charge(arithmetic, result, inputs.hi, inputs.lo,
    inputs.shi, inputs.slo, inputs.H, inputs.He, inputs.query) && result.valid,
    "actual common-support affine certificate acquired")) return;
  // Independent exact polynomial: p(z)=1/4+z/4, q(z)=1/2,
  // charge(z)=2*p+4*q=5/2+z/2, derivative=1/2.
  need(result.hydrogen_fraction.lower <= 0.5L && result.hydrogen_fraction.upper >= 0.5L &&
       result.helium_fraction.lower <= 0.5L && result.helium_fraction.upper >= 0.5L &&
       result.comoving_charge.lower <= 3 && result.comoving_charge.upper >= 3 &&
       result.charge_gradient_upper >= 0.5L && result.charge_gradient_upper < 0.625L &&
       result.separation_lower > 1.5L,
       "certificate encloses exact independent affine polynomial and its derivative");

  ref::Owner uniform_owner(assumptions);
  ref::Arithmetic ua(uniform_owner);
  Inputs uniform;
  if (!source(uniform_owner, ua, uniform, true)) return;
  ref::AffineCellCertificate range;
  if (!need(ref::certify_affine_charge(ua, range, uniform.hi, uniform.lo, uniform.shi,
    uniform.slo, uniform.H, uniform.He, uniform.query) && range.valid,
    "whole query interval is certified on the same common support")) return;
  need(range.hydrogen_fraction.lower <= 0.375L && range.hydrogen_fraction.upper >= 0.625L &&
       range.comoving_charge.lower <= 2.75L && range.comoving_charge.upper >= 3.25L,
       "uniform enclosure covers both exact affine extrema rather than only midpoint");
}
void source_and_support_refusals() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  Inputs inputs;
  if (!source(owner, arithmetic, inputs, false)) return;
  const ref::W clock_lo = 1.9375L, clock_hi = 2.0625L;
  if (!need(arithmetic.enclose(inputs.hi, clock_lo, clock_hi),
            "whole high-endpoint clock enclosure retained")) return;
  ref::AffineCellCertificate result;
  if (!need(ref::certify_affine_charge(arithmetic, result, inputs.hi, inputs.lo, inputs.shi,
    inputs.slo, inputs.H, inputs.He, inputs.query) && result.valid,
    "source clock uncertainty reaches represented affine charge")) return;
  need(result.hydrogen_fraction.lower < 0.5L && result.hydrogen_fraction.upper > 0.5L &&
       result.comoving_charge.lower < 3 && result.comoving_charge.upper > 3,
       "uncertain source coordinate is retained instead of reset to its nominal knot");

  ref::Owner gap_owner(assumptions);
  ref::Arithmetic ga(gap_owner);
  Inputs gap;
  if (!source(gap_owner, ga, gap, false)) return;
  ref::SourceNumbers numbers(ga);
  if (!need(numbers.integer(gap.query, 3), "query outside original cell support")) return;
  ref::AffineCellCertificate refused;
  need(!ref::certify_affine_charge(ga, refused, gap.hi, gap.lo, gap.shi, gap.slo,
    gap.H, gap.He, gap.query) && !refused.valid && gap_owner.failure().cause == ref::Cause::invalid_input,
    "missing common support refuses without clipping or inventing another cell");

  ref::Owner alias_owner(assumptions);
  ref::Arithmetic aa(alias_owner);
  Inputs alias;
  if (!source(alias_owner, aa, alias, false)) return;
  ref::SourceNumbers an(aa);
  ref::AffineCellCertificate acquired;
  if (!need(an.integer(acquired.comoving_charge, 2), "overlapping acquired source charge")) return;
  need(!ref::certify_affine_charge(aa, acquired, alias.hi, alias.lo, alias.shi, alias.slo,
    acquired.comoving_charge, alias.He, alias.query) &&
    alias_owner.failure().cause == ref::Cause::invalid_ownership &&
    acquired.comoving_charge.valid && acquired.comoving_charge.lower == 2,
    "output/input overlap refuses before damaging the owned nuclei source");
}
}
int main() {
  polynomial_and_uniform_controls(); source_and_support_refusals();
  std::printf("conditional_drag_reference_affine checks=%u failures=%u "
              "prefix_caller_qualified=false complete_arithmetic_qualified=false whole_payload_qualified=false\n",
              checks, failures);
  return failures ? 1 : 0;
}
