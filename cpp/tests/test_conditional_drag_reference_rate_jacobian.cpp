#include "conditional_drag_reference_rate_jacobian.hpp"
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
bool fixture(ref::Owner &o, ref::Arithmetic &a, ref::RateBox &box, bool uniform) {
  if (!need(o.admit_arithmetic(), "rate differential conditional arithmetic environment") ||
      !need(a.prepare_directions(), "rate once-owned direction source")) return false;
  ref::SourceNumbers n(a);
  if (!(n.rational(box.state[0], 1, 2) && n.rational(box.state[1], 1, 2) &&
      n.integer(box.state[2], 1) && n.integer(box.nH, 1) && n.integer(box.nHe, 1) &&
      n.integer(box.H, 1) && n.integer(box.Tr, 2) && n.integer(box.g, 1))) return false;
  for (unsigned i = 0; i != 2; ++i) {
    auto &r = box.ion[i];
    if (!(n.integer(r.alpha, 1) && n.integer(r.alpha_T, 0) && n.rational(r.ground, 1, 2) &&
        n.integer(r.ground_T, 0) && n.integer(r.beta, 1) && n.integer(r.beta_T, 0) &&
        n.integer(r.escape, 1) && n.integer(r.escape_T, 0) && n.integer(r.lambda, 1)))
      return false;
  }
  if (uniform) {
    const ref::W p_lo = 0.25L, p_hi = 0.75L;
    ref::Interval two_thirds, two;
    if (!a.enclose(box.state[0], p_lo, p_hi) || !n.rational(two_thirds, 2, 3) ||
        !n.integer(two, 2) || !a.enclose(box.ion[0].escape, two_thirds.lower, two.upper)) return false;
    // Whole escape law e_H=1/[2(1-p)] is covered, not a constant midpoint
    // escape pretending to have e_x=e/(1-p) on a nontrivial box.
  }
  return true;
}
void original_differential_algebra_control() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  ref::RateBox box{};
  if (!need(fixture(owner, arithmetic, box, false), "exact shared-ne synthetic rate fixture")) return;
  ref::RateDifferential result;
  if (!need(ref::certify_uniform_rate_differential(arithmetic, result, box) && result.valid,
            "actual original rate differential assembly acquired")) return;
  // Independent symbolic values at p=q=1/2, theta=1:
  // ne=1,N=2,C=2/3,r=1/4,C_x=2/9, f_p=f_q=-1/6;
  // own chemical derivative=-25/18, cross=-1/3; theta derivatives0.
  // Temperature f=-1, derivatives(0,0,-4/3). Dyadic brackets below
  // distinguish the actual C_x*r term from the incorrect -4/3 own derivative.
  need(result.value[0].lower > -0.16796875L && result.value[0].upper < -0.166015625L &&
       result.value[1].lower > -0.16796875L && result.value[1].upper < -0.166015625L &&
       result.jacobian[0].lower > -1.390625L && result.jacobian[0].upper < -1.38671875L &&
       result.jacobian[4].lower > -1.390625L && result.jacobian[4].upper < -1.38671875L &&
       result.jacobian[1].lower > -0.3359375L && result.jacobian[1].upper < -0.33203125L &&
       result.jacobian[3].lower > -0.3359375L && result.jacobian[3].upper < -0.33203125L,
       "true shared-charge/C derivative agrees with independent exact rational algebra");
  need(result.jacobian[2].lower == 0 && result.jacobian[2].upper == 0 &&
       result.jacobian[5].lower == 0 && result.jacobian[5].upper == 0 &&
       result.jacobian[6].lower == 0 && result.jacobian[7].upper == 0 &&
       result.value[2].lower <= -1 && result.value[2].upper >= -1 &&
       result.jacobian[8].lower > -1.3359375L && result.jacobian[8].upper < -1.33203125L,
       "constant thermal-rate limit retains exact-zero derivatives and original adiabatic/Compton term");
}
void whole_support_and_refusal_controls() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  ref::RateBox box{};
  if (!need(fixture(owner, arithmetic, box, true), "whole fraction/escape support fixture")) return;
  ref::RateDifferential result;
  if (!need(ref::certify_uniform_rate_differential(arithmetic, result, box) && result.valid,
            "whole same-support rate differential admitted")) return;
  // Exact-real p endpoints1/4,3/4: C=(3/2-p)/(5/2-2p),
  // r=p^2+p-1/2, so f_H is15/128 and-39/64 respectively.
  need(result.value[0].lower <= -0.609375L && result.value[0].upper >= 0.1171875L,
       "uniform source/escape image covers both exact extrema rather than central rate only");

  ref::Owner bad_owner(assumptions);
  ref::Arithmetic ba(bad_owner);
  ref::RateBox bad{};
  if (!need(fixture(bad_owner, ba, bad, false), "invalid independent nuclei fixture")) return;
  ref::SourceNumbers numbers(ba);
  if (!need(numbers.integer(bad.nH, 0), "actual zero supplied hydrogen source")) return;
  ref::RateDifferential refused;
  need(!ref::certify_uniform_rate_differential(ba, refused, bad) && !refused.valid &&
       bad_owner.failure().cause == ref::Cause::invalid_input,
       "missing independent positive nuclei refuses before pretending to have shared-ne derivatives");

  ref::Owner unsupported;
  ref::Arithmetic ua(unsupported);
  ref::RateBox unavailable{};
  ref::RateDifferential no_image;
  need(!ref::certify_uniform_rate_differential(ua, no_image, unavailable) && !no_image.valid &&
       unsupported.failure().cause == ref::Cause::unsupported_profile &&
       unsupported.served()[unsigned(ref::Counter::destinations)] == 0,
       "missing arithmetic profile refuses before acquired rate values or candidate writes");
}
void nonzero_T_chain_and_Compton_response() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  ref::RateBox box{};
  if (!need(fixture(owner, arithmetic, box, false), "nonzero-T synthetic point fixture")) return;
  ref::SourceNumbers numbers(arithmetic);
  for (unsigned i = 0; i != 2; ++i) {
    if (!need(numbers.integer(box.ion[i].alpha_T, 1) &&
              numbers.rational(box.ion[i].ground_T, 1, 2) &&
              numbers.integer(box.ion[i].beta_T, 1) &&
              numbers.integer(box.ion[i].escape_T, i == 0 ? 0 : 1),
              "actual nonzero thermal derivative images acquired")) return;
  }
  ref::RateDifferential derivative;
  if (!need(ref::certify_uniform_rate_differential(arithmetic, derivative, box) && derivative.valid,
            "actual nonzero-T chain assembled")) return;
  // Independent symbolic point algebra: dr/dT=1/4; C_H,T=-2/9,
  // C_He,T=-1/9, r=1/4, C=2/3, Tr=2. Thus df_H/dtheta=-2/9
  // and df_He/dtheta=-5/18. Omitting C_T gives-1/3; omitting Tr
  // halves both correct derivatives and fails these independent dyadic brackets.
  need(derivative.jacobian[2].lower > -0.224609375L && derivative.jacobian[2].upper < -0.220703125L &&
       derivative.jacobian[5].lower > -0.279296875L && derivative.jacobian[5].upper < -0.27734375L,
       "C_T and nonzero ground/alpha derivatives propagate through the actual Tr chain");

  ref::Owner thermal_owner(assumptions);
  ref::Arithmetic ta(thermal_owner);
  ref::RateBox thermal{};
  if (!need(fixture(thermal_owner, ta, thermal, false), "nonzero Compton charge-response fixture")) return;
  ref::SourceNumbers tn(ta);
  if (!need(tn.integer(thermal.state[2], 2) && tn.integer(thermal.nH, 2),
            "theta2 and distinct independent nuclei2:1 acquired")) return;
  ref::RateDifferential response;
  if (!need(ref::certify_uniform_rate_differential(ta, response, thermal) && response.valid,
            "nonzero same-source Compton response acquired")) return;
  // Exact N=3,ne=3/2,gamma=1/3 gives f_theta=-7/3,
  // derivatives(-8/27,-4/27,-4/3). Charge-role duplication or resetting
  // theta to1 cannot pass these strictly separated dyadic bounds.
  need(response.value[2].lower > -2.3359375L && response.value[2].upper < -2.33203125L &&
       response.jacobian[6].lower > -0.296875L && response.jacobian[6].upper < -0.294921875L &&
       response.jacobian[7].lower > -0.1484375L && response.jacobian[7].upper < -0.146484375L &&
       response.jacobian[8].lower > -1.3359375L && response.jacobian[8].upper < -1.33203125L,
       "Compton response retains nonzero theta shift and each independent nuclei derivative");
}
}
int main() {
  original_differential_algebra_control(); whole_support_and_refusal_controls();
  nonzero_T_chain_and_Compton_response();
  std::printf("conditional_drag_reference_rate_jacobian checks=%u failures=%u "
              "rate_inputs_qualified=false complete_arithmetic_qualified=false whole_payload_qualified=false\n",
              checks, failures);
  return failures ? 1 : 0;
}
