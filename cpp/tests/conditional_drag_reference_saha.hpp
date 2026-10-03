#pragma once
// G1's auxiliary certificate of the original scaled-cubic initial state.
// It neither replaces128 nominal bisections nor supplies a physical late tail.
#include "conditional_drag_reference_constants.hpp"

namespace conditional_drag_reference_certificate {
struct InitialSahaCertificate {
  Interval electron_ratio, hydrogen_fraction, helium_fraction;
  W electron_error, hydrogen_error, helium_error, temperature_ratio_error;
  bool valid;
  InitialSahaCertificate() noexcept : valid(false) {}
  InitialSahaCertificate(const InitialSahaCertificate &) = delete;
  InitialSahaCertificate &operator=(const InitialSahaCertificate &) = delete;
  InitialSahaCertificate(InitialSahaCertificate &&) = delete;
  InitialSahaCertificate &operator=(InitialSahaCertificate &&) = delete;
};

// nH/nHe are independently supplied positive nuclei m^-3; A/B are the
// original positive dimensionless Q*exp/N and4Q*exp/N enclosures, derived
// from that SAME exact N=nH+nHe source. Stored t/p/q/theta borrow the actual
// original cubic midpoint/returned fields. Full source/clock/helper ownership
// and original128-step instrumentation are still required in its caller.
inline bool certify_initial_saha(Arithmetic &a, InitialSahaCertificate &out,
                                 const Interval &nH, const Interval &nHe,
                                 const Interval &A, const Interval &B,
                                 const W &stored_t, const W &stored_p,
                                 const W &stored_q, const W &stored_theta) noexcept {
  Owner &o = a.owner();
  MacroScope certificate(o, Macro::saha);
  if (!certificate) return false;
  const std::array<const void *, 8> borrowed{&nH, &nHe, &A, &B,
                                            &stored_t, &stored_p, &stored_q, &stored_theta};
  const auto output_address = reinterpret_cast<std::uintptr_t>(&out);
  for (unsigned i = 0; i != borrowed.size(); ++i) {
    if (!o.guard() || !o.iteration()) return false;
    const auto input_address = reinterpret_cast<std::uintptr_t>(borrowed[i]);
    const U bytes = i < 4 ? sizeof(Interval) : sizeof(W);
    if (output_address <= input_address ? input_address - output_address < sizeof(out) :
                                         output_address - input_address < bytes)
      return o.fail(Cause::invalid_ownership);
  }
  out.valid = false;
  out.electron_ratio.valid = out.hydrogen_fraction.valid = out.helium_fraction.valid = false;
  Interval hn, hen, an, bn, t, p, q, theta;
  if (!a.copy(hn, nH) || !a.copy(hen, nHe) || !a.copy(an, A) || !a.copy(bn, B) ||
      !a.point(t, stored_t) || !a.point(p, stored_p) || !a.point(q, stored_q) ||
      !a.point(theta, stored_theta) || !o.guard()) return false;
  if (!(hn.lower > 0 && hen.lower > 0 && an.lower > 0 && bn.lower > 0 &&
        stored_t > 0 && stored_t < 1 && stored_p > 0 && stored_p < 1 &&
        stored_q > 0 && stored_q < 1 && stored_theta == 1))
    return o.fail(Cause::invalid_input);
  Interval N, hshare, heshare, aden, bden, ph, phe, hcharge, hecharge, total_charge, F;
  // F(t)=t-(nH/N)A/(A+t)-(nHe/N)B/(B+t). For every exact source
  // member, F'>=1, F(0)<0<F(1); hence |t_true-stored_t|<=|F(stored_t)|.
  // The point is the ACTUAL returned midpoint, not a Newton-size proxy.
  if (!a.add(N, hn, hen) || !a.divide(hshare, hn, N) || !a.divide(heshare, hen, N) ||
      !a.add(aden, an, t) || !a.add(bden, bn, t) ||
      !a.divide(ph, an, aden) || !a.divide(phe, bn, bden) ||
      !a.multiply(hcharge, hshare, ph) || !a.multiply(hecharge, heshare, phe) ||
      !a.add(total_charge, hcharge, hecharge) || !a.subtract(F, t, total_charge) ||
      !a.absolute_upper(out.electron_error, F)) return false;
  Interval error, root_lower, root_upper;
  if (!a.point(error, out.electron_error) || !a.subtract(root_lower, t, error) ||
      !a.add(root_upper, t, error) || !o.guard()) return false;
  if (!(root_lower.lower > 0 && root_upper.upper < 1))
    return o.fail(Cause::unresolved_positive);
  if (!a.enclose(out.electron_ratio, root_lower.lower, root_upper.upper)) return false;
  Interval true_aden, true_bden, pdiff, qdiff;
  if (!a.add(true_aden, an, out.electron_ratio) ||
      !a.add(true_bden, bn, out.electron_ratio) ||
      !a.divide(out.hydrogen_fraction, an, true_aden) ||
      !a.divide(out.helium_fraction, bn, true_bden) || !o.guard()) return false;
  if (!(out.hydrogen_fraction.lower > 0 && out.hydrogen_fraction.upper < 1 &&
        out.helium_fraction.lower > 0 && out.helium_fraction.upper < 1))
    return o.fail(Cause::unresolved_positive);
  if (!a.subtract(pdiff, out.hydrogen_fraction, p) ||
      !a.subtract(qdiff, out.helium_fraction, q) ||
      !a.absolute_upper(out.hydrogen_error, pdiff) ||
      !a.absolute_upper(out.helium_error, qdiff) ||
      !o.action(Action::floating(1, 0, 1))) return false;
  out.temperature_ratio_error = 0; // original literal theta=1 is exact.
  out.valid = true;
  return true;
}
} // namespace conditional_drag_reference_certificate
