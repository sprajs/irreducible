#pragma once
// G3's exact-real represented HII/HeII affine charge on common cell support.
// The unchanged nominal cell retains (1-f)*hi+f*lo grouping. This auxiliary
// proof uses the algebraically equal hi+f*(lo-hi) expression. It is not H/ne/R,
// a moving-boundary sliver, a Simpson prefix, or an earned continuous history.
#include "conditional_drag_reference_constants.hpp"

namespace conditional_drag_reference_certificate {
struct AffineCellCertificate {
  Interval hydrogen_fraction, helium_fraction, comoving_charge;
  W charge_gradient_upper, separation_lower;
  bool valid;
  AffineCellCertificate() noexcept : valid(false) {}
  AffineCellCertificate(const AffineCellCertificate &) = delete;
  AffineCellCertificate &operator=(const AffineCellCertificate &) = delete;
  AffineCellCertificate(AffineCellCertificate &&) = delete;
  AffineCellCertificate &operator=(AffineCellCertificate &&) = delete;
};
inline bool certify_affine_charge(Arithmetic &a, AffineCellCertificate &out,
    const Interval &z_hi, const Interval &z_lo,
    const std::array<Interval, 2> &state_hi, const std::array<Interval, 2> &state_lo,
    const Interval &nH0, const Interval &nHe0, const Interval &query_z) noexcept {
  Owner &o = a.owner();
  MacroScope cell(o, Macro::cell);
  if (!cell) return false;
  const std::array<const void *, 7> borrowed{
    &z_hi, &z_lo, &state_hi, &state_lo, &nH0, &nHe0, &query_z};
  const std::array<U, 7> bytes{sizeof(z_hi), sizeof(z_lo), sizeof(state_hi),
                             sizeof(state_lo), sizeof(nH0), sizeof(nHe0), sizeof(query_z)};
  const auto output_address = reinterpret_cast<std::uintptr_t>(&out);
  unsigned i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == borrowed.size()) break;
    if (!o.iteration() || !o.guard()) return false;
    const auto input_address = reinterpret_cast<std::uintptr_t>(borrowed[i]);
    if (output_address <= input_address ? input_address - output_address < sizeof(out) :
                                        output_address - input_address < bytes[i])
      return o.fail(Cause::invalid_ownership);
    ++i;
  }
  out.valid = false;
  Interval hi, lo, H, He, query;
  std::array<Interval, 2> shi, slo, differences;
  if (!a.copy(hi, z_hi) || !a.copy(lo, z_lo) || !a.copy(H, nH0) ||
      !a.copy(He, nHe0) || !a.copy(query, query_z) || !o.guard()) return false;
  // Every query member lies in EVERY endpoint member's retained cell.
  // This common-support gate precedes interpolation; gaps are never clipped.
  if (!(lo.lower >= 0 && hi.lower > lo.upper && query.lower >= lo.upper &&
        query.upper <= hi.lower && H.lower > 0 && He.lower > 0))
    return o.fail(Cause::invalid_input);
  i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 2) break;
    if (!o.iteration() || !a.copy(shi[i], state_hi[i]) || !a.copy(slo[i], state_lo[i]) ||
        !o.guard()) return false;
    if (!(shi[i].lower > 0 && shi[i].upper < 1 && slo[i].lower > 0 && slo[i].upper < 1))
      return o.fail(Cause::invalid_input);
    ++i;
  }
  Interval width, numerator, fraction, product;
  if (!a.subtract(width, hi, lo) || !a.subtract(numerator, hi, query) ||
      !a.divide(fraction, numerator, width) || !o.guard()) return false;
  if (!(width.lower > 0)) return o.fail(Cause::unresolved_positive);
  // Interval dependency may widen f beyond its exact [0,1] support. Retain
  // that conservative image; no endpoint or state interval is silently clipped.
  if (!a.subtract(differences[0], slo[0], shi[0]) ||
      !a.multiply(product, fraction, differences[0]) ||
      !a.add(out.hydrogen_fraction, shi[0], product) ||
      !a.subtract(differences[1], slo[1], shi[1]) ||
      !a.multiply(product, fraction, differences[1]) ||
      !a.add(out.helium_fraction, shi[1], product) || !o.guard()) return false;
  if (!(out.hydrogen_fraction.lower > 0 && out.hydrogen_fraction.upper < 1 &&
        out.helium_fraction.lower > 0 && out.helium_fraction.upper < 1))
    return o.fail(Cause::unresolved_positive);
  Interval hydrogen_charge, helium_charge, slope_H, slope_He, slope_sum, slope;
  // Densities are the same independently supplied retained nuclei m^-3 TODAY;
  // this output is comoving charge. u^3 and H/R factors remain in the caller.
  if (!a.multiply(hydrogen_charge, H, out.hydrogen_fraction) ||
      !a.multiply(helium_charge, He, out.helium_fraction) ||
      !a.add(out.comoving_charge, hydrogen_charge, helium_charge) ||
      !a.multiply(slope_H, H, differences[0]) ||
      !a.multiply(slope_He, He, differences[1]) ||
      !a.add(slope_sum, slope_H, slope_He) || !a.divide(slope, slope_sum, width) ||
      !a.absolute_upper(out.charge_gradient_upper, slope) || !o.guard()) return false;
  if (!(out.comoving_charge.lower > 0)) return o.fail(Cause::unresolved_positive);
  if (!o.action(Action::floating(1, 0, 1))) return false;
  out.separation_lower = width.lower;
  out.valid = true;
  return true;
}
} // namespace conditional_drag_reference_certificate
