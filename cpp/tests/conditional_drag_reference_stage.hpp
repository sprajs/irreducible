#pragma once
// G2's six-variable algebraic inclusion/transport, without a rate oracle.
// The caller still owes actual uniformly enclosed TRUE Radau J, actual-center
// residual/source/clock arithmetic, six original preconditioner solves, and
// original RHS callback accounting. This partial owner earns none of those.
#include "conditional_drag_reference_constants.hpp"

namespace conditional_drag_reference_certificate {
using StageVector = std::array<W, 6>;
using StageMatrix = std::array<Interval, 36>;
struct StageTransport;
class StageSupport {
  std::array<Interval, 6> center_, weight_, box_;
  bool prepared_ = false;
  friend bool prepare_stage_support(Arithmetic &, StageSupport &, const StageVector &) noexcept;
  friend bool certify_stage_transport(Arithmetic &, StageTransport &, const StageSupport &,
      const StageMatrix &, const StageMatrix &, const std::array<Interval, 6> &,
      const std::array<W, 3> &) noexcept;
public:
  // Automatic/runtime storage only. A USER-PROVIDED integer-only constructor
  // also prevents brace-value initialization from zeroing36 unassigned W slots.
  // Static storage's language-mandated initialization is outside this owner.
  StageSupport() noexcept {}
  StageSupport(const StageSupport &) = delete;
  StageSupport &operator=(const StageSupport &) = delete;
  StageSupport(StageSupport &&) = delete;
  StageSupport &operator=(StageSupport &&) = delete;
  bool prepared() const noexcept { return prepared_; }
  const Interval &center(unsigned i) const noexcept { return center_[i]; }
  const Interval &weight(unsigned i) const noexcept { return weight_[i]; }
  const Interval &box(unsigned i) const noexcept { return box_[i]; }
};
struct StageTransport {
  std::array<Interval, 6> radius;
  W contraction, scaled_displacement;
  bool valid;
  StageTransport() noexcept : valid(false) {}
  StageTransport(const StageTransport &) = delete;
  StageTransport &operator=(const StageTransport &) = delete;
  StageTransport(StageTransport &&) = delete;
  StageTransport &operator=(StageTransport &&) = delete;
};
inline bool stage_overlap(Owner &o, const void *out, U out_bytes,
                          const void *in, U in_bytes) noexcept {
  if (!o.guard()) return false;
  const auto x = reinterpret_cast<std::uintptr_t>(out);
  const auto y = reinterpret_cast<std::uintptr_t>(in);
  return (x <= y ? y - x >= out_bytes : x - y >= in_bytes) ||
         o.fail(Cause::invalid_ownership);
}
inline bool prepare_stage_support(Arithmetic &a, StageSupport &out,
                                   const StageVector &stored_v) noexcept {
  Owner &o = a.owner();
  MacroScope certificate(o, Macro::certificate);
  if (!certificate || !stage_overlap(o, &out, sizeof(out), &stored_v, sizeof(stored_v)) ||
      !o.guard()) return false;
  if (out.prepared_) return o.fail(Cause::invalid_input); // ONE box; no enlargement/retry.
  if (!o.arithmetic_ready() || !o.action(Action::floating(2, 0, 2))) return false;
  W one, sixteen; one = 1; sixteen = 16;
  unsigned i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 6) break;
    if (!o.iteration() || !a.point(out.center_[i], stored_v[i]) || !o.guard()) return false;
    if (!(stored_v[i] > 0 && (i % 3 == 2 || stored_v[i] < one)))
      return o.fail(Cause::invalid_input);
    W complement, weight;
    const W *small = &stored_v[i];
    if (!o.guard()) return false;
    if (i % 3 != 2) {
      if (!o.action(Action::floating(1, 1))) return false;
      complement = one - stored_v[i];
      if (!o.guard()) return false;
      if (!std::isnormal(complement) || !(complement > 0))
        return o.fail(Cause::unresolved_positive);
      if (!o.guard()) return false;
      if (complement < stored_v[i]) small = &complement;
    }
    if (!o.action(Action::floating(1, 1))) return false;
    weight = *small / sixteen;
    if (!o.guard()) return false;
    if (!std::isnormal(weight) || !(weight > 0)) return o.fail(Cause::unresolved_positive);
    Interval lo, hi;
    if (!a.point(out.weight_[i], weight) || !a.subtract(lo, out.center_[i], out.weight_[i]) ||
        !a.add(hi, out.center_[i], out.weight_[i]) || !o.guard()) return false;
    if (!(lo.lower > 0 && (i % 3 == 2 || hi.upper < one)))
      return o.fail(Cause::unresolved_positive);
    if (!a.enclose(out.box_[i], lo.lower, hi.upper)) return false;
    ++i;
  }
  out.prepared_ = true;
  return true;
}

// M is a fixed POINT enclosure of an arbitrary retained finite preconditioner;
// true_J encloses the actual Jacobian UNIFORMLY on the one prepared B.
// residual encloses r(v;yhat), excluding incoming radius. E*y=(y,y).
// Q=|I-MJ|, d=|M*r|+|M*E|*r_previous, q=max(sum Qij Dj/Di),
// t=max(di/Di)/(1-q). q<1 and t<=1 give the Banach inclusion on B;
// output is the tight Di*t radius, not the support weight Di.
inline bool certify_stage_transport(Arithmetic &a, StageTransport &out,
      const StageSupport &support, const StageMatrix &M, const StageMatrix &true_J,
      const std::array<Interval, 6> &residual, const std::array<W, 3> &r_previous) noexcept {
  Owner &o = a.owner();
  MacroScope certificate(o, Macro::certificate);
  if (!certificate ||
      !stage_overlap(o, &out, sizeof(out), &support, sizeof(support)) ||
      !stage_overlap(o, &out, sizeof(out), &M, sizeof(M)) ||
      !stage_overlap(o, &out, sizeof(out), &true_J, sizeof(true_J)) ||
      !stage_overlap(o, &out, sizeof(out), &residual, sizeof(residual)) ||
      !stage_overlap(o, &out, sizeof(out), &r_previous, sizeof(r_previous)) || !o.guard())
    return false;
  out.valid = false;
  if (!support.prepared_) return o.fail(Cause::invalid_input);
  SourceNumbers numbers(a);
  Interval zero, one, previous[3];
  if (!numbers.integer(zero, 0) || !numbers.integer(one, 1)) return false;
  unsigned i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 3) break;
    if (!o.iteration() || !a.point(previous[i], r_previous[i]) || !o.guard()) return false;
    if (!(r_previous[i] >= 0)) return o.fail(Cause::invalid_input);
    ++i;
  }
  i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 36) break;
    Interval checked;
    if (!o.iteration() || !a.copy(checked, M[i]) || !o.guard()) return false;
    if (checked.lower != checked.upper) return o.fail(Cause::invalid_input);
    ++i;
  }
  if (!o.action(Action::floating(2, 0, 2))) return false;
  W qmax, dmax; qmax = 0; dmax = 0;
  i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 6) break;
    if (!o.iteration()) return false;
    Interval rsum[2], product, e, incoming[2], qsum[2];
    unsigned current = 0, j = 0;
    if (!a.copy(rsum[0], zero)) return false;
    for (;;) {
      if (!o.guard()) return false;
      if (j == 6) break;
      if (!o.iteration() || !a.multiply(product, M[6 * i + j], residual[j]) ||
          !a.add(rsum[1 - current], rsum[current], product)) return false;
      current = 1 - current; ++j;
    }
    W emax;
    if (!a.absolute_upper(emax, rsum[current]) || !a.point(e, emax) ||
        !a.copy(incoming[0], zero)) return false;
    current = 0; j = 0;
    for (;;) {
      if (!o.guard()) return false;
      if (j == 3) break;
      Interval ME, absME;
      W bound;
      if (!o.iteration() || !a.add(ME, M[6 * i + j], M[6 * i + j + 3]) ||
          !a.absolute_upper(bound, ME) || !a.point(absME, bound) ||
          !a.multiply(product, absME, previous[j]) ||
          !a.add(incoming[1 - current], incoming[current], product)) return false;
      current = 1 - current; ++j;
    }
    Interval d, ratio;
    if (!a.add(d, e, incoming[current]) || !a.divide(ratio, d, support.weight_[i]) ||
        !o.guard() || !o.action(Action::floating(1, 0, 1))) return false;
    dmax = ratio.upper > dmax ? ratio.upper : dmax;
    if (!a.copy(qsum[0], zero)) return false;
    current = 0; j = 0;
    for (;;) {
      if (!o.guard()) return false;
      if (j == 6) break;
      if (!o.iteration()) return false;
      Interval sum[2], defect, abs_defect;
      unsigned k = 0, sum_current = 0;
      if (!a.copy(sum[0], zero)) return false;
      for (;;) {
        if (!o.guard()) return false;
        if (k == 6) break;
        if (!o.iteration() || !a.multiply(product, M[6 * i + k], true_J[6 * k + j]) ||
            !a.add(sum[1 - sum_current], sum[sum_current], product)) return false;
        sum_current = 1 - sum_current; ++k;
      }
      W bound;
      if (!o.guard() || !a.subtract(defect, i == j ? one : zero, sum[sum_current]) ||
          !a.absolute_upper(bound, defect) || !a.point(abs_defect, bound) ||
          !a.multiply(product, abs_defect, support.weight_[j]) ||
          !a.add(qsum[1 - current], qsum[current], product)) return false;
      current = 1 - current; ++j;
    }
    if (!a.divide(ratio, qsum[current], support.weight_[i]) ||
        !o.guard() || !o.action(Action::floating(1, 0, 1))) return false;
    qmax = ratio.upper > qmax ? ratio.upper : qmax;
    ++i;
  }
  if (!o.guard()) return false;
  if (!(qmax >= 0 && qmax < 1 && dmax >= 0)) return o.fail(Cause::unresolved_positive);
  Interval q, max_d, denominator, t;
  if (!a.point(q, qmax) || !a.point(max_d, dmax) || !a.subtract(denominator, one, q) ||
      !a.divide(t, max_d, denominator) || !o.guard()) return false;
  if (!(t.upper >= 0 && t.upper <= 1)) return o.fail(Cause::unresolved_positive);
  i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 6) break;
    if (!o.iteration() || !a.multiply(out.radius[i], support.weight_[i], t)) return false;
    ++i;
  }
  if (!o.action(Action::floating(2, 0, 2))) return false;
  out.contraction = qmax; out.scaled_displacement = t.upper;
  out.valid = true;
  return true;
}
} // namespace conditional_drag_reference_certificate
