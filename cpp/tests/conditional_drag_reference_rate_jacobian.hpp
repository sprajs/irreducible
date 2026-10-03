#pragma once
// Original singlet HII/HeII and Compton-ratio differential algebra (G2).
// These interval INPUTS must enclose the original functions and derivatives
// uniformly at T=theta*Tr on the SAME stage/source/clock support. This owner
// does not earn their power/exp/source ancestry or the original RHS40M charge.
#include "conditional_drag_reference_constants.hpp"

namespace conditional_drag_reference_certificate {
struct IonRateBox {
  Interval alpha, alpha_T, ground, ground_T, beta, beta_T;
  Interval escape, escape_T, lambda;
  // Automatic/runtime storage only. Static language-mandated zeroing is out
  // of scope. These user-provided constructors initialize no floating member.
  IonRateBox() noexcept {}
};
struct RateBox {
  std::array<Interval, 3> state; // p,q,theta; theta=T/Tr.
  Interval nH, nHe, H, Tr, g; // nuclei m^-3, H s^-1, Tr K, g dimensionless.
  std::array<IonRateBox, 2> ion;
  RateBox() noexcept {}
};
struct RateDifferential {
  std::array<Interval, 3> value;
  std::array<Interval, 9> jacobian; // rows f_p,f_q,f_theta; columns p,q,theta.
  bool valid;
  RateDifferential() noexcept : valid(false) {}
  RateDifferential(const RateDifferential &) = delete;
  RateDifferential &operator=(const RateDifferential &) = delete;
  RateDifferential(RateDifferential &&) = delete;
  RateDifferential &operator=(RateDifferential &&) = delete;
};
inline bool certify_uniform_rate_differential(Arithmetic &a, RateDifferential &out,
                                              const RateBox &in) noexcept {
  Owner &o = a.owner();
  MacroScope rate(o, Macro::rate);
  if (!rate || !o.guard()) return false;
  const auto x = reinterpret_cast<std::uintptr_t>(&out), y = reinterpret_cast<std::uintptr_t>(&in);
  if (x <= y ? y - x < sizeof(out) : x - y < sizeof(in))
    return o.fail(Cause::invalid_ownership);
  out.valid = false;
  SourceNumbers numbers(a);
  Interval zero, one, negative, N, nh, nhe, H, Tr, g;
  std::array<Interval, 3> state;
  if (!numbers.integer(zero, 0) || !numbers.integer(one, 1) || !numbers.integer(negative, -1) ||
      !a.copy(nh, in.nH) || !a.copy(nhe, in.nHe) || !a.copy(H, in.H) ||
      !a.copy(Tr, in.Tr) || !a.copy(g, in.g) || !o.guard()) return false;
  if (!(nh.lower > 0 && nhe.lower > 0 && H.lower > 0 && Tr.lower > 0 && g.lower > 0))
    return o.fail(Cause::invalid_input);
  unsigned i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 3) break;
    if (!o.iteration() || !a.copy(state[i], in.state[i]) || !o.guard()) return false;
    if (!(state[i].lower > 0 && (i == 2 || state[i].upper < 1)))
      return o.fail(Cause::invalid_input);
    ++i;
  }
  Interval charge_H, charge_He, ne, denominator, gamma_numerator, gamma;
  if (!a.add(N, nh, nhe) || !a.multiply(charge_H, nh, state[0]) ||
      !a.multiply(charge_He, nhe, state[1]) || !a.add(ne, charge_H, charge_He) ||
      !a.add(denominator, N, ne) || !a.multiply(gamma_numerator, g, ne) ||
      !a.divide(gamma, gamma_numerator, denominator)) return false;
  i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 2) break;
    if (!o.iteration()) return false;
    const auto &r = in.ion[i];
    // Acquire ALL nine rate/derivative images once; signed T derivatives are
    // permitted, while cold ground/beta lower0 retains absolute uncertainty.
    Interval alpha, alpha_T, ground, ground_T, beta, beta_T, escape, escape_T, lambda;
    if (!a.copy(alpha, r.alpha) || !a.copy(alpha_T, r.alpha_T) ||
        !a.copy(ground, r.ground) || !a.copy(ground_T, r.ground_T) ||
        !a.copy(beta, r.beta) || !a.copy(beta_T, r.beta_T) ||
        !a.copy(escape, r.escape) || !a.copy(escape_T, r.escape_T) ||
        !a.copy(lambda, r.lambda) || !o.guard()) return false;
    if (!(alpha.lower > 0 && ground.lower >= 0 && beta.lower >= 0 &&
          escape.lower > 0 && lambda.lower > 0)) return o.fail(Cause::invalid_input);
    Interval neutral, d0, d, d2, C, alpha_ne, forward, reverse, balance;
    Interval beta_escape, neutral_d2, C_x, beta_escape_T, d0_beta_T, CT_numerator, C_T;
    if (!a.subtract(neutral, one, state[i]) || !a.add(d0, lambda, escape) ||
        !a.add(d, d0, beta) || !a.multiply(d2, d, d) || !a.divide(C, d0, d) ||
        !a.multiply(alpha_ne, alpha, ne) || !a.multiply(forward, alpha_ne, state[i]) ||
        !a.multiply(reverse, ground, neutral) || !a.subtract(balance, forward, reverse) ||
        !a.multiply(beta_escape, beta, escape) || !a.multiply(neutral_d2, neutral, d2) ||
        !a.divide(C_x, beta_escape, neutral_d2) || !a.multiply(beta_escape_T, beta, escape_T) ||
        !a.multiply(d0_beta_T, d0, beta_T) || !a.subtract(CT_numerator, beta_escape_T, d0_beta_T) ||
        !a.divide(C_T, CT_numerator, d2)) return false;
    Interval Cr, negative_Cr;
    if (!a.multiply(Cr, C, balance) || !a.multiply(negative_Cr, negative, Cr) ||
        !a.divide(out.value[i], negative_Cr, H)) return false;
    unsigned j = 0;
    for (;;) {
      if (!o.guard()) return false;
      if (j == 2) break;
      if (!o.iteration()) return false;
      const Interval &nj = j == 0 ? nh : nhe;
      Interval alpha_nj, first, diagonal, dr, Cj_r, C_dr, total, minus_total;
      if (!a.multiply(alpha_nj, alpha, nj) || !a.multiply(first, alpha_nj, state[i]) ||
          !o.guard()) return false;
      if (j == i) {
        if (!a.add(diagonal, alpha_ne, ground) || !a.add(dr, first, diagonal) ||
            !a.multiply(Cj_r, C_x, balance)) return false;
      } else if (!a.copy(dr, first) || !a.copy(Cj_r, zero)) return false;
      if (!a.multiply(C_dr, C, dr) || !a.add(total, Cj_r, C_dr) ||
          !a.multiply(minus_total, negative, total) ||
          !a.divide(out.jacobian[3 * i + j], minus_total, H)) return false;
      ++j;
    }
    // T derivatives become theta derivatives through EXACT dT/dtheta=Tr.
    Interval alphaT_ne, forward_T, reverse_T, dr_T, CT_r, C_drT, total_T, theta_T, minus_T;
    if (!a.multiply(alphaT_ne, alpha_T, ne) || !a.multiply(forward_T, alphaT_ne, state[i]) ||
        !a.multiply(reverse_T, ground_T, neutral) || !a.subtract(dr_T, forward_T, reverse_T) ||
        !a.multiply(CT_r, C_T, balance) || !a.multiply(C_drT, C, dr_T) ||
        !a.add(total_T, CT_r, C_drT) || !a.multiply(theta_T, Tr, total_T) ||
        !a.multiply(minus_T, negative, theta_T) || !a.divide(out.jacobian[3 * i + 2], minus_T, H))
      return false;
    ++i;
  }
  Interval theta_minus_one, gamma_change, theta_total;
  if (!a.subtract(theta_minus_one, state[2], one) ||
      !a.multiply(gamma_change, gamma, theta_minus_one) ||
      !a.add(theta_total, state[2], gamma_change) ||
      !a.multiply(out.value[2], negative, theta_total)) return false;
  Interval denominator2, gN, theta_gN, prefactor;
  if (!a.multiply(denominator2, denominator, denominator) || !a.multiply(gN, g, N) ||
      !a.multiply(theta_gN, theta_minus_one, gN) || !a.divide(prefactor, theta_gN, denominator2))
    return false;
  i = 0;
  for (;;) {
    if (!o.guard()) return false;
    if (i == 2) break;
    if (!o.iteration()) return false;
    const Interval &ni = i == 0 ? nh : nhe;
    Interval derivative;
    if (!a.multiply(derivative, prefactor, ni) ||
        !a.multiply(out.jacobian[6 + i], negative, derivative)) return false;
    ++i;
  }
  Interval one_gamma;
  if (!a.add(one_gamma, one, gamma) || !a.multiply(out.jacobian[8], negative, one_gamma))
    return false;
  out.valid = true;
  return true;
}
} // namespace conditional_drag_reference_certificate
