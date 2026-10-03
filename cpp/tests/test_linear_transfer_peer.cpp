#include "irred/linear_transfer.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace irred::cosmology;
using S = irred::numerics::Status;
using W = long double;
void need(bool b, const char *m) {
  if (!b) {
    std::cerr << "FAIL " << m << '\n';
    std::exit(1);
  }
}
void close(W a, W b, W r, const char *m) {
  if (std::abs(a - b) > r * std::max(1.L, std::abs(b))) {
    std::cerr << (double)a << " vs " << (double)b << '\n';
    need(false, m);
  }
}
W tracer_integral(W s, unsigned panels) {
  W h = s / panels, sum = 0;
  for (unsigned i = 0; i <= panels; ++i) {
    W t = i * h, t2 = t * t;
    W f = t == 0               ? 0
          : std::abs(t) < .01L ? t / 6 - t * t2 / 120 + t * t2 * t2 / 5040
                               : (1 - std::sin(t) / t) / t;
    sum += (i == 0 || i == panels ? 1 : i % 2 ? 4 : 2) * f;
  }
  return sum * h / 3;
}
// Independent a-coordinate Einstein trace equation, unscaled conformal
// velocities and potential derivative: no production momentum/Delta RHS.
using Y = std::array<W, 6>; // phi,phi_eta,dc,theta_c,dr,theta_r
struct Model {
  W h, m, r, l;
};
Y rhs(W a, const Y &y, Model m, W k) {
  W e2 = m.r / (a * a * a * a) + m.m / (a * a * a) + m.l,
    H = a * m.h * std::sqrt(e2) / 299792.458L, fr = m.r / (a * a * a * a * e2),
    fm = m.m / (a * a * a * e2), fl = m.l / e2, L = -1 + fm / 2 + 2 * fl,
    Q = 1 / (a * H);
  return {Q * y[1],
          Q * (-3 * H * y[1] - (2 * L + 1) * H * H * y[0] +
               .5L * H * H * fr * y[4]),
          Q * (-y[3] + 3 * y[1]),
          Q * (-H * y[3] + k * k * y[0]),
          Q * (-4 * y[5] / 3 + 4 * y[1]),
          Q * k * k * (y[4] / 4 + y[0])};
}
Y reference(Model m, W k, W start, W target, W step) {
  W a = start, e2 = m.r / std::pow(a, 4) + m.m / std::pow(a, 3) + m.l,
    H = a * m.h * std::sqrt(e2) / 299792.458L, r = m.m * a / m.r,
    x2 = k * k / (H * H), p0 = -2.L / 3;
  W phi = p0 * (1 - r / 16 - x2 / 30), vc = p0 * (.5L + r / 16 - x2 / 120),
    dv = -p0 * x2 / 24, fr = m.r / std::pow(a, 4) / e2,
    fm = m.m / std::pow(a, 3) / e2, F = fm + 4 * fr / 3,
    delta = (-x2 * phi - 6 * fr * dv) / (1.5L * F), dc = delta - 3 * vc,
    dr = 4 * (delta - 3 * vc) / 3, dphi = -phi + 1.5L * F * vc + 2 * fr * dv;
  Y y{phi, H * dphi, dc, k * k / H * vc, dr, k * k / H * (vc + dv)};
  while (a < target) {
    e2 = m.r / std::pow(a, 4) + m.m / std::pow(a, 3) + m.l;
    H = a * m.h * std::sqrt(e2) / 299792.458L;
    W h = std::min({step * a, step * a * H / k, target - a});
    auto z1 = rhs(a, y, m, k);
    Y q;
    for (unsigned j = 0; j < 6; ++j)
      q[j] = y[j] + h * z1[j] / 2;
    auto z2 = rhs(a + h / 2, q, m, k);
    for (unsigned j = 0; j < 6; ++j)
      q[j] = y[j] + h * z2[j] / 2;
    auto z3 = rhs(a + h / 2, q, m, k);
    for (unsigned j = 0; j < 6; ++j)
      q[j] = y[j] + h * z3[j];
    auto z4 = rhs(a + h, q, m, k);
    for (unsigned j = 0; j < 6; ++j)
      y[j] += h * (z1[j] + 2 * z2[j] + 2 * z3[j] + z4[j]) / 6;
    a += h;
  }
  return y;
}
W top_hat(W x) {
  return std::abs(x) < .01L ? 1 - x * x / 10 + std::pow(x, 4) / 280
                            : 3 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
}
int main() {
  auto pure = prepare_perfect_fluid_transfer(
      prepare_thermal_background({70, 1, 0, 0, 0, {}}), 1e-9);
  const double k = .03;
  for (double a : {.02, .2}) {
    auto result = pure.evaluate(std::span(&k, 1), a, 3);
    need(result.rows[0].metric.status == S::ok &&
             result.rows[0].comoving_cdm.status == S::ok,
         "pure radiation native acceptance");
    W s = k * 299792.458L / 70 * a / std::sqrt(3.L), p0 = -2.L / 3,
      phi = 3 * p0 * (std::sin(s) - s * std::cos(s)) / (s * s * s),
      vc = 3 * p0 * (1 - std::sin(s) / s) / (s * s);
    W F = tracer_integral(s, 16384), refined = tracer_integral(s, 32768),
      delta = 3 * phi - 4.5L * p0 - 9 * p0 * refined + 3 * vc;
    close(F, refined, 1e-14L, "independent tracer integration refinement");
    close(*result.rows[0].metric.value, phi, 3e-6L,
          "radiation analytic potential");
    close(*result.rows[0].comoving_cdm.value, delta, 3e-6L,
          "radiation analytic CDM tracer");
  }
  Model m{70, .3, .000085, 1 - .3 - .000085};
  auto owner = prepare_perfect_fluid_transfer(
      prepare_thermal_background({70, (double)m.r, 0, 0, (double)m.m, {}}),
      1e-10);
  for (double wave : {.001, .01}) {
    auto native = owner.evaluate(std::span(&wave, 1), .1, 3);
    need(native.rows[0].comoving_cdm.status == S::ok &&
             native.rows[0].metric.status == S::ok,
         "mixed radiation matter native acceptance");
    auto coarse = reference(m, wave, 1e-10L, .1L, .002L),
         fine = reference(m, wave, 1e-10L, .1L, .001L);
    W H = .1L * m.h * std::sqrt(m.r / 1e-4L + m.m / .001L + m.l) / 299792.458L,
      dc = coarse[2] + 3 * H * coarse[3] / (wave * wave),
      df = fine[2] + 3 * H * fine[3] / (wave * wave);
    close(dc, df, 2e-7L, "independent trace ODE refinement");
    close(*native.rows[0].comoving_cdm.value, df, 1e-5L,
          "independent trace transfer");
    close(*native.rows[0].metric.value, fine[0], 1e-6L,
          "independent trace potential");
    need(native.rows[0].maximum_constraint_residual < 1e-6,
         "separate constraint gate");
    TransferPolicy tighter;
    tighter.maximum_log_step = .04;
    auto refined_native = owner.evaluate(std::span(&wave, 1), .1, 3, tighter);
    need(refined_native.rows[0].comoving_cdm.status == S::ok,
         "native time refinement accepted");
    close(*refined_native.rows[0].comoving_cdm.value,
          *native.rows[0].comoving_cdm.value, 2e-6L,
          "native time refinement agreement");
    need(refined_native.rows[0].comoving_cdm.time_refinement <
             native.rows[0].comoving_cdm.time_refinement,
         "native time diagnostic improves");
    std::cout << "mixed k=" << wave
              << " transfer=" << *native.rows[0].comoving_cdm.value
              << " reference=" << (double)df << " error="
              << native.rows[0].comoving_cdm.absolute_error_estimate
              << " constraint=" << native.rows[0].maximum_constraint_residual
              << " callbacks=" << native.callbacks << "\n";
  }
  auto dust = prepare_perfect_fluid_transfer(
      prepare_thermal_background({70, 0, 0, 0, 1, {}}), 1e-9);
  PrimordialBand law{2.1e-9, .965, .05, .001, .002};
  BandVariancePolicy p;
  p.base_panels = 16;
  auto band = dust.sigma8_band(law, 1, p);
  need(band.status == S::ok, "band accepted");
  auto integral = [&](unsigned panels) {
    W sum = 0, lo = std::log(W(law.minimum_mpc_inverse)),
      width = std::log(W(law.maximum_mpc_inverse) / law.minimum_mpc_inverse);
    for (unsigned i = 0; i <= panels; ++i) {
      W wave = std::exp(lo + width * i / panels),
        T = .4L * wave * wave * std::pow(299792.458L / 70, 2),
        window = top_hat(wave * 8 / .7L),
        v = law.amplitude *
            std::pow(wave / law.pivot_mpc_inverse, W(law.spectral_index) - 1) *
            T * T * window * window;
      sum += (i == 0 || i == panels ? 1 : i % 2 ? 4 : 2) * v;
    }
    return sum * width / (3 * panels);
  };
  W coarse = integral(4096), fine = integral(8192);
  close(coarse, fine, 1e-15L, "independent band window refinement");
  need(std::abs(*band.variance / fine - 1) < 1e-6L,
       "independent analytic transfer/window amplitude");
  std::cout << "linear_transfer_peer_contract passed analytic "
               "radiation/independent Einstein trace/band integration\n";
}
