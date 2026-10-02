#pragma once
// Original direct-SI, resolved fixed-step RK4 reference. No engine headers or
// production equations are read back. Source/reference contract frozen first.
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>
namespace recombination_reference {
using W = long double;
constexpr W kb = 1.380649e-23L, planck = 6.62607015e-34L, c = 299792458,
            ev = 1.602176634e-19L, me = 9.1093837139e-31L,
            mp = 1.67262192595e-27L, G = 6.67430e-11L,
            sigmaT = 6.6524587051e-29L, lya = 121.5682e-9L,
            B1 = 13.598434599702L * ev;
struct Model {
  W H0 = 67.4L, b = .02237L, cdm = .12L, T0 = 2.7255L, other = 1.7e-5L;
};
struct Coeff {
  W A, B, D, E, drag;
};
struct Ref {
  W z, x, tau;
};
struct Result {
  std::vector<Ref> rows;
  W root = 0, stability = 0;
};
W mpc() { return 1e6L * 648000 / std::numbers::pi_v<W> * 149597870700.L; }
W n0(const Model &m) {
  return 3 * std::pow(100000 / mpc(), 2) / (8 * std::numbers::pi_v<W> * G) *
         m.b / (mp + me - B1 / c / c);
}
W photon_physical(const Model &m) {
  W rho = std::numbers::pi_v<W> * std::numbers::pi_v<W> / 15 *
          std::pow(kb * m.T0, 4) /
          std::pow(planck * c / (2 * std::numbers::pi_v<W>), 3) / c / c;
  return rho /
         (3 * std::pow(100000 / mpc(), 2) / (8 * std::numbers::pi_v<W> * G));
}
Coeff coefficients(const Model &m, W z) {
  W u = 1 + z, T = m.T0 * u, n = n0(m) * u * u * u, h = m.H0 / 100,
    gamma = photon_physical(m), orad = (gamma + m.other) / (h * h),
    om = (m.b + m.cdm) / (h * h),
    rate = m.H0 * 1000 / mpc() *
           std::sqrt(orad * u * u * u * u + om * u * u * u + 1 - orad - om),
    alpha = 1e-19L * 4.309L * std::pow(T / 10000, -.6166L) /
            (1 + .6703L * std::pow(T / 10000, .53L)),
    excitation = planck * c / lya,
    beta = alpha *
           std::pow(2 * std::numbers::pi_v<W> * me * kb * T / (planck * planck),
                    1.5L) *
           std::exp(-(B1 - excitation) / (kb * T)),
    K = lya * lya * lya / (8 * std::numbers::pi_v<W> * rate),
    R = 3 * m.b / (4 * gamma * u);
  return {K * 8.22458L * n, K * (8.22458L + beta) * n, n * alpha / (rate * u),
          beta * std::exp(-excitation / (kb * T)) / (rate * u),
          c * sigmaT * n / (rate * u * R)};
}
W saha(const Model &m, W z) {
  W T = m.T0 * (1 + z),
    S = std::pow(2 * std::numbers::pi_v<W> * me * kb * T / (planck * planck),
                 1.5L) *
        std::exp(-B1 / (kb * T)) / (n0(m) * std::pow(1 + z, 3));
  return 2 / (1 + std::sqrt(1 + 4 / S));
}
std::array<W, 2> rhs(const Coeff &q, const std::array<W, 2> &y) {
  W neutral = 1 - y[0];
  return {(1 + q.A * neutral) / (1 + q.B * neutral) *
              (q.D * y[0] * y[0] - q.E * neutral),
          q.drag * y[0]};
}
W derivative(const Coeff &q, W x) {
  W t = 1 - x, den = 1 + q.B * t, cor = (1 + q.A * t) / den;
  return q.D * (2 * x * cor + x * x * (q.B - q.A) / (den * den)) +
         q.E * (1 + 2 * q.A * t + q.A * q.B * t * t) / (den * den);
}
Result integrate(Model m, W initial, W late, W step) {
  Result out;
  unsigned N = std::llround((initial - late) / step);
  step = (initial - late) / N;
  std::array<W, 2> y{saha(m, initial), 0};
  out.rows.push_back({initial, y[0], 0});
  auto add = [](auto y, auto k, W h) {
    for (unsigned j = 0; j < 2; ++j)
      y[j] += h * k[j];
    return y;
  };
  for (unsigned i = 0; i < N; ++i) {
    W z = initial - i * step;
    auto q1 = coefficients(m, z), q2 = coefficients(m, z - step / 2),
         q4 = coefficients(m, z - step);
    out.stability = std::max(out.stability, step * derivative(q1, y[0]));
    auto k1 = rhs(q1, y), a = add(y, k1, -step / 2), k2 = rhs(q2, a),
         b = add(y, k2, -step / 2), k3 = rhs(q2, b), d = add(y, k3, -step),
         k4 = rhs(q4, d);
    for (unsigned j = 0; j < 2; ++j)
      y[j] -= step * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
    if (!(y[0] > 0 && y[0] < 1 && a[0] > 0 && a[0] < 1 && b[0] > 0 &&
          b[0] < 1 && d[0] > 0 && d[0] < 1))
      throw std::runtime_error("RK stage positivity failure");
    out.rows.push_back({z - step, y[0], y[1]});
  }
  W end = out.rows.back().tau;
  for (auto &r : out.rows)
    r.tau -= end;
  for (unsigned i = 0; i < N; ++i)
    if (out.rows[i].tau >= 1 && out.rows[i + 1].tau <= 1) {
      auto a = out.rows[i], b = out.rows[i + 1];
      out.root = b.z + (1 - b.tau) * (a.z - b.z) / (a.tau - b.tau);
      break;
    }
  return out;
}
Ref sample(const Result &r, W z) {
  W step = r.rows[0].z - r.rows[1].z;
  return r.rows[std::llround((r.rows[0].z - z) / step)];
}
} // namespace recombination_reference
