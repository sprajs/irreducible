#pragma once
// Original direct-SI stiff comparison. No engine headers, helpers, outputs,
// copied external code or empirical history table. Atomic/rate closure ancestry
// is shared explicitly; polynomial initialization, log-expansion Radau IIA,
// finite-difference Newton and reference refinement are separate algorithms.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>
namespace hydrogen_helium_history_reference {
using W = long double;
using State = std::array<W, 3>;
constexpr W kb = 1.380649e-23L, h = 6.62607015e-34L, c = 299792458,
    ev = 1.602176634e-19L, me = 9.1093837139e-31L,
    G = 6.67430e-11L, sigma = 6.6524587051e-29L,
    chiH = 13.598434599702L * ev, chiHe = 24.587389011L * ev;
struct Model {
  double H0 = 67.4, baryon = .02237, cdm = .12, T0 = 2.7255,
      other = 1.7e-5, nH0 = .19, nHe0 = .015, initial = 2700, late = 300;
  std::vector<std::array<double, 2>> massless_species;
};
// Defining emitted working inputs, not reference expectations or propagated
// physical mapping uncertainty. Fractions are today's Omega, not omega=Omega*h^2.
struct EmittedThermalSource {
  double photon_fraction, baryon_fraction, cdm_fraction, other_massless_fraction;
};
struct Stats { std::size_t steps = 0, rhs = 0, newton = 0; W root_correction = 0; };
W mpc() { return 1e6L * 648000 / std::numbers::pi_v<W> * 149597870700.L; }
W radiation_constant() {
  return 8 * std::pow(std::numbers::pi_v<W>, 5) * std::pow(kb, 4) / (15 * h * h * h * c * c * c);
}
struct Physics {
  Model model;
  W photon, radiation, matter, lambda, H100;
  std::optional<EmittedThermalSource> emitted_source;
  std::optional<W> fixed_photon_energy_today;
  explicit Physics(Model m) : model(m) {
    H100 = 100000 / mpc();
    const W critical = 3 * H100 * H100 / (8 * std::numbers::pi_v<W> * G);
    photon = radiation_constant() * std::pow(W(m.T0), 4) / (c * c * critical);
    radiation = photon + m.other;
    for (const auto &s : m.massless_species)
      radiation += 7.L / 16 * s[1] * radiation_constant() * std::pow(W(s[0]), 4) / (c * c * critical);
    matter = W(m.baryon) + m.cdm;
    lambda = std::pow(W(m.H0) / 100, 2) - radiation - matter;
    if (!(lambda >= 0)) throw std::runtime_error("reference flat closure");
  }
  // Distinct exact-emitted source route. Original Model H0/T0/emitted nH0/nHe0
  // and boundaries remain retained. Physical baryon/CDM/other values are
  // provenance here and are deliberately not remapped into the working law.
  Physics(Model m, const EmittedThermalSource &source) : model(std::move(m)), emitted_source(source) {
    if (!model.massless_species.empty()) throw std::runtime_error("emitted reference requires empty species");
    for (double v : {model.H0, model.T0, model.nH0, model.nHe0, model.initial, model.late,
                     source.photon_fraction, source.baryon_fraction, source.cdm_fraction,
                     source.other_massless_fraction})
      if (!std::isfinite(v)) throw std::runtime_error("nonfinite emitted reference source");
    if (!(model.H0 > 0 && model.T0 > 0 && model.nH0 > 0 && model.nHe0 > 0 &&
          model.initial > model.late && model.late >= 0 && source.photon_fraction > 0 &&
          source.baryon_fraction > 0 && source.cdm_fraction > 0 && source.other_massless_fraction >= 0))
      throw std::runtime_error("emitted reference source domain");
    H100 = 100000 / mpc();
    const W h_ratio = W(model.H0) / 100, h2 = h_ratio * h_ratio,
        critical = 3 * H100 * H100 / (8 * std::numbers::pi_v<W> * G),
        radiation_fraction = W(source.photon_fraction) + source.other_massless_fraction,
        matter_fraction = W(source.baryon_fraction) + source.cdm_fraction,
        lambda_fraction = 1 - radiation_fraction - matter_fraction;
    if (!(lambda_fraction >= 0)) throw std::runtime_error("emitted reference flat closure");
    photon = W(source.photon_fraction) * h2;
    radiation = radiation_fraction * h2; matter = matter_fraction * h2;
    lambda = lambda_fraction * h2;
    fixed_photon_energy_today = photon * critical * c * c;
  }
  W hubble(W u) const { return H100 * std::sqrt(radiation * u * u * u * u + matter * u * u * u + lambda); }
  State initial() const {
    const W u = 1 + W(model.initial), T = W(model.T0) * u,
        Q = std::pow(2 * std::numbers::pi_v<W> * me * kb * T / (h * h), 1.5L),
        nH = W(model.nH0) * u * u * u, nHe = W(model.nHe0) * u * u * u, N = nH + nHe,
        A = Q * std::exp(-chiH / (kb * T)) / N,
        B = 4 * Q * std::exp(-chiHe / (kb * T)) / N;
    // Neutrality expanded as an original scaled cubic in t=ne/(nH+nHe).
    // Polynomial control does not call the production charge-root algorithm.
    auto polynomial = [&](W t) {
      return ((t / (A * B) + 1 / A + 1 / B) * t + 1 - nH / (N * B) - nHe / (N * A)) * t - 1;
    };
    W lo = 0, hi = 1;
    if (!(polynomial(lo) < 0 && polynomial(hi) > 0)) throw std::runtime_error("reference cubic bracket");
    for (unsigned k = 0; k < 128; ++k) {
      W t = (lo + hi) / 2;
      if (polynomial(t) < 0) lo = t; else hi = t;
    }
    const W e = (lo + hi) / 2;
    return {A / (A + e), B / (B + e), 1};
  }
  State operator()(W s, const State &state) const {
    const W u = (1 + W(model.initial)) * std::exp(-s),
        nH = W(model.nH0) * u * u * u, nHe = W(model.nHe0) * u * u * u,
        ne = nH * state[0] + nHe * state[1], H = hubble(u),
        Tr = W(model.T0) * u, T = state[2] * Tr,
        Q = std::pow(2 * std::numbers::pi_v<W> * me * kb * T / (h * h), 1.5L),
        t = T / 10000, alphaH = 1e-19L * 4.309L * std::pow(t, -.6166L) / (1 + .6703L * std::pow(t, .53L)),
        s0 = std::sqrt(T / std::pow(10.L, .477121L)), s1 = std::sqrt(T / std::pow(10.L, 5.114L)),
        alphaHe = std::pow(10.L, -16.744L) / (s0 * std::pow(1 + s0, .289L) * std::pow(1 + s1, 1.711L)),
        EH = h * c / 121.5682e-9L, EHe = h * c * 1.66277434e7L,
        delta = h * c * (1.71134891e7L - 1.66277434e7L),
        betaH = alphaH * Q * std::exp(-(chiH - EH) / (kb * T)),
        betaHe = 4 * alphaHe * Q * std::exp(-(chiHe - EHe) / (kb * T)),
        KH = std::pow(121.5682e-9L, 3) / (8 * std::numbers::pi_v<W> * H),
        KHe = 1 / (8 * std::numbers::pi_v<W> * H * std::pow(1.71134891e7L, 3)),
        neutralH = nH * (1 - state[0]), neutralHe = nHe * (1 - state[1]),
        escapeH = 1 / (KH * neutralH),
        escapeHe = std::exp(-delta / (kb * T)) / (KHe * neutralHe),
        CH = (8.22458L + escapeH) / (8.22458L + escapeH + betaH),
        CHe = (51.3L + escapeHe) / (51.3L + escapeHe + betaHe),
        groundH = alphaH * Q * std::exp(-chiH / (kb * T)),
        groundHe = 4 * alphaHe * Q * std::exp(-chiHe / (kb * T)),
        gamma_over_H = fixed_photon_energy_today
            ? 8 * sigma * *fixed_photon_energy_today * u * u * u * u /
                (3 * me * c * H) * ne / (nH + nHe + ne)
            : 8 * sigma * radiation_constant() * std::pow(Tr, 4) /
                (3 * me * c * H) * ne / (nH + nHe + ne);
    return {-CH * (alphaH * ne * state[0] - groundH * (1 - state[0])) / H,
            -CHe * (alphaHe * ne * state[1] - groundHe * (1 - state[1])) / H,
            -state[2] - gamma_over_H * (state[2] - 1)};
  }
  W opacity(W z, const State &state) const {
    const W u = 1 + z;
    return c * sigma * u * u * (W(model.nH0) * state[0] + W(model.nHe0) * state[1]) / hubble(u);
  }
};
bool positive(const std::array<W, 6> &y) {
  for (unsigned i = 0; i < 2; ++i)
    if (!(y[3 * i] > 0 && y[3 * i] < 1 && y[3 * i + 1] > 0 && y[3 * i + 1] < 1 &&
          y[3 * i + 2] > 0 && std::isfinite(y[3 * i + 2]))) return false;
  return true;
}
using Vector = std::array<W, 6>;
using Matrix = std::array<Vector, 6>;
Vector solve(Matrix a, Vector b) {
  for (unsigned k = 0; k < 6; ++k) {
    unsigned p = k;
    for (unsigned j = k + 1; j < 6; ++j) if (std::abs(a[j][k]) > std::abs(a[p][k])) p = j;
    if (!(std::abs(a[p][k]) > 1e-28L)) throw std::runtime_error("reference pivot");
    std::swap(a[p], a[k]); std::swap(b[p], b[k]);
    for (unsigned j = k + 1; j < 6; ++j) {
      const W ratio = a[j][k] / a[k][k];
      for (unsigned l = k; l < 6; ++l) a[j][l] -= ratio * a[k][l];
      b[j] -= ratio * b[k];
    }
  }
  Vector x{};
  for (int i = 5; i >= 0; --i) {
    W t = b[i]; for (unsigned j = i + 1; j < 6; ++j) t -= a[i][j] * x[j]; x[i] = t / a[i][i];
  }
  return x;
}
W norm(const Vector &d, const Vector &v) {
  W out = 0;
  for (unsigned j = 0; j < 6; ++j) {
    if (!std::isfinite(d[j]) || !(v[j] > 0)) return std::numeric_limits<W>::infinity();
    out = std::max(out, std::abs(d[j]) / v[j]);
  }
  return out;
}
State step(const Physics &f, W s, W ds, State y, Stats &stats) {
  constexpr W A[2][2]{{5.L / 12, -1.L / 12}, {3.L / 4, 1.L / 4}};
  Vector v{y[0], y[1], y[2], y[0], y[1], y[2]};
  auto residual = [&](const Vector &x) {
    State k[]{f(s + ds / 3, {x[0], x[1], x[2]}), f(s + ds, {x[3], x[4], x[5]})};
    stats.rhs += 2; Vector r{};
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 3; ++j) r[3 * i + j] = x[3 * i + j] - y[j] - ds * (A[i][0] * k[0][j] + A[i][1] * k[1][j]);
    return r;
  };
  for (unsigned iteration = 0; iteration < 48; ++iteration) {
    auto r = residual(v); ++stats.newton;
    Matrix J{};
    for (unsigned j = 0; j < 6; ++j) {
      const W scale = j % 3 == 2 ? v[j] : std::min(v[j], 1 - v[j]),
          d = 1e-5L * scale;
      auto a = v, b = v; a[j] += d; b[j] -= d;
      const auto ra = residual(a), rb = residual(b);
      for (unsigned i = 0; i < 6; ++i) J[i][j] = (ra[i] - rb[i]) / (2 * d);
    }
    for (W &a : r) a = -a;
    const auto direction = solve(J, r); const W n = norm(direction, v);
    if (!std::isfinite(n)) throw std::runtime_error("reference nonfinite correction");
    if (n < 2e-16L) {
      stats.root_correction = std::max(stats.root_correction, n); ++stats.steps;
      return {v[3], v[4], v[5]};
    }
    W damping = 1; bool accepted = false;
    for (unsigned k = 0; k < 64; ++k) {
      auto next = v; for (unsigned j = 0; j < 6; ++j) next[j] += damping * direction[j];
      if (positive(next)) {
        const auto rr = residual(next);
        W scaled = 0;
        // Compare residuals in the existing Newton metric; final acceptance
        // recomputes the finite-difference Jacobian and correction above.
        Vector negative = rr; for (W &a : negative) a = -a;
        const auto dd = solve(J, negative); scaled = norm(dd, next);
        if (scaled < n) { v = next; accepted = true; break; }
      }
      damping /= 2;
    }
    if (!accepted) throw std::runtime_error("reference line search");
  }
  throw std::runtime_error("reference Newton ceiling");
}
struct Row { W z, hydrogen, helium, temperature, electrons, opacity; };
struct Result { std::vector<Row> rows; Stats stats; };
Result integrate_physics(Physics f, const std::vector<W> &requested, unsigned refinement) {
  const Model &model = f.model;
  State state = f.initial(); W s = 0;
  std::vector<W> stops = requested; stops.push_back(model.initial); stops.push_back(model.late);
  stops.push_back(W(model.initial) - .1L);
  std::sort(stops.begin(), stops.end(), std::greater<W>());
  stops.erase(std::unique(stops.begin(), stops.end()), stops.end());
  if (stops.front() != model.initial || stops.back() != model.late || !refinement)
    throw std::runtime_error("reference domain");
  Result out;
  auto append = [&](W z) {
    const W u = 1 + z;
    out.rows.push_back({z, state[0], state[1], state[2] * W(model.T0) * u,
        u * u * u * (W(model.nH0) * state[0] + W(model.nHe0) * state[1]), f.opacity(z, state)});
  };
  if (std::find(requested.begin(), requested.end(), stops.front()) != requested.end()) append(stops.front());
  for (std::size_t j = 1; j < stops.size(); ++j) {
    const W end = std::log((1 + W(model.initial)) / (1 + stops[j])), length = end - s;
    std::size_t n = std::max<std::size_t>(1, std::ceil(length / 1e-4L));
    // The explicitly resolved initial Compton thermal boundary, while stiff
    // fractions themselves are advanced by the L-stable implicit method.
    if (stops[j] >= W(model.initial) - .1L) n = std::max<std::size_t>(n, 128);
    n *= refinement; const W ds = length / n, begin = s;
    for (std::size_t k = 0; k < n; ++k) {
      state = step(f, s, ds, state, out.stats); s = begin + (k + 1) * ds;
    }
    if (std::find(requested.begin(), requested.end(), stops[j]) != requested.end()) append(stops[j]);
  }
  return out;
}
// Existing actual consumer retains the physical-remapping entry point and
// sorted unique result convention. Solver, mesh, Newton and source facts are
// unchanged. The exact-emitted cohort has a separately identified entry point.
Result integrate(Model model, const std::vector<W> &requested, unsigned refinement) {
  return integrate_physics(Physics(std::move(model)), requested, refinement);
}
Result restore_original_rows(Result out, const std::vector<W> &requested) {
  // Own both buffers while restoring literal source order and duplicate
  // occurrences. The returned rows have exactly the original requested length;
  // caller row IDs remain attached by original index, never by a sorted join.
  std::vector<Row> restored;
  restored.reserve(requested.size());
  for (W z : requested) {
    const auto row = std::find_if(out.rows.begin(), out.rows.end(),
                                  [&](const Row &r) { return r.z == z; });
    if (row == out.rows.end()) throw std::runtime_error("missing emitted reference row");
    restored.push_back(*row);
  }
  out.rows.swap(restored);
  return out;
}
Result integrate_emitted(Model model, const EmittedThermalSource &source,
                         const std::vector<W> &requested, unsigned refinement) {
  for (W z : requested)
    if (!std::isfinite(z) || z < W(model.late) || z > W(model.initial))
      throw std::runtime_error("emitted reference query domain");
  auto out = integrate_physics(Physics(std::move(model), source), requested, refinement);
  return restore_original_rows(std::move(out), requested);
}
} // namespace hydrogen_helium_history_reference
