#pragma once
// Original explicit-FD direct-SI H provider for the retained reference thermal
// equations. No production kernels, readbacks, external solver or H table.
#include "hydrogen_thermal_reference.hpp"
namespace hydrogen_relic_reference {
namespace base = hydrogen_thermal_reference;
using base::State;
using base::W;
struct Species {
  W mass, temperature, weight;
};
struct Model {
  W H0 = 67.4, baryon = .02237, cdm = .12, T0 = 2.7255, other = 0;
  std::vector<Species> species;
};
template <unsigned N> struct MomentumRule {
  // Each original panel is bisected after the preserved GL16/32 y=.01
  // refinement failure. Orders, support and the tail envelope are unchanged.
  std::array<W, 14 * N> q{}, factor{};
  MomentumRule() {
    std::array<W, N> x{}, w{};
    for (unsigned i = 0; i < N; ++i) {
      W z = std::cos(std::numbers::pi_v<W> * (i + .75L) / (N + .5L));
      auto legendre = [](W v) {
        W a = 1, b = v;
        for (unsigned k = 2; k <= N; ++k) {
          W next = ((2 * k - 1) * v * b - (k - 1) * a) / k;
          a = b;
          b = next;
        }
        return std::array<W, 2>{b, N * (v * b - a) / (v * v - 1)};
      };
      bool resolved = false;
      for (unsigned j = 0; j < 64; ++j) {
        auto p = legendre(z);
        W delta = p[0] / p[1];
        z -= delta;
        if (std::abs(delta) < 4 * std::numeric_limits<W>::epsilon()) {
          resolved = true;
          break;
        }
      }
      if (!resolved)
        throw std::runtime_error("reference Legendre root ceiling");
      W derivative = legendre(z)[1];
      x[i] = z;
      w[i] = 2 / ((1 - z * z) * derivative * derivative);
    }
    constexpr W edge[]{0, 1, 2, 4, 8, 16, 32, 64};
    for (unsigned panel = 0; panel < 7; ++panel)
      for (unsigned part = 0; part < 2; ++part)
        for (unsigned i = 0; i < N; ++i) {
          W half = (edge[panel + 1] - edge[panel]) / 4;
          unsigned j = (2 * panel + part) * N + i;
          q[j] = edge[panel] + (2 * part + 1) * half + half * x[i];
          factor[j] = half * w[i] * q[j] * q[j] / (std::exp(q[j]) + 1);
        }
  }
  W moment(W y) const {
    W sum = 0;
    for (unsigned i = 0; i < q.size(); ++i)
      sum += factor[i] * std::sqrt(q[i] * q[i] + y * y);
    return sum;
  }
};
W tail(W y) {
  constexpr W L = 64;
  return std::exp(-L) *
         (L * L * L + 3 * L * L + 6 * L + 6 + y * (L * L + 2 * L + 2));
}
struct Background {
  Model model;
  MomentumRule<16> coarse;
  MomentumRule<32> fine;
  unsigned order;
  W H100, critical_energy, conversion, photon, lambda;
  mutable std::size_t h_queries = 0, momentum_terms = 0;
  Background(Model source, unsigned n) : model(std::move(source)), order(n) {
    if (order != 16 && order != 32)
      throw std::runtime_error("reference momentum order");
    H100 = 100000 / base::mpc();
    critical_energy = 3 * H100 * H100 * base::c * base::c /
                      (8 * std::numbers::pi_v<W> * base::G);
    W hbar = base::h / (2 * std::numbers::pi_v<W>);
    conversion = std::pow(base::ev, 4) / std::pow(hbar * base::c, 3);
    photon =
        base::radiation_constant() * std::pow(model.T0, 4) / critical_energy;
    W today = 0;
    for (auto a : model.species)
      today += omega(a, 1);
    lambda = std::pow(model.H0 / 100, 2) - photon - model.other - model.baryon -
             model.cdm - today;
    if (!(lambda > 0))
      throw std::runtime_error("reference massive flat closure");
  }
  W moment(W y) const {
    momentum_terms += 14 * order;
    return order == 16 ? coarse.moment(y) : fine.moment(y);
  }
  W omega(const Species &a, W u) const {
    W temperature = base::kb * a.temperature * u / base::ev;
    W integral = a.mass == 0 ? 7 * std::pow(std::numbers::pi_v<W>, 4) / 120
                             : moment(a.mass / temperature);
    return a.weight /
           (2 * std::numbers::pi_v<W> *
            std::numbers::pi_v<W>)*std::pow(temperature, 4) *
           integral * conversion / critical_energy;
  }
  W hubble(W z) const {
    ++h_queries;
    W u = 1 + z;
    W physical = (photon + model.other) * u * u * u * u +
                 (model.baryon + model.cdm) * u * u * u + lambda;
    for (auto a : model.species)
      physical += omega(a, u);
    if (!(physical > 0))
      throw std::runtime_error("reference positive H density");
    return H100 * std::sqrt(physical);
  }
};
struct Physics {
  Model model;
  W initial, late;
  bool evolved;
  base::Physics hydrogen;
  Background background;
  Physics(Model m, W zi, W zl, bool evolve, unsigned order)
      : model(m), initial(zi), late(zl), evolved(evolve),
        hydrogen({m.H0, m.baryon, m.cdm, m.T0, m.other, {}}, zi, zl),
        background(std::move(m), order) {}
  W hubble(W z) const { return background.hubble(z); }
  W density(W z) const { return hydrogen.density(z); }
  W gamma_ratio(W z, W x) const {
    return 8 * base::sigma * base::radiation_constant() *
           std::pow(model.T0 * (1 + z), 4) /
           (3 * base::me * base::c * hubble(z)) * x / (1 + x);
  }
  W initial_x() const { return hydrogen.initial_x(); }
  // Retained original reference equations; only supplied H changes.
  State rhs(W s, const State &y, W H) const {
    W u = (1 + initial) * std::exp(-s), z = u - 1, T = y[1] * model.T0 * u,
      n = density(z);
    W t = T / 10000, alpha = 1e-19L * 4.309L * std::pow(t, -.6166L) /
                             (1 + .6703L * std::pow(t, .53L));
    W beta = alpha *
             std::pow(2 * std::numbers::pi_v<W> * base::me * base::kb * T /
                          (base::h * base::h),
                      1.5L) *
             std::exp(-(base::binding - base::excitation) / (base::kb * T));
    W K = base::lya * base::lya * base::lya / (8 * std::numbers::pi_v<W> * H),
      neutral = 1 - y[0];
    W inhibition = (1 + K * 8.22458L * n * neutral) /
                   (1 + K * (8.22458L + beta) * n * neutral);
    W A = 8 * base::sigma * base::radiation_constant() *
          std::pow(model.T0 * u, 4) / (3 * base::me * base::c * H) * y[0] /
          (1 + y[0]);
    return {
        -inhibition *
            (n * alpha * y[0] * y[0] -
             beta * neutral * std::exp(-base::excitation / (base::kb * T))) /
            H,
        evolved ? -y[1] - A * (y[1] - 1) : 0};
  }
  State operator()(W s, const State &y) const {
    return rhs(s, y, hubble((1 + initial) * std::exp(-s) - 1));
  }
  std::array<W, 2> opacity_at(W s, W x, W H) const {
    W u = (1 + initial) * std::exp(-s),
      q = base::c * base::sigma * density(u - 1) * x / H;
    return {q, q / (3 * model.baryon / (4 * background.photon * u))};
  }
  std::array<W, 2> opacity_s(W s, W x) const {
    return opacity_at(s, x, hubble((1 + initial) * std::exp(-s) - 1));
  }
};
struct Result {
  Physics physics;
  base::Stats stats;
  std::vector<base::Node> nodes;
  std::vector<std::size_t> queries;
  W root = 0;
  base::Value query(std::size_t i) const {
    auto n = nodes.at(queries.at(i));
    W tau = nodes.back().thomson - n.thomson, drag = nodes.back().drag - n.drag;
    auto q = physics.opacity_s(n.s, n.state[0]);
    W survival = std::exp(-tau), opacity = q[0] / (1 + n.z);
    return {n.z,
            n.state[0],
            n.state[1] * physics.model.T0 * (1 + n.z),
            tau,
            drag,
            opacity,
            opacity * survival,
            survival};
  }
};
Result integrate(Model model, W initial, W late,
                 const std::vector<W> &requested, unsigned refinement,
                 unsigned momentum_order, bool evolved = true) {
  Result out{Physics(std::move(model), initial, late, evolved, momentum_order),
             {},
             {},
             {},
             0};
  std::vector<W> stops = requested;
  stops.push_back(initial);
  stops.push_back(late);
  stops.push_back(initial - .1L);
  std::sort(stops.begin(), stops.end(), std::greater<W>());
  stops.erase(std::unique(stops.begin(), stops.end()), stops.end());
  if (stops.front() != initial || stops.back() != late || !refinement)
    throw std::runtime_error("reference query domain");
  State state{out.physics.initial_x(), 1};
  W s = 0, at = 0, ad = 0;
  auto append = [&](W z) {
    auto d = out.physics(s, state);
    ++out.stats.rhs;
    out.nodes.push_back({s, z, state, d, at, ad});
  };
  append(initial);
  out.stats.initial_stiffness =
      evolved ? 1 + out.physics.gamma_ratio(initial, 1) : 1;
  for (std::size_t j = 1; j < stops.size(); ++j) {
    W end = std::log((1 + initial) / (1 + stops[j])), length = end - s;
    std::size_t n = std::max<std::size_t>(1, std::ceil(length / 4e-4L));
    if (stops[j] >= initial - .1L)
      n = std::max<std::size_t>(
          n, std::ceil(length * out.stats.initial_stiffness * 8));
    n *= refinement;
    W begin = s, ds = length / n;
    for (std::size_t k = 0; k < n; ++k) {
      if (stops[j] >= initial - .1L)
        out.stats.maximum_initial_step_stiffness =
            std::max(out.stats.maximum_initial_step_stiffness,
                     ds * out.stats.initial_stiffness);
      const W s1 = s + ds / 3, s2 = s + ds;
      // Freeze only the two exact stage coefficients during Newton. No H(z)
      // interpolation or provider table can enter the reference evolution.
      const W H1 = out.physics.hubble((1 + initial) * std::exp(-s1) - 1),
              H2 = out.physics.hubble((1 + initial) * std::exp(-s2) - 1);
      auto stage_equation = [&](W coordinate, const State &y) {
        if (coordinate != s1 && coordinate != s2)
          throw std::runtime_error("reference stage coefficient coordinate");
        return out.physics.rhs(coordinate, y, coordinate == s1 ? H1 : H2);
      };
      auto step = base::radau(stage_equation, s, ds, state, out.stats);
      auto a = out.physics.opacity_at(s1, step.first[0], H1),
           b = out.physics.opacity_at(s2, step.last[0], H2);
      at += ds * (3 * a[0] + b[0]) / 4;
      ad += ds * (3 * a[1] + b[1]) / 4;
      state = step.end;
      s = begin + (k + 1) * ds;
      append(k + 1 == n ? stops[j] : (1 + initial) * std::exp(-s) - 1);
    }
  }
  for (W z : requested) {
    auto i = std::find_if(out.nodes.begin(), out.nodes.end(),
                          [z](const base::Node &n) { return n.z == z; });
    if (i == out.nodes.end())
      throw std::runtime_error("reference exact stop absent");
    out.queries.push_back(i - out.nodes.begin());
  }
  for (std::size_t i = 1; i < out.nodes.size(); ++i) {
    W a = out.nodes.back().drag - out.nodes[i - 1].drag,
      b = out.nodes.back().drag - out.nodes[i].drag;
    if (a >= 1 && b <= 1) {
      out.root = out.nodes[i].z +
                 (1 - b) * (out.nodes[i - 1].z - out.nodes[i].z) / (a - b);
      break;
    }
  }
  return out;
}
base::Quadrature opacity_gl(const Result &r) {
  std::vector<W> cumulative(r.nodes.size());
  W drag = 0;
  for (std::size_t i = 1; i < r.nodes.size(); ++i) {
    auto a = r.nodes[i - 1], b = r.nodes[i];
    auto q = [&](W s) {
      return r.physics.opacity_s(s, base::dense_x(a, b, s));
    };
    cumulative[i] =
        cumulative[i - 1] + base::gl8([&](W s) { return q(s)[0]; }, a.s, b.s);
    drag += base::gl8([&](W s) { return q(s)[1]; }, a.s, b.s);
  }
  W total = cumulative.back(), mass = 0;
  for (std::size_t i = 1; i < r.nodes.size(); ++i) {
    auto a = r.nodes[i - 1], b = r.nodes[i];
    auto q = [&](W s) {
      return r.physics.opacity_s(s, base::dense_x(a, b, s))[0];
    };
    mass += base::gl8(
        [&](W s) {
          W part = base::gl8(q, a.s, s);
          return q(s) * std::exp(-(total - cumulative[i - 1] - part));
        },
        a.s, b.s);
  }
  return {total, drag, mass, std::exp(-total)};
}
} // namespace hydrogen_relic_reference
