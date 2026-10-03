#pragma once
// Original numerical reference for the exactly emitted working source.
// No production headers/helpers or numerical outputs occur in this file.
// Constants, singlet rates and the earlier independent Radau mathematics have
// declared shared ancestry. These empirical numerical controls qualify neither
// a late ionization law nor the physical recombination/drag endpoint.
#include "hydrogen_helium_history_reference.hpp"
#include <cstddef>
#include <functional>
#include <utility>

namespace conditional_hydrogen_helium_drag_reference {
namespace ancestor = hydrogen_helium_history_reference;
using W = long double;
using State = ancestor::State;
using Vector = ancestor::Vector;
using Matrix = ancestor::Matrix;
constexpr W epsilon = std::numeric_limits<W>::epsilon();
constexpr W endpoint_allowance_mpc = 5e-5L;
// First source checkpoint: all affecting categories are retained below, but the
// work-count/epsilon/Newton proxy has NOT earned sensitivity/libm coverage.
// A passing preliminary sum cannot earn root's complete reference gate.
constexpr bool complete_arithmetic_qualified = false;

struct Source {
  // Widen the actual binary64 source fields, never decimal replacements.
  double H0, T0, gamma, other, baryon, cdm, nH0, nHe0, initial, late;
};
struct Limits {
  // Fixed before execution; shared across all three histories and endpoints.
  std::size_t steps = 200000, rhs = 40000000, quadrature = 8000000;
  unsigned inverse_trials = 128;
  W log_step = 1e-4L;
  unsigned ruler_panels = 4096;
};
struct Work {
  std::size_t steps = 0, rhs = 0, quadrature = 0, inverse = 0, newton = 0;
  W accepted_relative_corrections = 0;
  static void charge(std::size_t &counter, std::size_t n, std::size_t cap) {
    if (counter > cap || n > cap - counter)
      throw std::runtime_error("independent reference work ceiling");
    counter += n;
  }
};
struct Physics {
  Source source;
  W H0_si, photon_energy_today, radiation, matter, lambda, loading;
  explicit Physics(Source s) : source(s) {
    for (double x : {s.H0, s.T0, s.gamma, s.other, s.baryon, s.cdm,
                     s.nH0, s.nHe0, s.initial, s.late})
      if (!std::isfinite(x)) throw std::runtime_error("reference nonfinite source");
    if (!(s.H0 > 0 && s.T0 > 0 && s.gamma > 0 && s.baryon > 0 &&
          s.nH0 > 0 && s.nHe0 > 0 && s.other >= 0 && s.cdm >= 0 &&
          s.initial > s.late && s.late >= 0))
      throw std::runtime_error("reference source domain");
    H0_si = W(s.H0) * 1000 / ancestor::mpc();
    const W critical = 3 * H0_si * H0_si * ancestor::c * ancestor::c /
                       (8 * std::numbers::pi_v<W> * ancestor::G);
    photon_energy_today = W(s.gamma) * critical;
    radiation = W(s.gamma) + s.other;
    matter = W(s.baryon) + s.cdm;
    // Independently grouped closure on fixed emitted fractions. Its arithmetic
    // is retained below; do not import native lambda, R0, H, roots or errors.
    lambda = 1 - radiation - matter;
    loading = 3 * W(s.baryon) / (4 * W(s.gamma));
    if (!(lambda >= 0 && loading > 0 && std::isnormal(H0_si)))
      throw std::runtime_error("reference flat closure/loading");
  }
  W polynomial(W a) const {
    return radiation + matter * a + lambda * a * a * a * a;
  }
  W hubble(W u) const {
    return H0_si * std::sqrt(radiation * u * u * u * u +
                            matter * u * u * u + lambda);
  }
  State initial() const {
    const W u = 1 + W(source.initial), T = W(source.T0) * u,
        Q = std::pow(2 * std::numbers::pi_v<W> * ancestor::me * ancestor::kb * T /
                         (ancestor::h * ancestor::h), 1.5L),
        nH = W(source.nH0) * u * u * u, nHe = W(source.nHe0) * u * u * u,
        N = nH + nHe, A = Q * std::exp(-ancestor::chiH / (ancestor::kb * T)) / N,
        B = 4 * Q * std::exp(-ancestor::chiHe / (ancestor::kb * T)) / N;
    auto cubic = [&](W x) {
      return ((x / (A * B) + 1 / A + 1 / B) * x +
              1 - nH / (N * B) - nHe / (N * A)) * x - 1;
    };
    W lo = 0, hi = 1;
    if (!(cubic(lo) < 0 && cubic(hi) > 0))
      throw std::runtime_error("reference independent Saha bracket");
    for (unsigned i = 0; i < 128; ++i) {
      const W x = (lo + hi) / 2;
      if (cubic(x) < 0) lo = x; else hi = x;
    }
    const W x = (lo + hi) / 2;
    return {A / (A + x), B / (B + x), 1};
  }
  State operator()(W s, const State &state) const {
    const W u = (1 + W(source.initial)) * std::exp(-s),
        nH = W(source.nH0) * u * u * u, nHe = W(source.nHe0) * u * u * u,
        ne = nH * state[0] + nHe * state[1], H = hubble(u),
        Tr = W(source.T0) * u, T = state[2] * Tr,
        Q = std::pow(2 * std::numbers::pi_v<W> * ancestor::me * ancestor::kb * T /
                         (ancestor::h * ancestor::h), 1.5L),
        t = T / 10000,
        alphaH = 1e-19L * 4.309L * std::pow(t, -.6166L) /
                 (1 + .6703L * std::pow(t, .53L)),
        s0 = std::sqrt(T / std::pow(10.L, .477121L)),
        s1 = std::sqrt(T / std::pow(10.L, 5.114L)),
        alphaHe = std::pow(10.L, -16.744L) /
                  (s0 * std::pow(1 + s0, .289L) * std::pow(1 + s1, 1.711L)),
        EH = ancestor::h * ancestor::c / 121.5682e-9L,
        EHe = ancestor::h * ancestor::c * 1.66277434e7L,
        delta = ancestor::h * ancestor::c * (1.71134891e7L - 1.66277434e7L),
        betaH = alphaH * Q * std::exp(-(ancestor::chiH - EH) / (ancestor::kb * T)),
        betaHe = 4 * alphaHe * Q * std::exp(-(ancestor::chiHe - EHe) / (ancestor::kb * T)),
        KH = std::pow(121.5682e-9L, 3) / (8 * std::numbers::pi_v<W> * H),
        KHe = 1 / (8 * std::numbers::pi_v<W> * H * std::pow(1.71134891e7L, 3)),
        escapeH = 1 / (KH * nH * (1 - state[0])),
        escapeHe = std::exp(-delta / (ancestor::kb * T)) /
                   (KHe * nHe * (1 - state[1])),
        CH = (8.22458L + escapeH) / (8.22458L + escapeH + betaH),
        CHe = (51.3L + escapeHe) / (51.3L + escapeHe + betaHe),
        groundH = alphaH * Q * std::exp(-ancestor::chiH / (ancestor::kb * T)),
        groundHe = 4 * alphaHe * Q * std::exp(-ancestor::chiHe / (ancestor::kb * T)),
        gamma = 8 * ancestor::sigma * photon_energy_today * u * u * u * u /
                (3 * ancestor::me * ancestor::c * H) * ne / (nH + nHe + ne);
    return {-CH * (alphaH * ne * state[0] - groundH * (1 - state[0])) / H,
            -CHe * (alphaHe * ne * state[1] - groundHe * (1 - state[1])) / H,
            -state[2] - gamma * (state[2] - 1)};
  }
  W drag(W z, W p, W q) const {
    const W u = 1 + z, ne = u * u * u * (W(source.nH0) * p + W(source.nHe0) * q),
        R = loading / u;
    return ancestor::c * ancestor::sigma * ne / (hubble(u) * u * R);
  }
};

inline State radau_step(const Physics &f, W s, W ds, State y,
                        Work &work, const Limits &limits) {
  constexpr W A[2][2]{{5.L / 12, -1.L / 12}, {3.L / 4, 1.L / 4}};
  Vector v{y[0], y[1], y[2], y[0], y[1], y[2]};
  auto residual = [&](const Vector &x) {
    Work::charge(work.rhs, 2, limits.rhs);
    const State k[]{f(s + ds / 3, {x[0], x[1], x[2]}),
                    f(s + ds, {x[3], x[4], x[5]})};
    Vector r{};
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 3; ++j)
        r[3 * i + j] = x[3 * i + j] - y[j] -
                       ds * (A[i][0] * k[0][j] + A[i][1] * k[1][j]);
    return r;
  };
  for (unsigned iteration = 0; iteration < 48; ++iteration) {
    auto r = residual(v); ++work.newton;
    Matrix J{};
    for (unsigned j = 0; j < 6; ++j) {
      const W d = 1e-5L * (j % 3 == 2 ? v[j] : std::min(v[j], 1 - v[j]));
      auto plus = v, minus = v; plus[j] += d; minus[j] -= d;
      const auto rp = residual(plus), rm = residual(minus);
      for (unsigned i = 0; i < 6; ++i) J[i][j] = (rp[i] - rm[i]) / (2 * d);
    }
    for (W &x : r) x = -x;
    const auto direction = ancestor::solve(J, r);
    const W n = ancestor::norm(direction, v);
    if (!std::isfinite(n)) throw std::runtime_error("reference nonfinite Newton correction");
    if (n < 2e-16L) {
      Work::charge(work.steps, 1, limits.steps);
      work.accepted_relative_corrections += n;
      return {v[3], v[4], v[5]};
    }
    W damping = 1; bool accepted = false;
    for (unsigned trial = 0; trial < 64; ++trial) {
      auto next = v;
      for (unsigned j = 0; j < 6; ++j) next[j] += damping * direction[j];
      if (ancestor::positive(next)) {
        auto rr = residual(next); for (W &x : rr) x = -x;
        if (ancestor::norm(ancestor::solve(J, rr), next) < n) {
          v = next; accepted = true; break;
        }
      }
      damping /= 2;
    }
    if (!accepted) throw std::runtime_error("reference Newton line search");
  }
  throw std::runtime_error("reference Newton iteration ceiling");
}

template<class F> W simpson(F &&f, W a, W b, unsigned panels) {
  if (a == b) return 0;
  if (!(a < b) || !panels || panels % 2)
    throw std::runtime_error("reference Simpson interval/panels");
  const W dx = (b - a) / panels;
  W sum = f(a) + f(b), correction = 0;
  for (unsigned i = 1; i < panels; ++i) {
    const W term = (i % 2 ? 4 : 2) * f(a + dx * i), adjusted = term - correction,
        next = sum + adjusted;
    correction = (next - sum) - adjusted; sum = next;
  }
  return dx * sum / 3;
}

inline W ruler(const Physics &physics, W z, unsigned panels, Work &work,
                const Limits &limits) {
  const W a = 1 / (1 + z);
  return ancestor::c / 1000 / W(physics.source.H0) / std::sqrt(3.L) *
      simpson([&](W x) {
        Work::charge(work.quadrature, 1, limits.quadrature);
        const W aa = x * x;
        return 2 * x / std::sqrt(physics.polynomial(aa) * (1 + physics.loading * aa));
      }, 0, std::sqrt(a), panels);
}

struct Node {
  W z; State state;
  W prefix = 0, depth_quadrature = 0, cell_quadrature = 0;
};
struct Endpoint {
  W z = 0, ruler = 0, root_radius = 0, depth_quadrature_z = 0,
      arithmetic_z = 0, slope = 0, ruler_quadrature = 0, ruler_arithmetic = 0;
};
struct UniformBracketDiagnostic { W root_radius, depth_radius, slope; };
template<class F>
UniformBracketDiagnostic uniform_bracket(F &&depth, W target, W lo, W hi,
                                         W bracket_lo, W bracket_hi,
                                         W depth_error, W minimum_rate, W slope) {
  for (W x : {target, lo, hi, bracket_lo, bracket_hi, depth_error, minimum_rate, slope})
    if (!std::isfinite(x)) throw std::runtime_error("reference nonfinite uniform bracket");
  if (!(bracket_lo <= lo && lo <= hi && hi <= bracket_hi && depth_error >= 0 &&
        minimum_rate > 0 && slope > 0))
    throw std::runtime_error("reference uniform bracket domain");
  const W rho = depth_error / minimum_rate;
  if (depth_error > 0 && !(std::isnormal(rho) && rho > 0))
    throw std::runtime_error("reference positive inverse correction unresolved");
  // One outward neighbour covers the rounding of each locator-shift addition.
  // These tests include the FULL locator [lo,hi], not only its midpoint.
  const W expanded_lo = std::nextafter(lo - rho, -std::numeric_limits<W>::infinity()),
      expanded_hi = std::nextafter(hi + rho, std::numeric_limits<W>::infinity());
  if (!(depth(lo) <= target && depth(hi) >= target &&
        expanded_lo >= bracket_lo && expanded_hi <= bracket_hi &&
        depth(bracket_lo) + depth_error <= target &&
        depth(bracket_hi) - depth_error >= target))
    throw std::runtime_error("reference expanded locator/sign bracket did not close");
  return {(hi - lo) / 2, rho, slope};
}

inline W loading_xi(W R, W ER, W amax) {
  // Stable exact-real identity demanded by the root admission. This independent
  // control never supplies a native loading/error value or reference oracle.
  if (!(std::isfinite(R) && std::isfinite(ER) && std::isfinite(amax) &&
        R > ER && ER >= 0 && amax > 0))
    throw std::runtime_error("reference loading correction domain");
  if (ER == 0) return 0;
  const W A = 1 + R * amax, B = 1 + (R - ER) * amax,
      xi = ER * amax / (std::sqrt(B) * (std::sqrt(A) + std::sqrt(B)));
  if (!(std::isnormal(xi) && xi > 0))
    throw std::runtime_error("reference positive loading correction unresolved");
  return xi;
}
class History {
  Physics physics_;
  std::vector<Node> nodes_; // decreasing z; positive prefix from the late end
  Work &work_;
  const Limits &limits_;
  W arithmetic_relative_ = 0;
  W cell(std::size_t i, W a, W b, unsigned panels) const {
    const auto &hi = nodes_[i], &lo = nodes_[i + 1];
    return simpson([&](W z) {
      Work::charge(work_.quadrature, 1, limits_.quadrature);
      const W f = (hi.z - z) / (hi.z - lo.z),
          p = (1 - f) * hi.state[0] + f * lo.state[0],
          q = (1 - f) * hi.state[1] + f * lo.state[1],
          w = physics_.drag(z, p, q);
      if (!(w > 0 && std::isfinite(w))) throw std::runtime_error("reference positive drag");
      return w;
    }, a, b, panels);
  }
  W cell_quadrature_bound(std::size_t i, unsigned panels) const {
    // Exact-real composite Simpson truncation bound, separately from floating
    // arithmetic. It bounds EVERY restricted part of this represented cell.
    // w=k*u^3*C(u)*Q(u)^(-1/2), C positive affine,
    // Q=r*u^4+m*u^3+lambda. Nonnegative background coefficients make Q monotone.
    // Leibniz and the chain rule bound |w''''| without a cancellation-prone
    // full-cell difference being used as a partial-cell error bound.
    const auto &hi = nodes_[i], &lo = nodes_[i + 1];
    const W length = hi.z - lo.z, u = 1 + hi.z, v = 1 + lo.z,
        chi = W(physics_.source.nH0) * hi.state[0] + W(physics_.source.nHe0) * hi.state[1],
        clo = W(physics_.source.nH0) * lo.state[0] + W(physics_.source.nHe0) * lo.state[1],
        charge = std::max(chi, clo), dcharge = std::abs(chi - clo) / length,
        r = physics_.radiation, m = physics_.matter,
        q = r * v * v * v * v + m * v * v * v + physics_.lambda,
        q1 = 4 * r * u * u * u + 3 * m * u * u,
        q2 = 12 * r * u * u + 6 * m * u,
        q3 = 24 * r * u + 6 * m, q4 = 24 * r;
    if (!(length > 0 && chi > 0 && clo > 0 && q > 0 && panels && panels % 2 == 0))
      throw std::runtime_error("reference cell quadrature domain");
    const W f0 = u * u * u * charge,
        f1 = 3 * u * u * charge + u * u * u * dcharge,
        f2 = 6 * u * charge + 6 * u * u * dcharge,
        f3 = 6 * charge + 18 * u * dcharge, f4 = 24 * dcharge,
        iq = 1 / q, g0 = 1 / std::sqrt(q),
        g1 = g0 * (.5L * q1 * iq),
        g2 = g0 * (.75L * q1 * q1 * iq * iq + .5L * q2 * iq),
        g3 = g0 * (1.875L * q1 * q1 * q1 * iq * iq * iq +
                    2.25L * q1 * q2 * iq * iq + .5L * q3 * iq),
        g4 = g0 * (6.5625L * q1 * q1 * q1 * q1 * iq * iq * iq * iq +
                    11.25L * q1 * q1 * q2 * iq * iq * iq +
                    (2.25L * q2 * q2 + 3 * q1 * q3) * iq * iq + .5L * q4 * iq),
        fourth = ancestor::c * ancestor::sigma / (physics_.H0_si * physics_.loading) *
                 (f4 * g0 + 4 * f3 * g1 + 6 * f2 * g2 + 4 * f1 * g3 + f0 * g4),
        n = panels,
        error = length * length * length * length * length * fourth / (180 * n * n * n * n);
    if (!(std::isfinite(error) && error > 0))
      throw std::runtime_error("reference positive cell quadrature allowance unresolved");
    return error;
  }
  std::size_t containing(W z) const {
    if (z < nodes_.back().z || z > nodes_.front().z)
      throw std::runtime_error("reference retained support");
    auto pos = std::lower_bound(nodes_.begin(), nodes_.end(), z,
                               [](const Node &n, W x) { return n.z > x; });
    auto j = std::size_t(pos - nodes_.begin());
    if (j == 0) j = 1;
    if (j == nodes_.size()) j = nodes_.size() - 1;
    return j - 1;
  }
  W depth(W z) const {
    if (z == nodes_.back().z) return 0;
    const auto i = containing(z);
    return nodes_[i + 1].prefix + cell(i, nodes_[i + 1].z, z, 4);
  }
  W depth_quadrature_bound(W z) const {
    if (z == nodes_.back().z) return 0;
    const auto i = containing(z);
    // Even if z is interior, the whole-cell fourth-derivative bound covers its
    // shorter Simpson4 primitive. Shared prefix terms are conservatively summed.
    return nodes_[i + 1].depth_quadrature + nodes_[i].cell_quadrature;
  }
public:
  History(Source source, unsigned refinement, Work &work, const Limits &limits)
      : physics_(source), work_(work), limits_(limits) {
    if (!refinement || !std::isnormal(limits.log_step))
      throw std::runtime_error("reference fixed mesh");
    const auto rhs_start = work.rhs;
    const W corrections_start = work.accepted_relative_corrections;
    State state = physics_.initial(); W s = 0;
    nodes_.push_back({W(source.initial), state});
    for (W z : {W(source.initial) - .1L, W(source.late)}) {
      const W end = std::log((1 + W(source.initial)) / (1 + z)), length = end - s;
      std::size_t count = std::max<std::size_t>(1, std::ceil(length / limits.log_step));
      if (z > W(source.initial) - 1) count = std::max<std::size_t>(count, 128);
      if (count > limits.steps / refinement) throw std::runtime_error("reference mesh ceiling");
      count *= refinement;
      if (work.steps > limits.steps || count > limits.steps - work.steps)
        throw std::runtime_error("reference aggregate mesh ceiling");
      const W ds = length / count, begin = s;
      nodes_.reserve(nodes_.size() + count);
      for (std::size_t k = 0; k < count; ++k) {
        state = radau_step(physics_, s, ds, state, work, limits);
        s = begin + (k + 1) * ds;
        const W node_z = k + 1 == count ? z : (1 + W(source.initial)) * std::exp(-s) - 1;
        nodes_.push_back({node_z, state});
      }
    }
    // Explicit arithmetic/iteration diagnostic, not a continuous proof. Includes
    // SI constants/pi/closure, all charged RHS/Newton work, and accumulated final
    // correction estimates. Whole-pipeline refinement below remains mandatory.
    arithmetic_relative_ = 4096 * epsilon + 512 * epsilon * (work.rhs - rhs_start) +
        (work.accepted_relative_corrections - corrections_start);
    W correction = 0;
    for (std::size_t i = nodes_.size() - 1; i-- > 0;) {
      const W coarse = cell(i, nodes_[i + 1].z, nodes_[i].z, 2),
          fine = cell(i, nodes_[i + 1].z, nodes_[i].z, 4),
          adjusted = fine - correction, next = nodes_[i + 1].prefix + adjusted;
      correction = (next - nodes_[i + 1].prefix) - adjusted;
      nodes_[i].prefix = next;
      nodes_[i].cell_quadrature = cell_quadrature_bound(i, 4);
      // Retain the independent 2/4 discrepancy in addition to the derivative
      // bound; neither discrepancy nor a Richardson divisor hides cancellation.
      nodes_[i].depth_quadrature = nodes_[i + 1].depth_quadrature +
                                 nodes_[i].cell_quadrature + std::abs(fine - coarse);
    }
  }
  Endpoint endpoint(double supplied_D, double supplied_late) const {
    const W D = supplied_D, late = supplied_late;
    if (!(std::isfinite(late) && late >= nodes_.back().z && late < nodes_.front().z))
      throw std::runtime_error("reference supplied late endpoint");
    const W base_depth = depth(late), target = base_depth + (1 - D);
    if (!(std::isfinite(D) && D >= 0 && D <= 1 && target <= nodes_.front().prefix))
      throw std::runtime_error("reference supplied depth support");
    W lo = late, hi = nodes_.front().z;
    if (D != 1) {
      auto pos = std::lower_bound(nodes_.begin(), nodes_.end(), target,
                                 [](const Node &n, W x) { return n.prefix > x; });
      auto j = std::size_t(pos - nodes_.begin());
      if (j == 0) j = 1;
      if (j >= nodes_.size()) throw std::runtime_error("reference inverse cell");
      lo = nodes_[j].z; hi = nodes_[j - 1].z;
      if (!(depth(lo) <= target && depth(hi) >= target))
        throw std::runtime_error("reference prefix/direct-primitive locator mismatch");
      bool done = false;
      for (unsigned trial = 0; trial < limits_.inverse_trials; ++trial) {
        ++work_.inverse;
        const W mid = lo + (hi - lo) / 2;
        if (depth(mid) < target) lo = mid; else hi = mid;
        if (hi - lo <= 256 * epsilon * (1 + std::abs(mid))) { done = true; break; }
      }
      if (!done) throw std::runtime_error("reference inverse ceiling");
    } else hi = lo;
    const W z = lo + (hi - lo) / 2;
    const auto i = containing(z);
    // The neighbouring reference cells give ONE uniform inverse/slope bracket.
    // Arithmetic/depth shifts must close inside it; a central slope is not a
    // propagation bound. No numerical outputs from the native locator enter it.
    const auto first = i == 0 ? i : i - 1,
               last = std::min(i + 1, nodes_.size() - 2);
    const W low_u = 1 + nodes_[last + 1].z, high_u = 1 + nodes_[first].z,
        amin = 1 / high_u, amax = 1 / low_u;
    W charge_min = std::numeric_limits<W>::infinity();
    for (std::size_t j = first; j <= last + 1; ++j)
      charge_min = std::min(charge_min, W(physics_.source.nH0) * nodes_[j].state[0] +
                                          W(physics_.source.nHe0) * nodes_[j].state[1]);
    const W
        minimum_rate = ancestor::c * ancestor::sigma * low_u * low_u * low_u * charge_min /
                       (physics_.hubble(high_u) * physics_.loading),
        slope = (ancestor::c / 1000 / W(physics_.source.H0)) * amax * amax /
                std::sqrt(3 * physics_.polynomial(amin) * (1 + physics_.loading * amin));
    if (!(minimum_rate > 0 && slope > 0)) throw std::runtime_error("reference positive inverse/slope");
    const W r0 = ruler(physics_, z, limits_.ruler_panels, work_, limits_),
        r1 = ruler(physics_, z, 2 * limits_.ruler_panels, work_, limits_),
        depth_arithmetic = D == 1 ? 0 : arithmetic_relative_ * nodes_.front().prefix +
                           512 * epsilon * nodes_.size() * nodes_.front().prefix +
                           16 * epsilon * (1 + std::abs(target) + base_depth),
        depth_quadrature = D == 1 ? 0 : depth_quadrature_bound(nodes_[first].z) +
                                        depth_quadrature_bound(late),
        depth_error = depth_arithmetic + depth_quadrature;
    const W bracket_lo = nodes_[last + 1].z, bracket_hi = nodes_[first].z,
        radius = (hi - lo) / 2;
    const auto diagnostic = D == 1 ? UniformBracketDiagnostic{radius, 0, slope} :
        uniform_bracket([&](W x) { return depth(x); }, target, lo, hi,
                        bracket_lo, bracket_hi, depth_error, minimum_rate, slope);
    return {z, r1, diagnostic.root_radius,
            depth_quadrature / minimum_rate,
            depth_arithmetic / minimum_rate, slope, std::abs(r1 - r0),
            1024 * epsilon * (1 + 3 * limits_.ruler_panels) * std::abs(r1)};
  }
};

struct CompleteError {
  // Complete category ledger with EMPIRICAL arithmetic proxies. The name does
  // not assert a rigorous bound or the withheld combined qualification above.
  W history_depth_refinement = 0, depth_quadrature = 0, ruler_quadrature = 0,
      root_bracket = 0, all_arithmetic = 0, root_redshift = 0;
  W total() const {
    return history_depth_refinement + depth_quadrature + ruler_quadrature +
           root_bracket + all_arithmetic;
  }
};
inline CompleteError complete_error(const Endpoint &coarse, const Endpoint &middle,
                                    const Endpoint &fine) {
  // Linear sum, with no presumed Richardson factor or independent-error RSS.
  return {std::abs(coarse.ruler - middle.ruler) + std::abs(middle.ruler - fine.ruler),
          fine.slope * fine.depth_quadrature_z, fine.ruler_quadrature,
          fine.slope * fine.root_radius,
          fine.slope * fine.arithmetic_z + fine.ruler_arithmetic +
              std::abs(W(double(fine.ruler)) - fine.ruler),
          std::abs(coarse.z - middle.z) + std::abs(middle.z - fine.z) +
              fine.root_radius + fine.depth_quadrature_z + fine.arithmetic_z +
              std::abs(W(double(fine.z)) - fine.z)};
}

// Independent exact-radical Gaussian controls. Stored constants are generated
// from these original expressions, not production coefficients or samples.
template<class F> W gauss4(F &&f, W a, W b) {
  const W r = std::sqrt(6.L / 5), x0 = std::sqrt((3 - 2 * r) / 7),
      x1 = std::sqrt((3 + 2 * r) / 7),
      w0 = (18 + std::sqrt(30.L)) / 36, w1 = (18 - std::sqrt(30.L)) / 36,
      mid = a + (b - a) / 2, half = (b - a) / 2;
  return half * (w0 * (f(mid - half * x0) + f(mid + half * x0)) +
                 w1 * (f(mid - half * x1) + f(mid + half * x1)));
}
template<class F> W gauss5(F &&f, W a, W b) {
  const W r = std::sqrt(10.L / 7), x0 = std::sqrt((5 - 2 * r) / 9),
      x1 = std::sqrt((5 + 2 * r) / 9),
      w0 = (322 + 13 * std::sqrt(70.L)) / 900,
      w1 = (322 - 13 * std::sqrt(70.L)) / 900,
      mid = a + (b - a) / 2, half = (b - a) / 2;
  return half * (128.L / 225 * f(mid) + w0 * (f(mid - half * x0) + f(mid + half * x0)) +
                 w1 * (f(mid - half * x1) + f(mid + half * x1)));
}
inline W rational_cell_mass(W a, W b) {
  W sum = 0;
  for (unsigned j = 0; j <= 6; ++j) sum += std::pow(1 + b, 6 - j) * std::pow(1 + a, j);
  return (b - a) * sum / 28;
}
} // namespace conditional_hydrogen_helium_drag_reference
