#pragma once
#include "irred/hydrogen_helium_history.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace irred::cosmology::detail {
// One retained within-cell law serves projected history queries and the drag
// primitive. A is the interpolated nodal c*sigma_T/(H*u), not a newly evaluated
// interior background quotient. No full second cell/history vector is owned.
struct HydrogenHeliumCell {
  using W = long double;
  W high = 0, low = 0;
  // Native order is high, low. Curvature stores the original second-slope estimate.
  std::array<W, 2> hydrogen{}, helium{}, temperature{}, coefficient{};
  std::array<W, 2> hydrogen_error{}, helium_error{}, temperature_error{}, coefficient_error{};
  std::array<W, 4> curvature{};
  struct Value {
    W hydrogen, helium, temperature, hydrogen_error, helium_error, temperature_error;
    W nH, nHe, ne, ne_error, opacity, opacity_error;
  };
  Value evaluate(W z, W hydrogen_today, W helium_today) const noexcept {
    const W dz = high - low, f = (high - z) / dz;
    auto curve = [&](unsigned k) { return f * (1 - f) * dz * dz * curvature[k] / 2; };
    auto blend = [&](const auto &v) { return (1 - f) * v[0] + f * v[1]; };
    const W p = blend(hydrogen), q = blend(helium),
        ep = blend(hydrogen_error) + curve(0),
        eq = blend(helium_error) + curve(1),
        T = blend(temperature), eT = blend(temperature_error) + curve(2),
        u = 1 + z, nH = hydrogen_today * u * u * u,
        nHe = helium_today * u * u * u, ne = nH * p + nHe * q,
        ene = nH * ep + nHe * eq + 128 * std::numeric_limits<double>::epsilon() * ne,
        A = blend(coefficient), eA = blend(coefficient_error) + curve(3);
    return {p, q, T, ep, eq, eT, nH, nHe, ne, ene, A * ne, A * ene + (ne + ene) * eA};
  }
};

struct HydrogenHeliumCellAccess {
  static std::size_t count(const HydrogenHeliumHistory &h) noexcept {
    return h.nodes_.empty() ? 0 : h.nodes_.size() - 1;
  }
  // Index is in the retained descending-z order. Only bounded callers with
  // status==ok and index<count enter here; no publicly supplied unchecked index.
  static HydrogenHeliumCell cell(const HydrogenHeliumHistory &h, std::size_t i) noexcept {
    using W = long double;
    using Node = HydrogenHeliumHistory::Node;
    const auto &hi = h.nodes_[i], &lo = h.nodes_[i + 1];
    const W dz = hi.z - lo.z;
    auto curvature = [&](auto member) {
      const W slope = (hi.*member - lo.*member) / dz;
      W second = 0;
      if (i > 0) {
        const auto &p = h.nodes_[i - 1];
        second = std::max(second, std::abs((p.*member - hi.*member) / (p.z - hi.z) - slope) / ((p.z - lo.z) / 2));
      }
      if (i + 2 < h.nodes_.size()) {
        const auto &n = h.nodes_[i + 2];
        second = std::max(second, std::abs(slope - (lo.*member - n.*member) / (lo.z - n.z)) / ((hi.z - n.z) / 2));
      }
      return second;
    };
    return {hi.z, lo.z, {hi.hydrogen, lo.hydrogen}, {hi.helium, lo.helium},
            {hi.temperature, lo.temperature}, {hi.opacity_coefficient, lo.opacity_coefficient},
            {hi.hydrogen_error, lo.hydrogen_error}, {hi.helium_error, lo.helium_error},
            {hi.temperature_error, lo.temperature_error},
            {hi.opacity_coefficient_error, lo.opacity_coefficient_error},
            {curvature(&Node::hydrogen), curvature(&Node::helium),
             curvature(&Node::temperature), curvature(&Node::opacity_coefficient)}};
  }
  static std::optional<std::size_t> retained_payload(const HydrogenHeliumHistory &h) noexcept {
    irred::detail::PayloadAccounting p(sizeof(h));
    p.vector(h.nodes_);
    p.vector(h.background_.source().species);
    if (h.source_) {
      p.string(h.source_->nuclei_origin);
      p.vector(h.source_->model.species);
    }
    return p.result();
  }
};

// Bounded nonnegative Bernstein algebra. Restrict before integration; small
// intervals never subtract nearby global antiderivatives. The arithmetic
// allowance is explicit empirical wide-operation diagnostics, not certification.
struct HHePositivePolynomial {
  using W = long double;
  std::array<W, 9> coefficient{};
  unsigned degree = 0;
  bool valid = true;
  static constexpr W arithmetic = 4096 * std::numeric_limits<W>::epsilon();
  static unsigned choose(unsigned n, unsigned k) noexcept {
    unsigned v = 1;
    for (unsigned i = 1; i <= k; ++i) v = v * (n - i + 1) / i;
    return v;
  }
  bool admissible() const noexcept {
    if (!valid || degree > 8) return false;
    for (unsigned i = 0; i <= degree; ++i)
      if (coefficient[i] < 0 || !std::isfinite(coefficient[i]) ||
          (coefficient[i] != 0 && !std::isnormal(coefficient[i]))) return false;
    return true;
  }
  HHePositivePolynomial scaled(W scale) const noexcept {
    auto out = *this;
    if (!(scale >= 0) || !std::isfinite(scale)) { out.valid = false; return out; }
    for (unsigned i = 0; i <= degree; ++i) out.coefficient[i] *= scale;
    return out;
  }
  HHePositivePolynomial elevated(unsigned n) const noexcept {
    auto out = *this;
    if (n < degree || n > 8) { out.valid = false; return out; }
    while (out.degree < n) {
      const auto old = out.coefficient;
      const unsigned d = ++out.degree;
      out.coefficient[0] = old[0]; out.coefficient[d] = old[d - 1];
      for (unsigned i = 1; i < d; ++i)
        out.coefficient[i] = (W(i) * old[i - 1] + W(d - i) * old[i]) / d;
    }
    return out;
  }
  HHePositivePolynomial plus(const HHePositivePolynomial &other) const noexcept {
    const auto n = std::max(degree, other.degree);
    auto out = elevated(n); const auto q = other.elevated(n);
    out.valid = out.valid && q.valid;
    for (unsigned i = 0; i <= n; ++i) out.coefficient[i] += q.coefficient[i];
    return out;
  }
  HHePositivePolynomial times(const HHePositivePolynomial &q) const noexcept {
    HHePositivePolynomial out;
    out.degree = degree + q.degree;
    out.valid = valid && q.valid && out.degree <= 8;
    if (!out.valid) return out;
    for (unsigned i = 0; i <= degree; ++i)
      for (unsigned j = 0; j <= q.degree; ++j)
        out.coefficient[i + j] += coefficient[i] * q.coefficient[j] *
            (W(choose(degree, i)) * choose(q.degree, j) / choose(out.degree, i + j));
    return out;
  }
  std::array<HHePositivePolynomial, 2> split(W t) const noexcept {
    auto left = *this, right = *this;
    auto triangle = coefficient;
    left.coefficient[0] = triangle[0]; right.coefficient[degree] = triangle[degree];
    for (unsigned level = 1; level <= degree; ++level) {
      for (unsigned i = 0; i <= degree - level; ++i)
        triangle[i] = (1 - t) * triangle[i] + t * triangle[i + 1];
      left.coefficient[level] = triangle[0];
      right.coefficient[degree - level] = triangle[degree - level];
    }
    left.valid = right.valid = valid && t >= 0 && t <= 1 && std::isfinite(t);
    return {left, right};
  }
  HHePositivePolynomial restricted(W low, W high) const noexcept {
    if (!(low >= 0) || !(high >= low) || high > 1) { auto out = *this; out.valid = false; return out; }
    if (low == 0) return split(high)[0];
    if (low == 1) return split(1)[1];
    return split(low)[1].split((high - low) / (1 - low))[0];
  }
  W integral(W width) const noexcept {
    W sum = 0, correction = 0;
    for (unsigned i = 0; i <= degree; ++i) {
      const W y = coefficient[i] - correction, next = sum + y;
      correction = (next - sum) - y; sum = next;
    }
    return width * (sum / (degree + 1));
  }
  W lower() const noexcept {
    W lo = coefficient[0], hi = coefficient[0];
    for (unsigned i = 1; i <= degree; ++i) { lo = std::min(lo, coefficient[i]); hi = std::max(hi, coefficient[i]); }
    return lo - arithmetic * hi;
  }
};
struct HHeDragCell {
  numerics::Status status = numerics::Status::invalid_input;
  long double low = 0, high = 0;
  HHePositivePolynomial rate, error;
  struct Integral { numerics::Status status; long double value, error; };
  Integral integrate(long double zlow, long double zhigh) const noexcept {
    using S = numerics::Status;
    if (status != S::ok) return {status, 0, 0};
    if (zlow < low || zhigh > high || zhigh < zlow) return {S::outside_domain, 0, 0};
    if (zlow == zhigh) return {S::ok, 0, 0};
    const long double width = high - low, l = (zlow - low) / width, h = (zhigh - low) / width;
    const auto r = rate.restricted(l, h), e = error.restricted(l, h);
    const long double v = r.integral(zhigh - zlow),
        diagnostic = e.integral(zhigh - zlow) + HHePositivePolynomial::arithmetic * v;
    if (!r.admissible() || !e.admissible() || !std::isnormal(v) || !std::isnormal(diagnostic))
      return {S::conditioning_budget_exceeded, 0, 0};
    return {S::ok, v, diagnostic};
  }
  long double lower(long double zlow, long double zhigh) const noexcept {
    const long double dz = high - low;
    return rate.restricted((zlow - low) / dz, (zhigh - low) / dz).lower();
  }
};
inline HHeDragCell hydrogen_helium_drag_cell(const HydrogenHeliumCell &c,
    long double nH, long double nHe, long double R, long double ER) noexcept {
  using W = long double;
  using P = HHePositivePolynomial;
  using S = numerics::Status;
  HHeDragCell out; out.low = c.low; out.high = c.high;
  if (!(R > ER) || !(ER >= 0) || !(c.high > c.low) || !(nH > 0) || !(nHe > 0)) return out;
  auto affine = [](const std::array<W, 2> &v) { P p; p.degree = 1; p.coefficient[0] = v[1]; p.coefficient[1] = v[0]; return p; };
  auto error = [&](const std::array<W, 2> &v, W curvature) {
    auto p = affine(v).elevated(2);
    const W dz = c.high - c.low;
    p.coefficient[1] += dz * dz * curvature / 4;
    return p;
  };
  P u; u.degree = 1; u.coefficient[0] = 1 + c.low; u.coefficient[1] = 1 + c.high;
  const auto u3 = u.times(u).times(u),
      C = affine(c.hydrogen).scaled(nH).plus(affine(c.helium).scaled(nHe)),
      eC = error(c.hydrogen_error, c.curvature[0]).scaled(nH).plus(error(c.helium_error, c.curvature[1]).scaled(nHe)),
      ne = u3.times(C),
      ene = u3.times(eC).plus(ne.scaled(128 * std::numeric_limits<double>::epsilon())),
      A = affine(c.coefficient), eA = error(c.coefficient_error, c.curvature[3]),
      qT = A.times(ne), eqT = A.times(ene).plus(ne.plus(ene).times(eA));
  out.rate = u.times(qT).scaled(1 / R);
  out.error = u.times(eqT.plus(qT.plus(eqT).scaled(ER / (R - ER)))).scaled(1 / R)
      .plus(out.rate.scaled(P::arithmetic));
  out.status = out.rate.admissible() && out.error.admissible() && out.rate.lower() > 0
      ? S::ok : S::conditioning_budget_exceeded;
  return out;
}
} // namespace irred::cosmology::detail
