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
    if (h.supplied_initial_) {
      p.string(h.supplied_initial_->origin);
      p.string(h.supplied_initial_->source_identity);
    }
    return p.result();
  }
};
} // namespace irred::cosmology::detail
