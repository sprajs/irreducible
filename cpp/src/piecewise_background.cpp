#include "irred/piecewise_background.hpp"
#include "flat_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
// Entire-function continuation, not a q threshold replacing the true integral.
long double exprel(long double t) {
  if (t == 0)
    return 1;
  if (std::abs(t) > .0001L)
    return std::expm1(t) / t;
  long double sum = 1, term = 1;
  for (int k = 1; k <= 12; ++k) {
    term *= t / (k + 1);
    sum += term;
  }
  return sum;
}
bool normal_or_zero(long double x) {
  if (!detail::physical_representable(x))
    return false;
  auto y = static_cast<double>(x);
  return x == 0 || (y != 0 && std::fpclassify(y) == FP_NORMAL);
}
numerics::Status range_cause(long double x) {
  return !std::isfinite(x) || std::abs(x) > std::numeric_limits<double>::max()
             ? numerics::Status::overflow
             : numerics::Status::outside_domain;
}
} // namespace
PiecewiseBackground prepare_piecewise_q(PiecewiseQParameters p) {
  PiecewiseBackground out;
  out.parameters_ = p;
  if (!std::isfinite(p.h0_km_s_mpc) || p.h0_km_s_mpc <= 0)
    return out;
  for (auto q : p.q) {
    if (!std::isfinite(q))
      return out;
    if (q < -3 || q > 2) {
      out.status_ = Status::unsupported_domain;
      return out;
    }
  }
  auto scale = detail::prepare_flat_scale(p.h0_km_s_mpc);
  if (scale.status != Status::ok) {
    out.status_ = scale.status;
    return out;
  }
  out.hubble_distance_mpc_ = scale.distance_mpc;
  out.hubble_time_seconds_ = scale.time_seconds;
  long double E = 1;
  for (std::size_t k = 0; k < 5; ++k) {
    out.start_E_[k] = E;
    auto lo = (long double)piecewise_q_edges[k],
         hi = (long double)piecewise_q_edges[k + 1];
    E *= std::exp((1 + (long double)p.q[k]) * std::log1p((hi - lo) / (1 + lo)));
  }
  out.status_ = Status::ok;
  return out;
}
PiecewiseBatch
PiecewiseBackground::evaluate_batch(std::span<const Query> queries,
                                    PiecewisePolicy policy) const {
  PiecewiseBatch out;
  out.status = status_;
  if (status_ != Status::ok)
    return out;
  if (queries.size() > policy.maximum_queries) {
    out.status = Status::work_limit;
    return out;
  }
  out.slots.resize(queries.size());
  std::size_t remaining = policy.maximum_segment_visits;
  for (std::size_t i = 0; i < queries.size(); ++i) {
    auto &s = out.slots[i];
    auto query = queries[i];
    s.source = query;
    if (query.convention != Convention::geometric_same_redshift &&
        query.convention != Convention::released_zhd_zhel) {
      s.status = Status::incompatible_convention;
      continue;
    }
    if (!std::isfinite(query.z_expansion) || !std::isfinite(query.z_observer)) {
      s.status = Status::invalid_input;
      continue;
    }
    if (query.z_expansion < 0 || query.z_expansion > 2.5) {
      s.status = Status::unsupported_domain;
      continue;
    }
    if (query.z_observer <= -1 ||
        (query.convention == Convention::geometric_same_redshift &&
         query.z_observer != query.z_expansion)) {
      s.status = Status::incompatible_convention;
      continue;
    }
    std::size_t k =
        std::upper_bound(piecewise_q_edges.begin() + 1, piecewise_q_edges.end(),
                         query.z_expansion) -
        piecewise_q_edges.begin() - 1;
    k = std::min(k, std::size_t(4));
    std::size_t visits = 0;
    for (std::size_t b = 0; b < 5; ++b)
      if (query.z_expansion > piecewise_q_edges[b])
        ++visits;
    if (visits > remaining) {
      s.status = Status::work_limit;
      s.numerical_status = numerics::Status::work_limit;
      continue;
    }
    remaining -= visits;
    s.segments_processed = visits;
    out.segments_processed += visits;
    long double radial = 0, clock = 0;
    for (std::size_t b = 0; b < visits; ++b) {
      auto lo = (long double)piecewise_q_edges[b];
      auto hi = std::min((long double)query.z_expansion,
                         (long double)piecewise_q_edges[b + 1]);
      auto L = std::log1p((hi - lo) / (1 + lo));
      auto q = (long double)parameters_.q[b];
      radial += (1 + lo) / start_E_[b] * L * exprel(-q * L);
      clock += 1 / start_E_[b] * L * exprel(-(1 + q) * L);
    }
    auto lo = (long double)piecewise_q_edges[k];
    auto L = std::log1p(((long double)query.z_expansion - lo) / (1 + lo));
    auto E = start_E_[k] * std::exp((1 + (long double)parameters_.q[k]) * L);
    auto H = (long double)parameters_.h0_km_s_mpc * E;
    bool valid = true;
    for (auto x : {E, H, radial, clock})
      if (!normal_or_zero(x)) {
        s.numerical_status = range_cause(x);
        valid = false;
        break;
      }
    if (query.z_expansion > 0 && (!(radial > 0) || !(clock > 0))) {
      valid = false;
      s.numerical_status = numerics::Status::outside_domain;
    }
    if (!valid) {
      s.status = Status::numerical_failure;
      continue;
    }
    auto g = detail::flat_geometry(query, E, (double)radial, (double)clock,
                                   hubble_distance_mpc_, hubble_time_seconds_);
    if (g.status != Status::ok) {
      s.status = g.status;
      s.numerical_status = g.numerical_status;
      continue;
    }
    for (auto x : {g.dc, g.da, g.shape, g.dl, g.lookback, g.volume})
      if (!normal_or_zero(x)) {
        s.numerical_status = range_cause(x);
        valid = false;
        break;
      }
    if (!valid) {
      s.status = Status::numerical_failure;
      continue;
    }
    s.status = Status::ok;
    s.numerical_status = numerics::Status::ok;
    s.bin = k;
    s.geometry = {
        static_cast<double>(E),          static_cast<double>(H),
        static_cast<double>(radial),     static_cast<double>(g.dc),
        static_cast<double>(g.dc),       static_cast<double>(g.da),
        static_cast<double>(g.dl),       static_cast<double>(g.shape),
        static_cast<double>(g.lookback), static_cast<double>(g.volume)};
    s.assigned_q = parameters_.q[k];
    s.q_convention = QConvention::interior_constant_bin;
    bool jump = k > 0 && query.z_expansion == piecewise_q_edges[k] &&
                parameters_.q[k] != parameters_.q[k - 1];
    if (jump) {
      s.q_convention = QConvention::right_limit_at_internal_jump;
      s.jerk_availability = JerkAvailability::unavailable_at_jump;
    } else {
      s.jerk = (double)((long double)parameters_.q[k] *
                        (2 * (long double)parameters_.q[k] + 1));
      s.jerk_availability = JerkAvailability::ordinary_within_bin;
    }
    if (query.z_expansion == 0) {
      s.q_convention = QConvention::right_limit_at_zero;
      s.jerk_availability = JerkAvailability::one_sided_endpoint;
      s.q0_within_piecewise_model = parameters_.q[0];
    }
    if (query.z_expansion == 2.5) {
      s.q_convention = QConvention::left_limit_at_final_endpoint;
      s.jerk_availability = JerkAvailability::one_sided_endpoint;
    }
  }
  out.status = Status::ok;
  return out;
}
} // namespace irred::cosmology
