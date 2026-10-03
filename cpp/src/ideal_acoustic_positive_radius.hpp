#pragma once
#include "ideal_acoustic_transport.hpp"

namespace irred::cosmology::detail {
using IdealAcousticRadiusVector = std::array<long double, 12>;
struct IdealAcousticRadiusDiagonal {
  std::array<long double, 6> mu_lower, mu_upper, positive;
  numerics::Status status = numerics::Status::invalid_input;
};

// Private conditional PC2 assembly. This is not a continuous upper-envelope
// certificate. Central Y and signed Z keep their separate classical RK4 law.
// Work retains the original owned-state/derivative destination semantics:
// charge each stored diagonal/output slot and one diagnostic per graph. Local
// Arithmetic operations are radius assembly, not a generic FP-write counter.
namespace ideal_acoustic_positive_internal {
using W = long double;
using S = numerics::Status;
using Arithmetic = ideal_acoustic_transport_internal::Arithmetic;
inline bool profile() noexcept {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
inline bool radius(W v) noexcept {
  return ideal_acoustic_transport_internal::radius(v);
}
inline bool valid(const IdealAcousticRadiusDiagonal &d) noexcept {
  if (d.status != S::ok)
    return false;
  for (std::size_t i = 0; i < 6; ++i)
    if (!radius(d.mu_lower[i]) || !radius(d.mu_upper[i]) ||
        !radius(d.positive[i]) || d.mu_lower[i] > d.mu_upper[i])
      return false;
  return true;
}
inline S begin(const IdealAcousticTransportAccounting &work) noexcept {
  const S status = ideal_acoustic_transport_internal::account(work);
  return status == S::ok && !profile() ? S::invalid_input : status;
}
inline bool store(W value, W &destination,
                  const IdealAcousticTransportAccounting &work) noexcept {
  if (!work.state_writes(work.context, 1))
    return false;
  destination = value;
  return true;
}
// Each exact-zero branch has an exact represented target. A positive product
// rounded to zero/subnormal is refused by Arithmetic, never replaced by zero.
inline W denominator_lower(W h, W mu, Arithmetic &a) noexcept {
  if (mu == 0)
    return 1;
  const W product = a.product_lower(h, mu);
  const W raw = a.add(1, product);
  const W result = a.checked(std::nextafter(raw, W(0)));
  if (!(result > 0))
    a.status = S::conditioning_budget_exceeded;
  return result;
}
inline W numerator_upper(W half_h, W mu_lower, W mu_upper,
                         Arithmetic &a) noexcept {
  if (mu_upper == 0)
    return 1;
  // The guard uses an UPPER damping product. The reported upper factor uses
  // LOWER damping. Reversing those two directions would overdamp the radius.
  (void)a.lower(1, a.times(half_h, mu_upper));
  if (a.status != S::ok)
    return 0;
  if (mu_lower == 0)
    return 1;
  const W factor = a.signed_upper(a.sub(1, a.product_lower(half_h, mu_lower)));
  if (!(factor > 0) || !radius(factor))
    a.status = S::conditioning_budget_exceeded;
  return factor;
}
inline bool positive_inputs(std::span<const W, 12> v) noexcept {
  for (W x : v)
    if (!radius(x))
      return false;
  return true;
}
inline bool positive_inputs(std::span<const W, 6> v) noexcept {
  for (W x : v)
    if (!radius(x))
      return false;
  return true;
}
} // namespace ideal_acoustic_positive_internal

inline numerics::Status ideal_acoustic_radius_diagonal(
    const IdealAcousticCoefficients &actual, IdealAcousticRadiusDiagonal &out,
    const IdealAcousticTransportAccounting &work) noexcept {
  namespace p = ideal_acoustic_positive_internal;
  using S = numerics::Status;
  using W = long double;
  out.status = p::begin(work);
  if (out.status != S::ok)
    return out.status;
  if (!ideal_acoustic_transport_internal::normal_or_zero(actual.L) ||
      !p::radius(actual.loading_over_one_plus_loading))
    return out.status = S::outside_domain;
  p::Arithmetic a;
  for (std::size_t i = 0; i < 6; ++i) {
    W low = 0, high = 0, positive = 0;
    if (i == 4) {
      low = high = 1; // phi diagonal is the exact represented constant -1.
    } else if (i == 2 || i == 3) {
      const W subtrahend = i == 2 ? actual.loading_over_one_plus_loading : 1;
      const W diagonal = a.sub(actual.L, subtrahend);
      const W allowance = a.loss(a.plus(a.magnitude(actual.L), subtrahend));
      const W upper = a.signed_upper(a.add(diagonal, allowance));
      const W raw_lower = a.sub(diagonal, allowance);
      const W lower = raw_lower == 0
                          ? 0
                          : a.checked(std::nextafter(
                                raw_lower, -std::numeric_limits<W>::infinity()));
      // Zero lower damping is a coefficient enclosure, not state clipping.
      low = upper < 0 ? -upper : 0;
      high = lower < 0 ? -lower : 0;
      positive = upper > 0 ? upper : 0;
    }
    if (a.status != S::ok || !p::radius(low) || !p::radius(high) ||
        !p::radius(positive) || low > high)
      return out.status = a.status == S::ok ? S::outside_domain : a.status;
    if (!p::store(low, out.mu_lower[i], work) ||
        !p::store(high, out.mu_upper[i], work) ||
        !p::store(positive, out.positive[i], work))
      return out.status = S::work_limit;
  }
  return out.status = S::ok;
}

inline numerics::Status ideal_acoustic_radius_predictor(
    std::span<const long double, 12> base,
    std::span<const long double, 12> p_start, long double h,
    const IdealAcousticRadiusDiagonal &diagonal,
    std::span<long double, 12> out,
    const IdealAcousticTransportAccounting &work) noexcept {
  namespace p = ideal_acoustic_positive_internal;
  using S = numerics::Status;
  const S begin = p::begin(work);
  if (begin != S::ok)
    return begin;
  if (!(h > 0) || !p::radius(h) || !p::valid(diagonal) ||
      !p::positive_inputs(base) || !p::positive_inputs(p_start))
    return S::invalid_input;
  p::Arithmetic a;
  for (std::size_t i = 0; i < 12; ++i) {
    const long double denominator =
        p::denominator_lower(h, diagonal.mu_lower[i % 6], a);
    const long double value =
        a.over(a.plus(base[i], a.times(h, p_start[i])), denominator);
    if (a.status != S::ok)
      return a.status;
    if (!p::store(value, out[i], work))
      return S::work_limit;
  }
  return S::ok;
}

inline numerics::Status ideal_acoustic_radius_corrector(
    std::span<const long double, 12> base,
    std::span<const long double, 12> p_start,
    std::span<const long double, 12> p_endpoint, long double h,
    const IdealAcousticRadiusDiagonal &start_diagonal,
    const IdealAcousticRadiusDiagonal &endpoint_diagonal,
    std::span<long double, 12> out,
    const IdealAcousticTransportAccounting &work) noexcept {
  namespace p = ideal_acoustic_positive_internal;
  using S = numerics::Status;
  const S begin = p::begin(work);
  if (begin != S::ok)
    return begin;
  if (!(h > 0) || !p::radius(h) || !p::valid(start_diagonal) ||
      !p::valid(endpoint_diagonal) || !p::positive_inputs(base) ||
      !p::positive_inputs(p_start) || !p::positive_inputs(p_endpoint))
    return S::invalid_input;
  p::Arithmetic a;
  const long double half_h = a.div(h, 2);
  if (a.status != S::ok)
    return a.status;
  for (std::size_t i = 0; i < 12; ++i) {
    const std::size_t coordinate = i % 6;
    const long double factor =
        p::numerator_upper(half_h, start_diagonal.mu_lower[coordinate],
                           start_diagonal.mu_upper[coordinate], a);
    const long double denominator = p::denominator_lower(
        half_h, endpoint_diagonal.mu_lower[coordinate], a);
    const long double numerator =
        a.plus(a.times(factor, base[i]),
               a.times(half_h, a.plus(p_start[i], p_endpoint[i])));
    const long double value = a.over(numerator, denominator);
    if (a.status != S::ok)
      return a.status;
    if (!p::store(value, out[i], work))
      return S::work_limit;
  }
  return S::ok;
}

inline numerics::Status ideal_acoustic_radius_stage(
    std::span<const long double, 12> force,
    std::span<const long double, 6> fresh_assembly, long double weight,
    std::span<long double, 12> out,
    const IdealAcousticTransportAccounting &work) noexcept {
  namespace p = ideal_acoustic_positive_internal;
  using S = numerics::Status;
  const S begin = p::begin(work);
  if (begin != S::ok)
    return begin;
  if (!p::radius(weight) || !p::positive_inputs(force) ||
      !p::positive_inputs(fresh_assembly))
    return S::invalid_input;
  p::Arithmetic a;
  for (std::size_t i = 0; i < 12; ++i) {
    const long double value =
        a.plus(a.times(weight, force[i]), i < 6 ? 0 : fresh_assembly[i - 6]);
    if (a.status != S::ok)
      return a.status;
    if (!p::store(value, out[i], work))
      return S::work_limit;
  }
  return S::ok;
}

inline numerics::Status ideal_acoustic_radius_endpoint(
    std::span<const long double, 12> flow,
    const std::array<const IdealAcousticRadiusVector *, 4> &forces,
    const std::array<long double, 4> &weights,
    std::span<const long double, 6> final_assembly,
    std::span<long double, 12> out,
    const IdealAcousticTransportAccounting &work) noexcept {
  namespace p = ideal_acoustic_positive_internal;
  using S = numerics::Status;
  const S begin = p::begin(work);
  if (begin != S::ok)
    return begin;
  if (!p::positive_inputs(flow) || !p::positive_inputs(final_assembly))
    return S::invalid_input;
  for (std::size_t s = 0; s < 4; ++s)
    if (!forces[s] || !p::radius(weights[s]) ||
        !p::positive_inputs(std::span<const long double, 12>(*forces[s])))
      return S::invalid_input;
  p::Arithmetic a;
  for (std::size_t i = 0; i < 12; ++i) {
    long double value = flow[i];
    for (std::size_t s = 0; s < 4; ++s)
      value = a.plus(value, a.times(weights[s], (*forces[s])[i]));
    value = a.plus(value, i < 6 ? 0 : final_assembly[i - 6]);
    if (a.status != S::ok)
      return a.status;
    if (!p::store(value, out[i], work))
      return S::work_limit;
  }
  return S::ok;
}
} // namespace irred::cosmology::detail
