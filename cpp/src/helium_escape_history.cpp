#include "irred/helium_escape_history.hpp"
#include "helium_escape_rates.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
namespace irred::cosmology {
namespace {
using W = long double;
using S = numerics::Status;
constexpr std::size_t hard_bytes = std::size_t(1) << 30;
constexpr W floor64 = 128 * std::numeric_limits<double>::epsilon();
struct DriverState { W ell, z, H, nH, T, xHI; };
struct MeshState { W q = 0, root_error = 0; };
bool arithmetic() noexcept {
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
  return false;
#else
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<double>::radix == 2 &&
         std::numeric_limits<double>::digits == 53 &&
         std::numeric_limits<double>::max_exponent == 1024 &&
         std::numeric_limits<W>::radix == 2 &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
#endif
}
bool budget(double a, double r) noexcept {
  return std::isfinite(a) && std::isfinite(r) && a >= 0 && r >= 0 &&
         (a > 0 || r > 0);
}
bool physical_fraction(W q) noexcept {
  return std::isfinite(q) && q > 0 && q < 1;
}
bool charge(std::size_t &part, const HeliumEscapeHistoryWork &work,
            std::size_t cap) noexcept {
  if (work.total() >= cap) return false;
  ++part; return true;
}
void range(std::optional<double> &lo, std::optional<double> &hi, W value) {
  if (!std::isfinite(value)) return;
  const double low = std::nextafter(static_cast<double>(value), -INFINITY),
      high = std::nextafter(static_cast<double>(value), INFINITY);
  if (!lo || low < *lo) lo = low;
  if (!hi || high > *hi) hi = high;
}
DriverState knot(const HeliumEscapeDriverKnot &s) {
  return {-std::log1p(W(s.redshift)), W(s.redshift), W(s.hubble_per_second),
          W(s.hydrogen_nuclei_per_cubic_metre),
          W(s.radiation_temperature_kelvin), W(s.hydrogen_neutral_fraction)};
}
DriverState between(const HeliumEscapeDriverKnot &a,
                    const HeliumEscapeDriverKnot &b, W ell, W z) {
  const W high = -std::log1p(W(a.redshift)), low = -std::log1p(W(b.redshift));
  if (ell == high) { auto c = knot(a); c.z = z; return c; }
  if (ell == low) { auto c = knot(b); c.z = z; return c; }
  const W t = (ell - high) / (low - high);
  auto blend = [&](double x, double y) {
    if (x == y) return W(x);
    return std::exp((1 - t) * std::log(W(x)) + t * std::log(W(y)));
  };
  return {ell, z, blend(a.hubble_per_second, b.hubble_per_second),
          blend(a.hydrogen_nuclei_per_cubic_metre, b.hydrogen_nuclei_per_cubic_metre),
          blend(a.radiation_temperature_kelvin, b.radiation_temperature_kelvin),
          blend(a.hydrogen_neutral_fraction, b.hydrogen_neutral_fraction)};
}
S admit(const HeliumEscapeDriverKnot &k) noexcept {
  for (double v : {k.redshift, k.hubble_per_second,
      k.hydrogen_nuclei_per_cubic_metre, k.radiation_temperature_kelvin,
      k.hydrogen_neutral_fraction}) if (!std::isfinite(v)) return S::nonfinite_input;
  if (k.redshift < 1900 || k.redshift > 2800 ||
      k.hubble_per_second < 1e-14 || k.hubble_per_second > 1e-12 ||
      k.hydrogen_nuclei_per_cubic_metre < 5e8 || k.hydrogen_nuclei_per_cubic_metre > 8e9 ||
      k.radiation_temperature_kelvin < 4300 || k.radiation_temperature_kelvin > 7800 ||
      k.hydrogen_neutral_fraction <= 0 || k.hydrogen_neutral_fraction >= 1)
    return S::outside_domain;
  return S::ok;
}
struct Evolution {
  HeliumEscapeHistoryWork &work;
  HeliumEscapeHistoryWitness &witness;
  std::size_t cap;
  W f;
  S status = S::ok;
  detail::HeliumEscapeRate rate(const DriverState &c, W q) {
    if (!charge(work.rate_evaluations, work, cap)) {
      status = S::work_limit; return {};
    }
    ++witness.attempted_rates;
    if (!std::isfinite(c.T)) {
      ++witness.nonfinite_temperature_attempts; status = S::nonfinite_input; return {};
    }
    range(witness.minimum_attempted_temperature_kelvin,
          witness.maximum_attempted_temperature_kelvin, c.T);
    const W tau = detail::helium_escape_optical_depth(c.H, c.nH, f, q);
    range(witness.minimum_attempted_tau, witness.maximum_attempted_tau, tau);
    if (!physical_fraction(q) || c.T < 4300 || c.T > 7800 ||
        !(tau >= 100) || !(f * q >= 1e-6L)) {
      status = S::outside_domain; return {};
    }
    const auto r = detail::helium_escape_rate(c.H, c.nH, c.T, c.xHI, f, q);
    for (W x : {r.tau, r.continuum_inverse_seconds, r.width_per_second,
                r.continuum_depth, r.enhancement, r.effective_escape,
                r.downward_per_second, r.saha, r.electron_density,
                r.absolute_rate_scale})
      if (!(x > 0) || !std::isfinite(x)) { status = S::conditioning_budget_exceeded; return {}; }
    if (!std::isfinite(r.value) || !std::isfinite(r.fraction_derivative))
      status = S::conditioning_budget_exceeded;
    return r;
  }
  MeshState step(const DriverState &c, MeshState previous, W h) {
    W q = previous.q;
    auto residual = [&](W x, W &correction, W &arithmetic_error) {
      const auto r = rate(c, x);
      if (status != S::ok) return false;
      const W J = 1 - h * r.fraction_derivative,
          R = x - previous.q - h * r.value,
          scale = std::max(W(1), std::abs(h * r.fraction_derivative));
      if (!(J > 64 * std::numeric_limits<W>::epsilon() * scale) || !std::isfinite(J)) {
        status = S::conditioning_budget_exceeded; return false;
      }
      correction = -R / J;
      // Empirical finite-driver/libm/residual sensitivity allowance. The
      // factor128 covers the bounded temperature exponent scales; it is not
      // a certified libm bound and faces the independent per-output gate.
      arithmetic_error = (256 * std::numeric_limits<W>::epsilon() + floor64) *
          (std::abs(x) + std::abs(previous.q) + 128 * h * r.absolute_rate_scale) / J;
      if (!std::isfinite(correction) || !std::isfinite(arithmetic_error)) {
        status = S::conditioning_budget_exceeded; return false;
      }
      return true;
    };
    for (unsigned trial = 0; trial < 80; ++trial) {
      W d = 0, e = 0;
      if (!residual(q, d, e)) return {};
      const W norm = std::abs(d) / q;
      if (norm <= 64 * std::numeric_limits<W>::epsilon())
        return {q, previous.root_error + std::abs(d) + e};
      W damping = 1; bool accepted = false;
      for (unsigned attempt = 0; attempt < 64; ++attempt) {
        const W next = q + damping * d;
        if (physical_fraction(next)) {
          W nd = 0, ne = 0;
          if (!residual(next, nd, ne)) return {};
          if (std::abs(nd) / next < norm) { q = next; accepted = true; break; }
        }
        damping /= 2;
      }
      if (!accepted) { status = S::conditioning_budget_exceeded; return {}; }
    }
    status = S::conditioning_budget_exceeded; return {};
  }
};
void project(HeliumEscapeHistoryValue &v, W x, W error, W absolute,
             W relative, bool fraction = false) {
  v.status = S::conditioning_budget_exceeded;
  if (!(x > 0) || !std::isfinite(x) || !(error >= 0) || !std::isfinite(error) ||
      (fraction && !(x < 1))) return;
  const double rounded = static_cast<double>(x);
  error += std::abs(x - W(rounded)) + floor64 * x;
  if (!(rounded > 0) || !std::isfinite(rounded) ||
      (fraction && !(rounded < 1)) || error > std::numeric_limits<double>::max()) return;
  v.absolute_error_estimate = std::nextafter(static_cast<double>(error), INFINITY);
  if (W(v.absolute_error_estimate) > absolute + relative * std::abs(x)) return;
  v.value = rounded; v.status = S::ok;
}
std::size_t byte_cap(const HeliumEscapeHistoryPolicy &p) noexcept {
  return std::min(p.maximum_native_bytes, hard_bytes);
}
void clone_preflight(const HeliumEscapeHistory &o, std::size_t extra) {
  const auto retained = o.retained_payload_bound();
  irred::detail::PayloadAccounting p(extra);
  if (!retained) throw std::length_error("helium history payload overflow");
  p.add(2, *retained);
  if (!p.result() || *p.result() > hard_bytes)
    throw std::length_error("helium history simultaneous copy payload");
}
} // namespace
std::optional<std::size_t> helium_escape_history_payload_bound(
    std::size_t fine, std::size_t knots, std::size_t identity, std::size_t outputs) noexcept {
  if (fine == std::numeric_limits<std::size_t>::max() ||
      identity == std::numeric_limits<std::size_t>::max()) return {};
  irred::detail::PayloadAccounting p(sizeof(HeliumEscapeHistory) +
      sizeof(HeliumEscapeHistoryBatch) + 4096);
  p.add(fine + 1, sizeof(DriverState) + sizeof(HeliumEscapeHistory::Node) +
                         4 * sizeof(MeshState));
  p.add(knots, sizeof(HeliumEscapeDriverKnot) + sizeof(std::size_t));
  p.add(identity + 1, 1); p.add(outputs, sizeof(HeliumEscapeHistoryRow));
  return p.result();
}
std::optional<std::size_t> HeliumEscapeHistory::retained_payload_bound() const noexcept {
  irred::detail::PayloadAccounting p(sizeof(*this)); p.vector(nodes_);
  if (source_) { p.vector(source_->knots); p.string(source_->producer_identity); }
  return p.result();
}
void HeliumEscapeHistory::discard_failed_payload() noexcept {
  // A refused trajectory is unavailable: release its capacity, not just size.
  std::vector<Node>{}.swap(nodes_);
  const auto exceeds = [&]() {
    const auto retained = retained_payload_bound();
    return !retained || *retained > byte_cap(policy_);
  };
  if (!source_ || !exceeds()) return;
  // Complete lawful captures fit at acquisition and survive later refusals.
  // An oversized acquisition must not retain its rejected reservation or claim
  // that discarded input is complete. Work and attempted witnesses stay intact.
  std::vector<HeliumEscapeDriverKnot>{}.swap(source_->knots);
  witness_.source_capture_complete = false;
  if (!exceeds()) return;
  std::string{}.swap(source_->producer_identity);
  if (exceeds()) source_.reset();
}
HeliumEscapeHistory::HeliumEscapeHistory(const HeliumEscapeHistory &o) {
  clone_preflight(o, 0);
  const auto old = o.retained_payload_bound();
  if (!old || *old > byte_cap(o.policy_) / 2)
    throw std::length_error("helium history simultaneous copy exceeds policy");
  source_ = o.source_; nodes_ = o.nodes_; status_ = o.status_;
  policy_ = o.policy_; work_ = o.work_; witness_ = o.witness_;
  const auto now = retained_payload_bound();
  if (!now || *now > byte_cap(policy_) - *old)
    throw std::length_error("helium history actual copy capacities exceed policy");
}
HeliumEscapeHistory &HeliumEscapeHistory::operator=(const HeliumEscapeHistory &o) {
  if (this == &o) return *this;
  const auto prior = retained_payload_bound(), other = o.retained_payload_bound();
  if (!prior || !other) throw std::length_error("helium history assignment payload overflow");
  const std::size_t cap = std::min(byte_cap(policy_), byte_cap(o.policy_));
  irred::detail::PayloadAccounting p(*prior); p.add(2, *other);
  if (!p.result() || *p.result() > cap)
    throw std::length_error("helium history replacement exceeds policy");
  HeliumEscapeHistory copy(o);
  const auto actual = copy.retained_payload_bound();
  irred::detail::PayloadAccounting complete(*prior); complete.add(1, *other);
  if (!actual) throw std::length_error("helium history replacement capacities overflow");
  complete.add(1, *actual);
  if (!complete.result() || *complete.result() > cap)
    throw std::length_error("helium history actual replacement capacities exceed policy");
  *this = std::move(copy); return *this;
}
HeliumEscapeHistory::HeliumEscapeHistory(HeliumEscapeHistory &&o) noexcept {
  *this = std::move(o);
}
HeliumEscapeHistory &HeliumEscapeHistory::operator=(HeliumEscapeHistory &&o) noexcept {
  if (this != &o) {
    status_ = o.status_; source_ = std::move(o.source_); policy_ = o.policy_;
    work_ = o.work_; witness_ = o.witness_; nodes_ = std::move(o.nodes_);
    o.source_.reset(); std::vector<Node>{}.swap(o.nodes_); o.status_ = S::invalid_input;
    o.work_ = {}; o.witness_ = {};
  }
  return *this;
}
HeliumEscapeHistory prepare_helium_escape_history(
    const HeliumEscapeHistoryRequest &request, HeliumEscapeHistoryPolicy policy) {
  HeliumEscapeHistory out; out.policy_ = policy;
  if (!arithmetic() || !budget(policy.absolute_fraction_tolerance, policy.relative_fraction_tolerance) ||
      !budget(policy.absolute_opacity_tolerance, policy.relative_opacity_tolerance) ||
      request.producer_identity.empty() || policy.base_intervals < 2) return out;
  if (request.knots.size() < 2) return out;
  if (request.knots.size() > 4097 || request.producer_identity.size() > 65536 ||
      policy.base_intervals > 16384 || policy.base_intervals < request.knots.size() - 1 ||
      4 * policy.base_intervals > policy.maximum_fine_intervals) {
    out.status_ = S::work_limit; return out;
  }
  const std::size_t N = 4 * policy.base_intervals, cap = byte_cap(policy);
  const auto planned = helium_escape_history_payload_bound(N, request.knots.size(),
      request.producer_identity.size(), 0);
  if (!planned || *planned > cap || request.knots.size() > policy.maximum_total_work) {
    out.status_ = S::work_limit; return out;
  }
  try {
    out.source_.emplace(); auto &source = *out.source_;
    source.helium_to_hydrogen_nuclei_ratio = request.helium_to_hydrogen_nuclei_ratio;
    source.initial_helium_singly_ionized_fraction = request.initial_helium_singly_ionized_fraction;
    source.role = request.role; source.producer_identity = request.producer_identity;
    source.knots.reserve(request.knots.size());
    auto acquired = out.retained_payload_bound();
    if (!acquired || *acquired > cap) {
      out.status_ = S::work_limit; out.discard_failed_payload(); return out;
    }
    for (const auto &row : request.knots) {
      if (!charge(out.work_.imported_rows, out.work_, policy.maximum_total_work)) {
        out.status_ = S::work_limit; return out;
      }
      source.knots.push_back(row);
    }
    out.witness_.source_capture_complete = true;
    // Complete acquisition remains available through every subsequent refusal.
    if (out.work_.total() >= policy.maximum_total_work) { out.status_ = S::work_limit; return out; }
    const W f = source.helium_to_hydrogen_nuclei_ratio,
        q0 = source.initial_helium_singly_ionized_fraction;
    if (!std::isfinite(f) || !std::isfinite(q0)) { out.status_ = S::nonfinite_input; return out; }
    if (f < .04L || f > .12L || !physical_fraction(q0)) { out.status_ = S::outside_domain; return out; }
    if (source.role != HeliumEscapeDriverRole::supplied_cosmological_history &&
        source.role != HeliumEscapeDriverRole::synthetic_prescribed_bath) return out;
    for (std::size_t i = 0; i < source.knots.size(); ++i) {
      const S s = admit(source.knots[i]);
      if (s != S::ok) { out.status_ = s; return out; }
      if (i && !(source.knots[i - 1].redshift > source.knots[i].redshift)) return out;
    }
    if (source.knots.front().redshift < 2600 || source.knots.front().redshift > 2800 ||
        source.knots.back().redshift < 1900 || source.knots.back().redshift > 2200) {
      out.status_ = S::outside_domain; return out;
    }
    std::vector<std::size_t> split(source.knots.size() - 1, 1);
    const W begin = -std::log1p(W(source.knots.front().redshift)),
        end = -std::log1p(W(source.knots.back().redshift)), span = end - begin;
    std::size_t left = policy.base_intervals - split.size();
    const std::size_t remaining = left;
    for (std::size_t i = 0; i < split.size(); ++i) {
      const W width = std::log1p(W(source.knots[i].redshift)) -
                      std::log1p(W(source.knots[i + 1].redshift));
      const auto extra = static_cast<std::size_t>(std::floor(W(remaining) * width / span));
      if (extra > left) { out.status_ = S::conditioning_budget_exceeded; return out; }
      split[i] += extra; left -= extra;
    }
    for (std::size_t i = 0; left; ++i, --left) ++split[i % split.size()];
    std::vector<DriverState> coefficients(N + 1);
    std::vector<MeshState> coarse, middle, fine;
    auto payload = [&]() {
      auto retained = out.retained_payload_bound();
      if (!retained) return false;
      irred::detail::PayloadAccounting p(*retained);
      p.add(sizeof(HeliumEscapeHistoryBatch) + 4096, 1);
      p.vector(split); p.vector(coefficients); p.vector(coarse); p.vector(middle); p.vector(fine);
      return p.result() && *p.result() <= cap;
    };
    if (!payload()) { out.status_ = S::work_limit; return out; }
    std::size_t j = 0;
    for (std::size_t cell = 0; cell < split.size(); ++cell) {
      const auto &high = source.knots[cell], &low = source.knots[cell + 1];
      const W a = -std::log1p(W(high.redshift)), b = -std::log1p(W(low.redshift));
      const std::size_t intervals = 4 * split[cell];
      for (std::size_t k = cell ? 1 : 0; k <= intervals; ++k) {
        if (!charge(out.work_.driver_evaluations, out.work_, policy.maximum_total_work)) {
          out.status_ = S::work_limit; return out;
        }
        const W ell = k == 0 ? a : k == intervals ? b : a + (b - a) * W(k) / intervals,
            z = k == 0 ? W(high.redshift) : k == intervals ? W(low.redshift) : std::expm1(-ell);
        coefficients[j] = between(high, low, ell, z);
        if (j && !(coefficients[j].ell > coefficients[j - 1].ell)) {
          out.status_ = S::conditioning_budget_exceeded; return out;
        }
        ++j;
      }
    }
    if (j != N + 1) { out.status_ = S::conditioning_budget_exceeded; return out; }
    Evolution evolution{out.work_, out.witness_, policy.maximum_total_work, f};
    auto mesh = [&](std::vector<MeshState> &states, std::size_t stride) {
      states.resize(N / stride + 1);
      if (!payload()) { evolution.status = S::work_limit; return; }
      states[0] = {q0, 0};
      // Boundary source rate is an actual checked attempt, with the same witnesses.
      (void)evolution.rate(coefficients.front(), q0);
      if (evolution.status != S::ok) return;
      for (std::size_t i = 1; i < states.size(); ++i) {
        const std::size_t k = i * stride;
        states[i] = evolution.step(coefficients[k], states[i - 1],
            coefficients[k].ell - coefficients[k - stride].ell);
        if (evolution.status != S::ok) return;
      }
    };
    mesh(coarse, 4); if (evolution.status != S::ok) { out.status_ = evolution.status; return out; }
    mesh(middle, 2); if (evolution.status != S::ok) { out.status_ = evolution.status; return out; }
    mesh(fine, 1); if (evolution.status != S::ok) { out.status_ = evolution.status; return out; }
    auto interpolate = [&](const std::vector<MeshState> &states, std::size_t stride,
                           std::size_t i) {
      const std::size_t k = i / stride;
      if (i % stride == 0) return states[k];
      const W a = coefficients[k * stride].ell, b = coefficients[(k + 1) * stride].ell,
          t = (coefficients[i].ell - a) / (b - a);
      return MeshState{(1 - t) * states[k].q + t * states[k + 1].q,
                       (1 - t) * states[k].root_error + t * states[k + 1].root_error};
    };
    auto interpolation_error = [&](const std::vector<MeshState> &states,
                                   std::size_t stride, std::size_t i) {
      if (i % stride == 0) return W(0);
      const std::size_t k = i / stride;
      const W a = coefficients[k * stride].ell, b = coefficients[(k + 1) * stride].ell,
          width = b - a, t = (coefficients[i].ell - a) / width,
          slope = (states[k + 1].q - states[k].q) / width;
      W second = 0;
      if (k) {
        const W previous = coefficients[(k - 1) * stride].ell;
        second = std::max(second, std::abs(slope - (states[k].q - states[k - 1].q) /
            (a - previous)) / ((b - previous) / 2));
      }
      if (k + 2 < states.size()) {
        const W next = coefficients[(k + 2) * stride].ell;
        second = std::max(second, std::abs((states[k + 2].q - states[k + 1].q) /
            (next - b) - slope) / ((next - a) / 2));
      }
      return t * (1 - t) * width * width * second / 2;
    };
    out.nodes_.resize(N + 1);
    if (!payload()) { out.status_ = S::work_limit; out.discard_failed_payload(); return out; }
    W max_activity = -INFINITY;
    for (std::size_t i = 0; i <= N; ++i) {
      const auto mid = interpolate(middle, 2, i), old = interpolate(coarse, 4, i);
      const W q = 2 * fine[i].q - mid.q,
          error = std::abs(q - (2 * mid.q - old.q)) / 3 +
              2 * interpolation_error(middle, 2, i) + interpolation_error(coarse, 4, i) / 3 +
              2 * (fine[i].root_error + mid.root_error + old.root_error) + floor64 * q;
      if (!physical_fraction(q) || !(error >= 0) || !std::isfinite(error)) {
        out.status_ = S::conditioning_budget_exceeded; out.discard_failed_payload(); return out;
      }
      const auto r = evolution.rate(coefficients[i], q);
      if (evolution.status != S::ok) { out.status_ = evolution.status; out.discard_failed_payload(); return out; }
      const W activity = detail::helium_escape_heiii_log_activity(coefficients[i].T,
                                                                 r.electron_density);
      if (!std::isfinite(activity)) { out.status_ = S::conditioning_budget_exceeded; out.discard_failed_payload(); return out; }
      max_activity = std::max(max_activity, activity);
      out.witness_.maximum_log_retained_heiii_activity =
          std::nextafter(static_cast<double>(max_activity), INFINITY);
      if (max_activity > std::log(1e-12L)) { out.status_ = S::outside_domain; out.discard_failed_payload(); return out; }
      out.nodes_[i] = {coefficients[i].ell, coefficients[i].z, q, error, 0};
    }
    for (std::size_t i = 0; i < N; ++i) {
      const W width = out.nodes_[i + 1].ell - out.nodes_[i].ell,
          slope = (out.nodes_[i + 1].fraction - out.nodes_[i].fraction) / width;
      W second = 0;
      if (i) second = std::max(second, std::abs(slope -
          (out.nodes_[i].fraction - out.nodes_[i - 1].fraction) /
          (out.nodes_[i].ell - out.nodes_[i - 1].ell)) /
          ((out.nodes_[i + 1].ell - out.nodes_[i - 1].ell) / 2));
      if (i + 2 <= N) second = std::max(second, std::abs(
          (out.nodes_[i + 2].fraction - out.nodes_[i + 1].fraction) /
          (out.nodes_[i + 2].ell - out.nodes_[i + 1].ell) - slope) /
          ((out.nodes_[i + 2].ell - out.nodes_[i].ell) / 2));
      out.nodes_[i].curvature = second;
    }
    out.status_ = S::ok;
  } catch (const std::bad_alloc &) { out.status_ = S::work_limit; out.discard_failed_payload(); }
  catch (const std::length_error &) { out.status_ = S::work_limit; out.discard_failed_payload(); }
  return out;
}
HeliumEscapeHistoryBatch HeliumEscapeHistory::evaluate(
    std::span<const double> redshifts, unsigned mask, std::size_t maximum_points,
    std::size_t maximum_bytes) const {
  HeliumEscapeHistoryBatch out; out.requested_outputs = mask;
  if (!mask || (mask & ~7u) || !arithmetic()) return out;
  if (status_ != S::ok) { out.status = status_; return out; }
  if (redshifts.size() > maximum_points || redshifts.size() > 4096) {
    out.status = S::work_limit; return out;
  }
  const std::size_t cap = std::min(maximum_bytes, byte_cap(policy_));
  auto retained = retained_payload_bound();
  irred::detail::PayloadAccounting planned(sizeof(out));
  if (!retained) { out.status = S::work_limit; return out; }
  planned.add(1, *retained); planned.add(redshifts.size(), sizeof(HeliumEscapeHistoryRow));
  if (!planned.result() || *planned.result() > cap ||
      redshifts.size() > policy_.maximum_total_work - work_.total()) {
    out.status = S::work_limit; return out;
  }
  try {
    out.rows.resize(redshifts.size());
    irred::detail::PayloadAccounting actual(sizeof(out));
    actual.add(1, *retained); actual.vector(out.rows);
    if (!actual.result() || *actual.result() > cap) {
      out.status = S::work_limit; std::vector<HeliumEscapeHistoryRow>{}.swap(out.rows); return out;
    }
    for (std::size_t i = 0; i < redshifts.size(); ++i) {
      auto &row = out.rows[i]; const double z = redshifts[i]; row.redshift = z;
      HeliumEscapeHistoryValue *values[]{&row.helium_singly_ionized_fraction,
          &row.electron_number_density_per_cubic_metre, &row.thomson_opacity_per_redshift};
      auto refuse = [&](S s) { for (unsigned k = 0; k < 3; ++k)
        if (mask & (1u << k)) values[k]->status = s; };
      if (!std::isfinite(z)) { refuse(S::nonfinite_input); continue; }
      if (z > source_->knots.front().redshift || z < source_->knots.back().redshift) {
        refuse(S::outside_domain); continue;
      }
      ++out.driver_evaluations;
      const W ell = -std::log1p(W(z));
      auto lower = std::lower_bound(nodes_.begin(), nodes_.end(), ell,
                                   [](const Node &n, W x) { return n.ell < x; });
      std::size_t j = lower - nodes_.begin();
      if (!j) j = 1;
      if (j == nodes_.size()) j = nodes_.size() - 1;
      const auto &a = nodes_[j - 1], &b = nodes_[j];
      const W width = b.ell - a.ell, t = (ell - a.ell) / width,
          q = (1 - t) * a.fraction + t * b.fraction,
          eq = (1 - t) * a.fraction_error + t * b.fraction_error +
               t * (1 - t) * width * width * a.curvature / 2 + floor64 * q;
      auto source_low = std::lower_bound(source_->knots.begin(), source_->knots.end(), z,
          [](const HeliumEscapeDriverKnot &n, double x) { return n.redshift > x; });
      std::size_t k = source_low - source_->knots.begin();
      if (!k) k = 1;
      if (k == source_->knots.size()) k = source_->knots.size() - 1;
      const auto c = between(source_->knots[k - 1], source_->knots[k], ell, W(z));
      const W f = source_->helium_to_hydrogen_nuclei_ratio,
          nHe = c.nH * f, ne = detail::hhe_electron_density(c.nH, 1 - c.xHI, nHe, q),
          ene = nHe * eq + floor64 * (ne + c.nH),
          A = detail::hhe_opacity_coefficient(c.H, 1 + W(z)),
          opacity = A * ne, eopacity = A * ene + floor64 * opacity,
          tau = detail::helium_escape_optical_depth(c.H, c.nH, f, q),
          activity = detail::helium_escape_heiii_log_activity(c.T, ne);
      if (!physical_fraction(q) || !(f * q >= 1e-6L) || !(tau >= 100) ||
          !std::isfinite(activity) || activity > std::log(1e-12L)) {
        refuse(S::outside_domain); continue;
      }
      if (mask & 1) project(*values[0], q, eq, policy_.absolute_fraction_tolerance,
                            policy_.relative_fraction_tolerance, true);
      if (mask & 2) project(*values[1], ne, ene, nHe * policy_.absolute_fraction_tolerance,
                            policy_.relative_fraction_tolerance);
      if (mask & 4) project(*values[2], opacity, eopacity, policy_.absolute_opacity_tolerance,
                            policy_.relative_opacity_tolerance);
    }
    out.status = S::ok;
  } catch (const std::bad_alloc &) { out.status = S::work_limit; std::vector<HeliumEscapeHistoryRow>{}.swap(out.rows); }
  catch (const std::length_error &) { out.status = S::work_limit; std::vector<HeliumEscapeHistoryRow>{}.swap(out.rows); }
  return out;
}
} // namespace irred::cosmology
