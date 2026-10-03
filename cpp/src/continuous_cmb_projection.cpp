#include "irred/continuous_cmb_projection.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace irred::projection {
namespace {
using S = numerics::Status;
using W = long double;
using V = std::array<W, 4>;
constexpr W eps = std::numeric_limits<W>::epsilon();
bool profile() {
  return std::numeric_limits<W>::digits >= 64 &&
         std::fegetround() == FE_TONEAREST;
}
std::optional<std::size_t> source_bytes(const ContinuousCmbSource &s) {
  irred::detail::PayloadAccounting b(sizeof(ContinuousCmbProjection) +
                                    sizeof(ContinuousCmbSource));
  for (const auto *v : {&s.k_mpc_inverse, &s.eta_mpc, &s.t0, &s.t1,
                        &s.t2, &s.polarization})
    b.add(v->capacity(), sizeof(double));
  for (const auto *id : {&s.producer_id, &s.signed_mode_id, &s.normalization_id})
    b.add(id->capacity(), 1);
  b.add(3, 1); // terminating bytes; allocator/control-block metadata excluded
  return b.result();
}
bool finite(W v) { return std::isfinite(v); }
bool usable(W v) {
  return finite(v) && (v == 0 || std::abs(v) >= std::numeric_limits<W>::min());
}
struct Kernel { S status = S::ok; V value{}; };
Kernel radial(unsigned l, W x) {
  Kernel out;
  const W spin = l < 2 ? 0 : std::sqrt(3.L * (l - 1) * l * (l + 1) *
                                      (l + 2) / 8);
  if (x == 0) {
    out.value = {l == 0 ? 1.L : 0.L, l == 1 ? 1.L / 3 : 0.L,
                 l == 2 ? 1.L / 5 : 0.L, l == 2 ? 1.L / 5 : 0.L};
    return out;
  }
  if (x <= .5L) {
    // Terms x^(l+2m)/(2l+1)!! * product[-x²/(2m(2l+2m+1))].
    // Derivative kernels are formed termwise, never dividing a rounded j by x².
    W leading = 1;
    for (unsigned n = 1; n <= l; ++n) leading *= x / (2 * n + 1);
    W term = leading;
    if (leading == 0) { out.status = S::outside_domain; return out; }
    W j = 0, first = 0, quadrupole = 0, spin_value = 0;
    if (l >= 2) quadrupole = 1.5L * l * (l - 1) * ((leading / x) / x);
    for (unsigned m = 0; m < 32; ++m) {
      const unsigned n = l + 2 * m;
      j += term;
      if (n) first += n * (term / x);
      // Align the next second-derivative term with this j term. In particular
      // the exact l=0,m=0 coefficient is zero, avoiding a false residual .5-.5
      // that would hide the regular x²/15 quadrupole at a tiny positive phase.
      quadrupole += term / 2 *
          (1 - 3 * W(n + 2) * (n + 1) / (2 * W(m + 1) * (2 * l + 2 * m + 3)));
      if (l >= 2) spin_value += spin * ((term / x) / x);
      const W next = term * (-x * x) / (2 * W(m + 1) * (2 * l + 2 * m + 3));
      if (m >= 3 && std::abs(next) <= eps * std::abs(leading) / 64) break;
      term = next;
    }
    out.value = {j, first, quadrupole, spin_value};
  } else {
    // Portable baseline. Independent amplitude/phase controls are required;
    // ODE or recurrence residuals alone are not a libm accuracy certificate.
    try {
      const W j = std::sph_bessel(l, x), next = std::sph_bessel(l + 1, x);
      const W first = l * j / x - next;
      const W second = (W(l) * (W(l) - 1) / (x * x) - 1) * j + 2 * next / x;
      out.value = {j, first, (3 * second + j) / 2, spin * j / (x * x)};
    } catch (const std::exception &) { out.status = S::outside_domain; }
  }
  for (W v : out.value) if (!usable(v)) out.status = S::outside_domain;
  return out;
}
struct Integral { S status = S::ok; V value{}, error{}, absolute{}; };
struct Node { W eta = 0; V value{}, absolute{}; };
struct Frame { Node a, m, b; V coarse{}; W absolute_tolerance = 0; };
struct Traversal {
  const ContinuousCmbSource &source;
  const ContinuousCmbPolicy &policy;
  std::size_t k_index, cell;
  unsigned ell;
  unsigned mask;
  std::size_t &all_calls, &row_calls;
  S status = S::ok;
  Node node(W eta) {
    Node out{eta, {}};
    if (status != S::ok) return out;
    if (all_calls >= policy.maximum_kernel_evaluations) {
      status = S::work_limit; return out;
    }
    ++all_calls; ++row_calls; // includes any subsequently refused radial query
    const W distance = W(source.observer_eta_mpc) - eta;
    const W x = W(source.k_mpc_inverse[k_index]) * distance;
    if (distance < 0 || (distance > 0 && (x == 0 || !usable(x)))) {
      status = S::outside_domain; return out;
    }
    const auto kernel = radial(ell, x);
    if (kernel.status != S::ok) { status = kernel.status; return out; }
    const W u = (eta - source.eta_mpc[cell]) /
                (W(source.eta_mpc[cell + 1]) - source.eta_mpc[cell]);
    if (!(u >= 0 && u <= 1)) { status = S::outside_domain; return out; }
    const std::size_t ia = cell * source.k_mpc_inverse.size() + k_index;
    const std::size_t ib = ia + source.k_mpc_inverse.size();
    const std::array<const std::vector<double> *, 4> channels{
        &source.t0, &source.t1, &source.t2, &source.polarization};
    for (unsigned n = 0; n < 4; ++n) {
      if ((n < 3 && !(mask & continuous_temperature)) ||
          (n == 3 && !(mask & continuous_e_mode))) continue;
      const W a = (*channels[n])[ia], b = (*channels[n])[ib];
      const W source_value = (1 - u) * a + u * b;
      const W source_scale = (1 - u) * std::abs(a) + u * std::abs(b);
      out.value[n] = source_value * kernel.value[n];
      out.absolute[n] = source_scale * std::abs(kernel.value[n]);
      if (!usable(out.value[n]) || !usable(out.absolute[n]) ||
          (source_value != 0 && kernel.value[n] != 0 && out.value[n] == 0) ||
          (source_scale != 0 && kernel.value[n] != 0 && out.absolute[n] == 0))
        status = S::outside_domain;
    }
    return out;
  }
  V simpson(const Node &a, const Node &m, const Node &b) {
    V out{};
    for (unsigned n = 0; n < 4; ++n) {
      const W weighted = a.value[n] + 4 * m.value[n] + b.value[n];
      out[n] = (b.eta - a.eta) * weighted / 6;
      if (!usable(out[n]) || (weighted != 0 && out[n] == 0)) status = S::outside_domain;
    }
    return out;
  }
  V assembly_scale(const Node &a, const Node &m, const Node &b) {
    V out{};
    for (unsigned n = 0; n < 4; ++n) {
      const W weighted = a.absolute[n] + 4 * m.absolute[n] + b.absolute[n];
      out[n] = (b.eta - a.eta) * weighted / 6;
      if (!usable(out[n]) || (weighted > 0 && out[n] == 0)) status = S::outside_domain;
    }
    return out;
  }
  Integral visit(const Frame &f, unsigned depth) {
    Integral out;
    if (status != S::ok) { out.status = status; return out; }
    const auto left = node((f.a.eta + f.m.eta) / 2);
    const auto right = node((f.m.eta + f.b.eta) / 2);
    if (status != S::ok) { out.status = status; return out; }
    if (!(left.eta > f.a.eta && left.eta < f.m.eta &&
          right.eta > f.m.eta && right.eta < f.b.eta) ||
        left.eta - f.a.eta != f.m.eta - left.eta ||
        right.eta - f.m.eta != f.b.eta - right.eta ||
        f.m.eta - f.a.eta != f.b.eta - f.m.eta) {
      out.status = S::conditioning_budget_exceeded; return out;
    }
    const V lo = simpson(f.a, left, f.m), hi = simpson(f.m, right, f.b);
    const V lo_scale = assembly_scale(f.a, left, f.m);
    const V hi_scale = assembly_scale(f.m, right, f.b);
    if (status != S::ok) { out.status = status; return out; }
    W et = 0, ee = 0, vt = 0;
    for (unsigned n = 0; n < 4; ++n) {
      out.value[n] = lo[n] + hi[n];
      out.error[n] = std::abs(out.value[n] - f.coarse[n]) / 15;
      out.absolute[n] = lo_scale[n] + hi_scale[n];
      if (n < 3) { et += out.error[n]; vt += out.value[n]; }
      else ee = out.error[n];
    }
    const bool phase = W(source.k_mpc_inverse[k_index]) * (f.b.eta - f.a.eta) <= .5L;
    const bool t = !(mask & continuous_temperature) ||
        et <= f.absolute_tolerance + policy.relative_tolerance * std::abs(vt);
    const bool e = !(mask & continuous_e_mode) ||
        ee <= f.absolute_tolerance + policy.relative_tolerance * std::abs(out.value[3]);
    if (phase && t && e) return out;
    if (depth >= policy.maximum_depth || left.eta == f.a.eta || right.eta == f.b.eta) {
      out.status = S::conditioning_budget_exceeded; return out;
    }
    const auto l = visit({f.a, left, f.m, lo, f.absolute_tolerance / 2}, depth + 1);
    if (l.status != S::ok) return l;
    const auto r = visit({f.m, right, f.b, hi, f.absolute_tolerance / 2}, depth + 1);
    if (r.status != S::ok) return r;
    for (unsigned n = 0; n < 4; ++n) {
      out.value[n] = l.value[n] + r.value[n];
      out.error[n] = l.error[n] + r.error[n];
      out.absolute[n] = l.absolute[n] + r.absolute[n];
    }
    return out;
  }
};
void add_compensated(W v, W &sum, W &correction) {
  const W t = sum + v;
  correction += std::abs(sum) >= std::abs(v) ? (sum - t) + v : (v - t) + sum;
  sum = t;
}
} // namespace
struct ContinuousCmbProjectionAccess {
  static auto source(const ContinuousCmbProjection &p) { return p.source_; }
};
ContinuousCmbProjection::ContinuousCmbProjection(ContinuousCmbProjection &&p) noexcept {
  *this = std::move(p);
}
ContinuousCmbProjection &ContinuousCmbProjection::operator=(ContinuousCmbProjection &&p) noexcept {
  if (this != &p) {
    status_ = p.status_; source_ = std::move(p.source_); p.status_ = S::invalid_input;
  }
  return *this;
}
std::optional<std::size_t> ContinuousCmbProjection::retained_payload_bound() const noexcept {
  return source_ ? source_bytes(*source_) : std::optional<std::size_t>{sizeof(*this)};
}
ContinuousCmbProjection prepare_continuous_cmb_projection(
    ContinuousCmbSource &&s, ContinuousCmbPreparationPolicy p) {
  ContinuousCmbProjection out;
  if (!profile()) return out;
  if (s.k_mpc_inverse.empty() || s.eta_mpc.size() < 2 || s.producer_id.empty() ||
      s.signed_mode_id.empty() || s.normalization_id.empty()) return out;
  if (s.k_mpc_inverse.size() > p.maximum_k || s.eta_mpc.size() > p.maximum_eta ||
      s.eta_mpc.size() > p.maximum_cells / s.k_mpc_inverse.size()) {
    out.status_ = S::work_limit; return out;
  }
  const std::size_t cells = s.eta_mpc.size() * s.k_mpc_inverse.size();
  for (const auto *v : {&s.t0, &s.t1, &s.t2, &s.polarization}) {
    if (v->size() != cells) return out;
    for (double x : *v) if (!std::isfinite(x)) { out.status_ = S::nonfinite_input; return out; }
  }
  if (!std::isfinite(s.observer_eta_mpc)) { out.status_ = S::nonfinite_input; return out; }
  for (const auto *axis : {&s.k_mpc_inverse, &s.eta_mpc})
    for (std::size_t i = 0; i < axis->size(); ++i) {
      if (!std::isfinite((*axis)[i])) { out.status_ = S::nonfinite_input; return out; }
      if ((i && !((*axis)[i] > (*axis)[i - 1])) || (*axis)[i] < 0 ||
          (axis == &s.k_mpc_inverse && (*axis)[i] == 0)) {
        out.status_ = S::outside_domain; return out;
      }
    }
  if (s.observer_eta_mpc < s.eta_mpc.back() ||
      W(s.k_mpc_inverse.back()) * (W(s.observer_eta_mpc) - s.eta_mpc.front()) > 512) {
    out.status_ = S::outside_domain; return out;
  }
  const auto bytes = source_bytes(s);
  if (!bytes || *bytes > p.maximum_payload_bytes) { out.status_ = S::work_limit; return out; }
  static_assert(std::is_nothrow_move_constructible_v<ContinuousCmbSource>);
  try { out.source_ = std::make_shared<const ContinuousCmbSource>(std::move(s)); }
  catch (const std::bad_alloc &) { out.status_ = S::work_limit; return out; }
  out.status_ = S::ok;
  return out;
}
ContinuousCmbResult project_continuous_cmb(
    const ContinuousCmbProjection &p, std::span<const unsigned> multipoles,
    unsigned mask, ContinuousCmbPolicy policy) {
  ContinuousCmbResult out;
  out.requested_outputs = mask;
  out.source_owner = ContinuousCmbProjectionAccess::source(p);
  if (!profile() || p.status() != S::ok || !out.source_owner || !mask ||
      (mask & ~(continuous_temperature | continuous_e_mode)) || multipoles.empty() ||
      !std::isfinite(policy.absolute_tolerance) || !std::isfinite(policy.relative_tolerance) ||
      policy.absolute_tolerance < 0 || policy.relative_tolerance < 0 ||
      (policy.absolute_tolerance == 0 && policy.relative_tolerance == 0) || policy.maximum_depth > 20)
    return out;
  const auto &s = *out.source_owner;
  if (multipoles.size() > policy.maximum_multipoles ||
      multipoles.size() > std::numeric_limits<std::size_t>::max() / s.k_mpc_inverse.size()) {
    out.status = S::work_limit; return out;
  }
  for (unsigned l : multipoles)
    if (l > 64 || ((mask & continuous_e_mode) && l < 2)) { out.status = S::outside_domain; return out; }
  const std::size_t count = multipoles.size() * s.k_mpc_inverse.size();
  const auto retained = p.retained_payload_bound();
  irred::detail::PayloadAccounting bytes(retained.value_or(std::numeric_limits<std::size_t>::max()));
  bytes.add(1, sizeof(ContinuousCmbResult)); bytes.add(count, sizeof(ContinuousCmbRow));
  // Fixed numerical objects, including recursive traversal frames and returned
  // child integrals. Allocator metadata/RSS/incidental call ABI are excluded.
  bytes.add(policy.maximum_depth + 2, sizeof(Frame) + 3 * sizeof(Integral) + 4 * sizeof(Node));
  bytes.add(1, sizeof(Traversal) + 8 * sizeof(V) + 2 * sizeof(ContinuousCmbRow));
  const auto peak = bytes.result();
  if (!retained || !peak || *peak > policy.maximum_payload_bytes) { out.status = S::work_limit; return out; }
  out.payload_bound = *peak;
  try { out.rows.resize(count); } catch (const std::bad_alloc &) { out.status = S::work_limit; return out; }
  if (out.rows.capacity() > count) {
    bytes.add(out.rows.capacity() - count, sizeof(ContinuousCmbRow));
    const auto actual = bytes.result();
    if (!actual || *actual > policy.maximum_payload_bytes) { out.status = S::work_limit; return out; }
    out.payload_bound = *actual;
  }
  out.status = S::ok;
  for (std::size_t li = 0; li < multipoles.size(); ++li)
    for (std::size_t ki = 0; ki < s.k_mpc_inverse.size(); ++ki) {
      auto &row = out.rows[li * s.k_mpc_inverse.size() + ki];
      row.multipole_index = li; row.k_index = ki; row.ell = multipoles[li]; row.k_mpc_inverse = s.k_mpc_inverse[ki];
      Integral total; V correction{};
      for (std::size_t cell = 0; cell + 1 < s.eta_mpc.size(); ++cell) {
        Traversal tr{s, policy, ki, cell, row.ell, mask, out.kernel_evaluations, row.kernel_evaluations};
        const auto a = tr.node(s.eta_mpc[cell]);
        const auto m = tr.node((W(s.eta_mpc[cell]) + s.eta_mpc[cell + 1]) / 2);
        const auto b = tr.node(s.eta_mpc[cell + 1]);
        const auto v = tr.visit({a, m, b, tr.simpson(a, m, b),
                                W(policy.absolute_tolerance) / (s.eta_mpc.size() - 1)}, 0);
        if (v.status != S::ok) { total.status = v.status; break; }
        ++row.completed_source_cells;
        for (unsigned n = 0; n < 4; ++n) {
          add_compensated(v.value[n], total.value[n], correction[n]);
          total.error[n] += v.error[n]; total.absolute[n] += v.absolute[n];
        }
      }
      row.status = total.status;
      W t = 0, te = 0, ta = 0;
      for (unsigned n = 0; n < 4; ++n) total.value[n] += correction[n];
      for (unsigned n = 0; n < 3; ++n) { t += total.value[n]; te += total.error[n]; ta += total.absolute[n]; }
      // On traversal refusal these are estimates for completed cells only;
      // completed_source_cells distinguishes that prefix from a full transfer.
      row.temperature_quadrature_estimate = static_cast<double>(te);
      row.e_quadrature_estimate = static_cast<double>(total.error[3]);
      auto store = [&](W value, W qe, W norm, std::optional<double> &dest,
                       double &quad, double &arith) {
        const double cast = static_cast<double>(value);
        // Positive assembly scale survives interpolation/Simpson cancellation.
        // Actual row calls enter the conservative accumulation allowance.
        const W accumulation = (128 + 16 * W(row.kernel_evaluations)) * eps * norm;
        const W ae = accumulation + std::abs(value - W(cast));
        quad = static_cast<double>(qe); arith = static_cast<double>(ae);
        if (!usable(value) || !std::isfinite(cast) || (value != 0 && cast == 0) ||
            !usable(qe) || !usable(ae) ||
            (norm > 0 && (accumulation == 0 || !usable(accumulation))) ||
            !std::isfinite(quad) || !std::isfinite(arith) ||
            (qe > 0 && quad == 0) || (ae > 0 && arith == 0)) { row.status = S::outside_domain; return; }
        if (qe + ae > policy.absolute_tolerance + policy.relative_tolerance * std::abs(value)) {
          row.status = S::conditioning_budget_exceeded; return;
        }
        dest = cast;
      };
      if (row.status == S::ok) {
        if (mask & continuous_temperature) store(t, te, ta, row.temperature, row.temperature_quadrature_estimate, row.temperature_arithmetic_estimate);
        if (mask & continuous_e_mode) store(total.value[3], total.error[3], total.absolute[3], row.e_mode, row.e_quadrature_estimate, row.e_arithmetic_estimate);
      }
      if (row.status != S::ok) {
        row.temperature.reset(); row.e_mode.reset();
        if (out.status == S::ok) out.status = row.status;
      }
    }
  return out;
}
} // namespace irred::projection
