#include "irred/hydrogen_helium_history.hpp"
#include "hydrogen_helium_rates.hpp"
#include "hydrogen_helium_cell.hpp"
#include "payload_accounting.hpp"
#include "thermal_constants.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <numbers>
namespace irred::cosmology {
namespace {
using W = long double;
using S = numerics::Status;
using State = std::array<W, 3>;
constexpr W arithmetic_floor = 128 * std::numeric_limits<double>::epsilon();
struct Coeff { W z, u, nH, nHe, H, Tr, compton, opacity, H_error; };
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
bool physical(const State &x) {
  return x[0] > 0 && x[0] < 1 && x[1] > 0 && x[1] < 1 &&
         x[2] > 0 && std::isfinite(x[2]);
}
bool budget(double a, double r) {
  return std::isfinite(a) && std::isfinite(r) && a >= 0 && r >= 0 && (a > 0 || r > 0);
}
std::size_t remaining_work(const HydrogenHeliumHistoryWork &work, std::size_t cap) {
  const auto spent = work.total();
  return spent < cap ? cap - spent : 0;
}
bool charge_momentum(HydrogenHeliumHistoryWork &work, std::size_t count,
                     std::size_t cap) {
  if (count > remaining_work(work, cap) ||
      count > SIZE_MAX - work.momentum_callbacks) return false;
  work.momentum_callbacks += count;
  return true;
}
void retain_range(std::optional<std::array<W, 2>> &range, W T) {
  if (!range) range = std::array<W, 2>{T, T};
  else {
    (*range)[0] = std::min((*range)[0], T);
    (*range)[1] = std::max((*range)[1], T);
  }
}
W heiii_log(W T, W ne) {
  const W kT = boltzmann_constant_joule_per_kelvin * T;
  return std::log(atomic::detail::electron_quantum_density_si(kT)) -
         atomic::helium_second_ionization_energy_ev * electron_volt_joule / kT -
         std::log(ne);
}
// The restricted two-stage equilibrium has its own model identity, but shares
// every Saha coefficient/asset with the ground-state atomic operator.
State initial(const Coeff &c, HydrogenHeliumHistoryWork &work,
              std::size_t cap, S &status) {
  const W kT = boltzmann_constant_joule_per_kelvin * c.Tr,
      Q = atomic::detail::electron_quantum_density_si(kT),
      a = Q * std::exp(-atomic::hydrogen_ionization_energy_ev * electron_volt_joule / kT),
      b = 4 * Q * std::exp(-atomic::helium_first_ionization_energy_ev * electron_volt_joule / kT);
  W low = 0, high = c.nH + c.nHe;
  for (unsigned iteration = 0; iteration < 128; ++iteration) {
    if (work.total() >= cap) { status = S::work_limit; return {}; }
    ++work.initial_charge_evaluations;
    const W ne = (low + high) / 2,
        charge = c.nH * a / (ne + a) + c.nHe * b / (ne + b);
    if (ne < charge) low = ne; else high = ne;
    if (high - low <= 8 * std::numeric_limits<W>::epsilon() * high) {
      const W e = (low + high) / 2;
      return {a / (e + a), b / (e + b), c.Tr};
    }
  }
  status = S::conditioning_budget_exceeded; return {};
}
struct Evolution {
  S status = S::ok;
  HydrogenHeliumHistoryWork &work;
  std::size_t cap;
  HydrogenHeliumRateDomainWitness &domain;
  State accumulated_root_error{};
  struct Residual {
    std::array<W, 2> r;
    std::array<std::array<W, 2>, 2> J;
    State state;
    std::array<W, 2> temperature_derivative{};
  };
  Residual residual(const Coeff &c, State previous, W h,
                    const std::array<W, 2> &x) {
    if (work.total() >= cap) { status = S::work_limit; return {}; }
    ++work.rhs_evaluations;
    const W ne = c.nH * x[0] + c.nHe * x[1], M = c.nH + c.nHe,
        A = c.compton * ne / (M + ne), den = 1 + h * A + 2 * h / c.u,
        T = (previous[2] + h * A * c.Tr) / den;
    if (!(T > 0) || !std::isfinite(T)) {
      domain.invalid_temperature_attempted = true;
      status = S::conditioning_budget_exceeded; return {};
    }
    // One observation site covers all meshes and accepted/failed Newton trials.
    retain_range(domain.attempted_kelvin_range, T);
    const auto f = detail::hhe_rhs(x[0], x[1], T, c.nH, c.nHe, c.H, c.u);
    const W n[]{c.nH, c.nHe};
    Residual out{{x[0] - previous[0] + h * f.value[0],
                  x[1] - previous[1] + h * f.value[1]}, {}, {x[0], x[1], T}};
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 2; ++j) {
        const W Tj = h * c.compton * M * n[j] * (c.Tr - T) /
                     ((M + ne) * (M + ne) * den);
        out.temperature_derivative[j] = Tj;
        out.J[i][j] = (i == j ? 1 : 0) +
            h * (f.fraction_derivative[i][j] + f.temperature_derivative[i] * Tj);
      }
    return out;
  }
  bool correction(const Residual &r, std::array<W, 2> &d) {
    const auto &j = r.J;
    const W det = j[0][0] * j[1][1] - j[0][1] * j[1][0],
        norm = std::max(std::abs(j[0][0]) + std::abs(j[0][1]),
                        std::abs(j[1][0]) + std::abs(j[1][1])),
        adj = std::max(std::abs(j[1][1]) + std::abs(j[0][1]),
                       std::abs(j[1][0]) + std::abs(j[0][0]));
    if (!(det > 0) || !std::isfinite(det) ||
        !(norm * adj / det < 1e18L)) return false;
    d = {(-j[1][1] * r.r[0] + j[0][1] * r.r[1]) / det,
         (j[1][0] * r.r[0] - j[0][0] * r.r[1]) / det};
    return std::isfinite(d[0]) && std::isfinite(d[1]);
  }
  State step(const Coeff &c, State previous, W h) {
    std::array<W, 2> x{previous[0], previous[1]};
    for (unsigned trial = 0; trial < 80; ++trial) {
      const auto r = residual(c, previous, h, x);
      if (status != S::ok) return {};
      std::array<W, 2> d{};
      if (!correction(r, d)) { status = S::conditioning_budget_exceeded; return {}; }
      // The resolvable root correction, rather than a raw stiff residual, is
      // the dimensionless convergence criterion. Both coupled rows are used.
      const W norm = std::max(std::abs(d[0]) / x[0], std::abs(d[1]) / x[1]);
      if (norm <= 64 * std::numeric_limits<W>::epsilon()) {
        accumulated_root_error[0] += std::abs(d[0]);
        accumulated_root_error[1] += std::abs(d[1]);
        accumulated_root_error[2] += std::abs(r.temperature_derivative[0] * d[0]) +
                                    std::abs(r.temperature_derivative[1] * d[1]);
        return r.state;
      }
      bool accepted = false;
      W damping = 1;
      for (unsigned k = 0; k < 64; ++k) {
        const std::array<W, 2> next{x[0] + damping * d[0], x[1] + damping * d[1]};
        if (next[0] > 0 && next[0] < 1 && next[1] > 0 && next[1] < 1) {
          const auto rr = residual(c, previous, h, next);
          if (status != S::ok) return {};
          std::array<W, 2> dd{};
          if (correction(rr, dd) &&
              std::max(std::abs(dd[0]) / next[0], std::abs(dd[1]) / next[1]) < norm) {
            x = next; accepted = true; break;
          }
        }
        damping /= 2;
      }
      if (!accepted) { status = S::conditioning_budget_exceeded; return {}; }
    }
    status = S::conditioning_budget_exceeded; return {};
  }
};
void project(HydrogenHeliumHistoryValue &v, W value, W error, W a, W r) {
  v.status = S::conditioning_budget_exceeded;
  if (!(value > 0) || !std::isfinite(value) || !(error >= 0) || !std::isfinite(error)) return;
  const double rounded = static_cast<double>(value);
  error += std::abs(value - W(rounded)) + arithmetic_floor * value;
  if (!(rounded > 0) || !std::isfinite(rounded) ||
      error > W(std::numeric_limits<double>::max())) return;
  v.absolute_error_estimate = std::nextafter(static_cast<double>(error), INFINITY);
  if (error > a + r * std::abs(value)) return;
  v.status = S::ok; v.value = rounded;
}
} // namespace
std::optional<std::size_t> hydrogen_helium_history_payload_bound(
    std::size_t fine, std::size_t output, std::size_t species,
    std::size_t origin, std::size_t supplied_origin,
    std::size_t supplied_identity) noexcept {
  irred::detail::PayloadAccounting p(sizeof(HydrogenHeliumHistory) +
                                    sizeof(HydrogenHeliumHistoryBatch) + 4096);
  if (fine == SIZE_MAX || origin == SIZE_MAX || supplied_origin == SIZE_MAX ||
      supplied_identity == SIZE_MAX) return {};
  // Simultaneous coefficient, output and three BE/refined state buffers.
  p.add(fine + 1, sizeof(Coeff) + sizeof(HydrogenHeliumHistory::Node) +
                         4 * sizeof(State));
  p.add(output, sizeof(HydrogenHeliumHistoryRow));
  p.add(species, 4 * (sizeof(ThermalSpecies) + sizeof(ThermalPhysicalSpecies)));
  p.add(origin + 1, 1);
  if (supplied_origin) p.add(supplied_origin + 1, 1);
  if (supplied_identity) p.add(supplied_identity + 1, 1);
  return p.result();
}
HydrogenHeliumHistory &HydrogenHeliumHistory::operator=(const HydrogenHeliumHistory &o) {
  if (this != &o) { HydrogenHeliumHistory copy(o); *this = std::move(copy); }
  return *this;
}
HydrogenHeliumHistory::HydrogenHeliumHistory(HydrogenHeliumHistory &&o) noexcept {
  *this = std::move(o);
}
HydrogenHeliumHistory &HydrogenHeliumHistory::operator=(HydrogenHeliumHistory &&o) noexcept {
  if (this != &o) {
    status_ = o.status_; source_ = std::move(o.source_);
    boundary_ = o.boundary_; supplied_initial_ = std::move(o.supplied_initial_);
    initial_import_ = o.initial_import_; rate_domain_ = o.rate_domain_;
    background_ = std::move(o.background_); policy_ = o.policy_; work_ = o.work_;
    mapping_witnesses_ = std::move(o.mapping_witnesses_);
    excluded_heiii_activity_ = o.excluded_heiii_activity_;
    maximum_log_heiii_activity_ = o.maximum_log_heiii_activity_;
    nodes_ = std::move(o.nodes_);
    o.source_.reset(); o.status_ = S::invalid_input; o.work_ = {};
    o.boundary_.reset(); o.supplied_initial_.reset(); o.initial_import_.reset();
    o.rate_domain_.reset();
    o.mapping_witnesses_.reset();
    o.excluded_heiii_activity_.reset(); o.maximum_log_heiii_activity_.reset();
    o.nodes_.clear();
  }
  return *this;
}
std::string_view HydrogenHeliumHistory::model_identity() const noexcept {
  if (!boundary_) return {};
  return *boundary_ == HydrogenHeliumHistoryBoundary::supplied_binary64
             ? hydrogen_helium_supplied_history_model_id
             : hydrogen_helium_history_model_id;
}
namespace detail {
struct HydrogenHeliumPreparationAccess {
static HydrogenHeliumHistory prepare(const HydrogenHeliumHistoryRequest &s,
    const HydrogenHeliumSuppliedInitialState *supplied, HydrogenHeliumHistoryPolicy p) {
  HydrogenHeliumHistory out;
  out.boundary_ = supplied ? HydrogenHeliumHistoryBoundary::supplied_binary64
                          : HydrogenHeliumHistoryBoundary::restricted_two_stage_saha;
  if (!arithmetic() || !budget(p.absolute_fraction_tolerance, p.relative_fraction_tolerance) ||
      !budget(p.absolute_temperature_tolerance_kelvin, p.relative_temperature_tolerance) ||
      !budget(p.absolute_opacity_tolerance, p.relative_opacity_tolerance) ||
      p.base_intervals < 2 || s.nuclei_origin.empty()) return out;
  if (supplied && (supplied->origin.empty() || supplied->source_identity.empty())) return out;
  if (p.base_intervals > 16384 || 4 * p.base_intervals > p.maximum_fine_intervals ||
      s.model.species.size() > 16 || s.nuclei_origin.size() > 65536) {
    out.status_ = S::work_limit; return out;
  }
  if (supplied && (supplied->origin.size() > 65536 - s.nuclei_origin.size() ||
      supplied->source_identity.size() >
          65536 - s.nuclei_origin.size() - supplied->origin.size())) {
    out.status_ = S::work_limit; return out;
  }
  const std::size_t N = 4 * p.base_intervals;
  auto bytes = hydrogen_helium_history_payload_bound(N, 0, s.model.species.size(),
      s.nuclei_origin.size(), supplied ? supplied->origin.size() : 0,
      supplied ? supplied->source_identity.size() : 0);
  if (!bytes || *bytes > p.maximum_native_bytes || *bytes > (1ull << 30)) {
    out.status_ = S::work_limit; return out;
  }
  for (double v : {s.initial_redshift, s.late_redshift,
       s.hydrogen_nuclei_today_per_cubic_metre, s.helium_nuclei_today_per_cubic_metre,
       s.model.h0_km_s_mpc, s.model.physical_baryon_density, s.model.physical_cdm_density,
       s.model.tcmb_kelvin, s.model.physical_massless_nonphoton_density})
    if (!std::isfinite(v)) { out.status_ = S::nonfinite_input; return out; }
  const double nH0 = s.hydrogen_nuclei_today_per_cubic_metre,
      nHe0 = s.helium_nuclei_today_per_cubic_metre;
  if (s.initial_redshift < 2600 || s.initial_redshift > 2800 ||
      s.late_redshift < 300 || s.late_redshift > 600 ||
      nH0 < .1 || nH0 > .3 || nHe0 < .005 || nHe0 > .035 ||
      nHe0 / nH0 < .04 || nHe0 / nH0 > .12 ||
      s.model.h0_km_s_mpc < 60 || s.model.h0_km_s_mpc > 80 ||
      s.model.physical_baryon_density < .015 || s.model.physical_baryon_density > .03 ||
      s.model.physical_cdm_density < .08 || s.model.physical_cdm_density > .15 ||
      s.model.tcmb_kelvin < 2.7 || s.model.tcmb_kelvin > 2.75 ||
      s.model.physical_massless_nonphoton_density < 0 ||
      s.model.physical_massless_nonphoton_density > 3e-5) {
    out.status_ = S::outside_domain; return out;
  }
  for (const auto &v : s.model.species) {
    if (!std::isfinite(v.mass_ev) || !std::isfinite(v.temperature_today_kelvin) ||
        !std::isfinite(v.statistical_weight)) { out.status_ = S::nonfinite_input; return out; }
    // This first new history profile qualifies only the massless thermal state.
    if (v.mass_ev != 0) { out.status_ = S::outside_domain; return out; }
  }
  if (supplied) {
    if (!remaining_work(out.work_, p.maximum_total_work)) {
      out.status_ = S::work_limit; return out;
    }
    ++out.work_.initial_boundary_evaluations;
    for (double value : {supplied->hydrogen_ionized_fraction,
                        supplied->helium_singly_ionized_fraction,
                        supplied->matter_temperature_kelvin})
      if (!std::isfinite(value)) { out.status_ = S::nonfinite_input; return out; }
    const W Tr = W(s.model.tcmb_kelvin) * (1 + W(s.initial_redshift));
    if (!(supplied->hydrogen_ionized_fraction > 0 &&
          supplied->hydrogen_ionized_fraction < 1 &&
          supplied->helium_singly_ionized_fraction > 0 &&
          supplied->helium_singly_ionized_fraction < 1 &&
          supplied->matter_temperature_kelvin >= 4000 &&
          W(supplied->matter_temperature_kelvin) <= Tr)) {
      out.status_ = S::outside_domain; return out;
    }
  }
  try {
    // Acquire this source once before any thermal/background work. A failed
    // emplace publishes no partial boundary or import witness.
    if (supplied) {
      out.supplied_initial_.emplace(*supplied);
      out.initial_import_ = HydrogenHeliumInitialImportWitness{
          {W(supplied->hydrogen_ionized_fraction),
           W(supplied->helium_singly_ionized_fraction),
           W(supplied->matter_temperature_kelvin)}, {0, 0, 0}};
    }
    std::size_t accounted_fine = N, mapped_species_capacity = s.model.species.size();
    auto payload_ok = [&]() {
      const auto species = std::max({mapped_species_capacity,
          out.background_.source().species.capacity(),
          out.source_ ? out.source_->model.species.capacity() : std::size_t(0)});
      const auto bound = hydrogen_helium_history_payload_bound(accounted_fine, 0,
          species, out.source_ ? out.source_->nuclei_origin.capacity() : s.nuclei_origin.size(),
          out.supplied_initial_ ? out.supplied_initial_->origin.capacity() : 0,
          out.supplied_initial_ ? out.supplied_initial_->source_identity.capacity() : 0);
      return bound && *bound <= p.maximum_native_bytes && *bound <= (1ull << 30);
    };
    if (!payload_ok()) { out.status_ = S::work_limit; return out; }
    if (supplied && !remaining_work(out.work_, p.maximum_total_work)) {
      out.status_ = S::work_limit; return out;
    }
    auto t = p.thermal;
    t.maximum_native_bytes = std::min(t.maximum_native_bytes, p.maximum_native_bytes);
    t.maximum_total_callbacks = std::min(t.maximum_total_callbacks,
                                         remaining_work(out.work_, p.maximum_total_work));
    const auto mapped = map_thermal_physical_model(s.model, t);
    if (!mapped.model) { out.status_ = mapped.status; return out; }
    out.mapping_witnesses_ = mapped.scalar_witnesses;
    mapped_species_capacity = mapped.model->species.capacity();
    if (!payload_ok()) { out.status_ = S::work_limit; return out; }
    out.background_ = prepare_thermal_background(*mapped.model, t);
    if (!charge_momentum(out.work_, out.background_.preparation_callbacks(), p.maximum_total_work)) {
      out.status_ = S::work_limit; return out;
    }
    if (out.background_.status() != S::ok) { out.status_ = out.background_.status(); return out; }
    out.source_ = s; out.policy_ = p;
    if (!payload_ok()) { out.status_ = S::work_limit; return out; }
    std::vector<Coeff> c(N + 1);
    accounted_fine = std::max(accounted_fine, c.capacity() - 1);
    if (!payload_ok()) { out.status_ = S::work_limit; return out; }
    W max_H_error = 0;
    const W span = W(s.initial_redshift) - s.late_redshift,
        photon = W(mapped.model->omega_gamma) *
                 detail::critical_energy_density_si(s.model.h0_km_s_mpc);
    for (std::size_t i = 0; i <= N; ++i) {
      if (out.work_.total() >= p.maximum_total_work) { out.status_ = S::work_limit; return out; }
      ++out.work_.background_evaluations;
      const W coordinate = W(i) / N,
          z = i == N ? W(s.late_redshift) : W(s.initial_redshift) - span * coordinate * coordinate,
          u = 1 + z, a = 1 / u;
      t.maximum_total_callbacks = std::min(p.thermal.maximum_total_callbacks,
                                           p.maximum_total_work - out.work_.total());
      const auto b = out.background_.scaled_expansion(a, t);
      if (!charge_momentum(out.work_, b.callbacks, p.maximum_total_work)) {
        out.status_ = S::work_limit; return out;
      }
      if (b.status != S::ok) { out.status_ = b.status; return out; }
      if (!(b.a4_e2 > b.error_estimate)) { out.status_ = S::conditioning_budget_exceeded; return out; }
      const W H = s.model.h0_km_s_mpc * 1000 / megaparsec_in_metres_wide() *
                  std::sqrt(b.a4_e2) / (a * a),
          he = b.error_estimate / (2 * (b.a4_e2 - b.error_estimate)) + arithmetic_floor;
      c[i] = {z, u, nH0 * u * u * u, nHe0 * u * u * u, H, s.model.tcmb_kelvin * u,
          8 * detail::hhe_thomson * photon * u * u * u * u /
              (3 * atomic::hydrogen_electron_mass_kg * speed_of_light_m_per_s * H * u),
          speed_of_light_m_per_s * detail::hhe_thomson / (H * u), he};
      max_H_error = std::max(max_H_error, he);
    }
    S initialization = S::ok;
    const State x0 = supplied ? out.initial_import_->promoted_values
                             : initial(c.front(), out.work_, p.maximum_total_work, initialization);
    if (initialization != S::ok || !physical(x0)) {
      out.status_ = initialization == S::ok ? S::conditioning_budget_exceeded : initialization;
      return out;
    }
    const W initial_log = heiii_log(x0[2], c.front().nH * x0[0] + c.front().nHe * x0[1]);
    out.excluded_heiii_activity_ = static_cast<double>(std::exp(initial_log));
    out.rate_domain_.emplace();
    Evolution evolution{S::ok, out.work_, p.maximum_total_work, *out.rate_domain_};
    auto mesh = [&](std::size_t stride) {
      std::vector<State> states(N / stride + 1);
      if (states.capacity() - 1 > SIZE_MAX / stride) {
        evolution.status = S::work_limit; return states;
      }
      accounted_fine = std::max(accounted_fine, (states.capacity() - 1) * stride);
      if (!payload_ok()) { evolution.status = S::work_limit; return states; }
      states[0] = x0;
      for (std::size_t i = 1; i < states.size(); ++i) {
        const std::size_t j = i * stride;
        states[i] = evolution.step(c[j], states[i - 1], c[j - stride].z - c[j].z);
        if (evolution.status != S::ok) break;
        if (!physical(states[i]) || states[i][2] > c[j].Tr) {
          evolution.status = S::conditioning_budget_exceeded; break;
        }
      }
      return states;
    };
    auto coarse = mesh(4);
    if (evolution.status != S::ok) { out.status_ = evolution.status; return out; }
    auto middle = mesh(2);
    if (evolution.status != S::ok) { out.status_ = evolution.status; return out; }
    auto fine = mesh(1);
    if (evolution.status != S::ok) { out.status_ = evolution.status; return out; }
    auto interpolate = [&](const std::vector<State> &v, std::size_t stride, std::size_t i) {
      const std::size_t j = i / stride;
      if (i % stride == 0) return v[j];
      const W f = (c[j * stride].z - c[i].z) /
                  (c[j * stride].z - c[(j + 1) * stride].z);
      State x{}; for (unsigned k = 0; k < 3; ++k) x[k] = (1 - f) * v[j][k] + f * v[j + 1][k];
      return x;
    };
    auto interpolation_error = [&](const std::vector<State> &v, std::size_t stride,
                                   std::size_t i, unsigned k) {
      if (i % stride == 0) return W(0);
      const std::size_t j = i / stride;
      const W high = c[j * stride].z, low = c[(j + 1) * stride].z,
          width = high - low, f = (high - c[i].z) / width,
          slope = (v[j][k] - v[j + 1][k]) / width;
      W second = 0;
      if (j > 0) {
        const W previous = c[(j - 1) * stride].z;
        second = std::max(second, std::abs((v[j - 1][k] - v[j][k]) / (previous - high) - slope) /
                                           ((previous - low) / 2));
      }
      if (j + 2 < v.size()) {
        const W next = c[(j + 2) * stride].z;
        second = std::max(second, std::abs(slope - (v[j + 1][k] - v[j + 2][k]) / (low - next)) /
                                           ((high - next) / 2));
      }
      return f * (1 - f) * width * width * second / 2;
    };
    out.nodes_.resize(N + 1);
    accounted_fine = std::max(accounted_fine, out.nodes_.capacity() - 1);
    if (!payload_ok()) { out.status_ = S::work_limit; return out; }
    W maximum_log = initial_log;
    std::optional<std::array<W, 2>> retained_temperature;
    for (std::size_t i = 0; i <= N; ++i) {
      const auto mid = interpolate(middle, 2, i), old_coarse = interpolate(coarse, 4, i);
      State state{}, error{};
      for (unsigned k = 0; k < 3; ++k) {
        state[k] = 2 * fine[i][k] - mid[k];
        // Richardson coefficients (2,-1) require at least twice the combined
        // mesh root contribution; adding the three unweighted totals would
        // undercount a correction occurring only in the fine mesh.
        // The interpolated-middle and interpolated-coarse errors can cancel
        // the three-mesh difference. For BE's leading two orders, actual
        // R2 error = -Delta/3 - 2*e_middle + e_coarse/3, so retain both
        // curvature terms with their algebraic weights.
        error[k] = std::abs(state[k] - (2 * mid[k] - old_coarse[k])) / 3 +
                   2 * interpolation_error(middle, 2, i, k) +
                   interpolation_error(coarse, 4, i, k) / 3 +
                   2 * evolution.accumulated_root_error[k] +
                   (arithmetic_floor + max_H_error * span) * std::abs(state[k]);
      }
      if (!physical(state) || state[2] > c[i].Tr) {
        out.status_ = S::conditioning_budget_exceeded; return out;
      }
      retain_range(retained_temperature, state[2]);
      const W activity = heiii_log(state[2], c[i].nH * state[0] + c[i].nHe * state[1]);
      maximum_log = std::max(maximum_log, activity);
      out.nodes_[i] = {c[i].z, state[0], state[1], state[2], error[0], error[1], error[2],
                       c[i].opacity, (c[i].H_error + arithmetic_floor) * c[i].opacity};
    }
    out.rate_domain_->retained_kelvin_range = retained_temperature;
    out.maximum_log_heiii_activity_ = static_cast<double>(maximum_log);
    if (maximum_log > std::log(1e-12L)) { out.status_ = S::outside_domain; return out; }
    out.status_ = S::ok;
  } catch (const std::bad_alloc &) { out.status_ = S::work_limit; }
  return out;
}
};
} // namespace detail
HydrogenHeliumHistory prepare_hydrogen_helium_history(
    const HydrogenHeliumHistoryRequest &s, HydrogenHeliumHistoryPolicy p) {
  return detail::HydrogenHeliumPreparationAccess::prepare(s, nullptr, p);
}
HydrogenHeliumHistory prepare_hydrogen_helium_supplied_history(
    const HydrogenHeliumSuppliedHistoryRequest &s, HydrogenHeliumHistoryPolicy p) {
  return detail::HydrogenHeliumPreparationAccess::prepare(s.history, &s.initial, p);
}
HydrogenHeliumHistoryBatch HydrogenHeliumHistory::evaluate(
    std::span<const double> redshifts, unsigned mask, std::size_t max_points,
    std::size_t max_bytes) const {
  HydrogenHeliumHistoryBatch out; out.requested_outputs = mask;
  if (!mask || (mask & ~31u) || !arithmetic()) return out;
  if (status_ != S::ok) { out.status = status_; return out; }
  if (redshifts.size() > max_points || redshifts.size() > 65536) { out.status = S::work_limit; return out; }
  const auto bytes = hydrogen_helium_history_payload_bound(
      nodes_.capacity() - 1, redshifts.size(),
      std::max(source_->model.species.capacity(), background_.source().species.capacity()),
      source_->nuclei_origin.capacity(),
      supplied_initial_ ? supplied_initial_->origin.capacity() : 0,
      supplied_initial_ ? supplied_initial_->source_identity.capacity() : 0);
  if (!bytes || *bytes > std::min(max_bytes, policy_.maximum_native_bytes) || *bytes > (1ull << 30)) {
    out.status = S::work_limit; return out;
  }
  try {
    out.rows.resize(redshifts.size());
    const auto actual_bytes = hydrogen_helium_history_payload_bound(
        nodes_.capacity() - 1, out.rows.capacity(),
        std::max(source_->model.species.capacity(), background_.source().species.capacity()),
        source_->nuclei_origin.capacity(),
        supplied_initial_ ? supplied_initial_->origin.capacity() : 0,
        supplied_initial_ ? supplied_initial_->source_identity.capacity() : 0);
    if (!actual_bytes || *actual_bytes > std::min(max_bytes, policy_.maximum_native_bytes) ||
        *actual_bytes > (1ull << 30)) {
      out.status = S::work_limit;
      out.rows = std::vector<HydrogenHeliumHistoryRow>{}; return out;
    }
    for (std::size_t row = 0; row < redshifts.size(); ++row) {
      const double z = redshifts[row]; auto &r = out.rows[row]; r.redshift = z;
      HydrogenHeliumHistoryValue *v[]{&r.hydrogen_ionized_fraction, &r.helium_singly_ionized_fraction,
          &r.electron_number_density_per_cubic_metre, &r.matter_temperature_kelvin, &r.thomson_opacity_per_redshift};
      if (!std::isfinite(z) || z < source_->late_redshift || z > source_->initial_redshift) {
        for (unsigned k = 0; k < 5; ++k) if (mask & (1u << k))
          v[k]->status = std::isfinite(z) ? S::outside_domain : S::nonfinite_input;
        continue;
      }
      auto lower = std::lower_bound(nodes_.begin(), nodes_.end(), W(z),
                                  [](const Node &n, W x) { return n.z > x; });
      std::size_t j = lower - nodes_.begin();
      if (j == 0) j = 1;
      if (j == nodes_.size()) j = nodes_.size() - 1;
      const auto cell = detail::HydrogenHeliumCellAccess::cell(*this, j - 1);
      const auto vcell = cell.evaluate(W(z), source_->hydrogen_nuclei_today_per_cubic_metre,
                                      source_->helium_nuclei_today_per_cubic_metre);
      const W p = vcell.hydrogen, q = vcell.helium, ep = vcell.hydrogen_error,
          eq = vcell.helium_error, T = vcell.temperature, eT = vcell.temperature_error,
          nH = vcell.nH, nHe = vcell.nHe, ne = vcell.ne, ene = vcell.ne_error,
          opacity = vcell.opacity, eopacity = vcell.opacity_error;
      if (mask & 1) project(*v[0], p, ep, policy_.absolute_fraction_tolerance, policy_.relative_fraction_tolerance);
      if (mask & 2) project(*v[1], q, eq, policy_.absolute_fraction_tolerance, policy_.relative_fraction_tolerance);
      if (mask & 4) project(*v[2], ne, ene, nH * policy_.absolute_fraction_tolerance + nHe * policy_.absolute_fraction_tolerance,
                            policy_.relative_fraction_tolerance);
      if (mask & 8) project(*v[3], T, eT, policy_.absolute_temperature_tolerance_kelvin, policy_.relative_temperature_tolerance);
      if (mask & 16) project(*v[4], opacity, eopacity, policy_.absolute_opacity_tolerance, policy_.relative_opacity_tolerance);
    }
    out.status = S::ok;
  } catch (const std::bad_alloc &) { out.status = S::work_limit; out.rows.clear(); }
  return out;
}
} // namespace irred::cosmology
