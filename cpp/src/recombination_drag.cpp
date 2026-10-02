#include "irred/recombination_drag.hpp"
#include "irred/quantities.hpp"
#include "payload_accounting.hpp"
#include "thermal_constants.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <numbers>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
constexpr W pi = std::numbers::pi_v<W>,
            floor = 128 * std::numeric_limits<double>::epsilon();
constexpr W proton_mass = 1.67262192595e-27L, thomson = 6.6524587051e-29L,
            lya = 121.5682e-9L, two_photon = 8.22458L;
struct Coeff {
  W A = 0, B = 0, D = 0, E = 0, drag = 0, h_relative = 0;
};
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
bool budget(double a, double r) {
  return std::isfinite(a) && std::isfinite(r) && a >= 0 && r >= 0 &&
         (a > 0 || r > 0);
}
void project(HydrogenHistoryValue &out, W v, W e, double absolute,
             double relative) {
  double rounded = (double)v;
  e += std::abs(v - rounded) + floor * std::abs(v);
  if (!std::isfinite(v) || !std::isfinite(e) || v < 0 || e < 0 ||
      (v != 0 && !std::isnormal(rounded)) || !std::isfinite((double)e)) {
    out.status = S::outside_domain;
    return;
  }
  if (e > absolute + relative * std::abs(v)) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out = {S::ok, rounded, (double)e};
}
struct Evolution {
  S status = S::ok;
  const PureHydrogenPolicy &policy;
  HydrogenHistoryWork &work;
  W rhs(const Coeff &q, W x, W *derivative = nullptr) {
    if (work.total() >= policy.maximum_total_work) {
      status = S::work_limit;
      return 0;
    }
    ++work.rhs_evaluations;
    W neutral = 1 - x, den = 1 + q.B * neutral, C = (1 + q.A * neutral) / den;
    if (derivative)
      *derivative =
          q.D * (2 * x * C + x * x * (q.B - q.A) / (den * den)) +
          q.E * (1 + 2 * q.A * neutral + q.A * q.B * neutral * neutral) /
              (den * den);
    return C * (q.D * x * x - q.E * neutral);
  }
  W implicit(const Coeff &q, W previous, W h) {
    W lower = 0, upper = 1, x = previous;
    for (unsigned iteration = 0; iteration < 80; ++iteration) {
      W derivative = 0, f = rhs(q, x, &derivative),
        residual = x - previous + h * f;
      if (status != S::ok)
        return 0;
      if (std::abs(residual) <
          32 * std::numeric_limits<W>::epsilon() * std::max(x, 1e-8L))
        return x;
      if (residual > 0)
        upper = x;
      else
        lower = x;
      W next = x - residual / (1 + h * derivative);
      if (!(next > lower && next < upper))
        next = (lower + upper) / 2;
      if (next == x)
        return x;
      x = next;
    }
    status = S::conditioning_budget_exceeded;
    return 0;
  }
};
} // namespace
std::optional<size_t> pure_hydrogen_payload_bound(size_t fine, size_t outputs,
                                                  size_t species) noexcept {
  irred::detail::PayloadAccounting b(sizeof(PureHydrogenHistory) +
                                     sizeof(HydrogenHistoryBatch) + 4096);
  if (fine == SIZE_MAX)
    return {};
  b.add(fine + 1,
        sizeof(PureHydrogenHistory::Node) + sizeof(Coeff) + 4 * sizeof(W));
  b.add(outputs, sizeof(HydrogenHistoryRow));
  b.add(species, 4 * (sizeof(ThermalPhysicalSpecies) + sizeof(ThermalSpecies)));
  return b.result();
}
PureHydrogenHistory::PureHydrogenHistory(PureHydrogenHistory &&o) noexcept {
  *this = std::move(o);
}
PureHydrogenHistory &
PureHydrogenHistory::operator=(PureHydrogenHistory &&o) noexcept {
  if (this != &o) {
    status_ = o.status_;
    source_ = std::move(o.source_);
    background_ = std::move(o.background_);
    policy_ = o.policy_;
    nodes_ = std::move(o.nodes_);
    work_ = o.work_;
    root_ = o.root_;
    o.source_.reset();
    o.nodes_.clear();
    o.status_ = S::invalid_input;
    o.work_ = {};
    o.root_ = {};
  }
  return *this;
}
PureHydrogenHistory
prepare_pure_hydrogen_history(const PureHydrogenRequest &source,
                              PureHydrogenPolicy p) {
  PureHydrogenHistory out;
  if (!arithmetic() ||
      !budget(p.absolute_x_tolerance, p.relative_x_tolerance) ||
      !budget(p.absolute_tau_tolerance, p.relative_tau_tolerance) ||
      !std::isfinite(p.absolute_root_tolerance) ||
      !(p.absolute_root_tolerance > 0) || !p.base_intervals)
    return out;
  if (p.base_intervals > 16384 ||
      4 * p.base_intervals > p.maximum_fine_intervals ||
      source.model.species.size() > 16) {
    out.status_ = S::work_limit;
    return out;
  }
  size_t N = 4 * p.base_intervals;
  auto bytes = pure_hydrogen_payload_bound(N, 0, source.model.species.size());
  if (!bytes || *bytes > p.maximum_native_bytes ||
      *bytes > size_t(1024) * 1024 * 1024) {
    out.status_ = S::work_limit;
    return out;
  }
  for (double v :
       {source.initial_redshift, source.late_redshift, source.model.h0_km_s_mpc,
        source.model.physical_baryon_density, source.model.physical_cdm_density,
        source.model.tcmb_kelvin,
        source.model.physical_massless_nonphoton_density})
    if (!std::isfinite(v)) {
      out.status_ = S::nonfinite_input;
      return out;
    }
  if (source.initial_redshift < 1550 || source.initial_redshift > 1650 ||
      source.late_redshift < 300 || source.late_redshift > 600 ||
      source.model.h0_km_s_mpc < 60 || source.model.h0_km_s_mpc > 80 ||
      source.model.physical_baryon_density < .015 ||
      source.model.physical_baryon_density > .03 ||
      source.model.physical_cdm_density < .08 ||
      source.model.physical_cdm_density > .15 ||
      source.model.tcmb_kelvin < 2.7 || source.model.tcmb_kelvin > 2.75 ||
      source.model.physical_massless_nonphoton_density < 0 ||
      source.model.physical_massless_nonphoton_density > 3e-5) {
    out.status_ = S::outside_domain;
    return out;
  }
  for (const auto &s : source.model.species)
    if (s.mass_ev != 0) {
      out.status_ = S::outside_domain;
      return out;
    }
  auto thermal = p.thermal;
  thermal.maximum_native_bytes =
      std::min(thermal.maximum_native_bytes, p.maximum_native_bytes);
  thermal.maximum_total_callbacks =
      std::min(thermal.maximum_total_callbacks, p.maximum_total_work);
  try {
    auto mapped = map_thermal_physical_model(source.model, thermal);
    if (!mapped.model) {
      out.status_ = mapped.status;
      return out;
    }
    out.background_ = prepare_thermal_background(*mapped.model, thermal);
    out.work_.momentum_callbacks = out.background_.preparation_callbacks();
    if (out.background_.status() != S::ok) {
      out.status_ = out.background_.status();
      return out;
    }
    out.source_ = source;
    out.policy_ = p;
    W h = (source.initial_redshift - source.late_redshift) / N,
      B1 = atomic::hydrogen_ionization_energy_ev * electron_volt_joule,
      E21 = planck_constant_joule_second * speed_of_light_m_per_s / lya,
      mass = proton_mass + atomic::hydrogen_electron_mass_kg -
             B1 / (speed_of_light_m_per_s * speed_of_light_m_per_s),
      n0 = detail::critical_mass_density_si(100) *
           source.model.physical_baryon_density / mass;
    atomic::HydrogenState initial{
        source.model.tcmb_kelvin * (1 + source.initial_redshift),
        (double)(n0 * std::pow(1 + W(source.initial_redshift), 3))};
    atomic::HydrogenPolicy equilibrium_policy;
    equilibrium_policy.requested_outputs = atomic::ionized_fraction;
    equilibrium_policy.maximum_solves =
        std::min(size_t(1), p.maximum_total_work - out.work_.total());
    auto equilibrium = atomic::evaluate_hydrogen_equilibrium(
        std::span(&initial, 1), equilibrium_policy);
    out.work_.equilibrium_solves = equilibrium.solves;
    if (equilibrium.status != S::ok || equilibrium.rows.size() != 1 ||
        equilibrium.rows[0].admission_status != S::ok ||
        equilibrium.rows[0].ionized.status != S::ok ||
        !equilibrium.rows[0].ionized.value) {
      out.status_ = equilibrium.rows.empty()
                        ? equilibrium.status
                        : equilibrium.rows[0].ionized.status;
      return out;
    }
    W x0 = *equilibrium.rows[0].ionized.value;
    std::vector<Coeff> q(N + 1);
    W maximum_h_error = 0;
    for (size_t i = 0; i <= N; ++i) {
      if (out.work_.total() >= p.maximum_total_work) {
        out.status_ = S::work_limit;
        return out;
      }
      ++out.work_.background_evaluations;
      W z = source.initial_redshift - i * h, u = 1 + z, a = 1 / u,
        T = source.model.tcmb_kelvin * u, n = n0 * u * u * u;
      thermal.maximum_total_callbacks =
          std::min(p.thermal.maximum_total_callbacks,
                   p.maximum_total_work - out.work_.total());
      auto state = out.background_.scaled_expansion(a, thermal);
      out.work_.momentum_callbacks += state.callbacks;
      if (state.status != S::ok) {
        out.status_ = state.status;
        return out;
      }
      if (!(state.a4_e2 > state.error_estimate)) {
        out.status_ = S::conditioning_budget_exceeded;
        return out;
      }
      W H = source.model.h0_km_s_mpc * 1000.L / megaparsec_in_metres_wide() *
            std::sqrt(state.a4_e2) * u * u,
        alpha = 1e-19L * 4.309L * std::pow(T / 10000, -.6166L) /
                (1 + .6703L * std::pow(T / 10000, .5300L)),
        beta =
            alpha *
            std::pow(2 * pi * atomic::hydrogen_electron_mass_kg *
                         boltzmann_constant_joule_per_kelvin * T /
                         (planck_constant_joule_second *
                          planck_constant_joule_second),
                     1.5L) *
            std::exp(-(B1 - E21) / (boltzmann_constant_joule_per_kelvin * T)),
        K = lya * lya * lya / (8 * pi * H),
        R = 3 * W(mapped.model->omega_b) /
            (4 * W(mapped.model->omega_gamma) * u),
        h_error =
            state.error_estimate / (2 * (state.a4_e2 - state.error_estimate)) +
            floor;
      q[i] = {K * two_photon * n,
              K * (two_photon + beta) * n,
              n * alpha / (H * u),
              beta *
                  std::exp(-E21 / (boltzmann_constant_joule_per_kelvin * T)) /
                  (H * u),
              speed_of_light_m_per_s * thomson * n / (H * u * R),
              h_error};
      maximum_h_error = std::max(maximum_h_error, h_error);
    }
    Evolution evolution{S::ok, p, out.work_};
    auto mesh = [&](size_t stride) {
      std::vector<W> xs(N / stride + 1);
      xs[0] = x0;
      for (size_t i = 1; i < xs.size(); ++i) {
        xs[i] = evolution.implicit(q[i * stride], xs[i - 1], h * stride);
        if (evolution.status != S::ok)
          break;
        if (!(xs[i] > 0 && xs[i] < 1)) {
          evolution.status = S::outside_domain;
          break;
        }
      }
      return xs;
    };
    auto coarse = mesh(4), middle = mesh(2), fine = mesh(1);
    if (evolution.status != S::ok) {
      out.status_ = evolution.status;
      return out;
    }
    out.nodes_.resize(N + 1);
    std::vector<W> r1(N / 2 + 1), tau_coarse(N / 2 + 1);
    for (size_t i = 0; i < r1.size(); ++i) {
      W c = i % 2 ? (coarse[i / 2] + coarse[i / 2 + 1]) / 2 : coarse[i / 2];
      r1[i] = 2 * middle[i] - c;
    }
    for (size_t i = 0; i <= N; ++i) {
      W m = i % 2 ? (middle[i / 2] + middle[i / 2 + 1]) / 2 : middle[i / 2],
        older = i % 2 ? (r1[i / 2] + r1[i / 2 + 1]) / 2 : r1[i / 2];
      auto &node = out.nodes_[i];
      node.x = 2 * fine[i] - m;
      if (!(node.x > 0 && node.x < 1)) {
        out.status_ = S::outside_domain;
        return out;
      }
      node.x_error = std::abs(node.x - older) / 3 + floor * node.x +
                     maximum_h_error *
                         (source.initial_redshift - source.late_redshift) *
                         std::max(node.x, .001L);
    }
    for (size_t i = N; i > 0; --i) {
      auto &low = out.nodes_[i];
      auto &high = out.nodes_[i - 1];
      high.tau = low.tau + h / 2 * (q[i].drag * low.x + q[i - 1].drag * high.x);
      high.tau_error =
          low.tau_error +
          h / 2 * (q[i].drag * low.x_error + q[i - 1].drag * high.x_error) +
          floor * std::abs(high.tau - low.tau) +
          maximum_h_error * std::abs(high.tau - low.tau);
    }
    for (size_t i = N / 2; i > 0; --i)
      tau_coarse[i - 1] = tau_coarse[i] + h * (q[2 * i].drag * r1[i] +
                                               q[2 * (i - 1)].drag * r1[i - 1]);
    for (size_t i = 0; i < N; ++i) {
      W older = i % 2 ? (tau_coarse[i / 2] + tau_coarse[i / 2 + 1]) / 2
                      : tau_coarse[i / 2];
      out.nodes_[i].tau_error +=
          std::abs(out.nodes_[i].tau - older) / 3 + floor * out.nodes_[i].tau;
    }
    out.status_ = S::ok;
    out.root_.status = S::outside_domain;
    for (size_t i = 0; i < N; ++i) {
      const auto &high = out.nodes_[i];
      const auto &low = out.nodes_[i + 1];
      if (high.tau >= 1 && low.tau <= 1) {
        W slope = (high.tau - low.tau) / h,
          z = source.initial_redshift - (i + 1) * h + (1 - low.tau) / slope,
          error = std::max(high.tau_error, low.tau_error) / slope +
                  h * h / 8 *
                      std::abs(q[i].drag * high.x - q[i + 1].drag * low.x) /
                      (h * slope) +
                  floor * std::abs(z);
        if (!(slope > 0 && std::isfinite(error)))
          out.root_.status = S::conditioning_budget_exceeded;
        else
          project(out.root_, z, error, p.absolute_root_tolerance, 0);
        break;
      }
    }
    return out;
  } catch (const std::bad_alloc &) {
    out.status_ = S::work_limit;
    return out;
  }
}
HydrogenHistoryBatch
PureHydrogenHistory::evaluate(std::span<const double> redshifts,
                              unsigned outputs, size_t maximum_points,
                              size_t maximum_bytes) const {
  HydrogenHistoryBatch out;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  if (!arithmetic() || !outputs || (outputs & ~3u))
    return out;
  auto bytes = pure_hydrogen_payload_bound(nodes_.size() - 1, redshifts.size(),
                                           source_->model.species.size());
  if (redshifts.size() > 65536 ||
      redshifts.size() >
          std::min(maximum_points, policy_.maximum_output_points) ||
      !bytes || *bytes > maximum_bytes || *bytes > size_t(1024) * 1024 * 1024) {
    out.status = S::work_limit;
    return out;
  }
  try {
    out.rows.reserve(redshifts.size());
  } catch (const std::bad_alloc &) {
    out.status = S::work_limit;
    return out;
  }
  out.status = S::ok;
  out.requested = outputs;
  W start = source_->initial_redshift, end = source_->late_redshift,
    h = (start - end) / (nodes_.size() - 1);
  for (double z : redshifts) {
    out.rows.push_back({z, {}, {}});
    auto &row = out.rows.back();
    auto fail = [&](S cause) {
      if (outputs & 1)
        row.electron_fraction.status = cause;
      if (outputs & 2)
        row.drag_depth.status = cause;
    };
    if (!std::isfinite(z)) {
      fail(S::nonfinite_input);
      continue;
    }
    if (z < end || z > start) {
      fail(S::outside_domain);
      continue;
    }
    W coordinate = (start - z) / h;
    size_t i = std::min(size_t(coordinate), nodes_.size() - 2);
    W t = coordinate - i;
    const auto &a = nodes_[i];
    const auto &b = nodes_[i + 1];
    auto curvature = [&](bool optical) {
      auto v = [&](size_t n) { return optical ? nodes_[n].tau : nodes_[n].x; };
      W c = 0;
      if (i > 0)
        c = std::max(c, std::abs(v(i - 1) - 2 * v(i) + v(i + 1)));
      if (i + 2 < nodes_.size())
        c = std::max(c, std::abs(v(i) - 2 * v(i + 1) + v(i + 2)));
      return c * t * (1 - t) / 2;
    };
    if (outputs & 1)
      project(row.electron_fraction, (1 - t) * a.x + t * b.x,
              (1 - t) * a.x_error + t * b.x_error + curvature(false),
              policy_.absolute_x_tolerance, policy_.relative_x_tolerance);
    if (outputs & 2)
      project(row.drag_depth, (1 - t) * a.tau + t * b.tau,
              (1 - t) * a.tau_error + t * b.tau_error + curvature(true),
              policy_.absolute_tau_tolerance, policy_.relative_tau_tolerance);
  }
  return out;
}
} // namespace irred::cosmology
