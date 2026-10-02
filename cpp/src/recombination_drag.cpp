#include "irred/recombination_drag.hpp"
#include "hydrogen_rates.hpp"
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

struct ThermalCoeff {
  W z, radiation_temperature, H, u, n, K, compton, opacity, drag, h_relative;
};
struct ThermalState {
  W x, temperature;
};
struct ThermalEvolution {
  S status = S::ok;
  const PureHydrogenPolicy &policy;
  HydrogenHistoryWork &work;
  W B1, E21;
  ThermalState implicit(const ThermalCoeff &q, ThermalState previous, W h) {
    const W B = 1 + 2 * h / q.u;
    W lower = 0, upper = 1, x = previous.x;
    for (unsigned iteration = 0; iteration < 80; ++iteration) {
      if (work.total() > policy.maximum_total_work ||
          policy.maximum_total_work - work.total() < 2) {
        status = S::work_limit;
        return {};
      }
      ++work.rhs_evaluations;
      ++work.temperature_rhs_evaluations;
      const W coupling = q.compton * x / (1 + x), denT = B + h * coupling,
              T = (previous.temperature +
                   h * coupling * q.radiation_temperature) /
                  denT,
              Tx = h * q.compton *
                   (q.radiation_temperature * B - previous.temperature) /
                   std::pow(B + (B + h * q.compton) * x, 2);
      if (!(T > 0 && T <= q.radiation_temperature) || !std::isfinite(T)) {
        status = S::outside_domain;
        return {};
      }
      const auto rates = detail::hydrogen_rates(T, B1, E21);
      const W A = q.K * two_photon * q.n, beta_factor = q.K * q.n,
              rateB = beta_factor * (two_photon + rates.beta),
              D = q.n * rates.alpha / (q.H * q.u),
              boltzmann =
                  std::exp(-E21 / (boltzmann_constant_joule_per_kelvin * T)),
              E = rates.beta * boltzmann / (q.H * q.u),
              DT = q.n * rates.alpha_temperature_derivative / (q.H * q.u),
              ET = (rates.beta_temperature_derivative +
                    rates.beta * E21 /
                        (boltzmann_constant_joule_per_kelvin * T * T)) *
                   boltzmann / (q.H * q.u),
              neutral = 1 - x, den = 1 + rateB * neutral,
              C = (1 + A * neutral) / den, f = C * (D * x * x - E * neutral),
              fx = D * (2 * x * C + x * x * (rateB - A) / (den * den)) +
                   E * (1 + 2 * A * neutral + A * rateB * neutral * neutral) /
                       (den * den),
              CT = -(1 + A * neutral) * beta_factor *
                   rates.beta_temperature_derivative * neutral / (den * den),
              fT = C * (DT * x * x - ET * neutral) +
                   CT * (D * x * x - E * neutral),
              jacobian = 1 + h * (fx + fT * Tx),
              residual = x - previous.x + h * f,
              temperature_residual =
                  T - previous.temperature +
                  h * (coupling * (T - q.radiation_temperature) + 2 * T / q.u),
              eps = 64 * std::numeric_limits<W>::epsilon();
      // Trial-point conditioning, not a global monotonicity proof.
      if (!(jacobian > eps && std::isfinite(jacobian)) ||
          !std::isfinite(residual) || !std::isfinite(temperature_residual)) {
        status = S::conditioning_budget_exceeded;
        return {};
      }
      if (std::abs(residual) <= eps * std::max(x, 1e-8L) &&
          std::abs(temperature_residual) <= eps * std::max(T, 1.L))
        return {x, T};
      if (residual > 0)
        upper = x;
      else
        lower = x;
      W next = x - residual / jacobian;
      if (!(next > lower && next < upper))
        next = (lower + upper) / 2;
      if (next == x) {
        status = S::conditioning_budget_exceeded;
        return {};
      }
      x = next;
    }
    status = S::conditioning_budget_exceeded;
    return {};
  }
};
} // namespace
std::optional<size_t>
pure_hydrogen_payload_bound(size_t fine, size_t outputs, size_t species,
                            HydrogenTemperatureModel model) noexcept {
  irred::detail::PayloadAccounting b(sizeof(PureHydrogenHistory) +
                                     sizeof(HydrogenHistoryBatch) + 4096);
  if (fine == SIZE_MAX)
    return {};
  b.add(fine + 1,
        sizeof(PureHydrogenHistory::Node) + sizeof(Coeff) + 4 * sizeof(W));
  if (model == HydrogenTemperatureModel::evolved_compton_adiabatic)
    b.add(fine + 1, sizeof(PureHydrogenHistory::ThermalNode) + 10 * sizeof(W));
  b.add(outputs, sizeof(HydrogenHistoryRow));
  b.add(species, 4 * (sizeof(ThermalPhysicalSpecies) + sizeof(ThermalSpecies)));
  return b.result();
}
PureHydrogenHistory &
PureHydrogenHistory::operator=(const PureHydrogenHistory &o) {
  if (this != &o) {
    PureHydrogenHistory copy(o);
    *this = std::move(copy);
  }
  return *this;
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
    thermal_nodes_ = std::move(o.thermal_nodes_);
    work_ = o.work_;
    root_ = o.root_;
    o.source_.reset();
    o.nodes_.clear();
    o.thermal_nodes_.clear();
    o.status_ = S::invalid_input;
    o.work_ = {};
    o.root_ = {};
  }
  return *this;
}
std::string_view PureHydrogenHistory::model_identity() const noexcept {
  if (!source_)
    return {};
  return source_->temperature_model ==
                 HydrogenTemperatureModel::evolved_compton_adiabatic
             ? evolved_hydrogen_history_id
             : pure_hydrogen_history_id;
}
std::string_view PureHydrogenHistory::method_identity() const noexcept {
  if (!source_)
    return {};
  return source_->temperature_model ==
                 HydrogenTemperatureModel::evolved_compton_adiabatic
             ? "coupled-BE-analytic-temperature-elimination-quadratic-z-"
               "Richardson3/v1"
             : "scalar-BE-uniform-z-Richardson3/v1";
}
PureHydrogenHistory
prepare_pure_hydrogen_history(const PureHydrogenRequest &source,
                              PureHydrogenPolicy p) {
  PureHydrogenHistory out;
  const bool evolved = source.temperature_model ==
                       HydrogenTemperatureModel::evolved_compton_adiabatic;
  if ((!evolved && source.temperature_model !=
                       HydrogenTemperatureModel::prescribed_radiation) ||
      !arithmetic() ||
      !budget(p.absolute_x_tolerance, p.relative_x_tolerance) ||
      !budget(p.absolute_tau_tolerance, p.relative_tau_tolerance) ||
      !budget(p.absolute_temperature_tolerance_kelvin,
              p.relative_temperature_tolerance) ||
      !budget(p.absolute_thomson_tau_tolerance,
              p.relative_thomson_tau_tolerance) ||
      !budget(p.absolute_opacity_tolerance, p.relative_opacity_tolerance) ||
      !budget(p.absolute_visibility_tolerance,
              p.relative_visibility_tolerance) ||
      !budget(p.absolute_survival_tolerance, p.relative_survival_tolerance) ||
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
  auto bytes = pure_hydrogen_payload_bound(N, 0, source.model.species.size(),
                                           source.temperature_model);
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
  if (source.initial_redshift < 1550 ||
      source.initial_redshift > (evolved ? 1600 : 1650) ||
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
    W h = (W(source.initial_redshift) - source.late_redshift) / N,
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
    if (evolved) {
      std::vector<ThermalCoeff> q(N + 1);
      W maximum_h_error = 0;
      const W span = W(source.initial_redshift) - source.late_redshift,
              photon_today =
                  W(mapped.model->omega_gamma) *
                  detail::critical_energy_density_si(source.model.h0_km_s_mpc),
              photon_relative = floor; // mapped binary64 cast and SI arithmetic
      for (size_t i = 0; i <= N; ++i) {
        if (out.work_.total() >= p.maximum_total_work) {
          out.status_ = S::work_limit;
          return out;
        }
        ++out.work_.background_evaluations;
        W coordinate = W(i) / N,
          z = i == N ? W(source.late_redshift)
                     : source.initial_redshift - span * coordinate * coordinate,
          u = 1 + z, a = 1 / u, n = n0 * u * u * u;
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
          h_error = state.error_estimate /
                        (2 * (state.a4_e2 - state.error_estimate)) +
                    floor,
          opacity = speed_of_light_m_per_s * thomson * n / (H * u),
          R = 3 * W(mapped.model->omega_b) /
              (4 * W(mapped.model->omega_gamma) * u),
          energy = photon_today * u * u * u * u,
          compton = 8 * thomson * energy /
                    (3 * atomic::hydrogen_electron_mass_kg *
                     speed_of_light_m_per_s * H * u);
        q[i] = {z,           source.model.tcmb_kelvin * u,
                H,           u,
                n,           lya * lya * lya / (8 * pi * H),
                compton,     opacity,
                opacity / R, h_error};
        maximum_h_error = std::max(maximum_h_error, h_error);
      }
      ThermalEvolution evolution{S::ok, p, out.work_, B1, E21};
      const W T0 = q.front().radiation_temperature;
      auto mesh = [&](size_t stride) {
        std::vector<ThermalState> states(N / stride + 1);
        states[0] = {x0, T0};
        for (size_t i = 1; i < states.size(); ++i) {
          size_t j = i * stride;
          states[i] =
              evolution.implicit(q[j], states[i - 1], q[j - stride].z - q[j].z);
          if (evolution.status != S::ok)
            break;
          if (!(states[i].x > 0 && states[i].x < 1 &&
                states[i].temperature > 0 &&
                states[i].temperature <= q[j].radiation_temperature)) {
            evolution.status = S::outside_domain;
            break;
          }
        }
        return states;
      };
      auto coarse = mesh(4), middle = mesh(2), fine = mesh(1);
      if (evolution.status != S::ok) {
        out.status_ = evolution.status;
        return out;
      }
      auto interpolate = [&](const std::vector<ThermalState> &v, size_t stride,
                             size_t i) {
        if (i % stride == 0)
          return v[i / stride];
        size_t j = i / stride;
        W t = (q[j * stride].z - q[i].z) /
              (q[j * stride].z - q[(j + 1) * stride].z);
        return ThermalState{(1 - t) * v[j].x + t * v[j + 1].x,
                            (1 - t) * v[j].temperature +
                                t * v[j + 1].temperature};
      };
      std::vector<ThermalState> older(N / 2 + 1);
      for (size_t i = 0; i < older.size(); ++i) {
        auto c = interpolate(coarse, 4, 2 * i);
        older[i] = {2 * middle[i].x - c.x,
                    2 * middle[i].temperature - c.temperature};
      }
      out.nodes_.resize(N + 1);
      out.thermal_nodes_.resize(N + 1);
      const W initial_error =
          equilibrium.rows[0].ionized.relative_arithmetic_estimate * x0;
      for (size_t i = 0; i <= N; ++i) {
        auto m = interpolate(middle, 2, i), o = interpolate(older, 2, i);
        auto &node = out.nodes_[i];
        auto &tn = out.thermal_nodes_[i];
        node.x = 2 * fine[i].x - m.x;
        tn.temperature = 2 * fine[i].temperature - m.temperature;
        if (!(node.x > 0 && node.x < 1 && tn.temperature > 0 &&
              tn.temperature <= q[i].radiation_temperature)) {
          out.status_ = S::outside_domain;
          return out;
        }
        node.x_error = std::abs(node.x - o.x) / 3 + floor * node.x +
                       initial_error +
                       maximum_h_error * span * std::max(node.x, .001L);
        tn.temperature_error =
            std::abs(tn.temperature - o.temperature) / 3 +
            floor * tn.temperature +
            (maximum_h_error + photon_relative) * span * tn.temperature +
            initial_error * tn.temperature;
        tn.opacity = q[i].opacity * node.x;
        tn.opacity_error = q[i].opacity * node.x_error +
                           (q[i].h_relative + floor) * tn.opacity;
      }
      for (size_t i = N; i > 0; --i) {
        auto &low = out.nodes_[i];
        auto &high = out.nodes_[i - 1];
        auto &tl = out.thermal_nodes_[i];
        auto &th = out.thermal_nodes_[i - 1];
        W dz = q[i - 1].z - q[i].z;
        high.tau =
            low.tau + dz / 2 * (q[i].drag * low.x + q[i - 1].drag * high.x);
        high.tau_error =
            low.tau_error +
            dz / 2 * (q[i].drag * low.x_error + q[i - 1].drag * high.x_error) +
            (maximum_h_error + floor) * (high.tau - low.tau);
        th.thomson_tau = tl.thomson_tau + dz / 2 * (tl.opacity + th.opacity);
        th.thomson_tau_error = tl.thomson_tau_error +
                               dz / 2 * (tl.opacity_error + th.opacity_error) +
                               floor * (th.thomson_tau - tl.thomson_tau);
      }
      std::vector<W> old_drag(N / 2 + 1), old_thomson(N / 2 + 1);
      for (size_t i = N / 2; i > 0; --i) {
        W dz = q[2 * (i - 1)].z - q[2 * i].z;
        old_drag[i - 1] =
            old_drag[i] + dz / 2 *
                              (q[2 * i].drag * older[i].x +
                               q[2 * (i - 1)].drag * older[i - 1].x);
        old_thomson[i - 1] =
            old_thomson[i] + dz / 2 *
                                 (q[2 * i].opacity * older[i].x +
                                  q[2 * (i - 1)].opacity * older[i - 1].x);
      }
      for (size_t i = 0; i < N; ++i) {
        size_t j = i / 2;
        W oldd = old_drag[j], oldt = old_thomson[j];
        if (i % 2) {
          W dz = q[2 * j].z - q[2 * j + 2].z, l = q[i].z - q[2 * j + 2].z,
            t = l / dz, qlo = q[2 * j + 2].opacity * older[j + 1].x,
            qhi = q[2 * j].opacity * older[j].x;
          oldt = old_thomson[j + 1] + l * (qlo + (qhi - qlo) * t / 2);
          oldd = (1 - t) * old_drag[j + 1] + t * old_drag[j];
        }
        out.nodes_[i].tau_error +=
            std::abs(out.nodes_[i].tau - oldd) / 3 + floor * out.nodes_[i].tau;
        out.thermal_nodes_[i].thomson_tau_error +=
            std::abs(out.thermal_nodes_[i].thomson_tau - oldt) / 3 +
            floor * out.thermal_nodes_[i].thomson_tau;
      }
      out.status_ = S::ok;
      out.root_.status = S::outside_domain;
      for (size_t i = 0; i < N; ++i) {
        const auto &high = out.nodes_[i];
        const auto &low = out.nodes_[i + 1];
        if (high.tau >= 1 && low.tau <= 1) {
          W dz = q[i].z - q[i + 1].z, slope = (high.tau - low.tau) / dz,
            z = q[i + 1].z + (1 - low.tau) / slope,
            error = std::max(high.tau_error, low.tau_error) / slope +
                    dz / 8 *
                        std::abs(q[i].drag * high.x - q[i + 1].drag * low.x) /
                        slope +
                    floor * std::abs(z);
          if (!(slope > 0 && std::isfinite(error)))
            out.root_.status = S::conditioning_budget_exceeded;
          else
            project(out.root_, z, error, p.absolute_root_tolerance, 0);
          break;
        }
      }
      return out;
    }
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
      const auto rates = detail::hydrogen_rates(T, B1, E21);
      W H = source.model.h0_km_s_mpc * 1000.L / megaparsec_in_metres_wide() *
            std::sqrt(state.a4_e2) * u * u,
        alpha = rates.alpha, beta = rates.beta,
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
  if (!arithmetic() || !outputs || (outputs & ~127u))
    return out;
  auto bytes = pure_hydrogen_payload_bound(nodes_.size() - 1, redshifts.size(),
                                           source_->model.species.size(),
                                           source_->temperature_model);
  if (redshifts.size() > 65536 ||
      redshifts.size() >
          std::min(maximum_points, policy_.maximum_output_points) ||
      !bytes ||
      *bytes > std::min(maximum_bytes, policy_.maximum_native_bytes) ||
      *bytes > size_t(1024) * 1024 * 1024) {
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
  const bool evolved = source_->temperature_model ==
                       HydrogenTemperatureModel::evolved_compton_adiabatic;
  W start = source_->initial_redshift, end = source_->late_redshift,
    h = (start - end) / (nodes_.size() - 1);
  for (double z : redshifts) {
    out.rows.push_back({z, {}, {}, {}, {}, {}, {}, {}});
    auto &row = out.rows.back();
    auto fail = [&](S cause) {
      if (outputs & 1)
        row.electron_fraction.status = cause;
      if (outputs & 2)
        row.drag_depth.status = cause;
      if (outputs & 4)
        row.matter_temperature_kelvin.status = cause;
      if (outputs & 8)
        row.thomson_depth.status = cause;
      if (outputs & 16)
        row.thomson_opacity_per_redshift.status = cause;
      if (outputs & 32)
        row.visibility_per_redshift.status = cause;
      if (outputs & 64)
        row.finite_endpoint_survival.status = cause;
    };
    if (!std::isfinite(z)) {
      fail(S::nonfinite_input);
      continue;
    }
    if (z < end || z > start) {
      fail(S::outside_domain);
      continue;
    }
    // Exact declared endpoints select exact retained nodes, including meshes
    // whose interval spacing is not binary-representable.
    W coordinate =
        z == end     ? W(nodes_.size() - 1)
        : z == start ? W(0)
        : evolved ? std::sqrt((start - z) / (start - end)) * (nodes_.size() - 1)
                  : (start - z) / h;
    size_t i = std::min(size_t(coordinate), nodes_.size() - 2);
    auto grid_z = [&](size_t n) {
      if (!evolved)
        return start - n * h;
      if (n == nodes_.size() - 1)
        return end;
      W c = W(n) / (nodes_.size() - 1);
      return start - (start - end) * c * c;
    };
    W dz = grid_z(i) - grid_z(i + 1),
      t = evolved ? (grid_z(i) - z) / dz : coordinate - i;
    const auto &a = nodes_[i];
    const auto &b = nodes_[i + 1];
    auto curvature = [&](bool optical) {
      auto v = [&](size_t n) { return optical ? nodes_[n].tau : nodes_[n].x; };
      W c = 0;
      if (i > 0)
        c = std::max(c, std::abs(v(i - 1) - 2 * v(i) + v(i + 1)));
      if (i + 2 < nodes_.size())
        c = std::max(c, std::abs(v(i) - 2 * v(i + 1) + v(i + 2)));
      if (evolved) {
        c = 0;
        if (i > 0)
          c = std::max(
              c, std::abs((v(i - 1) - v(i)) / (grid_z(i - 1) - grid_z(i)) -
                          (v(i) - v(i + 1)) / dz) *
                     2 / (grid_z(i - 1) - grid_z(i + 1)));
        if (i + 2 < nodes_.size())
          c = std::max(c, std::abs((v(i) - v(i + 1)) / dz -
                                   (v(i + 1) - v(i + 2)) /
                                       (grid_z(i + 1) - grid_z(i + 2))) *
                              2 / (grid_z(i) - grid_z(i + 2)));
        return c * dz * dz * t * (1 - t) / 2;
      }
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
    auto strict_project = [&](HydrogenHistoryValue &group, W value, W error,
                              double ab, double rel) {
      project(group, value, error, ab, rel);
      if (group.status == S::ok && value > 0 &&
          !(group.absolute_error_estimate > 0)) {
        group = {};
        group.status = S::outside_domain;
      }
    };
    if (!evolved) {
      if (outputs & 4) {
        W T = W(source_->model.tcmb_kelvin) * (1 + W(z));
        strict_project(row.matter_temperature_kelvin, T, floor * T,
                       policy_.absolute_temperature_tolerance_kelvin,
                       policy_.relative_temperature_tolerance);
      }
      if (outputs & 8)
        row.thomson_depth.status = S::outside_domain;
      if (outputs & 16)
        row.thomson_opacity_per_redshift.status = S::outside_domain;
      if (outputs & 32)
        row.visibility_per_redshift.status = S::outside_domain;
      if (outputs & 64)
        row.finite_endpoint_survival.status = S::outside_domain;
      continue;
    }
    const auto &ta = thermal_nodes_[i];
    const auto &tb = thermal_nodes_[i + 1];
    auto thermal_curvature = [&](bool opacity) {
      auto v = [&](size_t n) {
        return opacity ? thermal_nodes_[n].opacity
                       : thermal_nodes_[n].temperature;
      };
      W c = 0;
      if (i > 0)
        c = std::max(c,
                     std::abs((v(i - 1) - v(i)) / (grid_z(i - 1) - grid_z(i)) -
                              (v(i) - v(i + 1)) / dz) *
                         2 / (grid_z(i - 1) - grid_z(i + 1)));
      if (i + 2 < nodes_.size())
        c = std::max(c, std::abs((v(i) - v(i + 1)) / dz -
                                 (v(i + 1) - v(i + 2)) /
                                     (grid_z(i + 1) - grid_z(i + 2))) *
                            2 / (grid_z(i) - grid_z(i + 2)));
      return c * dz * dz * t * (1 - t) / 2;
    };
    W T = (1 - t) * ta.temperature + t * tb.temperature,
      Te = (1 - t) * ta.temperature_error + t * tb.temperature_error +
           thermal_curvature(false),
      opacity = (1 - t) * ta.opacity + t * tb.opacity,
      oe = (1 - t) * ta.opacity_error + t * tb.opacity_error +
           thermal_curvature(true),
      length = (1 - t) * dz,
      tau = tb.thomson_tau + detail::hydrogen_opacity_cell_integral(
                                 tb.opacity, ta.opacity, dz, 1 - t),
      tau_error =
          std::max((1 - t) * ta.thomson_tau_error + t * tb.thomson_tau_error,
                   tb.thomson_tau_error +
                       detail::hydrogen_opacity_cell_integral(
                           tb.opacity_error, ta.opacity_error, dz, 1 - t)) +
          length * thermal_curvature(true) + floor * std::abs(tau),
      survival = std::exp(-tau),
      se = survival * std::expm1(tau_error) + floor * survival,
      visibility = opacity * survival,
      ve = survival * (oe + (std::abs(opacity) + oe) * std::expm1(tau_error)) +
           floor * visibility;
    if (!(T > 0 && T <= W(source_->model.tcmb_kelvin) * (1 + W(z)) &&
          opacity > 0 && tau >= 0)) {
      if (outputs & 4)
        row.matter_temperature_kelvin.status = S::outside_domain;
      if (outputs & 8)
        row.thomson_depth.status = S::outside_domain;
      if (outputs & 16)
        row.thomson_opacity_per_redshift.status = S::outside_domain;
      if (outputs & 32)
        row.visibility_per_redshift.status = S::outside_domain;
      if (outputs & 64)
        row.finite_endpoint_survival.status = S::outside_domain;
      continue;
    }
    if (outputs & 4)
      strict_project(row.matter_temperature_kelvin, T, Te,
                     policy_.absolute_temperature_tolerance_kelvin,
                     policy_.relative_temperature_tolerance);
    if (outputs & 8)
      strict_project(row.thomson_depth, tau, tau_error,
                     policy_.absolute_thomson_tau_tolerance,
                     policy_.relative_thomson_tau_tolerance);
    if (outputs & 16)
      strict_project(row.thomson_opacity_per_redshift, opacity, oe,
                     policy_.absolute_opacity_tolerance,
                     policy_.relative_opacity_tolerance);
    if (outputs & 32) {
      if (!(survival > 0 && visibility > 0 && ve > 0))
        row.visibility_per_redshift.status = S::outside_domain;
      else
        strict_project(row.visibility_per_redshift, visibility, ve,
                       policy_.absolute_visibility_tolerance,
                       policy_.relative_visibility_tolerance);
    }
    if (outputs & 64) {
      if (!(survival > 0 && se > 0))
        row.finite_endpoint_survival.status = S::outside_domain;
      else
        strict_project(row.finite_endpoint_survival, survival, se,
                       policy_.absolute_survival_tolerance,
                       policy_.relative_survival_tolerance);
    }
  }
  return out;
}
} // namespace irred::cosmology
