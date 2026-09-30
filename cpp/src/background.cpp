#include "irred/background.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
long double expansion(const Parameters &p, long double z) {
  const auto u = 1 + z;
  if (p.model == Model::flat_lcdm_late_v1)
    return std::sqrt(p.omega_m * u * u * u +
                     (1 - static_cast<long double>(p.omega_m)));
  return std::exp((1 + static_cast<long double>(p.constant_q)) * std::log1p(z));
}
struct Context {
  const Parameters *parameters;
  bool clock;
};
double integrand(double z, const void *opaque) {
  const auto &c = *static_cast<const Context *>(opaque);
  auto value = 1 / expansion(*c.parameters, z);
  if (c.clock)
    value /= 1 + static_cast<long double>(z);
  return static_cast<double>(value);
}
bool physical_representable(long double v) {
  return std::isfinite(v) &&
         std::abs(v) <= std::numeric_limits<double>::max() &&
         (v == 0 || std::abs(v) >= std::numeric_limits<double>::min());
}
bool representable(long double v) {
  return std::isfinite(v) && std::abs(v) <= std::numeric_limits<double>::max();
}
} // namespace
Background prepare(Parameters p) {
  Background out;
  out.parameters_ = p;
  if (p.model != Model::flat_lcdm_late_v1 &&
      p.model != Model::constant_q_flat_v1)
    return out;
  if (!std::isfinite(p.h0_km_s_mpc) || p.h0_km_s_mpc <= 0 ||
      !std::isfinite(p.omega_m) || !std::isfinite(p.constant_q))
    return out;
  if (p.model == Model::flat_lcdm_late_v1) {
    if (p.constant_q != 0)
      return out;
    if (p.omega_m < 0 || p.omega_m > 1) {
      out.status_ = Status::unsupported_domain;
      return out;
    }
  } else {
    if (p.omega_m != 0)
      return out;
    if (p.constant_q < -2 || p.constant_q > 2) {
      out.status_ = Status::unsupported_domain;
      return out;
    }
  }
  const auto h0 =
      convert({p.h0_km_s_mpc, Unit::km_per_s_per_mpc, Role::expansion_rate},
              {Unit::inverse_second, Role::expansion_rate});
  if (h0.status != QuantityStatus::ok) {
    out.status_ = Status::numerical_failure;
    return out;
  }
  const auto time = 1.L / h0.target.value;
  const auto length = static_cast<long double>(speed_of_light_m_per_s) * time;
  if (!representable(time) || !representable(length)) {
    out.status_ = Status::numerical_failure;
    return out;
  }
  const auto dh =
      convert({static_cast<double>(length), Unit::metre, Role::physical_length,
               Frame::none, LengthConvention::physical},
              {Unit::megaparsec, Role::physical_length, Frame::none,
               LengthConvention::physical});
  if (dh.status != QuantityStatus::ok) {
    out.status_ = Status::numerical_failure;
    return out;
  }
  out.hubble_time_seconds_ = static_cast<double>(time);
  out.hubble_distance_mpc_ = dh.target.value;
  out.status_ = Status::ok;
  return out;
}
std::string_view Background::model_id() const noexcept {
  if (parameters_.model == Model::flat_lcdm_late_v1)
    return "P01/flat-lcdm-radiation-free/v1";
  if (parameters_.model == Model::constant_q_flat_v1)
    return "P01/constant-q-flat-kinematic/v1";
  return {};
}
BatchResult Background::evaluate_batch(std::span<const Query> queries,
                                       Policy policy) const {
  BatchResult batch;
  batch.status = status_;
  if (status_ != Status::ok)
    return batch;
  // Reject an oversized batch before result allocation or input dereference.
  if (queries.size() > policy.maximum_queries) {
    batch.status = Status::work_limit;
    return batch;
  }
  const auto &ip = policy.integration;
  if (!std::isfinite(ip.absolute_tolerance) ||
      !std::isfinite(ip.relative_tolerance) || ip.absolute_tolerance < 0 ||
      ip.relative_tolerance < 0 ||
      (ip.absolute_tolerance == 0 && ip.relative_tolerance == 0) ||
      ip.max_evaluations < 3 || ip.max_depth == 0) {
    batch.status = Status::invalid_input;
    return batch;
  }
  batch.slots.resize(queries.size());
  std::size_t remaining = policy.maximum_total_evaluations;
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto &query = queries[i];
    auto &slot = batch.slots[i];
    slot.source = query;
    if (query.convention == Convention::geometric_same_redshift) {
      slot.luminosity_equation_id = "P01/flat-geometric-luminosity-distance/v1";
      slot.shape_equation_id =
          "P01/geometric-same-redshift-dimensionless-shape/v1";
    } else if (query.convention == Convention::released_zhd_zhel) {
      slot.luminosity_equation_id =
          "P01/released-zhd-zhel-luminosity-prefactor/v1";
      slot.shape_equation_id = "P01/released-zhd-zhel-intercept-free-shape/v1";
    } else {
      slot.status = Status::incompatible_convention;
      continue;
    }
    if (!std::isfinite(query.z_expansion) || !std::isfinite(query.z_observer)) {
      slot.status = Status::invalid_input;
      continue;
    }
    if (query.z_expansion < 0 || query.z_expansion > 5) {
      slot.status = Status::unsupported_domain;
      continue;
    }
    if (query.z_observer <= -1 ||
        (query.convention == Convention::geometric_same_redshift &&
         query.z_observer != query.z_expansion)) {
      slot.status = Status::incompatible_convention;
      continue;
    }
    const auto e = expansion(parameters_, query.z_expansion);
    const auto h = static_cast<long double>(parameters_.h0_km_s_mpc) * e;
    if (!representable(e) || !(e > 0) || !representable(h) || !(h > 0)) {
      slot.status = Status::numerical_failure;
      continue;
    }
    const auto u = 1 + static_cast<long double>(query.z_expansion);
    const auto q = parameters_.model == Model::flat_lcdm_late_v1
                       ? 1.5L * parameters_.omega_m * u * u * u / (e * e) - 1
                       : static_cast<long double>(parameters_.constant_q);
    const auto j =
        parameters_.model == Model::flat_lcdm_late_v1 ? 1.L : q * (2 * q + 1);
    numerics::ScalarResult radial{}, clock{};
    radial.status = clock.status = numerics::Status::ok;
    if (query.z_expansion != 0) {
      Context context{&parameters_, false};
      auto integrate = [&](bool is_clock) {
        numerics::ScalarResult value;
        value.status = numerics::Status::work_limit;
        if (remaining < 3)
          return value;
        auto local = ip;
        local.max_evaluations = std::min(local.max_evaluations, remaining);
        context.clock = is_clock;
        value = numerics::integrate(integrand, &context, 0, query.z_expansion,
                                    local);
        // The F02 evaluator counts actual callbacks on success and failure.
        if (value.evaluations > remaining) {
          remaining = 0;
          value.status = numerics::Status::work_limit;
        } else
          remaining -= value.evaluations;
        return value;
      };
      radial = integrate(false);
      slot.evaluations = radial.evaluations;
      if (radial.status != numerics::Status::ok) {
        slot.numerical_status = radial.status;
        slot.status = radial.status == numerics::Status::work_limit
                          ? Status::work_limit
                          : Status::numerical_failure;
        continue;
      }
      clock = integrate(true);
      slot.evaluations += clock.evaluations;
      if (clock.status != numerics::Status::ok) {
        slot.numerical_status = clock.status;
        slot.status = clock.status == numerics::Status::work_limit
                          ? Status::work_limit
                          : Status::numerical_failure;
        continue;
      }
    }
    if (query.z_expansion > 0 && (!(radial.value > 0) || !(clock.value > 0))) {
      slot.status = Status::numerical_failure;
      continue;
    }
    const auto dc =
        static_cast<long double>(hubble_distance_mpc_) * radial.value;
    const auto da = dc / u;
    const auto shape =
        (1 + static_cast<long double>(query.z_observer)) * radial.value;
    const auto dl = static_cast<long double>(hubble_distance_mpc_) * shape;
    const auto lookback =
        static_cast<long double>(hubble_time_seconds_) * clock.value;
    const auto volume = static_cast<long double>(hubble_distance_mpc_) * dc *
                        dc / e; // dVc/(dz dOmega)
    if (query.z_expansion > 0 &&
        (!(dc > 0) || !(da > 0) || !(shape > 0) || !(dl > 0) ||
         !(lookback > 0) || !(volume > 0))) {
      slot.status = Status::numerical_failure;
      continue;
    }
    if (!physical_representable(dc) || !physical_representable(da) ||
        !representable(shape) || !physical_representable(dl) ||
        !physical_representable(lookback) || !physical_representable(volume)) {
      slot.status = Status::numerical_failure;
      continue;
    }
    slot.status = Status::ok;
    slot.numerical_status = numerics::Status::ok;
    slot.expansion_E = static_cast<double>(e);
    slot.h_km_s_mpc = static_cast<double>(h);
    slot.deceleration_q = static_cast<double>(q);
    slot.jerk = static_cast<double>(j);
    slot.radial_integral = radial.value;
    slot.radial_mpc = slot.transverse_mpc = static_cast<double>(dc);
    slot.angular_diameter_mpc = static_cast<double>(da);
    slot.luminosity_mpc = static_cast<double>(dl);
    slot.dimensionless_luminosity_shape = static_cast<double>(shape);
    slot.lookback_seconds = static_cast<double>(lookback);
    slot.volume_mpc3_per_sr_per_redshift = static_cast<double>(volume);
    slot.radial_integral_error = radial.error_estimate;
    slot.lookback_integral_error = clock.error_estimate;
  }
  batch.status = Status::ok;
  return batch;
}
} // namespace irred::cosmology
