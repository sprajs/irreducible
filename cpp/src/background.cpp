#include "irred/background.hpp"
#include "flat_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
long double cpl_density(const Parameters &p, long double z) {
  return std::exp(3 * (1 + static_cast<long double>(p.w0) + p.wa) *
                      std::log1p(z) -
                  3 * static_cast<long double>(p.wa) * z / (1 + z));
}
long double expansion(const Parameters &p, long double z) {
  const auto u = 1 + z;
  if (p.model == Model::flat_lcdm_late_v1)
    return std::sqrt(p.omega_m * u * u * u +
                     (1 - static_cast<long double>(p.omega_m)));
  if (p.model == Model::flat_cpl_late_v1) {
    const auto matter = static_cast<long double>(p.omega_m) * u * u * u;
    const auto dark =
        p.omega_m == 1
            ? 0.L
            : (1 - static_cast<long double>(p.omega_m)) * cpl_density(p, z);
    return std::sqrt(matter + dark);
  }
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
using detail::representable;
} // namespace
Background prepare(Parameters p) {
  Background out;
  out.parameters_ = p;
  if (p.model != Model::flat_lcdm_late_v1 &&
      p.model != Model::constant_q_flat_v1)
    return out;
  if (p.w0 != -1 || p.wa != 0)
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
  const auto scale = detail::prepare_flat_scale(p.h0_km_s_mpc);
  if (scale.status != Status::ok) {
    out.status_ = scale.status;
    return out;
  }
  out.hubble_time_seconds_ = scale.time_seconds;
  out.hubble_distance_mpc_ = scale.distance_mpc;
  out.status_ = Status::ok;
  return out;
}
Background prepare_cpl(CplParameters p) {
  const Parameters attempted{
      Model::flat_cpl_late_v1, p.h0_km_s_mpc, p.omega_m, 0, p.w0, p.wa};
  Background out;
  out.parameters_ = attempted;
  if (!std::isfinite(p.w0) || !std::isfinite(p.wa))
    return out;
  if (p.w0 < -2 || p.w0 > 0 || p.wa < -2 || p.wa > 2) {
    out.status_ = Status::unsupported_domain;
    return out;
  }
  // Reuse the existing owned H0/unit/flat matter-domain preparation, then
  // retain the attempted CPL identity even when that preparation reports
  // failure.
  out = prepare({Model::flat_lcdm_late_v1, p.h0_km_s_mpc, p.omega_m, 0});
  out.parameters_ = attempted;
  return out;
}
std::string_view Background::model_id() const noexcept {
  if (parameters_.model == Model::flat_lcdm_late_v1)
    return "P01/flat-lcdm-radiation-free/v1";
  if (parameters_.model == Model::constant_q_flat_v1)
    return "P01/constant-q-flat-kinematic/v1";
  if (parameters_.model == Model::flat_cpl_late_v1)
    return "P01/flat-cpl-radiation-free/v1";
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
    auto q = parameters_.model == Model::flat_lcdm_late_v1
                 ? 1.5L * parameters_.omega_m * u * u * u / (e * e) - 1
                 : static_cast<long double>(parameters_.constant_q);
    auto j =
        parameters_.model == Model::flat_lcdm_late_v1 ? 1.L : q * (2 * q + 1);
    if (parameters_.model == Model::flat_cpl_late_v1) {
      const auto dark =
          parameters_.omega_m == 1
              ? 0.L
              : (1 - static_cast<long double>(parameters_.omega_m)) *
                    cpl_density(parameters_, query.z_expansion);
      const auto f = dark / (e * e);
      const auto w =
          static_cast<long double>(parameters_.w0) +
          parameters_.wa * static_cast<long double>(query.z_expansion) / u;
      q = .5L + 1.5L * w * f;
      j = 1 + 4.5L * f * w * (1 + w) + 1.5L * f * parameters_.wa / u;
    }
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
    const auto geometry =
        detail::flat_geometry(query, e, radial.value, clock.value,
                              hubble_distance_mpc_, hubble_time_seconds_);
    if (geometry.status != Status::ok) {
      slot.status = geometry.status;
      continue;
    }
    slot.status = Status::ok;
    slot.numerical_status = numerics::Status::ok;
    slot.expansion_E = static_cast<double>(e);
    slot.h_km_s_mpc = static_cast<double>(h);
    slot.deceleration_q = static_cast<double>(q);
    slot.jerk = static_cast<double>(j);
    slot.radial_integral = radial.value;
    slot.radial_mpc = slot.transverse_mpc = static_cast<double>(geometry.dc);
    slot.angular_diameter_mpc = static_cast<double>(geometry.da);
    slot.luminosity_mpc = static_cast<double>(geometry.dl);
    slot.dimensionless_luminosity_shape = static_cast<double>(geometry.shape);
    slot.lookback_seconds = static_cast<double>(geometry.lookback);
    slot.volume_mpc3_per_sr_per_redshift = static_cast<double>(geometry.volume);
    slot.radial_integral_error = radial.error_estimate;
    slot.lookback_integral_error = clock.error_estimate;
  }
  batch.status = Status::ok;
  return batch;
}
} // namespace irred::cosmology
