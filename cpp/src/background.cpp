#include "irred/background.hpp"
#include "flat_geometry.hpp"
#include "lcdm_state.hpp"
#include "irred/piecewise_background.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <unordered_map>
namespace irred::cosmology {
namespace {
constexpr auto radial_bit = (std::uint32_t)Observable::radial,
               shape_bit = (std::uint32_t)Observable::luminosity_shape,
               clock_bit = (std::uint32_t)Observable::clock,
               physical_bit = (std::uint32_t)Observable::flat_distances_volume,
               kin_bit = (std::uint32_t)Observable::kinematics,
               expansion_bit = (std::uint32_t)Observable::expansion;
long double dark(const CPL &p, long double z) {
  return std::exp(3 * (1 + (long double)p.w0 + p.wa) * std::log1p(z) -
                  3 * (long double)p.wa * z / (1 + z));
}
long double expansion(const ExpansionSpec &s, long double z) {
  const auto u = 1 + z;
  return std::visit(
      [&](const auto &p) -> long double {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, LCDM>)
          return detail::lcdm_expansion(p.omega_m, z);
        else if constexpr (std::is_same_v<T, CPL>)
          return std::sqrt(p.omega_m * u * u * u +
                           (p.omega_m == 1
                                ? 0.L
                                : (1 - (long double)p.omega_m) * dark(p, z)));
        else if constexpr (std::is_same_v<T, ConstantQ>)
          return std::exp((1 + (long double)p.q) * std::log1p(z));
        else
          return 0;
      },
      s);
}
struct Context {
  const ExpansionSpec *s;
  bool clock;
};
double integrand(double z, const void *p) {
  const auto &c = *static_cast<const Context *>(p);
  const auto e = expansion(*c.s, z);
  return (double)(1 / (e * (c.clock ? 1 + (long double)z : 1.L)));
}
template <class T> void fail(Outcome<T> &o, Status s, numerics::Status n) {
  o.availability = Availability::failed;
  o.status = s;
  o.numerical_status = n;
  o.value.reset();
}
template <class T> void success(Outcome<T> &o, T v) {
  o.availability = Availability::available;
  o.status = Status::ok;
  o.numerical_status = numerics::Status::ok;
  o.value = std::move(v);
}
bool normal(long double v) {
  return detail::physical_representable(v) &&
         (v == 0 || std::fpclassify((double)v) == FP_NORMAL);
}
numerics::Status cause(long double v) {
  return !std::isfinite(v) || std::abs(v) > std::numeric_limits<double>::max()
             ? numerics::Status::overflow
             : numerics::Status::outside_domain;
}
struct Node {
  explicit Node(double value) : z(value) {}
  double z;
  std::uint32_t mask = 0;
  Outcome<double> expansion;
  long double precise_E = 0;
  Outcome<Radial> radial;
  Outcome<Clock> clock;
  Outcome<Kinematics> kin;
  Work work;
};
} // namespace
Expansion prepare(ExpansionSpec spec, FlatFLRW) {
  Expansion out(std::move(spec));
  out.status_ = std::visit(
      [&](const auto &p) {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, LCDM>)
          return !std::isfinite(p.omega_m)        ? Status::invalid_input
                 : p.omega_m < 0 || p.omega_m > 1 ? Status::unsupported_domain
                                                  : Status::ok;
        else if constexpr (std::is_same_v<T, ConstantQ>)
          return !std::isfinite(p.q)   ? Status::invalid_input
                 : p.q < -2 || p.q > 2 ? Status::unsupported_domain
                                       : Status::ok;
        else if constexpr (std::is_same_v<T, CPL>) {
          if (!std::isfinite(p.omega_m) || !std::isfinite(p.w0) ||
              !std::isfinite(p.wa))
            return Status::invalid_input;
          return p.omega_m < 0 || p.omega_m > 1 || p.w0 < -2 || p.w0 > 0 ||
                         p.wa < -2 || p.wa > 2
                     ? Status::unsupported_domain
                     : Status::ok;
        } else {
          for (auto q : p.q)
            if (!std::isfinite(q))
              return Status::invalid_input;
          for (auto q : p.q)
            if (q < -3 || q > 2)
              return Status::unsupported_domain;
          out.start_E_[0] = 1;
          for (size_t i = 1; i < 5; ++i)
            out.start_E_[i] =
                out.start_E_[i - 1] *
                std::exp(
                    (1 + (long double)p.q[i - 1]) *
                    std::log1p(((long double)piecewise_q_edges[i] -
                                piecewise_q_edges[i - 1]) /
                               (1 + (long double)piecewise_q_edges[i - 1])));
          return Status::ok;
        }
      },
      out.spec_);
  return out;
}
std::string_view Expansion::model_id() const noexcept {
  switch (spec_.index()) {
  case 0:
    return "P01/flat-lcdm-radiation-free/v1";
  case 1:
    return "P01/constant-q-flat-kinematic/v1";
  case 2:
    return "P01/flat-cpl-radiation-free/v1";
  default:
    return "P01/fixed-five-bin-q-flat-kinematic/v1";
  }
}
std::optional<size_t>
Expansion::workspace_payload_bound(size_t count) noexcept {
  // Supported std::unordered_map layout: key/index value plus at most four
  // link/hash/alignment words per node, and at most ten bucket pointers per
  // key. This is 128 bytes/key on the tested x64 libstdc++ implementation.
  // Bucket reserve is checked below; allocator bookkeeping is explicitly
  // excluded.
  constexpr size_t map_payload =
      sizeof(std::pair<const std::uint64_t, size_t>) + 14 * sizeof(void *);
  constexpr size_t width =
      sizeof(Slot) + sizeof(Node) + sizeof(NodeWork) + map_payload;
  if (count > std::numeric_limits<size_t>::max() / width)
    return {};
  return count * width;
}
BatchResult Expansion::evaluate(std::span<const Request> requests,
                                EvaluationPolicy p) const {
  BatchResult out;
  if (status_ != Status::ok) {
    out.status = status_;
    return out;
  }
  const bool analytic = std::holds_alternative<FixedFiveBinQ>(spec_);
  const auto payload = workspace_payload_bound(requests.size());
  if (requests.size() > p.maximum_queries || !payload ||
      *payload > p.maximum_native_bytes) {
    out.status = Status::work_limit;
    return out;
  }
  if (requests.empty()) {
    out.status = Status::ok;
    return out;
  }
  std::vector<Node> nodes;
  nodes.reserve(requests.size());
  out.slots.reserve(requests.size());
  out.nodes.reserve(requests.size());
  std::unordered_map<std::uint64_t, size_t> keys;
  keys.reserve(requests.size());
  if (keys.bucket_count() > requests.size() * 10) {
    out.status = Status::work_limit;
    return out;
  }
  // Request masks are structural: reject the whole batch before any kernel.
  for (const auto &r : requests)
    if (!r.requested || (r.requested & ~63u))
      return out;
  bool integration_needed = false;
  for (const auto &r : requests)
    integration_needed |=
        bool(r.requested & (radial_bit | shape_bit | clock_bit | physical_bit));
  if (!analytic && integration_needed) {
    if (!p.integration)
      return out;
    const auto &ip = *p.integration;
    if (!std::isfinite(ip.absolute_tolerance) ||
        !std::isfinite(ip.relative_tolerance) || ip.absolute_tolerance < 0 ||
        ip.relative_tolerance < 0 ||
        (ip.absolute_tolerance == 0 && ip.relative_tolerance == 0) ||
        !ip.max_depth || ip.max_evaluations < 3)
      return out;
  }
  for (const auto &r : requests) {
    out.slots.emplace_back(r);
    auto &s = out.slots.back();
    Status error = Status::ok;
    if (!std::isfinite(r.z_expansion))
      error = Status::invalid_input;
    else if (r.z_expansion < 0 || r.z_expansion > (analytic ? 2.5 : 5.))
      error = Status::unsupported_domain;
    s.admission_status = error;
    if (error != Status::ok) {
      if (r.requested & expansion_bit)
        fail(s.expansion, error, numerics::Status::invalid_input);
      if (r.requested & radial_bit)
        fail(s.radial, error, numerics::Status::invalid_input);
      if (r.requested & shape_bit)
        fail(s.luminosity_shape, error, numerics::Status::invalid_input);
      if (r.requested & clock_bit)
        fail(s.clock, error, numerics::Status::invalid_input);
      if (r.requested & physical_bit)
        fail(s.physical, error, numerics::Status::invalid_input);
      if (r.requested & kin_bit)
        fail(s.kinematics, error, numerics::Status::invalid_input);
      continue;
    }
    const auto key = std::bit_cast<std::uint64_t>(r.z_expansion);
    auto [it, fresh] = keys.emplace(key, nodes.size());
    if (fresh)
      nodes.emplace_back(r.z_expansion);
    s.node_index = it->second;
    nodes[it->second].mask |= r.requested;
  }
  size_t callbacks = p.maximum_callbacks, segments = p.maximum_segment_visits;
  for (auto &node : nodes) {
    const bool need_r = node.mask & (radial_bit | shape_bit | physical_bit);
    const bool need_c = node.mask & clock_bit;
    long double E = 0, I = 0, C = 0;
    double rerr = 0, cerr = 0;
    numerics::Status rs = numerics::Status::ok, cs = numerics::Status::ok;
    size_t bin = 0;
    if (analytic) {
      const auto visits =
          (need_r || need_c) ? detail::piecewise_visits(node.z) : 0;
      if (visits > segments) {
        rs = cs = numerics::Status::work_limit;
        const auto v = detail::piecewise_integral(
            std::get<FixedFiveBinQ>(spec_), start_E_, node.z, false, false);
        E = v.E;
        bin = v.bin;
      } else {
        segments -= visits;
        node.work.segment_visits = visits;
        auto v = detail::piecewise_integral(std::get<FixedFiveBinQ>(spec_),
                                            start_E_, node.z, need_r, need_c);
        E = v.E;
        I = v.radial;
        C = v.clock;
        bin = v.bin;
      }
    } else {
      E = expansion(spec_, node.z);
      Context context{&spec_, false};
      auto integrate = [&](bool clock) {
        numerics::ScalarResult v{};
        v.status = numerics::Status::ok;
        if (node.z == 0)
          return v;
        if (callbacks < 3) {
          v.status = numerics::Status::work_limit;
          return v;
        }
        auto local = *p.integration;
        local.max_evaluations = std::min(local.max_evaluations, callbacks);
        context.clock = clock;
        v = numerics::integrate(integrand, &context, 0, node.z, local);
        if (v.evaluations > callbacks) {
          v.status = numerics::Status::work_limit;
          callbacks = 0;
        } else
          callbacks -= v.evaluations;
        node.work.callbacks += v.evaluations;
        return v;
      };
      if (need_r) {
        auto v = integrate(false);
        rs = v.status;
        I = v.value;
        rerr = v.error_estimate;
      }
      if (need_c) {
        auto v = integrate(true);
        cs = v.status;
        C = v.value;
        cerr = v.error_estimate;
      }
    }
    node.precise_E = E;
    if (node.mask & expansion_bit) {
      if (!normal(E) || !(E > 0))
        fail(node.expansion, Status::numerical_failure, cause(E));
      else
        success(node.expansion, (double)E);
    }
    if (need_r) {
      if (rs != numerics::Status::ok)
        fail(node.radial,
             rs == numerics::Status::work_limit ? Status::work_limit
                                                : Status::numerical_failure,
             rs);
      else if (!normal(E) || !(E > 0) || !normal(I) || (node.z > 0 && !(I > 0)))
        fail(node.radial, Status::numerical_failure, cause(!normal(E) ? E : I));
      else {
        Radial r;
        r.expansion_E = (double)E;
        r.integral = (double)I;
        r.error_estimate = rerr;
        r.precise_E_ = E;
        r.precise_integral_ = I;
        success(node.radial, r);
      }
    }
    if (need_c) {
      if (cs != numerics::Status::ok)
        fail(node.clock,
             cs == numerics::Status::work_limit ? Status::work_limit
                                                : Status::numerical_failure,
             cs);
      else if (!normal(C) || (node.z > 0 && !(C > 0)))
        fail(node.clock, Status::numerical_failure, cause(C));
      else
        success(node.clock, Clock{(double)C, cerr, {}});
    }
    if (node.mask & kin_bit) {
      Kinematics k;
      if (analytic) {
        const auto &v = std::get<FixedFiveBinQ>(spec_);
        k.q = v.q[bin];
        k.bin = bin;
        k.q_convention = QConvention::interior_constant_bin;
        bool jump = bin > 0 && node.z == piecewise_q_edges[bin] &&
                    v.q[bin] != v.q[bin - 1];
        if (jump) {
          k.q_convention = QConvention::right_limit_at_internal_jump;
          k.jerk_availability = JerkAvailability::unavailable_at_jump;
        } else {
          k.jerk = (double)((long double)k.q * (2 * (long double)k.q + 1));
          k.jerk_availability = JerkAvailability::ordinary_within_bin;
        }
        if (node.z == 0) {
          k.q0_within_piecewise_model = v.q[0];
          k.q_convention = QConvention::right_limit_at_zero;
          k.jerk_availability = JerkAvailability::one_sided_endpoint;
        }
        if (node.z == 2.5) {
          k.q_convention = QConvention::left_limit_at_final_endpoint;
          k.jerk_availability = JerkAvailability::one_sided_endpoint;
        }
      } else {
        const auto u = 1 + (long double)node.z;
        std::visit(
            [&](const auto &v) {
              using T = std::decay_t<decltype(v)>;
              if constexpr (std::is_same_v<T, LCDM>) {
                k.q = (double)(1.5L * v.omega_m * u * u * u / (E * E) - 1);
                k.jerk = 1;
              } else if constexpr (std::is_same_v<T, ConstantQ>) {
                k.q = v.q;
                k.jerk =
                    (double)((long double)v.q * (2 * (long double)v.q + 1));
              } else if constexpr (std::is_same_v<T, CPL>) {
                const auto f = v.omega_m == 1 ? 0.L
                                              : (1 - (long double)v.omega_m) *
                                                    dark(v, node.z) / (E * E);
                const auto w = (long double)v.w0 + v.wa * node.z / u;
                k.q = (double)(.5L + 1.5L * w * f);
                k.jerk =
                    (double)(1 + 4.5L * f * w * (1 + w) + 1.5L * f * v.wa / u);
              }
            },
            spec_);
        k.q_convention = QConvention::ordinary_smooth_model;
        k.jerk_availability = JerkAvailability::ordinary_within_bin;
      }
      if (!std::isfinite(k.q) || (k.jerk && !std::isfinite(*k.jerk)))
        fail(node.kin, Status::numerical_failure, numerics::Status::overflow);
      else
        success(node.kin, k);
    }
    out.work.callbacks += node.work.callbacks;
    out.work.segment_visits += node.work.segment_visits;
    out.nodes.push_back({node.z, node.work});
  }
  for (auto &s : out.slots) {
    if (!s.node_index)
      continue;
    const auto &node = nodes[*s.node_index];
    const auto &r = s.source;
    if (r.requested & expansion_bit) {
      if (!node.expansion.value)
        fail(s.expansion, node.expansion.status,
             node.expansion.numerical_status);
      else {
        ExpansionValue value;
        value.expansion_E = *node.expansion.value;
        value.precise_E_ = node.precise_E;
        if (r.physical_scale) {
          const auto h = r.physical_scale->h0_km_s_mpc;
          if (!std::isfinite(h) || !(h > 0))
            fail(value.h_km_s_mpc, Status::invalid_input,
                 numerics::Status::invalid_input);
          else {
            const auto hz = (long double)h * value.expansion_E;
            if (!normal(hz) || !(hz > 0))
              fail(value.h_km_s_mpc, Status::numerical_failure, cause(hz));
            else
              success(value.h_km_s_mpc, (double)hz);
          }
        }
        success(s.expansion, value);
      }
    }
    if (r.requested & radial_bit)
      s.radial = node.radial;
    if (r.requested & kin_bit)
      s.kinematics = node.kin;
    if (r.requested & clock_bit) {
      s.clock = node.clock;
      if (s.clock.value && r.physical_scale) {
        const auto h = r.physical_scale->h0_km_s_mpc;
        if (!std::isfinite(h) || !(h > 0)) {
          fail(s.clock.value->lookback_seconds, Status::invalid_input,
               numerics::Status::invalid_input);
        } else {
          const auto scale = detail::prepare_flat_scale(h);
          const auto t =
              (long double)scale.time_seconds * s.clock.value->integral;
          if (scale.status != Status::ok || !normal(t) ||
              (r.z_expansion > 0 && !(t > 0)))
            fail(s.clock.value->lookback_seconds, Status::numerical_failure,
                 cause(t));
          else
            success(s.clock.value->lookback_seconds, (double)t);
        }
      }
    }
    if (r.requested & (shape_bit | physical_bit)) {
      bool valid_o =
          r.observer && std::isfinite(r.observer->redshift) &&
          r.observer->redshift > -1 &&
          (r.observer->convention == Convention::released_zhd_zhel ||
           (r.observer->convention == Convention::geometric_same_redshift &&
            r.observer->redshift == r.z_expansion));
      if (!valid_o) {
        if (r.requested & shape_bit)
          fail(s.luminosity_shape, Status::incompatible_convention,
               numerics::Status::invalid_input);
        if (r.requested & physical_bit)
          fail(s.physical, Status::incompatible_convention,
               numerics::Status::invalid_input);
        continue;
      }
      if (!node.radial.value) {
        if (r.requested & shape_bit)
          fail(s.luminosity_shape, node.radial.status,
               node.radial.numerical_status);
        if (r.requested & physical_bit)
          fail(s.physical, node.radial.status, node.radial.numerical_status);
        continue;
      }
      const auto shape = (1 + (long double)r.observer->redshift) *
                         node.radial.value->precise_integral_;
      if (r.requested & shape_bit) {
        if (!normal(shape) || (r.z_expansion > 0 && !(shape > 0)))
          fail(s.luminosity_shape, Status::numerical_failure, cause(shape));
        else
          success(s.luminosity_shape, (double)shape);
      }
      if (r.requested & physical_bit) {
        if (!r.physical_scale ||
            !std::isfinite(r.physical_scale->h0_km_s_mpc) ||
            !(r.physical_scale->h0_km_s_mpc > 0)) {
          fail(s.physical, Status::invalid_input,
               numerics::Status::invalid_input);
          continue;
        }
        const auto scale =
            detail::prepare_flat_scale(r.physical_scale->h0_km_s_mpc);
        if (scale.status != Status::ok) {
          fail(s.physical, scale.status, numerics::Status::outside_domain);
          continue;
        }
        const auto dc = (long double)scale.distance_mpc *
                        node.radial.value->precise_integral_,
                   da = dc / (1 + (long double)r.z_expansion),
                   dl = (long double)scale.distance_mpc * shape,
                   volume = (long double)scale.distance_mpc * dc * dc /
                            node.radial.value->precise_E_;
        bool valid = true;
        for (auto x : {dc, da, dl, volume})
          if (!normal(x) || (r.z_expansion > 0 && !(x > 0))) {
            fail(s.physical, Status::numerical_failure, cause(x));
            valid = false;
            break;
          }
        if (valid)
          success(s.physical, Distances{(double)dc, (double)dc, (double)da,
                                        (double)dl, (double)volume});
      }
    }
  }
  out.status = Status::ok;
  return out;
}
} // namespace irred::cosmology
