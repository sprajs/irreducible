// Original rational EdS/cofactor controls; no Cholesky reference algorithm.
#include "irred/rsd.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred;
namespace {
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void near(long double x, long double y, const char *s) {
  need(std::abs(x - y) <= 1e-8L, s);
}
rsd::DensityInput input() {
  rsd::DensityInput x;
  auto &s = x.source;
  s.scale_factors = {.25, .5, 1};
  s.observed = {.3, .4, 1.1};
  s.ordered_ids = {"r0", "r1", "r2"};
  s.event_ids = {"e0", "e1", "e2"};
  s.covariance_axis_ids = s.ordered_ids;
  s.role = rsd::RowRole::synthetic_control;
  s.covariance_unit = rsd::CovarianceUnit::dimensionless_f_sigma8_squared;
  s.amplitude_convention = cosmology::Sigma8Convention::
      linear_pressureless_total_matter_top_hat_8_over_h_mpc;
  s.table_identity = "synthetic-eds-table";
  s.covariance_identity = "rational-full-covariance";
  s.ordering_provenance = "r0,r1,r2";
  s.calibration_provenance = "fixed supplied amplitude";
  s.dependence_provenance =
      "full supplied covariance; no cross-probe independence claim";
  x.covariance = {.01, .002, -.001, .002, .015, .0005, -.001, .0005, .02};
  return x;
}
rsd::PreparationPolicy preparation() {
  return {3,           9,    4096,
          1024 * 1024, 1e-8, numerics::Arithmetic::longdouble_cpu_v1};
}
rsd::DensityPolicy policy(unsigned requested = 7) {
  rsd::DensityPolicy p;
  p.maximum_models = 8;
  p.maximum_queries = 3;
  p.maximum_string_bytes = 4096;
  p.maximum_native_bytes = 1024 * 1024;
  p.maximum_total_callbacks = 4000000;
  p.maximum_forward_sensitivity = 1e-8;
  p.requested = requested;
  return p;
}
rsd::ModelPoint model(double matter = 1, double amplitude = 1) {
  return {cosmology::LCDM(matter),
          {},
          {amplitude, 1,
           cosmology::Sigma8Convention::
               linear_pressureless_total_matter_top_hat_8_over_h_mpc,
           cosmology::AmplitudeTreatment::fixed_supplied,
           "fixed synthetic sigma8", "exact supplied control"}};
}
std::array<long double, 4> cofactor(const rsd::DensityInput &x) {
  const auto &c = x.covariance;
  const long double a = c[0], b = c[1], d = c[2], e = c[4], f = c[5], g = c[8];
  const std::array<long double, 9> adj{
      e * g - f * f, d * f - b * g, b * f - d * e, d * f - b * g, a * g - d * d,
      b * d - a * f, b * f - d * e, b * d - a * f, a * e - b * b};
  const long double det = a * adj[0] + b * adj[3] + d * adj[6];
  std::array<long double, 3> r{};
  for (size_t j = 0; j < 3; ++j)
    r[j] = (long double)x.source.observed[j] - x.source.scale_factors[j];
  long double q = 0;
  for (size_t j = 0; j < 3; ++j)
    for (size_t k = 0; k < 3; ++k)
      q += r[j] * adj[j * 3 + k] * r[k] / det;
  const long double ld = std::log(det),
                    norm = 3 * std::log(2 * std::numbers::pi_v<long double>);
  return {q, ld, norm, -(q + ld + norm) / 2};
}
} // namespace
int main() {
  try {
    const auto original = input();
    auto owner = rsd::prepare_density(original, preparation());
    need(owner.status() == statistics::DensityStatus::finite &&
             owner.retained_payload_bound(),
         "prepared ordered covariance");
    need(owner.covariance().size() == 9 && owner.covariance()[2] == -.001 &&
             owner.source().event_ids[1] == "e1",
         "source and full matrix retained");
    const auto expected = cofactor(original);
    const std::array points{model(), model(.315, .811)};
    auto result = owner.evaluate(points, policy());
    need(result.status == statistics::DensityStatus::finite &&
             result.slots[0].result && result.slots[1].result,
         "EdS and ordinary density accepted");
    const auto &s = result.slots[0];
    near(s.result->quadratic, expected[0], "original cofactor quadratic");
    near(s.result->log_determinant, expected[1],
         "original cofactor determinant");
    near(s.result->normalization, expected[2],
         "original cofactor normalization");
    near(s.result->density.log_value, expected[3],
         "original cofactor normalized density");
    need(s.predictions.size() == 3 && s.residuals.size() == 3 &&
             s.predictions[0].value == .25 && s.callbacks == 0,
         "ordered requested prediction/residual arrays");
    for (unsigned mask : {1u, 2u, 4u}) {
      auto omitted = owner.evaluate(std::span(points).first(1), policy(mask));
      need(bool(omitted.slots[0].result) == bool(mask & 1) &&
               omitted.slots[0].predictions.empty() == !bool(mask & 2) &&
               omitted.slots[0].residuals.empty() == !bool(mask & 4),
           "omitted groups remain absent");
    }
    auto strict = policy();
    strict.maximum_projection_log_density_error = 1e-30;
    auto refused = owner.evaluate(points, strict);
    need(!refused.slots[0].result && refused.slots[0].predictions[0].value &&
             refused.slots[0].residuals[0].value &&
             refused.slots[0].projection_log_density_error_estimate &&
             refused.slots[0].density_state.numerical_status ==
                 numerics::Status::conditioning_budget_exceeded,
         "density projection refusal preserves successful predictions and "
         "residuals");
    auto varied = input();
    varied.source.ordering_provenance =
        "explicit covariance/table permutation 2,0,1";
    const std::array<size_t, 3> perm{2, 0, 1};
    for (size_t j = 0; j < 3; ++j) {
      varied.source.scale_factors[j] = original.source.scale_factors[perm[j]];
      varied.source.observed[j] = original.source.observed[perm[j]];
      varied.source.ordered_ids[j] = original.source.ordered_ids[perm[j]];
      varied.source.event_ids[j] = original.source.event_ids[perm[j]];
      varied.source.covariance_axis_ids[j] = varied.source.ordered_ids[j];
      for (size_t k = 0; k < 3; ++k)
        varied.covariance[j * 3 + k] =
            original.covariance[perm[j] * 3 + perm[k]];
    }
    auto permuted = rsd::prepare_density(varied, preparation())
                        .evaluate(std::span(points).first(1), policy());
    near(permuted.slots[0].result->density.log_value,
         s.result->density.log_value, "full covariance row permutation");
    need(permuted.slots[0].predictions[0].value == 1,
         "prediction follows permuted source order");
    std::array mixed{model(.315, .811), model(), model(.315, 0), model()};
    mixed[1].expansion = cosmology::CPL(.3, -1, 0);
    mixed[3].amplitude.reference_scale_factor = 0;
    auto failures = owner.evaluate(mixed, policy());
    need(failures.slots[0].result && !failures.slots[1].result &&
             failures.slots[2].result && !failures.slots[3].result,
         "mixed physical/source failures independent");
    auto cap = policy();
    cap.maximum_total_callbacks = 17;
    auto work = owner.evaluate(mixed, cap);
    need(work.callbacks <= 17 && !work.slots[0].result && work.slots[2].result,
         "joint callback cap includes failed reference");
    auto wrong = policy();
    wrong.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    need(owner.evaluate(points, wrong).slots.empty(),
         "factor arithmetic mismatch");
    for (unsigned fault = 0; fault < 7; ++fault) {
      auto x = input();
      if (fault == 0)
        x.source.role = rsd::RowRole::unknown;
      if (fault == 1)
        x.source.covariance_unit = rsd::CovarianceUnit::unknown;
      if (fault == 2)
        x.source.amplitude_convention = cosmology::Sigma8Convention::unknown;
      if (fault == 3)
        x.source.covariance_axis_ids[0] = "wrong";
      if (fault == 4)
        x.source.ordered_ids[1] = x.source.ordered_ids[0];
      if (fault == 5)
        x.source.scale_factors[0] = 0;
      if (fault == 6)
        x.source.dependence_provenance.clear();
      need(rsd::prepare_density(std::move(x), preparation()).status() !=
               statistics::DensityStatus::finite,
           "source semantic or axis failure");
    }
    auto not_spd = input();
    not_spd.covariance[0] = -1;
    need(rsd::prepare_density(std::move(not_spd), preparation())
                 .numerical_status() == numerics::Status::not_positive_definite,
         "indefinite covariance no repair");
    auto prep = preparation();
    auto boundary_input = input();
    const auto peak =
        rsd::preparation_payload_bound(boundary_input, prep.arithmetic);
    need(peak.has_value(), "preparation payload bound");
    prep.maximum_native_bytes = *peak - 1;
    need(rsd::prepare_density(std::move(boundary_input), prep)
                 .numerical_status() == numerics::Status::work_limit,
         "preparation payload refusal");
    boundary_input = input();
    prep.maximum_native_bytes =
        *rsd::preparation_payload_bound(boundary_input, prep.arithmetic);
    need(rsd::prepare_density(std::move(boundary_input), prep).status() ==
             statistics::DensityStatus::finite,
         "preparation exact payload boundary");
    auto payload = policy();
    size_t strings = 0;
    for (const auto &p : points)
      for (const auto *t :
           {&p.amplitude.amplitude_identity, &p.amplitude.amplitude_provenance})
        strings += std::max(size_t(32), t->size() + 1);
    const auto eval_peak = rsd::density_payload_bound(2, 3, strings, 7);
    need(eval_peak.has_value(), "evaluation payload bound");
    payload.maximum_native_bytes = *eval_peak - 1;
    need(owner.evaluate(points, payload).slots.empty(),
         "evaluation preallocation refusal");
    payload.maximum_native_bytes = *eval_peak;
    need(owner.evaluate(points, payload).slots.size() == 2,
         "evaluation exact payload boundary");
    need(!rsd::density_payload_bound(SIZE_MAX, 3, 0, 7),
         "payload overflow explicit");
    auto moved = std::move(owner);
    auto *self = &moved;
    moved = std::move(*self);
    need(owner.evaluate(points, policy()).slots.empty() &&
             moved.evaluate(points, policy()).slots[0].result,
         "move invalidates source and self move preserves state");
    const auto rounding = std::fegetround();
    need(std::fesetround(FE_DOWNWARD) == 0, "set rounding");
    auto unsupported = moved.evaluate(points, policy());
    need(std::fesetround(rounding) == 0, "restore rounding");
    need(unsupported.slots.empty(), "unsupported rounding rejected");
    std::cout << "PASS " << checks << " RSD owner controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
