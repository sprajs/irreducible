// Independent direct-z and sqrt(a) composite 8-point Gauss-Legendre references.
// Frozen ratio budget 1e-11+5e-11|ref|; density components absolute1e-8.
// Reference refinement consumes <=5% of each allocation. Synthetic only.
#include "irred/bao_conditional.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace irred;
namespace {
using W = long double;
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void near(W x, W y, W budget, const char *s) {
  need(std::abs(x - y) <= budget, s);
}
template <class F> W gl(F f, W end, unsigned panels) {
  constexpr W x[] = {
      .183434642495649804939476142360184L, .525532409916328985817739049189246L,
      .796666477413626739591553936475831L, .960289856497536231683560868569473L};
  constexpr W w[] = {
      .362683783378361982965150449277195L, .313706645877887287337962201986601L,
      .222381034453374470544355994426240L, .101228536290376259152531354309962L};
  W s = 0, h = end / (2 * panels);
  for (unsigned j = 0; j < panels; ++j) {
    W mid = (2 * j + 1) * h;
    for (unsigned k = 0; k < 4; ++k)
      s += h * w[k] * (f(mid - h * x[k]) + f(mid + h * x[k]));
  }
  return s;
}
std::vector<W> reference(cosmology::SoundHorizonRequest q, unsigned n) {
  auto m = q.model;
  W l = 1 - W(m.omega_m) - m.omega_r,
    b = 3 * W(m.omega_b) / (4 * m.omega_gamma);
  auto ez = [&](W z) {
    W y = 1 + z;
    return std::sqrt(W(m.omega_r) * y * y * y * y + W(m.omega_m) * y * y * y +
                     l);
  };
  W rs = 299792.458L / m.h0_km_s_mpc / std::sqrt(3.L) *
         gl(
             [&](W u) {
               W a = u * u;
               return 2 * u /
                      std::sqrt(W(m.omega_r) + W(m.omega_m) * a +
                                l * a * a * a * a) /
                      std::sqrt(1 + b * a);
             },
             1 / std::sqrt(1 + W(q.z_drag)), n);
  W z = .7L,
    dm = 299792.458L / m.h0_km_s_mpc * gl([&](W t) { return 1 / ez(t); }, z, n),
    dh = 299792.458L / m.h0_km_s_mpc / ez(z);
  return {dm / rs, dh / rs, std::cbrt(dm * dm * z * dh) / rs};
}
bao::ConditionalDensityPolicy policy() {
  bao::ConditionalDensityPolicy p;
  p.predictions = {1e-11,
                   2e-13,
                   1e-13,
                   2e-13,
                   1000000,
                   4000000,
                   40,
                   64,
                   4 * 1024 * 1024,
                   {1e-11, 2e-13, 1000000, 40, 64, 4000000, 4 * 1024 * 1024}};
  p.maximum_models = 64;
  p.maximum_queries = 64;
  p.maximum_string_bytes = 4096;
  p.maximum_native_bytes = 8 * 1024 * 1024;
  p.maximum_total_callbacks = 8000000;
  p.maximum_forward_sensitivity = 1e-8;
  p.requested = 7;
  return p;
}
bao::PreparedDensity prepared(const std::vector<W> &ref, bool reorder = false) {
  bao::DensityInput d;
  d.queries = {{.7, bao::Observable::transverse_over_ruler},
               {.7, bao::Observable::hubble_over_ruler},
               {.7, bao::Observable::volume_over_ruler}};
  d.observed = {double(ref[0] + .2L), double(ref[1] - .3L),
                double(ref[2] + .1L)};
  d.covariance = {1, .2, .1, .2, 2, .3, .1, .3, 1.5};
  d.ordered_ids = {"dm", "dh", "dv"};
  d.role = bao::RowRole::synthetic_control;
  d.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  d.table_identity = "synthetic-v1";
  d.covariance_identity = "full-SPD-v1";
  d.ordering_provenance = "dm,dh,dv";
  d.calibration_provenance = "synthetic";
  d.dependence_provenance = "full covariance declared";
  if (reorder) {
    std::swap(d.queries[0], d.queries[2]);
    std::swap(d.observed[0], d.observed[2]);
    std::swap(d.ordered_ids[0], d.ordered_ids[2]);
    auto c = d.covariance;
    unsigned ids[] = {2, 1, 0};
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        d.covariance[j * 3 + k] = c[ids[j] * 3 + ids[k]];
  }
  return bao::prepare_density(std::move(d),
                              {64, 4096, 4096, 8 * 1024 * 1024, 1e-8,
                               numerics::Arithmetic::longdouble_cpu_v1});
}
} // namespace
int main() {
  try {
    cosmology::SoundHorizonRequest point{
        {70, .3, 9e-5, .05, 5e-5}, 1059, "supplied synthetic drag"};
    auto ref = reference(point, 256), coarse = reference(point, 128);
    for (unsigned j = 0; j < 3; ++j)
      near(ref[j], coarse[j], .05L * (1e-11L + 5e-11L * std::abs(ref[j])),
           "reference ratio refinement");
    auto d = prepared(ref);
    need(d.status() == statistics::DensityStatus::finite, "prepare");
    auto p = policy();
    auto batch = d.evaluate_conditional(std::span(&point, 1), p);
    need(batch.slots.size() == 1, "one point");
    auto &s = batch.slots[0];
    need(s.result &&
             s.result->density.status == statistics::DensityStatus::finite,
         "density finite");
    for (unsigned j = 0; j < 3; ++j)
      near(s.predictions[j], ref[j], 1e-11L + 5e-11L * std::abs(ref[j]),
           "independent ratio");
    // Explicit symmetric 3x3 cofactors: independent from production Cholesky.
    W a = 1, b = .2L, c = .1L, e = 2, f = .3L, g = 1.5L,
      det = a * (e * g - f * f) - b * (b * g - c * f) + c * (b * f - c * e);
    W inv[] = {e * g - f * f, c * f - b * g, b * f - c * e,
               c * f - b * g, a * g - c * c, b * c - a * f,
               b * f - c * e, b * c - a * f, a * e - b * b};
    W r[] = {W(d.source().observed[0]) - ref[0],
             W(d.source().observed[1]) - ref[1],
             W(d.source().observed[2]) - ref[2]},
      q = 0;
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        q += r[j] * inv[j * 3 + k] * r[k] / det;
    W coarse_r[] = {W(d.source().observed[0]) - coarse[0],
                    W(d.source().observed[1]) - coarse[1],
                    W(d.source().observed[2]) - coarse[2]},
      coarse_q = 0;
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        coarse_q += coarse_r[j] * inv[j * 3 + k] * coarse_r[k] / det;
    near(q, coarse_q, 5e-10L, "independent quadratic refinement allocation");
    near(q / 2, coarse_q / 2, 5e-10L,
         "independent logdensity refinement allocation");
    W logdet = std::log(det), norm = 3 * std::log(2 * std::acos(-1.L));
    near(s.result->quadratic, q, 1e-8, "q independent");
    near(s.result->log_determinant, logdet, 1e-8, "logdet");
    near(s.result->normalization, norm, 1e-8, "normalization");
    near(s.result->density.log_value, -(q + logdet + norm) / 2, 1e-8,
         "logdensity");
    auto h = point;
    h.model.h0_km_s_mpc = 140;
    auto hs = d.evaluate_conditional(std::span(&h, 1), p);
    near(hs.slots[0].result->density.log_value, s.result->density.log_value,
         1e-8, "H0 density cancellation");
    for (unsigned j = 0; j < 3; ++j)
      near(hs.slots[0].predictions[j], s.predictions[j], 1e-11,
           "H0 ratio cancellation");
    auto perm =
        prepared(ref, true).evaluate_conditional(std::span(&point, 1), p);
    near(perm.slots[0].result->density.log_value, s.result->density.log_value,
         1e-8, "permutation");
    p.maximum_projection_log_density_error = 1e-30;
    auto tight = d.evaluate_conditional(std::span(&point, 1), p);
    need(!tight.slots[0].result && tight.slots[0].predictions.size() == 3,
         "projection refusal preserves predictions");
    need(tight.slots[0].numerical_status ==
             numerics::Status::conditioning_budget_exceeded,
         "projection status");
    p = policy();
    p.requested = 2;
    auto only = d.evaluate_conditional(std::span(&point, 1), p);
    need(!only.slots[0].result && only.slots[0].residuals.empty(),
         "omitted density residuals");
    auto changed = point;
    changed.z_drag = 1200;
    auto drag = d.evaluate_conditional(std::span(&changed, 1), p);
    need(drag.slots[0].predictions[0] > s.predictions[0],
         "supplied drag dependence");
    changed = point;
    changed.model.omega_b = 0;
    auto baryon = d.evaluate_conditional(std::span(&changed, 1), p);
    need(baryon.slots[0].predictions[0] < s.predictions[0],
         "baryon loading convention");
    auto invalid = point;
    invalid.model.omega_gamma = 2 * invalid.model.omega_r;
    cosmology::SoundHorizonRequest mixed[] = {point, invalid, point};
    auto mix = d.evaluate_conditional(mixed, p);
    need(mix.slots.size() == 3 && mix.slots[0].predictions.size() == 3 &&
             mix.slots[1].predictions.empty() &&
             mix.slots[2].predictions.size() == 3,
         "mixed independent states");
    p.maximum_total_callbacks = only.callbacks;
    auto global = d.evaluate_conditional(mixed, p);
    need(global.callbacks <= p.maximum_total_callbacks &&
             global.slots[2].predictions.empty(),
         "global callback quota");
    p = policy();
    p.maximum_native_bytes = 1;
    need(d.evaluate_conditional(std::span(&point, 1), p).slots.empty(),
         "payload quota");
    p = policy();
    p.maximum_string_bytes = 1;
    need(d.evaluate_conditional(std::span(&point, 1), p).slots.empty(),
         "string quota");
    cosmology::SoundHorizonRequest radiation{
        {70, 0, 1, 0, 1}, 9, "radiation analytic control"};
    auto radiation_ref = reference(radiation, 128);
    auto rd = prepared(radiation_ref);
    auto rad = rd.evaluate_conditional(std::span(&radiation, 1), policy());
    W rad_rs = 299792.458L / 70 / std::sqrt(3.L) / 10,
      rad_dm = 299792.458L / 70 * (.7L / 1.7L),
      rad_dh = 299792.458L / 70 / (1.7L * 1.7L);
    near(rad.slots[0].predictions[0], rad_dm / rad_rs,
         1e-11L + 5e-11L * rad_dm / rad_rs, "radiation analytic dm ratio");
    near(rad.slots[0].predictions[1], rad_dh / rad_rs,
         1e-11L + 5e-11L * rad_dh / rad_rs, "radiation analytic dh ratio");
    auto tiny_input = d.source();
    tiny_input.queries = {{0, bao::Observable::transverse_over_ruler}};
    tiny_input.observed = {1e-300};
    tiny_input.covariance = {1};
    tiny_input.ordered_ids = {"zero-dm-tiny-residual"};
    auto tiny = bao::prepare_density(std::move(tiny_input),
                                     {64, 4096, 4096, 8 * 1024 * 1024, 1e-8,
                                      numerics::Arithmetic::longdouble_cpu_v1});
    auto cast_failure =
        tiny.evaluate_conditional(std::span(&point, 1), policy());
    need(cast_failure.slots[0].predictions.size() == 1 &&
             !cast_failure.slots[0].result &&
             !cast_failure.slots[0].projection_log_density_error_estimate &&
             cast_failure.slots[0].numerical_status ==
                 numerics::Status::outside_domain,
         "positive diagnostic underflow refused honestly");
    auto moved = std::move(d);
    auto &self = moved;
    moved = std::move(self);
    need(moved.evaluate_conditional(std::span(&point, 1), policy())
             .slots[0]
             .result.has_value(),
         "self move");
    need(d.evaluate_conditional(std::span(&point, 1), policy()).slots.empty(),
         "moved source invalid");
    std::cout << checks << " conditional BAO checks passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
