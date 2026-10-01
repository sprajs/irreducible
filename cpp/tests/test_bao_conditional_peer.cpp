// Independent direct-z and u=sqrt(a) GL4 reference, then scalar dense LDLT.
// Synthetic ordered ratio-coordinate Gaussian only. Frozen
// ratios1e-11+5e-11rel, density components absolute1e-8; reference
// refinement<=5% of allocation.
#include "irred/bao_conditional.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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
void near(W x, W y, W b, const char *s) { need(std::abs(x - y) <= b, s); }
template <class F> W gl4(F f, W end, unsigned n) {
  constexpr W x[]{.339981043584856264802665759103245L,
                  .861136311594052575223946488892810L},
      w[]{.652145154862546142626936050778001L,
          .347854845137453857373063949221999L};
  W h = end / (2 * n), sum = 0;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < 2; ++j) {
      W mid = (2 * i + 1) * h;
      sum += h * w[j] * (f(mid - h * x[j]) + f(mid + h * x[j]));
    }
  return sum;
}
std::vector<W> predictions(cosmology::SoundHorizonRequest s,
                           std::span<const bao::Query> queries, unsigned n) {
  auto m = s.model;
  W l = 1 - W(m.omega_m) - m.omega_r,
    b = 3 * W(m.omega_b) / (4 * m.omega_gamma),
    scale = 299792.458L / m.h0_km_s_mpc;
  W ruler = scale / std::sqrt(3.L) *
            gl4(
                [&](W u) {
                  W a = u * u;
                  return 2 * u /
                         std::sqrt(W(m.omega_r) + W(m.omega_m) * a +
                                   l * a * a * a * a) /
                         std::sqrt(1 + b * a);
                },
                1 / std::sqrt(1 + W(s.z_drag)), n);
  auto ez = [&](W z) {
    W y = 1 + z;
    return std::sqrt(W(m.omega_r) * y * y * y * y + W(m.omega_m) * y * y * y +
                     l);
  };
  std::vector<W> out;
  for (auto q : queries) {
    W dm = scale * gl4([&](W z) { return 1 / ez(z); }, q.z, n),
      dh = scale / ez(q.z);
    switch (q.observable) {
    case bao::Observable::transverse_over_ruler:
      out.push_back(dm / ruler);
      break;
    case bao::Observable::hubble_over_ruler:
      out.push_back(dh / ruler);
      break;
    case bao::Observable::volume_over_ruler:
      out.push_back(std::cbrt(dm * dm * q.z * dh) / ruler);
      break;
    }
  }
  return out;
}
struct Density {
  W q, ld, norm, logp;
};
Density density(std::span<const double> c, std::span<const double> y,
                std::span<const W> mu) {
  size_t n = y.size();
  std::vector<W> l(n * n), d(n), r(n);
  W ld = 0;
  for (size_t i = 0; i < n; ++i) {
    l[i * n + i] = 1;
    for (size_t j = 0; j < i; ++j) {
      W s = c[i * n + j];
      for (size_t k = 0; k < j; ++k)
        s -= l[i * n + k] * d[k] * l[j * n + k];
      l[i * n + j] = s / d[j];
    }
    W s = c[i * n + i];
    for (size_t k = 0; k < i; ++k)
      s -= l[i * n + k] * l[i * n + k] * d[k];
    need(s > 0, "independent LDLT SPD");
    d[i] = s;
    ld += std::log(s);
    r[i] = W(y[i]) - mu[i];
    for (size_t k = 0; k < i; ++k)
      r[i] -= l[i * n + k] * r[k];
  }
  W q = 0;
  for (size_t i = 0; i < n; ++i)
    q += r[i] * r[i] / d[i];
  W norm = n * std::log(2 * std::acos(-1.L));
  return {q, ld, norm, -.5L * (q + ld + norm)};
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
  p.maximum_models = 8;
  p.maximum_queries = 16;
  p.maximum_string_bytes = 4096;
  p.maximum_native_bytes = 8 * 1024 * 1024;
  p.maximum_total_callbacks = 8000000;
  p.maximum_forward_sensitivity = 1e-8;
  p.requested = 7;
  return p;
}
} // namespace
int main() {
  try {
    cosmology::SoundHorizonRequest point{
        {67, .27, 8e-5, .043, 4e-5}, 1027, "synthetic supplied drag"};
    bao::DensityInput in;
    in.queries = {{.1, bao::Observable::hubble_over_ruler},
                  {.6, bao::Observable::transverse_over_ruler},
                  {1.2, bao::Observable::volume_over_ruler},
                  {2.5, bao::Observable::transverse_over_ruler}};
    in.covariance = {1,  .2, .1,  .05, .2,  2,   .3,  .07,
                     .1, .3, 1.5, .11, .05, .07, .11, .8};
    in.ordered_ids = {"DH01", "DM06", "DV12", "DM25"};
    in.role = bao::RowRole::synthetic_control;
    in.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
    in.table_identity = "peer synthetic";
    in.covariance_identity = "peer SPD4";
    in.ordering_provenance = "explicit query order";
    in.calibration_provenance = "conditional drag";
    in.dependence_provenance = "full covariance";
    auto ref = predictions(point, in.queries, 512),
         coarse = predictions(point, in.queries, 256);
    for (size_t i = 0; i < 4; ++i) {
      near(ref[i], coarse[i], .05L * (1e-11L + 5e-11L * std::abs(ref[i])),
           "ratio refinement");
      in.observed.push_back(double(ref[i] + (i % 2 ? .4L : -.2L)));
    }
    auto expected = density(in.covariance, in.observed, ref),
         less = density(in.covariance, in.observed, coarse);
    near(expected.q, less.q, 5e-10L, "quadratic reference refinement");
    near(expected.logp, less.logp, 5e-10L, "logdensity reference refinement");
    auto retained = bao::prepare_density(
        std::move(in), {16, 256, 4096, 8 * 1024 * 1024, 1e-8,
                        numerics::Arithmetic::longdouble_cpu_v1});
    need(retained.status() == statistics::DensityStatus::finite,
         "prepared source");
    auto p = policy();
    auto got = retained.evaluate_conditional({&point, 1}, p);
    need(got.slots.size() == 1 && got.slots[0].result, "conditional density");
    auto &slot = got.slots[0];
    for (size_t i = 0; i < 4; ++i)
      near(slot.predictions[i], ref[i], 1e-11L + 5e-11L * std::abs(ref[i]),
           "independent ratio");
    near(slot.result->quadratic, expected.q, 1e-8L, "LDLT quadratic");
    near(slot.result->log_determinant, expected.ld, 1e-8L, "LDLT logdet");
    near(slot.result->normalization, expected.norm, 1e-8L,
         "coordinate normalization");
    near(slot.result->density.log_value, expected.logp, 1e-8L,
         "normalized logdensity");
    point.model.h0_km_s_mpc *= 2;
    auto doubled = retained.evaluate_conditional({&point, 1}, p);
    need(doubled.slots[0].result.has_value(), "scaled H0 density");
    near(doubled.slots[0].result->density.log_value,
         slot.result->density.log_value, 1e-8L, "H0 cancels density");
    p.maximum_projection_log_density_error = 1e-30;
    auto refused = retained.evaluate_conditional({&point, 1}, p);
    need(!refused.slots[0].result && refused.slots[0].predictions.size() == 4,
         "projection failure preserves predictions");
    p = policy();
    p.requested = 2;
    auto predictions_only = retained.evaluate_conditional({&point, 1}, p);
    need(!predictions_only.slots[0].result &&
             predictions_only.slots[0].residuals.empty() &&
             predictions_only.slots[0].predictions.size() == 4,
         "omission preserved");
    p = policy();
    p.maximum_total_callbacks = 3;
    cosmology::SoundHorizonRequest models[]{point, point};
    auto capped = retained.evaluate_conditional(models, p);
    need(capped.callbacks <= 3, "global callback quota");
    p = policy();
    p.maximum_models = 0;
    need(retained.evaluate_conditional({&point, 1}, p).slots.empty(),
         "preallocation model quota");
    auto huge_source = retained.source();
    std::fill(huge_source.covariance.begin(), huge_source.covariance.end(), 0.);
    for (size_t i = 0; i < 4; ++i)
      huge_source.covariance[i * 4 + i] = 1e300;
    auto huge = bao::prepare_density(std::move(huge_source),
                                     {16, 256, 4096, 8 * 1024 * 1024, 1e-8,
                                      numerics::Arithmetic::longdouble_cpu_v1});
    need(huge.status() == statistics::DensityStatus::finite,
         "huge diagonal source admitted");
    auto tiny_diagnostic = huge.evaluate_conditional({&point, 1}, policy());
    need(tiny_diagnostic.slots.size() == 1 &&
             !tiny_diagnostic.slots[0].result &&
             tiny_diagnostic.slots[0].predictions.size() == 4 &&
             tiny_diagnostic.slots[0].density_state.numerical_status ==
                 numerics::Status::outside_domain &&
             !tiny_diagnostic.slots[0].projection_log_density_error_estimate,
         "positive diagnostic underflow refused without hiding predictions");
    point.model.omega_r = 0;
    auto bad = retained.evaluate_conditional({&point, 1}, policy());
    need(!bad.slots[0].result && bad.slots[0].predictions.empty(),
         "invalid shared radiation refused");
    std::cout << "PASS " << checks
              << " conditional BAO independent peer controls\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
