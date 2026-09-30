// Independent retained-consumer controls migrated from obsolete adapters.
// SN: scale-factor Simpson plus exact diagonal weighted intercept profiling.
// BAO: split scale-factor GL8 plus exact determinant544/cofactor inverse.
// Same conservation equations/libm; different quadrature and no production
// solve in references. Budgets and synthetic source values are unchanged.
#include "irred/bao.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
using namespace irred;
namespace {
unsigned checks = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
void near(long double a, long double b) {
  check(std::abs(a - b) < 1e-10L, "SN assembled strict1e-10 budget");
}
void observable(double a, long double b) {
  check(std::abs((long double)a - b) <= 2e-12L + 2e-10L * std::abs(b),
        "BAO observable budget");
}
struct ReferencePoint {
  std::array<double, 5> q;
  bao::Ruler ruler;
  ReferencePoint(std::array<double, 5> x, bao::Ruler r) : q(x), ruler(r) {}
};
cosmology::EvaluationPolicy background() {
  return {{numerics::IntegrationPolicy{1e-14, 1e-13, 100000, 30}},
          16,
          20000000,
          2000,
          1048576};
}
long double integral(double z, double om, double w0, double wa, int panels) {
  const auto lower = 1 / (1 + (long double)z), h = (1 - lower) / panels;
  long double sum = 0;
  for (int i = 0; i <= panels; ++i) {
    auto a = lower + i * h;
    auto e2 = om / (a * a * a) +
              (1 - om) * std::pow(a, -3 * (1 + (long double)w0 + wa)) *
                  std::exp(-3 * (long double)wa * (1 - a));
    sum += (i == 0 || i == panels ? 1
            : i % 2               ? 4
                                  : 2) /
           (a * a * std::sqrt(e2));
  }
  return sum * h / 3;
}
void sn_control() {
  observations::Input in{};
  in.profile = observations::Profile::typed_magnitude_covariance;
  in.role = observations::Role::synthetic_control;
  in.unit = observations::Unit::magnitude;
  in.calibration = observations::Calibration::not_applicable;
  in.uncertainty = observations::Uncertainty::covariance;
  in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  in.component = observations::Component::total;
  in.table_sha256 = std::string(64, 'c');
  in.uncertainty_sha256 = std::string(64, 'd');
  in.ordering_provenance = "independent synthetic diagonal row order";
  in.measurement_ids = {"first", "second", "third", "excluded"};
  in.event_ids = in.measurement_ids;
  in.uncertainty_axis_ids = in.measurement_ids;
  in.values = {12, 15, 18, 99};
  in.uncertainty_matrix = {4, 0, 0, 8, 0, 9, 0, 7, 0, 0, 16, 6, 1, 2, 3, -1};
  in.missing = {0, 0, 0, 0};
  in.quality = {0, 0, 0, 0};
  in.source_selection = {1, 1, 1, 1};
  auto original = in;
  auto raw = observations::prepare(std::move(in), {4, 16, 4096});
  check(raw.status() == observations::Status::ok, "SN source structural");
  auto shared = std::make_shared<const observations::Prepared>(std::move(raw));
  // Explicit selection replaces obsolete implicit z>.01 preparation, preserving
  // precisely the same three-row marginal and untouched excluded asymmetric
  // data.
  std::vector<supernova::MagnitudeCoordinate> queries{
      {.1, {.11, cosmology::Convention::released_zhd_zhel}},
      {.7, {.71, cosmology::Convention::released_zhd_zhel}},
      {2, {2.2, cosmology::Convention::released_zhd_zhel}}};
  auto consumer = supernova::prepare(
      {shared, {0, 1, 2}, {"first", "second", "third"}, queries},
      {numerics::Arithmetic::longdouble_cpu_v1, 4, 16, 1e-10, 1048576});
  check(consumer.status() == supernova::Status::ok,
        "SN selected diagonal prepared");
  check(consumer.selected_source().source_indices ==
                std::vector<std::size_t>{0, 1, 2} &&
            consumer.selected_source().ordered_ids ==
                std::vector<std::string>{"first", "second", "third"},
        "exact selected source order");
  check(consumer.selected_source().source.get() == shared.get(),
        "same shared source");
  supernova::Policy p;
  p.background = background();
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.maximum_models = 4;
  p.maximum_native_bytes = 1048576;
  p.maximum_forward_sensitivity = 1e-10;
  p.requested = 63;
  std::array<cosmology::CPL, 4> cpls{
      {{0, -2. / 3, 0}, {.3, -.9, .4}, {.3, -1.1, -.4}, {1, -2, 2}}};
  std::vector<supernova::ModelPoint> models;
  for (auto c : cpls)
    models.emplace_back(c, cosmology::FlatFLRW{},
                        supernova::NoMagnitudeEffect{});
  auto batch = consumer.evaluate_batch(models, p);
  check(batch.status == supernova::Status::ok && batch.slots.size() == 4,
        "SN model batch");
  for (size_t m = 0; m < 4; ++m) {
    auto &s = batch.slots[m];
    check(s.status == supernova::Status::ok &&
              s.profile_status == statistics::DensityStatus::finite &&
              s.score.has_value(),
          "SN finite before payload");
    check(s.geometry.availability == cosmology::Availability::available &&
              s.geometry.numerical_status == numerics::Status::ok &&
              s.profile.availability == cosmology::Availability::available &&
              s.profile.numerical_status == numerics::Status::ok,
          "SN requested states successful");
    check(s.geometric_shape.size() == 3 && s.profiled_residuals.size() == 3,
          "SN array sizes");
    long double residual[3], weighted = 0, gram = 0;
    constexpr long double variance[]{4, 9, 16};
    for (size_t i = 0; i < 3; ++i) {
      auto f = integral(queries[i].z_expansion, cpls[m].omega_m, cpls[m].w0,
                        cpls[m].wa, 4096),
           c = integral(queries[i].z_expansion, cpls[m].omega_m, cpls[m].w0,
                        cpls[m].wa, 2048);
      check(std::abs(f - c) < 2e-13L, "SN independent Simpson refinement");
      auto pred = 5 * std::log10((1 + queries[i].observer.redshift) * f);
      near(s.geometric_shape[i], pred);
      residual[i] = original.values[i] - pred;
      weighted += residual[i] / variance[i];
      gram += 1 / variance[i];
    }
    auto coefficient = weighted / gram;
    long double q = 0;
    for (size_t i = 0; i < 3; ++i) {
      q += (residual[i] - coefficient) * (residual[i] - coefficient) /
           variance[i];
      near(s.profiled_residuals[i], residual[i] - coefficient);
    }
    near(s.score->offset_coefficient, coefficient);
    near(s.score->quadratic, q);
    near(s.score->relative_profile_score, -q / 2);
    check(consumer.score_id ==
              std::string_view("F03/unconstrained-offset-relative-profile/v1"),
          "relative profile not normalized evidence");
  }
  auto tight = p;
  tight.maximum_forward_sensitivity = 1e-30;
  auto fail = consumer.evaluate_batch(models, tight);
  check(!fail.slots.empty() &&
            fail.slots[0].numerical_status ==
                numerics::Status::conditioning_budget_exceeded &&
            !fail.slots[0].score && fail.slots[0].profiled_residuals.empty(),
        "SN cached sensitivity failure no profile payload");
  auto cap = p;
  cap.maximum_models = 0;
  check(consumer.evaluate_batch(models, cap).slots.empty(), "SN model cap");
  check(shared->source().uncertainty_matrix == original.uncertainty_matrix,
        "SN raw excluded matrix untouched");
}
constexpr std::array<double, 6> edges{0, .1, .3, .6, 1, 2.5};
const std::array<std::array<double, 5>, 8> q_grid{{{0, 0, 0, 0, 0},
                                                   {-1, -1, -1, -1, -1},
                                                   {.5, .5, .5, .5, .5},
                                                   {-.4, -.4, -.2, .1, .3},
                                                   {-1, 0, -1, 0, -1},
                                                   {-3, 2, -3, 2, -3},
                                                   {-3, -3, -3, -3, -3},
                                                   {2, 2, 2, 2, 2}}};
std::vector<ReferencePoint> grid() {
  std::vector<ReferencePoint> result;
  for (auto q : q_grid)
    for (double ruler : {5000., 10000., 15000.})
      result.emplace_back(q, bao::Ruler(ruler));
  return result;
}
long double reference_E(long double z, const std::array<double, 5> &q) {
  long double E = 1;
  for (unsigned k = 0; k < 5 && z > (long double)edges[k]; ++k) {
    const auto hi = std::min(z, (long double)edges[k + 1]);
    E *=
        std::pow((1 + hi) / (1 + (long double)edges[k]), 1 + (long double)q[k]);
  }
  return E;
}
// Independent fixed GL8 scale-factor quadrature, split at model edges.
// This avoids the production exprel antiderivative and its removable limits.
long double reference_I(double z, const std::array<double, 5> &q,
                        unsigned panels) {
  constexpr long double nodes[]{
      .183434642495649804939476142360L, .525532409916328985817739049189L,
      .796666477413626739591553936476L, .960289856497536231683560868569L};
  constexpr long double weights[]{
      .362683783378361982965150449277L, .313706645877887287337962201987L,
      .222381034453374470544355994426L, .101228536290376259152531354310L};
  long double sum = 0;
  for (unsigned k = 0; k < 5 && z > edges[k]; ++k) {
    const auto hi = std::min((long double)z, (long double)edges[k + 1]);
    const auto a0 = 1 / (1 + hi), a1 = 1 / (1 + (long double)edges[k]);
    const auto h = (a1 - a0) / panels;
    for (unsigned panel = 0; panel < panels; ++panel)
      for (unsigned j = 0; j < 4; ++j)
        for (int sign : {-1, 1}) {
          const auto a = a0 + (panel + .5L) * h + sign * nodes[j] * h / 2;
          sum += h / 2 * weights[j] / (a * a * reference_E(1 / a - 1, q));
        }
  }
  return sum;
}
long double reference_prediction(bao::Query query, const ReferencePoint &p,
                                 unsigned panels) {
  const auto scale = 299792.458L / (long double)p.ruler.h0_rd_km_s;
  const auto dm = scale * reference_I(query.z, p.q, panels);
  const auto dh = scale / reference_E(query.z, p.q);
  switch (query.observable) {
  case bao::Observable::transverse_over_ruler:
    return dm;
  case bao::Observable::hubble_over_ruler:
    return dh;
  case bao::Observable::volume_over_ruler:
    return std::cbrt((long double)query.z * dm * dm * dh);
  }
  throw std::runtime_error("unknown reference observable");
}
bao::DensityInput tiny_source() {
  bao::DensityInput in;
  in.queries = {{.295, bao::Observable::volume_over_ruler},
                {.51, bao::Observable::transverse_over_ruler},
                {2.33, bao::Observable::hubble_over_ruler}};
  in.observed = {8, 13, 9};
  in.covariance = {4, 1, 0, 1, 9, 2, 0, 2, 16};
  in.ordered_ids = {"tiny-dv", "tiny-dm", "tiny-dh"};
  in.role = bao::RowRole::synthetic_control;
  in.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  in.table_identity = std::string(64, 'a');
  in.covariance_identity = std::string(64, 'b');
  in.ordering_provenance = "synthetic exact correlated3x3 source order";
  in.calibration_provenance = "synthetic free ruler; no physical calibration";
  in.dependence_provenance = "explicit synthetic full covariance";
  return in;
}

void bao_control() {
  auto input = tiny_source();
  auto source = input;
  bao::PreparationPolicy prep{
      3, 9, 4096, 1048576, 1e-10, numerics::Arithmetic::longdouble_cpu_v1};
  auto owner = bao::prepare_density(std::move(input), prep);
  check(owner.status() == statistics::DensityStatus::finite,
        "BAO prepared finite");
  auto refs = grid();
  std::vector<bao::ModelPoint> models;
  for (auto &r : refs)
    models.emplace_back(cosmology::FixedFiveBinQ(r.q), cosmology::FlatFLRW{},
                        r.ruler);
  bao::DensityPolicy p;
  p.maximum_models = 24;
  p.maximum_native_bytes = 1048576;
  p.maximum_forward_sensitivity = 1e-10;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.requested = 7;
  p.observables = {background(), 3, 1048576};
  auto batch = owner.evaluate(models, p);
  check(batch.status == statistics::DensityStatus::finite &&
            batch.slots.size() == 24,
        "BAO24 owns attempts");
  constexpr long double adj[]{140, -16, 2, -16, 64, -8, 2, -8, 35};
  for (size_t m = 0; m < 24; ++m) {
    auto &s = batch.slots[m];
    check(s.result &&
              s.result->density.status == statistics::DensityStatus::finite &&
              s.background_status == cosmology::Status::ok,
          "BAO finite before values");
    check(s.predictions_state.availability ==
                  cosmology::Availability::available &&
              s.predictions_state.numerical_status == numerics::Status::ok &&
              s.residuals_state.availability ==
                  cosmology::Availability::available &&
              s.residuals_state.numerical_status == numerics::Status::ok &&
              s.density_state.availability ==
                  cosmology::Availability::available &&
              s.density_state.numerical_status == numerics::Status::ok,
          "BAO requested states successful");
    check(s.predictions.size() == 3 && s.residuals.size() == 3,
          "BAO ordered vectors");
    std::array<long double, 3> residual{};
    for (size_t i = 0; i < 3; ++i) {
      auto a = reference_prediction(source.queries[i], refs[m], 32),
           b = reference_prediction(source.queries[i], refs[m], 64);
      check(std::abs(a - b) <= .1L * (2e-12L + 2e-10L * std::abs(b)),
            "BAO splitGL8 refinement");
      observable(s.predictions[i], b);
      residual[i] = (long double)source.observed[i] - b;
      check(s.residuals[i] == source.observed[i] - s.predictions[i],
            "BAO canonical residual");
    }
    long double q = 0;
    for (size_t i = 0; i < 3; ++i)
      for (size_t k = 0; k < 3; ++k)
        q += residual[i] * adj[3 * i + k] * residual[k] / 544;
    auto expected =
        -.5L * (q + std::log(544.L) + 3 * std::log(2 * std::acos(-1.L)));
    check(std::abs((long double)s.result->density.log_value - expected) <=
              1e-8L,
          "BAO full normalized cofactor density");
  }
  // Covariance and both ordered axes travel with a row permutation.
  auto reversed = source;
  constexpr size_t order[]{2, 0, 1};
  for (size_t i = 0; i < 3; ++i) {
    reversed.queries[i] = source.queries[order[i]];
    reversed.observed[i] = source.observed[order[i]];
    reversed.ordered_ids[i] = source.ordered_ids[order[i]];
    for (size_t k = 0; k < 3; ++k)
      reversed.covariance[3 * i + k] =
          source.covariance[3 * order[i] + order[k]];
  }
  auto perm = bao::prepare_density(std::move(reversed), prep);
  check(perm.status() == statistics::DensityStatus::finite,
        "permuted covariance valid");
  auto pb = perm.evaluate(models, p);
  check(pb.slots.size() == 24, "permuted rows retained");
  for (size_t m = 0; m < 24; ++m) {
    check(pb.slots[m].result && pb.slots[m].result->density.status ==
                                    statistics::DensityStatus::finite,
          "permuted density finite");
    check(std::abs(pb.slots[m].result->density.log_value -
                   batch.slots[m].result->density.log_value) <= 1e-12,
          "normalized permutation budget");
    for (size_t i = 0; i < 3; ++i)
      check(pb.slots[m].predictions[i] == batch.slots[m].predictions[order[i]],
            "inverse permutation predictions");
  }
  auto mismatch = p;
  mismatch.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
  check(owner.evaluate(models, mismatch).slots.empty(),
        "BAO precision mismatch no fallback");
  auto cap = p;
  cap.maximum_native_bytes = 0;
  check(owner.evaluate(models, cap).numerical_status ==
            numerics::Status::work_limit,
        "BAO byte cap");
}
} // namespace
int main() {
  try {
    sn_control();
    bao_control();
    std::printf("consumer independent controls %u PASS\n", checks);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL %s after %u\n", e.what(), checks);
    return 1;
  }
}
