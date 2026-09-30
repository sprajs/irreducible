// Independent analytic BAO ratio fixtures (flat de Sitter, EdS, constant-q)
// derive from direct antiderivatives; no astronomy/reference runtime
// dependency.
#include "fixtures/bao_reference.hpp"
#include "fixtures/ldlt_reference.hpp"
#include "irred/bao.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace irred;
namespace {
int checks = 0;
double worst = 0;
void check(bool ok, const char *why) {
  ++checks;
  if (!ok)
    throw std::runtime_error(why);
}
void close(double actual, long double expected) {
  const auto difference = std::abs(static_cast<long double>(actual) - expected);
  const auto budget = 2e-12L + 2e-10L * std::abs(expected);
  worst = std::max(worst, static_cast<double>(difference / budget));
  check(difference <= budget, "analytic observable budget");
}
cosmology::Background cq(double q, double h0 = 70) {
  return cosmology::prepare({cosmology::Model::constant_q_flat_v1, h0, 0, q});
}
} // namespace
long double reference_integral(long double z, const cosmology::Parameters &p,
                               int panels) {
  constexpr long double nodes[]{
      .183434642495649804939476142360L, .525532409916328985817739049189L,
      .796666477413626739591553936476L, .960289856497536231683560868569L};
  constexpr long double weights[]{
      .362683783378361982965150449277L, .313706645877887287337962201987L,
      .222381034453374470544355994426L, .101228536290376259152531354310L};
  const auto a0 = 1 / (1 + z), h = (1 - a0) / panels;
  long double sum = 0;
  for (int panel = 0; panel < panels; ++panel)
    for (int k = 0; k < 4; ++k)
      for (int sign : {-1, 1}) {
        const auto a = a0 + (panel + .5L) * h + sign * nodes[k] * h / 2;
        const auto dark = (1 - (long double)p.omega_m) *
                          std::pow(a, -3 * (1 + (long double)p.w0 + p.wa)) *
                          std::exp(-3 * (long double)p.wa * (1 - a));
        const auto E = std::sqrt(p.omega_m / (a * a * a) + dark);
        sum += h / 2 * weights[k] / (a * a * E);
      }
  return sum;
}
void original_reference(const bao::DensityInput &input,
                        const std::vector<bao::ModelQuery> &models) {
  const auto n = input.queries.size();
  auto f = test_reference::factor(input.covariance, n);
  if (!f.valid)
    throw std::runtime_error("independent original LDLT failed");
  std::vector<long double> inverse(n * n), row_norm(n);
  long double bnorm = 0, rnorm = 0, enorm = 0;
  const auto epsilon = std::numeric_limits<long double>::epsilon();
  const auto gamma = (3 * n + 3) * epsilon / (1 - (3 * n + 3) * epsilon);
  for (std::size_t j = 0; j < n; ++j) {
    std::vector<long double> e(n);
    e[j] = 1;
    auto x = test_reference::solve_longdouble(f, e);
    for (std::size_t i = 0; i < n; ++i)
      inverse[i * n + j] = x[i];
  }
  for (std::size_t i = 0; i < n; ++i) {
    long double sum = 0, rsum = 0, esum = 0;
    for (std::size_t j = 0; j < n; ++j) {
      sum += std::abs(inverse[i * n + j]);
      long double product = 0, absproduct = 0, reconstructed = 0,
                  absreconstructed = 0;
      for (std::size_t k = 0; k < n; ++k) {
        auto term =
            (long double)input.covariance[i * n + k] * inverse[k * n + j];
        product += term;
        absproduct += std::abs(term);
        auto term2 = f.lower[i * n + k] * f.diagonal[k] * f.lower[j * n + k];
        reconstructed += term2;
        absreconstructed += std::abs(term2);
      }
      rsum +=
          std::abs((i == j ? 1.L : 0.L) - product) + gamma * (absproduct + 1);
      esum +=
          std::abs((long double)input.covariance[i * n + j] - reconstructed) +
          gamma * (absreconstructed + std::abs(input.covariance[i * n + j]));
    }
    bnorm = std::max(bnorm, sum);
    rnorm = std::max(rnorm, rsum);
    enorm = std::max(enorm, esum);
  }
  if (!(rnorm < 1))
    throw std::runtime_error("independent inverse residual bound fails");
  const auto inverse_bound = bnorm / (1 - rnorm),
             perturbation = inverse_bound * enorm;
  if (!(perturbation < 1))
    throw std::runtime_error("independent determinant perturbation fails");
  long double abslog = 0;
  for (auto d : f.diagonal)
    abslog += std::abs(std::log(d));
  const auto logdet_allowance =
      n * (-std::log1p(-perturbation)) + 64 * n * epsilon * (1 + abslog);
  for (int panels : {16, 32, 64})
    for (std::size_t m = 0; m < models.size(); ++m) {
      std::vector<long double> predictions(n), residual(n);
      const auto &p = models[m].background.parameters();
      const auto scale = 299792.458L / models[m].ruler.h0_rd_km_s;
      for (std::size_t i = 0; i < n; ++i) {
        const auto z = (long double)input.queries[i].z, u = 1 + z;
        const auto E =
            std::sqrt(p.omega_m * u * u * u +
                      (1 - (long double)p.omega_m) *
                          std::pow(u, 3 * (1 + (long double)p.w0 + p.wa)) *
                          std::exp(-3 * (long double)p.wa * z / u));
        const auto dm = scale * reference_integral(z, p, panels),
                   dh = scale / E;
        switch (input.queries[i].observable) {
        case bao::Observable::transverse_over_ruler:
          predictions[i] = dm;
          break;
        case bao::Observable::hubble_over_ruler:
          predictions[i] = dh;
          break;
        case bao::Observable::volume_over_ruler:
          predictions[i] = std::cbrt(z * dm * dm * dh);
          break;
        }
        residual[i] = (long double)input.observed[i] - predictions[i];
      }
      auto x = test_reference::solve_longdouble(f, residual);
      long double q = 0, absq = 0, residual_bound = 0, rl1 = 0;
      for (std::size_t i = 0; i < n; ++i) {
        auto term = residual[i] * x[i];
        q += term;
        absq += std::abs(term);
        rl1 += std::abs(residual[i]);
        long double cx = 0, abscx = 0;
        for (std::size_t j = 0; j < n; ++j) {
          auto t = (long double)input.covariance[i * n + j] * x[j];
          cx += t;
          abscx += std::abs(t);
        }
        residual_bound = std::max(residual_bound,
                                  std::abs(residual[i] - cx) +
                                      gamma * (abscx + std::abs(residual[i])));
      }
      const auto allowance =
          .5L * (rl1 * inverse_bound * residual_bound + gamma * absq +
                 logdet_allowance + 64 * n * epsilon);
      const auto normalization = n * std::log(2 * std::acos(-1.L));
      std::printf("BAO_REFERENCE %d %zu %.21Lg %.21Lg %.21Lg %.21Lg %.21Lg "
                  "%.21Lg %.21Lg",
                  panels, m, -.5L * (q + f.logdet + normalization), q, f.logdet,
                  normalization, allowance, inverse_bound, rnorm);
      for (auto prediction : predictions)
        std::printf(" %.21Lg", prediction);
      std::printf("\n");
    }
}
int original(const char *mean_path, const char *covariance_path) {
  bao::DensityInput input;
  std::ifstream table(mean_path), covariance(covariance_path);
  if (!table || !covariance)
    throw std::runtime_error("original file unavailable");
  std::string line;
  while (std::getline(table, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream row(line);
    double z, value;
    std::string kind, extra;
    if (!(row >> z >> value >> kind) || (row >> extra))
      throw std::runtime_error("mean row shape");
    bao::Observable observable;
    if (kind == "DM_over_rs")
      observable = bao::Observable::transverse_over_ruler;
    else if (kind == "DH_over_rs")
      observable = bao::Observable::hubble_over_ruler;
    else if (kind == "DV_over_rs")
      observable = bao::Observable::volume_over_ruler;
    else
      throw std::runtime_error("unknown original observable");
    input.queries.push_back({z, observable});
    input.observed.push_back(value);
    input.ordered_ids.push_back("DESI-DR2-ALL-GCcomb:row:" +
                                std::to_string(input.queries.size() - 1));
  }
  double element;
  while (covariance >> element)
    input.covariance.push_back(element);
  if (!covariance.eof() || input.queries.size() != 13 ||
      input.covariance.size() != 169)
    throw std::runtime_error("original 13x13 shape");
  input.role = bao::RowRole::released_fitted_distance_summary;
  input.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  input.table_identity =
      "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585";
  input.covariance_identity =
      "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509";
  input.ordering_provenance =
      "declared original DESI13 release row order; hash-gated caller";
  input.calibration_provenance =
      "released free ruler distance summaries; no early ruler calibration";
  input.dependence_provenance =
      "supplied full BAO covariance; crossprobe dependence unknown";
  bao::Policy p;
  bao::DensityPolicy dp{p,
                        20,
                        169,
                        10000,
                        1000000,
                        1e-10,
                        numerics::Arithmetic::longdouble_cpu_v1};
  auto prepared = bao::prepare_density(input, dp);
  if (prepared.status() != statistics::DensityStatus::finite)
    throw std::runtime_error("original Gaussian preparation failed");
  std::vector<bao::ModelQuery> models;
  for (double om : {0., .3, 1.})
    models.push_back(
        {cosmology::prepare_cpl({70, om, -1, 0}), bao::Ruler(10000)});
  models.push_back(
      {cosmology::prepare_cpl({70, .3, -.8, -.5}), bao::Ruler(10000)});
  models.push_back(
      {cosmology::prepare_cpl({70, .353, -.42, -1.75}), bao::Ruler(10000)});
  for (double hrd : {5000., 15000.})
    for (double om : {0., .3, 1.})
      models.push_back(
          {cosmology::prepare_cpl({70, om, -1, 0}), bao::Ruler(hrd)});
  for (double tolerance : {1e-8, 1e-10, 1e-12, 1e-14}) {
    dp.observables.background.integration = {tolerance, tolerance, 100000, 30};
    dp.observables.background.maximum_total_evaluations = 2000000;
    auto batch = prepared.evaluate(models, dp);
    if (batch.slots.size() != 11)
      throw std::runtime_error("original model count");
    for (std::size_t i = 0; i < batch.slots.size(); ++i) {
      const auto &slot = batch.slots[i];
      std::printf(
          "BAO_ORIGINAL %.17g %zu %u %u %.17g %.17g %.17g %.17g %.17g %zu",
          tolerance, i, static_cast<unsigned>(slot.result.density.status),
          static_cast<unsigned>(slot.numerical_status),
          slot.result.density.log_value, slot.result.quadratic,
          slot.result.log_determinant, slot.result.normalization,
          slot.result.estimated_forward_sensitivity, slot.evaluations);
      for (double prediction : slot.predictions)
        std::printf(" %.17g", prediction);
      std::printf("\n");
      if (tolerance == 1e-14) {
        const auto &expected = test_reference::bao_points[i];
        if (std::abs(slot.result.density.log_value - expected.log_density) >
            1e-8)
          throw std::runtime_error("frozen original normalized density budget");
        if (slot.predictions.size() != expected.prediction.size())
          throw std::runtime_error("original frozen prediction shape");
        for (std::size_t j = 0; j < expected.prediction.size(); ++j)
          if (std::abs(slot.predictions[j] - expected.prediction[j]) >
              2e-12 + 2e-10 * std::abs(expected.prediction[j]))
            throw std::runtime_error("frozen original observable budget");
      }
      if (slot.result.density.status != statistics::DensityStatus::finite)
        throw std::runtime_error("original fixedpoint failure retained");
    }
  }
  original_reference(input, models);
  return 0;
}
int main(int argc, char **argv) {
  try {
    if (argc == 4 && std::string(argv[1]) == "--verified-original-assets")
      return original(argv[2], argv[3]);
    if (argc != 1)
      throw std::runtime_error("strict optional original arguments");
    bao::Policy p;
    p.background.integration = {1e-14, 1e-13, 100000, 30};
    for (double q : {-1., -.5, -1e-10, 0., 1e-10, .5})
      for (double z : {0., .0001, .295, .51, 1.484, 2.33, 5.})
        for (double hrd : {5000., 10000., 15000.}) {
          std::vector<bao::Query> queries{
              {z, bao::Observable::transverse_over_ruler},
              {z, bao::Observable::hubble_over_ruler},
              {z, bao::Observable::volume_over_ruler}};
          auto b = bao::evaluate(cq(q), bao::Ruler(hrd), queries, p);
          check(b.status == cosmology::Status::ok && b.slots.size() == 3,
                "analytic batch");
          const auto l = std::log1p(static_cast<long double>(z));
          const auto integral = q == 0 ? l : -std::expm1(-q * l) / q;
          const auto dm = 299792.458L / hrd * integral;
          const auto dh = 299792.458L / hrd * std::exp(-(1 + q) * l);
          const long double expected[]{dm, dh, std::cbrt(z * dm * dm * dh)};
          for (std::size_t i = 0; i < 3; ++i) {
            check(b.slots[i].status == cosmology::Status::ok,
                  "analytic finite");
            check(b.slots[i].numerical_status == numerics::Status::ok,
                  "underlying finite");
            close(b.slots[i].dimensionless_value, expected[i]);
          }
        }
    std::vector<bao::Query> qs{{.1, bao::Observable::transverse_over_ruler},
                               {1., bao::Observable::hubble_over_ruler},
                               {2.33, bao::Observable::volume_over_ruler}};
    auto ref = bao::evaluate(cq(-.5, 70), bao::Ruler(10000), qs, p);
    for (double h0 : {40., 100.}) {
      auto got = bao::evaluate(cq(-.5, h0), bao::Ruler(10000), qs, p);
      for (std::size_t i = 0; i < 3; ++i)
        check(got.slots[i].dimensionless_value ==
                  ref.slots[i].dimensionless_value,
              "computational H0 cancellation exact");
    }
    auto doubled = bao::evaluate(cq(-.5), bao::Ruler(5000), qs, p);
    for (std::size_t i = 0; i < 3; ++i)
      close(doubled.slots[i].dimensionless_value,
            2.L * ref.slots[i].dimensionless_value);
    auto lambda = cosmology::prepare_cpl({70, .3, -1, 0});
    auto lcdm =
        cosmology::prepare({cosmology::Model::flat_lcdm_late_v1, 70, .3, 0});
    auto l = bao::evaluate(lambda, bao::Ruler(10000), qs, p);
    auto m = bao::evaluate(lcdm, bao::Ruler(10000), qs, p);
    for (std::size_t i = 0; i < 3; ++i)
      check(l.slots[i].dimensionless_value == m.slots[i].dimensionless_value,
            "Lambda exact limit");
    qs.push_back({1, static_cast<bao::Observable>(99)});
    qs.push_back({-1, bao::Observable::transverse_over_ruler});
    qs.push_back({std::numeric_limits<double>::quiet_NaN(),
                  bao::Observable::hubble_over_ruler});
    auto bad = bao::evaluate(cq(0), bao::Ruler(10000), qs, p);
    check(bad.slots[3].status == cosmology::Status::invalid_input,
          "unknown tag rejected");
    check(bad.slots[4].status == cosmology::Status::unsupported_domain,
          "negative z rejected");
    check(bad.slots[5].status == cosmology::Status::invalid_input,
          "NaN rejected");
    check(bad.slots[3].dimensionless_value == 0 &&
              bad.slots[4].dimensionless_value == 0,
          "no failed payload");
    for (double hrd :
         {0., 4999., 15001., std::numeric_limits<double>::infinity()}) {
      auto b = bao::evaluate(cq(0), bao::Ruler(hrd), qs, p);
      check(b.status != cosmology::Status::ok && b.slots.empty(),
            "ruler invalid empty");
    }
    auto small = p;
    small.maximum_queries = 1;
    check(bao::evaluate(cq(0), bao::Ruler(10000), qs, small).slots.empty(),
          "count cap before copy");
    small = p;
    small.maximum_native_bytes = sizeof(bao::Slot) - 1;
    check(bao::evaluate(cq(0), bao::Ruler(10000), qs, small).status ==
              cosmology::Status::work_limit,
          "byte cap");
    small = p;
    small.background.maximum_total_evaluations = 3;
    auto capped = bao::evaluate(cq(0), bao::Ruler(10000), qs, small);
    check(capped.evaluations <= 3, "global callback cap");
    for (std::size_t i = 0; i < 3; ++i)
      check(capped.slots[i].status == cosmology::Status::work_limit,
            "budget not reset per slot");
    auto empty = bao::evaluate(cq(0), bao::Ruler(10000), {}, p);
    check(empty.status == cosmology::Status::ok && empty.slots.empty() &&
              empty.evaluations == 0,
          "empty batch");
    bao::DensityInput input;
    input.queries = {{.1, bao::Observable::transverse_over_ruler},
                     {.2, bao::Observable::hubble_over_ruler}};
    input.observed = {299792.458 / 10000 * .1 + 1, 299792.458 / 10000 + 2};
    input.covariance = {4, 1, 1, 9};
    input.ordered_ids = {"synthetic:DM:.1", "synthetic:DH:.2"};
    input.role = bao::RowRole::synthetic_control;
    input.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
    input.table_identity = "analytic deSitter synthetic ratios";
    input.covariance_identity = "rational covariance determinant35";
    input.ordering_provenance = "explicit fixture row order";
    input.calibration_provenance = "synthetic no calibration";
    input.dependence_provenance =
        "two correlated synthetic coordinates; no crossprobe claim";
    bao::DensityPolicy dp{p,
                          10,
                          4,
                          10000,
                          100000,
                          1e-10,
                          numerics::Arithmetic::longdouble_cpu_v1};
    auto prepared = bao::prepare_density(input, dp);
    check(prepared.status() == statistics::DensityStatus::finite,
          "density prepares");
    std::vector<bao::ModelQuery> models{{cq(-1), bao::Ruler(10000)}};
    auto density = prepared.evaluate(models, dp);
    check(density.slots.size() == 1 && density.slots[0].result.density.status ==
                                           statistics::DensityStatus::finite,
          "density finite");
    close(density.slots[0].result.quadratic, .6L);
    const auto expected_log =
        -.5L * (.6L + std::log(35.L) + 2 * std::log(2 * std::acos(-1.L)));
    check(std::abs(density.slots[0].result.density.log_value - expected_log) <
              1e-8,
          "normalized determinant density budget");
    check(density.slots[0].predictions.size() == 2 &&
              density.slots[0].residuals.size() == 2,
          "retained ordered arrays");
    auto tightdp = dp;
    tightdp.maximum_forward_sensitivity = 1e-30;
    auto failed = prepared.evaluate(models, tightdp);
    check(failed.slots[0].numerical_status ==
              numerics::Status::conditioning_budget_exceeded,
          "tight exact failure cause");
    check(failed.slots[0].predictions.empty() &&
              failed.slots[0].residuals.empty(),
          "failed density arrays cleared");
    models.push_back({cq(-1), bao::Ruler(0)});
    auto mixed = prepared.evaluate(models, dp);
    check(mixed.slots[1].attempted_h0_rd_km_s == 0 &&
              mixed.slots[1].attempted_model.constant_q == -1,
          "attempted failure identity");
    check(mixed.slots[1].result.density.status !=
                  statistics::DensityStatus::finite &&
              mixed.slots[1].predictions.empty(),
          "invalid ruler has no finite density");
    auto original = input;
    input.ordered_ids[1] = input.ordered_ids[0];
    check(bao::prepare_density(input, dp).status() !=
              statistics::DensityStatus::finite,
          "duplicate row identity rejects");
    input = original;
    input.covariance[1] = 2;
    check(bao::prepare_density(input, dp).status() !=
              statistics::DensityStatus::finite,
          "source asymmetry no repair");
    input = original;
    auto capdp = dp;
    capdp.maximum_native_bytes = 1;
    auto capped_density = bao::prepare_density(input, capdp);
    check(capped_density.status() != statistics::DensityStatus::finite &&
              capped_density.numerical_status() ==
                  numerics::Status::work_limit &&
              capped_density.source().queries.empty(),
          "source cap before ownedcopy");
    capdp = dp;
    capdp.maximum_native_bytes = 1;
    auto capped_batch = prepared.evaluate(models, capdp);
    check(capped_batch.slots.empty() &&
              capped_batch.numerical_status == numerics::Status::work_limit,
          "density output cap before allocate");
    const auto perquery = sizeof(bao::Slot) + sizeof(cosmology::Query) +
                          sizeof(std::size_t) + sizeof(cosmology::Slot);
    auto peak = p;
    peak.maximum_native_bytes = perquery - 1;
    check(bao::evaluate(cq(-1), bao::Ruler(10000),
                        std::span<const bao::Query>(qs.data(), 1), peak)
              .slots.empty(),
          "peak scratch boundary rejects");
    peak.maximum_native_bytes = perquery;
    check(bao::evaluate(cq(-1), bao::Ruler(10000),
                        std::span<const bao::Query>(qs.data(), 1), peak)
                  .slots.size() == 1,
          "exact payload cap boundary admits");
    std::printf("BAO owner checks %d PASS max_budget_fraction %.17g\n", checks,
                worst);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL %s after %d\n", e.what(), checks);
    return 1;
  }
}
