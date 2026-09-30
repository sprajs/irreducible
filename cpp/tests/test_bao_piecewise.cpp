// Native piecewise BAO owner evidence. Shared kernels/libm are not independent
// ancestry. The legacy transcript captures named observable/density bits before
// extraction; original-data comparisons are optional and SHA-guarded
// externally.
#include "fixtures/bao_reference.hpp"
#include "fixtures/ldlt_reference.hpp"
#include "irred/bao.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
using namespace irred;
namespace {
unsigned checks = 0;
long double worst_observable_fraction = 0, worst_density_error = 0;
void check(bool ok, const char *why) {
  ++checks;
  if (!ok)
    throw std::runtime_error(why);
}
void observable(double actual, long double expected) {
  const auto error = std::abs((long double)actual - expected);
  const auto budget = 2e-12L + 2e-10L * std::abs(expected);
  worst_observable_fraction =
      std::max(worst_observable_fraction, error / budget);
  check(error <= budget, "observable budget");
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
std::vector<bao::PiecewiseModelPoint> grid() {
  std::vector<bao::PiecewiseModelPoint> result;
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
long double reference_prediction(bao::Query query,
                                 const bao::PiecewiseModelPoint &p,
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
bao::DensityPolicy density_policy() {
  bao::DensityPolicy p;
  p.maximum_models = 24;
  p.maximum_matrix_elements = 9;
  p.maximum_string_bytes = 4096;
  p.maximum_native_bytes = 1048576;
  p.maximum_forward_sensitivity = 1e-10;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.observables.background.integration = {1e-12, 1e-12, 100000, 30};
  p.observables.background.maximum_total_evaluations = 20000000;
  return p;
}
bao::PiecewiseDensityPolicy analytic_policy() {
  bao::PiecewiseDensityPolicy p;
  p.maximum_models = 24;
  p.maximum_native_bytes = 1048576;
  p.maximum_forward_sensitivity = 1e-10;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.observables.background.maximum_segment_visits = 2000;
  return p;
}
void ordinary() {
  auto in = tiny_source();
  auto prepared = bao::prepare_density(in, density_policy());
  check(prepared.status() == statistics::DensityStatus::finite,
        "synthetic preparation");
  auto models = grid();
  auto policy = analytic_policy();
  auto batch = prepared.evaluate_piecewise(models, policy);
  check(batch.slots.size() == 24 &&
            batch.status == statistics::DensityStatus::finite,
        "24 synthetic models");
  std::size_t used = 0;
  // Exact rational C determinant544 and adjugate; entirely separate from
  // factorization.
  constexpr long double adj[9]{140, -16, 2, -16, 64, -8, 2, -8, 35};
  for (unsigned m = 0; m < 24; ++m) {
    const auto &row = batch.slots[m];
    check(row.source.q == models[m].q &&
              row.source.ruler.h0_rd_km_s == models[m].ruler.h0_rd_km_s,
          "full attempted identity");
    check(row.background_status == cosmology::Status::ok &&
              row.result.density.status == statistics::DensityStatus::finite,
          "finite synthetic density");
    check(row.predictions.size() == 3 && row.residuals.size() == 3,
          "ordered arrays");
    std::array<long double, 3> residual{};
    long double q = 0;
    for (unsigned i = 0; i < 3; ++i) {
      auto a = reference_prediction(in.queries[i], models[m], 32),
           b = reference_prediction(in.queries[i], models[m], 64);
      check(std::abs(a - b) <= .1L * (2e-12L + 2e-10L * std::abs(b)),
            "reference refinement");
      observable(row.predictions[i], b);
      residual[i] = (long double)in.observed[i] - b;
      check(row.residuals[i] == in.observed[i] - row.predictions[i],
            "canonical residual");
    }
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        q += residual[i] * adj[3 * i + j] * residual[j] / 544;
    const auto expected =
        -.5L * (q + std::log(544.L) + 3 * std::log(2 * std::acos(-1.L)));
    auto error = std::abs((long double)row.result.density.log_value - expected);
    worst_density_error = std::max(worst_density_error, error);
    check(error <= 1e-8L, "normalized synthetic density budget");
    used += row.segment_visits;
  }
  check(used == batch.segment_visits && used == 240, "total analytic visits");
  // All-equal q analytic endpoint controls: no reference quadrature required.
  for (double q : {-3., -1., 0., .5, 2.}) {
    auto bg = cosmology::prepare_piecewise_q({70, {q, q, q, q, q}});
    std::array<bao::Query, 3> queries{
        {{1, bao::Observable::transverse_over_ruler},
         {1, bao::Observable::hubble_over_ruler},
         {1, bao::Observable::volume_over_ruler}}};
    auto values = bao::evaluate_piecewise(bg, bao::Ruler(10000), queries,
                                          policy.observables);
    long double I = q == -3   ? 7.L / 3
                    : q == -1 ? 1
                    : q == 0  ? std::log(2.L)
                    : q == .5 ? 2 * (1 - 1 / std::sqrt(2.L))
                              : 3.L / 8;
    const auto E = std::pow(2.L, 1 + (long double)q), dm = 29.9792458L * I,
               dh = 29.9792458L / E;
    observable(values.slots[0].dimensionless_value, dm);
    observable(values.slots[1].dimensionless_value, dh);
    observable(values.slots[2].dimensionless_value, std::cbrt(dm * dm * dh));
  }
  auto limited = policy;
  limited.observables.background.maximum_segment_visits = 10;
  auto capped =
      prepared.evaluate_piecewise(std::span(models).first(2), limited);
  check(capped.slots[0].result.density.status ==
            statistics::DensityStatus::finite,
        "first admitted model");
  check(capped.slots[1].numerical_status == numerics::Status::work_limit &&
            capped.segment_visits == 10,
        "global segments failure");
  check(capped.slots[1].predictions.empty() &&
            capped.slots[1].residuals.empty(),
        "failed payload absent");
  bao::PiecewiseDensityPolicy unset;
  check(prepared.evaluate_piecewise(models, unset).slots.empty(),
        "unset precision rejected");
  auto bad = policy;
  bad.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
  check(prepared.evaluate_piecewise(models, bad).slots.empty(),
        "precision mismatch no fallback");
  bad = policy;
  bad.maximum_models = 0;
  auto resource = prepared.evaluate_piecewise(models, bad);
  check(resource.slots.empty() &&
            resource.numerical_status == numerics::Status::work_limit,
        "count cap");
  bad = policy;
  bad.maximum_native_bytes = 0;
  check(prepared.evaluate_piecewise(models, bad).numerical_status ==
            numerics::Status::work_limit,
        "byte cap");
  std::array<bao::PiecewiseModelPoint, 3> invalid{
      {{{3, 0, 0, 0, 0}, bao::Ruler(10000)},
       {{0, 0, 0, 0, 0}, bao::Ruler(4999)},
       {{std::numeric_limits<double>::quiet_NaN(), 0, 0, 0, 0},
        bao::Ruler(10000)}}};
  auto failures = prepared.evaluate_piecewise(invalid, policy);
  for (const auto &row : failures.slots)
    check(row.result.density.status != statistics::DensityStatus::finite &&
              row.predictions.empty() && row.segment_visits == 0,
          "invalid row no science/work");
  std::array<bao::Query, 2> out{{{2.6, bao::Observable::transverse_over_ruler},
                                 {.1, static_cast<bao::Observable>(99)}}};
  auto bg = cosmology::prepare_piecewise_q({70, {0, 0, 0, 0, 0}});
  auto outside =
      bao::evaluate_piecewise(bg, bao::Ruler(10000), out, policy.observables);
  check(outside.slots[0].status == cosmology::Status::unsupported_domain &&
            outside.slots[1].status == cosmology::Status::invalid_input,
        "domain/enum causes");
}

bao::DensityInput original_source(const char *mean, const char *cov) {
  bao::DensityInput in;
  std::ifstream table(mean), matrix(cov);
  check(bool(table) && bool(matrix), "original files available");
  std::string line;
  while (std::getline(table, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream row(line);
    double z, value;
    std::string kind, extra;
    check(bool(row >> z >> value >> kind) && !(row >> extra),
          "original mean row");
    bao::Observable o;
    if (kind == "DM_over_rs")
      o = bao::Observable::transverse_over_ruler;
    else if (kind == "DH_over_rs")
      o = bao::Observable::hubble_over_ruler;
    else if (kind == "DV_over_rs")
      o = bao::Observable::volume_over_ruler;
    else
      throw std::runtime_error("unknown original role");
    in.queries.push_back({z, o});
    in.observed.push_back(value);
    in.ordered_ids.push_back("DESI-DR2-ALL-GCcomb:row:" +
                             std::to_string(in.queries.size() - 1));
  }
  double value;
  while (matrix >> value)
    in.covariance.push_back(value);
  check(matrix.eof() && in.queries.size() == 13 && in.covariance.size() == 169,
        "original shape");
  check(in.queries[11].observable == bao::Observable::hubble_over_ruler &&
            in.queries[12].observable ==
                bao::Observable::transverse_over_ruler &&
            in.queries[11].z == 2.33 && in.queries[12].z == 2.33,
        "released last DH/DM order");
  in.role = bao::RowRole::released_fitted_distance_summary;
  in.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  in.table_identity = std::string(test_reference::bao_mean_sha256);
  in.covariance_identity = std::string(test_reference::bao_cov_sha256);
  in.ordering_provenance =
      "declared original13 source order; external SHA guard assertion";
  in.calibration_provenance =
      "released free ruler distance summaries; no early ruler calibration";
  in.dependence_provenance =
      "supplied full BAO covariance; other probe dependence unknown";
  return in;
}
long double gamma_count(size_t count) {
  const auto v = count * std::numeric_limits<long double>::epsilon();
  return v / (1 - v);
}
long double ulp(double x) {
  return std::max(
      std::abs((long double)std::nextafter(
                   x, std::numeric_limits<double>::infinity()) -
               x),
      std::abs((long double)x -
               std::nextafter(x, -std::numeric_limits<double>::infinity())));
}
long double l1(std::span<const long double> a) {
  long double sum = 0;
  for (auto x : a)
    sum += std::abs(x);
  return sum;
}
long double linf(std::span<const long double> a) {
  long double norm = 0;
  for (auto x : a)
    norm = std::max(norm, std::abs(x));
  return norm;
}
struct ReferenceContext {
  const bao::DensityInput &input;
  test_reference::LDLT factor;
  long double inverse_norm = 0, rho = 0, logdet_allowance = 0,
              normalization = 0;
  explicit ReferenceContext(const bao::DensityInput &in)
      : input(in),
        factor(test_reference::factor(in.covariance, in.queries.size())) {
    check(factor.valid, "canonical C reference factor");
    const auto n = factor.n;
    std::vector<long double> B(n * n), basis(n);
    for (size_t j = 0; j < n; ++j) {
      basis[j] = 1;
      auto x = test_reference::solve_longdouble(factor, basis);
      basis[j] = 0;
      for (size_t i = 0; i < n; ++i)
        B[i * n + j] = x[i];
    }
    long double bnorm = 0, enorm = 0;
    const auto g = gamma_count(3 * n + 3);
    for (size_t i = 0; i < n; ++i) {
      long double row = 0, rrow = 0, erow = 0;
      for (size_t j = 0; j < n; ++j) {
        row += std::abs(B[i * n + j]);
        long double cb = 0, acb = 0, recon = 0, arecon = 0;
        for (size_t k = 0; k < n; ++k) {
          auto t = (long double)input.covariance[i * n + k] * B[k * n + j];
          cb += t;
          acb += std::abs(t);
          auto v = factor.lower[i * n + k] * factor.diagonal[k] *
                   factor.lower[j * n + k];
          recon += v;
          arecon += std::abs(v);
        }
        rrow += std::abs((i == j ? 1.L : 0.L) - cb) + g * (acb + 1);
        erow += std::abs((long double)input.covariance[i * n + j] - recon) +
                g * (arecon + std::abs(input.covariance[i * n + j]));
      }
      bnorm = std::max(bnorm, row);
      rho = std::max(rho, rrow);
      enorm = std::max(enorm, erow);
    }
    check(rho < 1, "inverse residual finite estimate");
    inverse_norm = bnorm / (1 - rho);
    const auto perturbation = inverse_norm * enorm;
    check(perturbation < 1, "determinant perturbation estimate");
    long double abslog = 0;
    for (auto d : factor.diagonal)
      abslog += std::abs(std::log(d));
    logdet_allowance =
        n * (-std::log1p(-perturbation)) +
        64 * n * std::numeric_limits<long double>::epsilon() * (1 + abslog);
    normalization = n * std::log(2 * std::acos(-1.L));
  }
  struct Solve {
    std::vector<long double> x;
    long double q = 0, density = 0, dx = 0, absdot = 0;
  };
  Solve solve(std::span<const long double> rhs) const {
    Solve result;
    result.x = test_reference::solve_longdouble(factor, rhs);
    const auto n = factor.n;
    long double residual_inf = 0;
    const auto g = gamma_count(3 * n + 3);
    for (size_t i = 0; i < n; ++i) {
      const auto t = rhs[i] * result.x[i];
      result.q += t;
      result.absdot += std::abs(t);
      long double cx = 0, acx = 0;
      for (size_t j = 0; j < n; ++j) {
        const auto p = (long double)input.covariance[i * n + j] * result.x[j];
        cx += p;
        acx += std::abs(p);
      }
      residual_inf = std::max(residual_inf, std::abs(rhs[i] - cx) +
                                                g * (acx + std::abs(rhs[i])));
    }
    result.dx = inverse_norm * residual_inf;
    result.density = -.5L * (result.q + factor.logdet + normalization);
    return result;
  }
};
// Transport-only native transcript. Scientific acceptance is the independently
// directed interval gate; the older failed conservative ledger remains intact.
void direct_native(const char *mean, const char *cov) {
  auto input = original_source(mean, cov);
  // Transport adapter alone adopts the structural reader's canonical identity.
  // The scientific original_source reader and its accepted evidence are
  // unchanged.
  for (size_t i = 0; i < input.queries.size(); ++i) {
    const auto tag = input.queries[i].observable;
    const char *label =
        tag == bao::Observable::transverse_over_ruler ? "DM_over_rs"
        : tag == bao::Observable::hubble_over_ruler   ? "DH_over_rs"
                                                      : "DV_over_rs";
    input.ordered_ids[i] =
        input.table_identity + ":row:" + std::to_string(i) + ":" + label;
  }
  input.ordering_provenance =
      "supplied released original row order; not independently verified from "
      "unlabeled covariance";
  input.calibration_provenance =
      "released fitted-distance Gaussian; free empirical ruler";
  input.dependence_provenance =
      "internal covariance supplied; external probe dependence not assessed";

  auto prep = density_policy();
  prep.maximum_matrix_elements = 169;
  auto prepared = bao::prepare_density(input, prep);
  check(prepared.status() == statistics::DensityStatus::finite,
        "direct prepare");
  auto points = grid();
  auto batch = prepared.evaluate_piecewise(points, analytic_policy());
  check(batch.slots.size() == 24, "direct24 rows");
  std::printf(
      "{\"suite\":\"piecewise_BAO_direct_native_transport\",\"asset_identity_"
      "provenance\":\"external Rust SHA guard assertion; C++ does not "
      "hash\",\"mean_sha256\":\"%s\",\"covariance_sha256\":\"%s\",\"status\":%"
      "u,\"numerical_status\":%u,\"segments\":%zu,\"source\":{",
      test_reference::bao_mean_sha256.data(),
      test_reference::bao_cov_sha256.data(), (unsigned)batch.status,
      (unsigned)batch.numerical_status, batch.segment_visits);
  std::printf(
      "\"table_identity\":\"%s\",\"covariance_identity\":\"%s\","
      "\"ordering_provenance\":\"%s\",\"calibration_provenance\":\"%s\","
      "\"dependence_provenance\":\"%s\",\"role\":%u,\"covariance_unit\":%u,"
      "\"ordered_ids\":[",
      input.table_identity.c_str(), input.covariance_identity.c_str(),
      input.ordering_provenance.c_str(), input.calibration_provenance.c_str(),
      input.dependence_provenance.c_str(), (unsigned)input.role,
      (unsigned)input.covariance_unit);
  for (size_t i = 0; i < input.ordered_ids.size(); ++i)
    std::printf("%s\"%s\"", i ? "," : "", input.ordered_ids[i].c_str());
  std::printf("],\"queries\":[");
  for (size_t i = 0; i < input.queries.size(); ++i)
    std::printf("%s{\"z\":%.17g,\"observable\":%u}", i ? "," : "",
                input.queries[i].z, (unsigned)input.queries[i].observable);
  std::printf("],\"observed\":[");
  for (size_t i = 0; i < input.observed.size(); ++i)
    std::printf("%s%.17g", i ? "," : "", input.observed[i]);
  std::printf("],\"covariance\":[");
  for (size_t i = 0; i < input.covariance.size(); ++i)
    std::printf("%s%.17g", i ? "," : "", input.covariance[i]);
  std::printf("]},\"rows\":[");
  for (size_t j = 0; j < batch.slots.size(); ++j) {
    auto &s = batch.slots[j];
    check(s.result.density.status == statistics::DensityStatus::finite,
          "direct finite");
    if (j)
      std::printf(",");
    std::printf("{\"point\":%zu,\"q\":[", j);
    for (size_t k = 0; k < 5; ++k)
      std::printf("%s%.17g", k ? "," : "", s.source.q[k]);
    std::printf("],\"h0rd\":%.17g,\"background_status\":%u,\"numerical_"
                "status\":%u,\"density_status\":%u,\"segments\":%zu,"
                "\"density\":%.17g,\"quadratic\":%.17g,\"logdet\":%.17g,"
                "\"normalization\":%.17g,\"backward_residual\":%.17g,\"forward_"
                "sensitivity\":%.17g,\"prediction\":[",
                s.source.ruler.h0_rd_km_s, (unsigned)s.background_status,
                (unsigned)s.numerical_status, (unsigned)s.result.density.status,
                s.segment_visits, s.result.density.log_value,
                s.result.quadratic, s.result.log_determinant,
                s.result.normalization, s.result.backward_residual,
                s.result.estimated_forward_sensitivity);
    for (size_t k = 0; k < s.predictions.size(); ++k)
      std::printf("%s%.17g", k ? "," : "", s.predictions[k]);
    std::printf("],\"residual\":[");
    for (size_t k = 0; k < s.residuals.size(); ++k)
      std::printf("%s%.17g", k ? "," : "", s.residuals[k]);
    std::printf("]}");
  }
  std::printf("]}\n");
}
void original(const char *mean, const char *cov) {
  check(std::numeric_limits<long double>::digits >= 64 &&
            std::fegetround() == FE_TONEAREST,
        "reference arithmetic domain");
  auto input = original_source(mean, cov);
  auto prep = density_policy();
  prep.maximum_matrix_elements = 169;
  auto prepared = bao::prepare_density(input, prep);
  check(prepared.status() == statistics::DensityStatus::finite,
        "original preparation");
  auto models = grid();
  auto policy = analytic_policy();
  auto batch = prepared.evaluate_piecewise(models, policy);
  check(batch.slots.size() == 24, "original24 materialized attempts");
  ReferenceContext ref(input);
  const auto n = input.queries.size();
  const auto g = gamma_count(2 * n + 2),
             eps = std::numeric_limits<long double>::epsilon();
  unsigned failures = 0;
  size_t used = 0;
  for (size_t m = 0; m < models.size(); ++m) {
    const auto &s = batch.slots[m];
    used += s.segment_visits;
    if (s.result.density.status != statistics::DensityStatus::finite) {
      ++failures;
      std::printf("{\"point\":%zu,\"finite\":false,\"numerical_status\":%u,"
                  "\"segments\":%zu}\n",
                  m, (unsigned)s.numerical_status, s.segment_visits);
      continue;
    }
    std::vector<long double> r32(n), r64(n), r128(n), rc(n), p128(n), delta(n);
    long double observable_fraction = 0, refinement_fraction = 0;
    for (size_t i = 0; i < n; ++i) {
      auto a = reference_prediction(input.queries[i], models[m], 32),
           b = reference_prediction(input.queries[i], models[m], 64),
           c = reference_prediction(input.queries[i], models[m], 128);
      p128[i] = c;
      r32[i] = (long double)input.observed[i] - a;
      r64[i] = (long double)input.observed[i] - b;
      r128[i] = (long double)input.observed[i] - c;
      rc[i] = s.residuals[i];
      delta[i] = std::abs(rc[i] - r128[i]) + std::abs(c - b) +
                 .5L * ulp((double)c) +
                 eps * (std::abs(input.observed[i]) + std::abs(c));
      const auto budget = 2e-12L + 2e-10L * std::abs(c);
      observable_fraction =
          std::max(observable_fraction,
                   std::abs((long double)s.predictions[i] - c) / budget);
      refinement_fraction =
          std::max(refinement_fraction, std::abs(c - b) / (.1L * budget));
    }
    const auto a = ref.solve(r32), b = ref.solve(r64), c = ref.solve(r128),
               canonical = ref.solve(rc);
    const auto reference_refinement = std::abs(b.density - c.density);
    const auto reference_error = .5L * (l1(r128) * c.dx + g * c.absdot +
                                        ref.logdet_allowance + 64 * n * eps) +
                                 reference_refinement;
    const auto eta = (long double)s.result.estimated_forward_sensitivity;
    check(eta >= 0 && eta < 1 && eta <= 1e-10L, "runtime quality guard");
    const auto xnorm = linf(canonical.x) + canonical.dx;
    const auto factor_bound_base =
        .5L * l1(rc) * xnorm * eta / (1 - eta) +
        .5L * g * l1(rc) * xnorm / (1 - eta) + .25L * ulp(s.result.quadratic) +
        .5L * ulp(s.result.density.log_value) +
        .5L * std::abs((long double)s.result.log_determinant -
                       ref.factor.logdet) +
        .5L *
            std::abs((long double)s.result.normalization - ref.normalization) +
        .5L * ulp(s.result.normalization);
    const auto density_reduction =
        .5L * gamma_count(3) *
        (std::abs(canonical.q) +
         std::abs((long double)s.result.log_determinant) +
         std::abs(ref.normalization));
    const auto factor_bound = factor_bound_base + density_reduction;
    const auto factor_direct =
        std::abs((long double)s.result.density.log_value - canonical.density);
    long double weighted = 0;
    for (size_t i = 0; i < n; ++i)
      weighted += std::abs(c.x[i]) * delta[i];
    const auto background_bound =
        weighted + c.dx * l1(delta) +
        .5L * ref.inverse_norm * l1(delta) * linf(delta);
    const auto background_direct = std::abs(canonical.q - c.q) / 2;
    const auto assembled_direct =
        std::abs((long double)s.result.density.log_value - c.density);
    const bool reference_pass =
        reference_error <= 2e-9L && refinement_fraction <= 1;
    const bool factor_pass = factor_bound <= 3e-9L && factor_direct <= 3e-9L;
    const bool background_pass =
        background_bound <= 5e-9L && background_direct <= 5e-9L;
    const bool pass = reference_pass && factor_pass && background_pass &&
                      assembled_direct <= 1e-8L && observable_fraction <= 1;
    if (!pass)
      ++failures;
    std::printf(
        "{\"ledger_point\":%zu,\"solve_term\":%.21Lg,\"LD_dot_term\":%.21Lg,"
        "\"q_report_cast\":%.21Lg,\"density_report_cast\":%.21Lg,\"LD_density_"
        "reduction\":%.21Lg,\"logdet_difference_half\":%.21Lg,\"normalization_"
        "difference_half\":%.21Lg,\"normalization_report_cast_half\":%.21Lg,"
        "\"gradient_solve_correction\":%.21Lg,\"inverse_norm_estimate\":%.21Lg}"
        "\n",
        m, .5L * l1(rc) * xnorm * eta / (1 - eta),
        .5L * g * l1(rc) * xnorm / (1 - eta), .25L * ulp(s.result.quadratic),
        .5L * ulp(s.result.density.log_value), density_reduction,
        .5L *
            std::abs((long double)s.result.log_determinant - ref.factor.logdet),
        .5L * std::abs((long double)s.result.normalization - ref.normalization),
        .5L * ulp(s.result.normalization), c.dx * l1(delta), ref.inverse_norm);
    std::printf(
        "{\"point\":%zu,\"q\":[%.17g,%.17g,%.17g,%.17g,%.17g],\"h0rd\":%.17g,"
        "\"density\":%.17g,\"quadratic\":%.17g,\"logdet\":%.17g,"
        "\"normalization\":%.17g,\"segments\":%zu,\"eta\":%.17g,\"reference_"
        "error\":%.21Lg,\"factor_bound\":%.21Lg,\"factor_direct\":%.21Lg,"
        "\"background_bound\":%.21Lg,\"background_direct\":%.21Lg,\"assembled_"
        "direct\":%.21Lg,\"reference_refinement\":%.21Lg,\"coarse_refinement\":"
        "%.21Lg,\"reference_solve_dx\":%.21Lg,\"canonical_solve_dx\":%.21Lg,"
        "\"observable_fraction\":%.21Lg,\"refinement_fraction\":%.21Lg,"
        "\"reference_pass\":%s,\"factor_pass\":%s,\"background_pass\":%s,"
        "\"pass\":%s,\"prediction\":[",
        m, models[m].q[0], models[m].q[1], models[m].q[2], models[m].q[3],
        models[m].q[4], models[m].ruler.h0_rd_km_s, s.result.density.log_value,
        s.result.quadratic, s.result.log_determinant, s.result.normalization,
        s.segment_visits, s.result.estimated_forward_sensitivity,
        reference_error, factor_bound, factor_direct, background_bound,
        background_direct, assembled_direct, reference_refinement,
        std::abs(a.density - b.density), c.dx, canonical.dx,
        observable_fraction, refinement_fraction,
        reference_pass ? "true" : "false", factor_pass ? "true" : "false",
        background_pass ? "true" : "false", pass ? "true" : "false");
    for (size_t i = 0; i < n; ++i)
      std::printf("%s%.17g", i ? "," : "", s.predictions[i]);
    std::printf("],\"residual\":[");
    for (size_t i = 0; i < n; ++i)
      std::printf("%s%.17g", i ? "," : "", s.residuals[i]);
    std::printf("],\"reference_prediction\":[");
    for (size_t i = 0; i < n; ++i)
      std::printf("%s%.21Lg", i ? "," : "", p128[i]);
    std::printf("]}\n");
  }
  check(used == batch.segment_visits && used <= 2000, "actual total work");
  std::printf("{\"suite\":\"piecewise_BAO_original24\",\"failures\":%u,"
              "\"segments\":%zu,\"inverse_norm_estimate\":%.21Lg,\"inverse_"
              "residual_estimate\":%.21Lg,\"logdet_allowance\":%.21Lg}\n",
              failures, used, ref.inverse_norm, ref.rho, ref.logdet_allowance);
  check(failures == 0, "original24 unchanged allocated gates");
}
void legacy_transcript() {
  auto in = tiny_source();
  auto p = density_policy();
  auto prepared = bao::prepare_density(in, p);
  std::array<bao::ModelQuery, 9> models{
      {{cosmology::prepare({cosmology::Model::flat_lcdm_late_v1, 70, 0, 0}),
        bao::Ruler(5000)},
       {cosmology::prepare({cosmology::Model::flat_lcdm_late_v1, 70, .3, 0}),
        bao::Ruler(10000)},
       {cosmology::prepare({cosmology::Model::flat_lcdm_late_v1, 70, 1, 0}),
        bao::Ruler(15000)},
       {cosmology::prepare({cosmology::Model::constant_q_flat_v1, 70, 0, -1}),
        bao::Ruler(5000)},
       {cosmology::prepare({cosmology::Model::constant_q_flat_v1, 70, 0, 0}),
        bao::Ruler(10000)},
       {cosmology::prepare({cosmology::Model::constant_q_flat_v1, 70, 0, .5}),
        bao::Ruler(15000)},
       {cosmology::prepare_cpl({70, .3, -1, 0}), bao::Ruler(5000)},
       {cosmology::prepare_cpl({70, .3, -.8, .2}), bao::Ruler(10000)},
       {cosmology::prepare_cpl({70, .3, -1.2, -.3}), bao::Ruler(15000)}}};
  auto batch = prepared.evaluate(models, p);
  std::printf("%u %u %zu %zu\n", (unsigned)batch.status,
              (unsigned)batch.numerical_status, batch.evaluations,
              batch.slots.size());
  for (const auto &s : batch.slots) {
    std::printf("%u %u %u %u %zu\n", (unsigned)s.background_status,
                (unsigned)s.numerical_status, (unsigned)s.result.density.status,
                (unsigned)s.result.density.numerical_status, s.evaluations);
    for (double v :
         {s.result.density.log_value, s.result.quadratic,
          s.result.log_determinant, s.result.normalization,
          s.result.backward_residual, s.result.estimated_forward_sensitivity})
      std::printf("%a ", v);
    for (double v : s.predictions)
      std::printf("%a ", v);
    for (double v : s.residuals)
      std::printf("%a ", v);
    std::printf("\n");
  }
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--legacy-transcript") {
      legacy_transcript();
      return 0;
    }
    if (argc == 5 && std::string(argv[1]) == "--verified-original-assets" &&
        std::string(argv[4]) == "--direct-native-transcript") {
      direct_native(argv[2], argv[3]);
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "--verified-original-assets" &&
        std::getenv("IRRED_BAO_DIRECT_NATIVE_TRANSCRIPT")) {
      direct_native(argv[2], argv[3]);
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "--verified-original-assets") {
      original(argv[2], argv[3]);
      return 0;
    }
    if (argc != 1)
      throw std::runtime_error("unknown test mode");
    ordinary();
    std::printf("piecewise BAO owner %u PASS observable_fraction %.18Lg "
                "density_error %.18Lg\n",
                checks, worst_observable_fraction, worst_density_error);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
