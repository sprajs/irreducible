// Owned piecewise SN validation. Reference conservation equations/libm are
// shared; GL8 fixed redshift panels and LDLT are separate numerical routes, not
// an independent author or interval certificate. Synthetic C=[[4,1],[1,9]]
// gives profiled q=(r0-r1)^2/11. Score/offset budget1e-10 on this tiny case.
// Optional original mode must be invoked by the Rust SHA256 guard; no CI data.
#include "fixtures/ldlt_reference.hpp"
#include "fixtures/w01_historical.hpp"
#include "fixtures/w01_piecewise.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
using namespace irred;
namespace {
unsigned checks = 0;
long double max_error = 0;
void check(bool v, const char *why) {
  ++checks;
  if (!v)
    throw std::runtime_error(why);
}
void near(long double got, long double ref, long double budget,
          const char *why) {
  max_error = std::max(max_error, std::abs(got - ref));
  check(std::abs(got - ref) <= budget, why);
}
observations::Prepared synthetic() {
  observations::Input in{};
  in.profile = observations::Profile::gaussian_fixture_v1;
  in.role = observations::Role::synthetic_control;
  in.unit = observations::Unit::magnitude;
  in.calibration = observations::Calibration::not_applicable;
  in.uncertainty = observations::Uncertainty::covariance;
  in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  in.component = observations::Component::total;
  in.ordering_provenance = "generated exact source order";
  in.table_sha256 = std::string(64, 'a');
  in.uncertainty_sha256 = std::string(64, 'b');
  in.measurement_ids = {"generated:A", "generated:B", "generated:excluded"};
  in.event_ids = in.measurement_ids;
  in.uncertainty_axis_ids = in.measurement_ids;
  in.values = {static_cast<double>(5 * std::log10(1.9L) + 12),
               static_cast<double>(5 * std::log10(6.2L) + 7), 0};
  in.uncertainty_matrix = {4, 1, 8, 1, 9, 5, 7, 5, -100};
  in.missing = {0, 0, 0};
  in.source_selection = {1, 1, 1};
  in.quality = {0, 0, 0};
  return observations::prepare(std::move(in), {3, 9, 4096});
}
std::vector<std::string> words(const std::string &line) {
  std::istringstream in(line);
  std::vector<std::string> out;
  for (std::string x; in >> x;)
    out.push_back(x);
  return out;
}
observations::Prepared original(const char *cov_path, const char *table_path) {
  observations::Input in{};
  in.profile = observations::Profile::pantheon_plus_released_v1;
  in.role = observations::Role::released_fitted_summary;
  in.unit = observations::Unit::magnitude;
  in.calibration = observations::Calibration::released_corrected;
  in.uncertainty = observations::Uncertainty::covariance;
  in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  in.component = observations::Component::total;
  in.ordering_provenance =
      "exact original release order, SHA256 guard verified before/after";
  in.quality_dictionary =
      "absent decoded quality; zeros are structural placeholders";
  in.table_sha256 = std::string(test_reference::w01_table_sha256);
  in.uncertainty_sha256 = std::string(test_reference::w01_cov_sha256);
  std::ifstream table(table_path), cov(cov_path);
  std::string line;
  std::getline(table, line);
  auto header = words(line);
  auto column = [&](const char *name) {
    auto it = std::find(header.begin(), header.end(), name);
    if (it == header.end())
      throw std::runtime_error("missing release column");
    return static_cast<std::size_t>(it - header.begin());
  };
  auto zm = column("zHD"), zh = column("zHEL"), zc = column("zCMB"),
       mag = column("m_b_corr"), cid = column("CID"),
       survey = column("IDSURVEY");
  while (std::getline(table, line)) {
    if (line.empty())
      continue;
    auto row = words(line);
    if (row.size() != header.size())
      throw std::runtime_error("table shape");
    auto i = in.values.size();
    in.values.push_back(std::stod(row[mag]));
    in.zhd.push_back(std::stod(row[zm]));
    in.zhel.push_back(std::stod(row[zh]));
    in.zcmb.push_back(std::stod(row[zc]));
    in.measurement_ids.push_back(row[cid] + "/" + row[survey] +
                                 "/original-row:" + std::to_string(i));
    in.event_ids.push_back(row[cid]);
  }
  std::size_t n;
  if (!(cov >> n) || n != 1701 || in.values.size() != n)
    throw std::runtime_error("release shape");
  in.uncertainty_matrix.resize(n * n);
  for (auto &x : in.uncertainty_matrix)
    if (!(cov >> x) || !std::isfinite(x))
      throw std::runtime_error("covariance entry");
  if (cov >> line)
    throw std::runtime_error("covariance trailing values");
  in.uncertainty_axis_ids = in.measurement_ids;
  in.missing.assign(n, 0);
  in.zhd_missing.assign(n, 0);
  in.zhel_missing.assign(n, 0);
  in.zcmb_missing.assign(n, 0);
  in.source_selection.assign(n, 1);
  in.quality.assign(n, 0);
  return observations::prepare(std::move(in), {1701, 1701 * 1701, 1048576});
}
struct ReferenceScore {
  long double q, a;
};
ReferenceScore score(const test_reference::LDLT &f,
                     std::span<const long double> residual, long double gram) {
  auto wr = test_reference::solve_longdouble(f, residual);
  long double numerator = 0;
  for (auto x : wr)
    numerator += x;
  auto a = numerator / gram;
  std::vector<long double> u(residual.size());
  for (std::size_t i = 0; i < u.size(); ++i)
    u[i] = residual[i] - a;
  auto wu = test_reference::solve_longdouble(f, u);
  long double q = 0;
  for (std::size_t i = 0; i < u.size(); ++i)
    q += u[i] * wu[i];
  return {q, a};
}
// Independent split scale-factor GL8; multiplicative E, same canonical q/edges.
long double radial(long double z, const std::array<double, 5> &q,
                   unsigned panels) {
  constexpr long double x[]{
      .183434642495649804939476142360L, .525532409916328985817739049189L,
      .796666477413626739591553936476L, .960289856497536231683560868569L};
  constexpr long double w[]{
      .362683783378361982965150449277L, .313706645877887287337962201987L,
      .222381034453374470544355994426L, .101228536290376259152531354310L};
  long double start = 1, sum = 0;
  for (unsigned k = 0; k < 5; ++k) {
    long double lo = cosmology::piecewise_q_edges[k],
                hi = std::min(z,
                              (long double)cosmology::piecewise_q_edges[k + 1]);
    if (hi > lo) {
      auto a0 = 1 / (1 + hi), a1 = 1 / (1 + lo), h = (a1 - a0) / panels;
      for (unsigned j = 0; j < panels; ++j)
        for (unsigned t = 0; t < 4; ++t)
          for (int sign : {-1, 1}) {
            auto a = a0 + (j + .5L) * h + sign * x[t] * h / 2;
            auto E =
                start * std::pow(1 / (a * (1 + lo)), 1 + (long double)q[k]);
            sum += h * w[t] / (2 * a * a * E);
          }
    }
    if (z <= cosmology::piecewise_q_edges[k + 1])
      break;
    start *= std::pow((1 + (long double)cosmology::piecewise_q_edges[k + 1]) /
                          (1 + lo),
                      1 + (long double)q[k]);
  }
  return sum;
}
void actual(const char *cov, const char *table) {
  auto source = original(cov, table);
  check(source.status() == observations::Status::ok, "original typed source");
  supernova::Policy policy;
  policy.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  policy.maximum_forward_sensitivity = 1e-10;
  policy.background.maximum_total_evaluations = 20000000;
  auto consumer = supernova::prepare(std::move(source), policy);
  check(consumer.status() == supernova::Status::ok,
        "original wide preparation");
  auto indices = consumer.selected_source_indices();
  check(indices.size() == 1590, "original selection");
  auto n = indices.size();
  const auto &in = consumer.observations().source();
  std::vector<double> C(n * n);
  long double Cnorm = 0;
  for (std::size_t i = 0; i < n; ++i) {
    long double row = 0;
    for (std::size_t j = 0; j < n; ++j) {
      C[i * n + j] = in.uncertainty_matrix[indices[i] * 1701 + indices[j]];
      row += std::abs((long double)C[i * n + j]);
    }
    Cnorm = std::max(Cnorm, row);
  }
  auto ref = test_reference::factor(C, n);
  check(ref.valid, "reference LDLT");
  std::vector<long double> ones(n, 1), inverse_rows(n), basis(n);
  auto wone = test_reference::solve_longdouble(ref, ones);
  long double gram = 0;
  for (auto x : wone)
    gram += x;
  for (std::size_t j = 0; j < n; ++j) {
    basis[j] = 1;
    auto x = test_reference::solve_longdouble(ref, basis);
    basis[j] = 0;
    for (std::size_t i = 0; i < n; ++i)
      inverse_rows[i] += std::abs(x[i]);
  }
  auto inverse_norm =
      *std::max_element(inverse_rows.begin(), inverse_rows.end());
  std::array<supernova::PiecewiseModelPoint, 6> grid{
      supernova::PiecewiseModelPoint({0, 0, 0, 0, 0}),
      supernova::PiecewiseModelPoint({-1, -1, -1, -1, -1}),
      supernova::PiecewiseModelPoint({.5, .5, .5, .5, .5}),
      supernova::PiecewiseModelPoint({-.4, -.4, -.2, .1, .3}),
      supernova::PiecewiseModelPoint({-1, 0, -1, 0, -1}),
      supernova::PiecewiseModelPoint({-3, 2, -3, 2, -3})};
  supernova::PiecewiseEvaluationPolicy evaluation;
  evaluation.arithmetic = policy.arithmetic;
  evaluation.maximum_forward_sensitivity = 1e-10;
  evaluation.maximum_models = 6;
  evaluation.background.maximum_segment_visits = 50000;
  auto batch = consumer.evaluate_piecewise_batch(grid, evaluation);
  check(batch.slots.size() == 6, "original six points");
  long double gamma = (2 * n + 2) * std::numeric_limits<long double>::epsilon();
  gamma /= 1 - gamma;
  long double eps64 = std::numeric_limits<double>::epsilon();
  for (unsigned m = 0; m < 6; ++m) {
    const auto &p = grid[m];
    const auto &s = batch.slots[m];
    check(s.status == supernova::Status::ok, "original finite point");
    std::vector<long double> rcoarse(n), r32(n), r64(n), canonical(n);
    long double mu_error = 0;
    for (std::size_t i = 0; i < n; ++i) {
      auto z = in.zhd[indices[i]];
      auto pref = 1 + (long double)in.zhel[indices[i]];
      auto mu32 = 5 * std::log10(pref * radial(z, p.q, 64)),
           mu64 = 5 * std::log10(pref * radial(z, p.q, 128));
      rcoarse[i] = (long double)in.values[indices[i]] -
                   5 * std::log10(pref * radial(z, p.q, 32));
      r32[i] = (long double)in.values[indices[i]] - mu32;
      r64[i] = (long double)in.values[indices[i]] - mu64;
      canonical[i] = s.base_residuals[i];
      mu_error = std::max(mu_error,
                          std::abs((long double)s.shape_magnitudes[i] - mu64));
    }
    auto coarse = score(ref, rcoarse, gram);
    auto a = score(ref, r32, gram), b = score(ref, r64, gram),
         same = score(ref, canonical, gram);
    auto reference_refinement = std::abs(a.q - b.q) / 2,
         factor_error = std::abs((long double)s.quadratic - same.q) / 2,
         background_error = std::abs(same.q - b.q) / 2;
    const auto &d = s.solve_diagnostics;
    auto wx = consumer.profile_operator().cached_response_solution();
    long double wxnorm = 0, h = 0, hdot = 0, du1 = 0, duinf = 0;
    for (std::size_t i = 0; i < n; ++i) {
      wxnorm = std::max(wxnorm, std::abs((long double)wx[i]));
      h += (long double)wx[i] * s.profiled_residuals[i];
      hdot += std::abs((long double)wx[i] * s.profiled_residuals[i]);
      auto du =
          std::abs((long double)s.profiled_residuals[i] -
                   ((long double)s.base_residuals[i] - s.offset_coefficient));
      du1 += du;
      duinf = std::max(duinf, du);
    }
    long double
        ex = consumer.profile_operator().cached_response_forward_sensitivity() *
             wxnorm,
        er = d.coefficient_solve_forward_sensitivity * d.solution_norm_inf,
        eu = d.estimated_forward_sensitivity * d.adjusted_solution_norm_inf;
    auto G = consumer.profile_operator().gram();
    auto dg = n * ex + gamma * n * wxnorm,
         db = n * er + gamma * n * d.solution_norm_inf;
    check(G > dg, "positive empirical Gram bound");
    auto da =
        (db + std::abs((long double)s.offset_coefficient) * dg) / (G - dg) +
        eps64 / (1 - eps64) * std::abs((long double)s.offset_coefficient);
    auto dh = d.adjusted_residual_l1 * ex + gamma * hdot;
    auto allowance =
        .5L * d.adjusted_residual_l1 * eu +
        gamma * d.adjusted_residual_l1 * d.adjusted_solution_norm_inf +
        (std::abs(h) + dh) * da + .5L * (G + dg) * da * da +
        du1 * d.adjusted_solution_norm_inf + .5L * du1 * inverse_norm * duinf +
        .5L * eps64 / (1 - eps64) * std::abs((long double)s.quadratic);
    auto score_error =
        std::abs((long double)s.relative_profile_score + b.q / 2);
    bool pass = reference_refinement <= 2e-7L && factor_error <= 3e-7L &&
                allowance <= 3e-7L && background_error <= 5e-7L &&
                mu_error <= 1e-8L && score_error <= 1e-6L;
    std::printf("point %u score %.17g reference %.21Lg factor %.21Lg allowance "
                "%.21Lg background %.21Lg mu %.21Lg total %.21Lg pass %d\n",
                m, s.relative_profile_score, reference_refinement, factor_error,
                allowance, background_error, mu_error, score_error, pass);
    std::printf("point %u reference32_to64 %.21Lg\n", m,
                std::abs(coarse.q - a.q) / 2);
    near(s.relative_profile_score, test_reference::w01_piecewise_scores[m],
         1e-6L, "migrated frozen score");
    near(s.relative_profile_score,
         test_reference::w01_piecewise_historical_stable[m], 2e-7L,
         "historical stable migration");
    near(s.relative_profile_score,
         test_reference::w01_piecewise_historical_full_inverse[m], 2e-7L,
         "historical full migration");
    auto compressed_error =
        std::abs((long double)s.relative_profile_score -
                 test_reference::w01_piecewise_historical_compressed[m]);
    check((compressed_error <= 2e-7L) == (m == 1),
          "historical compressed precise failure retained");
    check(p.q == test_reference::w01_piecewise_q[m],
          "frozen parameter identity");
    check(pass, "original fixedpoint budgets");
  }
  std::printf(
      "{\"suite\":\"supernova_piecewise_original_fixedpoints\",\"n\":%zu,"
      "\"reference_condition_estimate\":%.21Lg,\"passed\":true}\n",
      n, Cnorm * inverse_norm);
}

// Recover baseline supernova.cpp from the pre-refactor public commit and use
// IRRED_LEGACY_TRANSCRIPT_ONLY + function sections/GC to compare this
// transcript.
void legacy_transcript() {
  auto input = synthetic();
  std::array<cosmology::Query, 3> queries{
      {{1, .9, cosmology::Convention::released_zhd_zhel},
       {2, 2.1, cosmology::Convention::released_zhd_zhel},
       {.01, .01, cosmology::Convention::released_zhd_zhel}}};
  supernova::Policy p;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.maximum_forward_sensitivity = 1e-10;
  p.maximum_matrix_elements = 9;
  auto c = supernova::prepare_synthetic(input, input.source().measurement_ids,
                                        queries, p);
  auto emit = [](const supernova::Slot &s) {
    std::printf("%u %u %u %u %zu\n", (unsigned)s.status,
                (unsigned)s.background_status, (unsigned)s.numerical_status,
                (unsigned)s.profile_status, s.background_evaluations);
    const auto d = s.solve_diagnostics;
    for (double v : {s.offset_coefficient, s.quadratic,
                     s.relative_profile_score, d.coefficient, d.quadratic,
                     d.backward_residual, d.estimated_forward_sensitivity,
                     d.coefficient_solve_backward_residual,
                     d.coefficient_solve_forward_sensitivity, d.residual_l1,
                     d.solution_norm_inf, d.adjusted_residual_l1,
                     d.adjusted_solution_norm_inf})
      std::printf("%a\n", v);
    for (const auto *v : {&s.shape_magnitudes, &s.base_residuals,
                          &s.profiled_residuals, &d.adjusted_residuals})
      for (double x : *v)
        std::printf("%a\n", x);
  };
  std::array<supernova::ModelPoint, 4> old{
      {{cosmology::Model::flat_lcdm_late_v1, .3, 0},
       {cosmology::Model::constant_q_flat_v1, 0, 0},
       {cosmology::Model::constant_q_flat_v1, 0, -1},
       {cosmology::Model::flat_cpl_late_v1, .3, 0}}};
  for (const auto &row : c.evaluate_batch(old, p).slots)
    emit(row);
  std::array<supernova::ModelPointV2, 4> v2{
      {{cosmology::Model::flat_cpl_late_v1, .3, 0, -.9, .4},
       {cosmology::Model::flat_lcdm_late_v1, .3, 0, -1, 0},
       {cosmology::Model::constant_q_flat_v1, 0, -1, -1, 0},
       {cosmology::Model::flat_cpl_late_v1, .3, 1, -1, 0}}};
  for (const auto &row : c.evaluate_batch_v2(v2, p).slots)
    emit(row.calculation);
}
void ordinary() {
  static_assert(
      !std::is_default_constructible_v<supernova::PiecewiseModelPoint>);
  auto source = synthetic();
  std::array<cosmology::Query, 3> queries{
      {{1, .9, cosmology::Convention::released_zhd_zhel},
       {2, 2.1, cosmology::Convention::released_zhd_zhel},
       {.01, .01, cosmology::Convention::released_zhd_zhel}}};
  supernova::Policy prep;
  prep.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  prep.maximum_forward_sensitivity = 1e-10;
  prep.maximum_matrix_elements = 9;
  auto consumer = supernova::prepare_synthetic(
      source, source.source().measurement_ids, queries, prep);
  check(consumer.status() == supernova::Status::ok,
        "prepared selected control");
  supernova::PiecewiseEvaluationPolicy policy;
  policy.arithmetic = prep.arithmetic;
  policy.maximum_forward_sensitivity = 1e-10;
  policy.background.maximum_segment_visits = 50000;
  for (double q : {-3., -1., 0., .5, 2.}) {
    std::array<supernova::PiecewiseModelPoint, 1> point{
        supernova::PiecewiseModelPoint({q, q, q, q, q})};
    auto batch = consumer.evaluate_piecewise_batch(point, policy);
    auto &s = batch.slots[0];
    check(s.status == supernova::Status::ok, "analytic constant q profile");
    auto I = [&](long double z) {
      return q == 0 ? std::log1p(z) : -std::expm1(-q * std::log1p(z)) / q;
    };
    long double r0 = source.source().values[0] - 5 * std::log10(1.9L * I(1));
    long double r1 = source.source().values[1] - 5 * std::log10(3.1L * I(2));
    near(s.quadratic, (r0 - r1) * (r0 - r1) / 11, 1e-10,
         "exact adjugate projected quadratic");
    near(s.offset_coefficient, (8 * r0 + 3 * r1) / 11, 1e-10, "exact offset");
    check(s.segment_visits == 9, "analytic visited segments");
  }
  std::array<supernova::PiecewiseModelPoint, 2> points{
      supernova::PiecewiseModelPoint({-1, 0, -1, 0, -1}),
      supernova::PiecewiseModelPoint({0, 0, 0, 0, 0})};
  auto cap = policy;
  cap.background.maximum_segment_visits = 9;
  auto limited = consumer.evaluate_piecewise_batch(points, cap);
  check(limited.slots[0].status == supernova::Status::ok,
        "first point consumes shared cap");
  check(limited.slots[1].status == supernova::Status::work_limit,
        "cap not reset between models");
  check(limited.slots[1].shape_magnitudes.empty(), "failed payload absent");
  points[0].q[2] = 3;
  auto bad = consumer.evaluate_piecewise_batch(points, policy);
  check(bad.slots[0].status != supernova::Status::ok &&
            bad.slots[0].source.q[2] == 3,
        "attempt retained rejected");
  check(bad.slots[1].status == supernova::Status::ok, "valid later row");
  check(limited.segment_visits == 9, "total includes consumed work");
  auto unknown = policy;
  unknown.arithmetic = static_cast<numerics::Arithmetic>(UINT32_MAX);
  check(consumer.evaluate_piecewise_batch(points, unknown).slots.empty(),
        "invalid precision rejected");
  supernova::PiecewiseEvaluationPolicy defaults;
  defaults.maximum_forward_sensitivity = 1e-10;
  check(consumer.evaluate_piecewise_batch(points, defaults).status ==
            supernova::Status::incompatible_metadata,
        "precision required");
  cap = policy;
  cap.background.maximum_queries = 0;
  auto querycap = consumer.evaluate_piecewise_batch(points, cap);
  check(querycap.slots[1].status == supernova::Status::work_limit &&
            querycap.slots[1].numerical_status == numerics::Status::work_limit,
        "query cap precise cause");
  cap = policy;
  cap.maximum_models = 1;
  check(consumer.evaluate_piecewise_batch(points, cap).slots.empty(),
        "model cap before allocation");
}
} // namespace
int main(int argc, char **argv) {
  try {
#ifdef IRRED_LEGACY_TRANSCRIPT_ONLY
    legacy_transcript();
    return 0;
#else
    if (argc == 2 && std::string(argv[1]) == "--legacy-transcript") {
      legacy_transcript();
      return 0;
    }
    if (argc == 4 && std::string(argv[3]) == "--verified-original-assets")
      actual(argv[1], argv[2]);
    else if (argc != 1)
      throw std::runtime_error("original mode requires SHA guard");
    ordinary();
    std::printf("piecewise SN owner PASS %u max_error %.18Lg\n", checks,
                max_error);
#endif
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
