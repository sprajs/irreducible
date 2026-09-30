// Owned CPL SN validation. Reference conservation equations/libm are shared;
// GL8 fixed redshift panels and LDLT are separate numerical routes, not an
// independent author or interval certificate. Synthetic C=[[4,1],[1,9]] gives
// profiled q=(r0-r1)^2/11. Score/offset budget1e-10 on this tiny case.
// Optional original mode must be invoked by the Rust SHA256 guard; no CI data.
#include "fixtures/ldlt_reference.hpp"
#include "fixtures/w01_historical.hpp"
#include "irred/supernova.hpp"
#include <array>
#include <algorithm>
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
long double radial(long double z, double om, double w0, double wa,
                   unsigned panels) {
  constexpr long double x[]{
      .183434642495649804939476142360L, .525532409916328985817739049189L,
      .796666477413626739591553936476L, .960289856497536231683560868569L};
  constexpr long double weight[]{
      .362683783378361982965150449277L, .313706645877887287337962201987L,
      .222381034453374470544355994426L, .101228536290376259152531354310L};
  auto h = z / panels;
  long double sum = 0;
  for (unsigned p = 0; p < panels; ++p)
    for (unsigned k = 0; k < 4; ++k)
      for (int sign : {-1, 1}) {
        auto t = (p + .5L) * h + sign * x[k] * h / 2, u = 1 + t;
        auto dark = std::pow(u, 3 * (1 + (long double)w0 + wa)) *
                    std::exp(-3 * (long double)wa * t / u);
        sum += h / 2 * weight[k] / std::sqrt(om * u * u * u + (1 - om) * dark);
      }
  return sum;
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
  std::array<supernova::ModelPointV2, 4> grid{
      {{cosmology::Model::flat_cpl_late_v1, .3, 0, -1, 0},
       {cosmology::Model::flat_cpl_late_v1, 0, 0, -2. / 3, 0},
       {cosmology::Model::flat_cpl_late_v1, .3, 0, -.9, .4},
       {cosmology::Model::flat_cpl_late_v1, .3, 0, -1.1, -.4}}};
  std::array<supernova::BatchResultV2, 3> batches;
  for (unsigned level = 0; level < 3; ++level) {
    auto eps = level == 0 ? 1e-8 : level == 1 ? 1e-10 : 1e-12;
    policy.background.integration = {eps, eps, 100000, 30};
    batches[level] = consumer.evaluate_batch_v2(grid, policy);
    check(batches[level].slots.size() == 4, "original model batch");
  }
  long double gamma = (2 * n + 2) * std::numeric_limits<long double>::epsilon();
  gamma /= 1 - gamma;
  long double eps64 = std::numeric_limits<double>::epsilon();
  for (unsigned m = 0; m < 4; ++m) {
    const auto &p = grid[m];
    const auto &s = batches[2].slots[m].calculation;
    check(s.status == supernova::Status::ok, "original finite point");
    std::vector<long double> r32(n), r64(n), canonical(n);
    long double mu_error = 0;
    for (std::size_t i = 0; i < n; ++i) {
      auto z = in.zhd[indices[i]];
      auto pref = 1 + (long double)in.zhel[indices[i]];
      auto mu32 = 5 * std::log10(pref * radial(z, p.omega_m, p.w0, p.wa, 32)),
           mu64 = 5 * std::log10(pref * radial(z, p.omega_m, p.w0, p.wa, 64));
      r32[i] = (long double)in.values[indices[i]] - mu32;
      r64[i] = (long double)in.values[indices[i]] - mu64;
      canonical[i] = s.base_residuals[i];
      mu_error = std::max(mu_error,
                          std::abs((long double)s.shape_magnitudes[i] - mu64));
    }
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
    long double ex =
             consumer.profile_operator().cached_response_forward_sensitivity() *
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
    auto refinement = std::abs(
        (long double)batches[1].slots[m].calculation.relative_profile_score -
        s.relative_profile_score);
    auto score_error =
        std::abs((long double)s.relative_profile_score + b.q / 2);
    bool pass = reference_refinement <= 2e-7L && factor_error <= 3e-7L &&
                allowance <= 3e-7L && background_error <= 5e-7L &&
                refinement <= 5e-7L && mu_error <= 1e-8L &&
                score_error <= 1e-6L;
    if (m == 0)
      near(s.relative_profile_score,
           test_reference::w01_historical[1].stable_score, 1e-6,
           "original Lambda frozen score");
    if (m == 1)
      near(s.relative_profile_score,
           test_reference::w01_historical[4].stable_score, 1e-6,
           "original constantw frozen score");
    std::printf(
        "{\"point\":%u,\"omega_m\":%.17g,\"w0\":%.17g,\"wa\":%.17g,\"relative_"
        "score\":%.17g,\"reference_refinement\":%.21Lg,\"factor_error\":%.21Lg,"
        "\"estimated_factor_allowance\":%.21Lg,\"background_error\":%.21Lg,"
        "\"production_refinement\":%.21Lg,\"magnitude_error\":%.21Lg,\"total_"
        "score_error\":%.21Lg,\"passed\":%s}\n",
        m, p.omega_m, p.w0, p.wa, s.relative_profile_score,
        reference_refinement, factor_error, allowance, background_error,
        refinement, mu_error, score_error, pass ? "true" : "false");
    check(pass, "original fixedpoint budgets");
  }
  std::printf("{\"suite\":\"supernova_cpl_original_fixedpoints\",\"n\":%zu,"
              "\"reference_condition_estimate\":%.21Lg,\"passed\":true}\n",
              n, Cnorm * inverse_norm);
}

void ordinary() {
  static_assert(!std::is_default_constructible_v<supernova::ModelPointV2>);
  auto source = synthetic();
  std::array<cosmology::Query, 3> queries{
      {{1, .9, cosmology::Convention::released_zhd_zhel},
       {2, 2.1, cosmology::Convention::released_zhd_zhel},
       {.01, .01, cosmology::Convention::released_zhd_zhel}}};
  supernova::Policy policy;
  policy.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  policy.maximum_forward_sensitivity = 1e-10;
  policy.maximum_matrix_elements = 9;
  auto consumer = supernova::prepare_synthetic(
      source, source.source().measurement_ids, queries, policy);
  check(consumer.status() == supernova::Status::ok, "synthetic preparation");
  check(consumer.selected_source_indices().size() == 2,
        "exact threshold exclusion");
  std::array<supernova::ModelPointV2, 4> points{
      {{cosmology::Model::flat_lcdm_late_v1, 0, 0, -1, 0},
       {cosmology::Model::constant_q_flat_v1, 0, -1, -1, 0},
       {cosmology::Model::flat_cpl_late_v1, 0, 0, -1, 0},
       {cosmology::Model::flat_cpl_late_v1, .3, 0, -.9, .4}}};
  auto result = consumer.evaluate_batch_v2(points, policy);
  check(result.status == supernova::Status::ok && result.slots.size() == 4,
        "mixed batch");
  for (unsigned i = 0; i < result.slots.size(); ++i) {
    const auto &s = result.slots[i];
    check(s.calculation.status == supernova::Status::ok, "finite mixed row");
    check(s.source.w0 == points[i].w0 && s.source.wa == points[i].wa,
          "owned explicit source");
    const auto &r = s.calculation.base_residuals;
    near(s.calculation.quadratic,
         (static_cast<long double>(r[0]) - r[1]) *
             (static_cast<long double>(r[0]) - r[1]) / 11,
         1e-10, "adjugate score");
  }
  for (unsigned i = 0; i < 3; ++i) {
    near(result.slots[i].calculation.quadratic, 25.L / 11, 1e-10,
         "Lambda analytic score");
    near(result.slots[i].calculation.offset_coefficient, 10.L + 7.L / 11, 1e-10,
         "Lambda offset");
  }
  for (double om : {0., .3, 1.}) {
    std::array<supernova::ModelPoint, 1> old{
        {{cosmology::Model::flat_lcdm_late_v1, om, 0}}};
    std::array<supernova::ModelPointV2, 1> lambda{
        {{cosmology::Model::flat_cpl_late_v1, om, 0, -1, 0}}};
    auto a = consumer.evaluate_batch(old, policy).slots[0];
    auto b = consumer.evaluate_batch_v2(lambda, policy).slots[0].calculation;
    near(a.relative_profile_score, b.relative_profile_score, 1e-10,
         "Lambda limit legacy");
    check(a.shape_magnitudes == b.shape_magnitudes, "Lambda shape bits");
  }
  for (double w : {-1., -2. / 3, -1. / 3, 0.}) {
    std::array<supernova::ModelPointV2, 1> p{
        {{cosmology::Model::flat_cpl_late_v1, 0, 0, w, 0}}};
    auto s = consumer.evaluate_batch_v2(p, policy).slots[0].calculation;
    check(s.status == supernova::Status::ok, "constantw");
    auto power = 1.5L * (1 + (long double)w), a = 1 - power;
    for (unsigned i = 0; i < 2; ++i) {
      auto L = std::log1p((long double)queries[i].z_expansion);
      auto I = a == 0 ? L : std::expm1(a * L) / a;
      near(s.shape_magnitudes[i],
           5 * std::log10((1 + queries[i].z_observer) * I), 1e-10,
           "constantw exact radial");
    }
  }
  for (auto p : {supernova::ModelPointV2{cosmology::Model::flat_cpl_late_v1, .3,
                                         0, -.9, .4},
                 supernova::ModelPointV2{cosmology::Model::flat_cpl_late_v1, .3,
                                         0, -1.1, -.4}}) {
    auto s = consumer.evaluate_batch_v2(std::span(&p, 1), policy)
                 .slots[0]
                 .calculation;
    for (unsigned i = 0; i < 2; ++i) {
      auto r32 = radial(queries[i].z_expansion, p.omega_m, p.w0, p.wa, 32),
           r64 = radial(queries[i].z_expansion, p.omega_m, p.w0, p.wa, 64);
      near(r32, r64, 2e-13, "GL refinement");
      near(s.shape_magnitudes[i],
           5 * std::log10((1 + queries[i].z_observer) * r64), 1e-10,
           "CPL independent radial");
    }
  }
  std::array<supernova::ModelPoint, 1> forbidden{
      {{cosmology::Model::flat_cpl_late_v1, .3, 0}}};
  check(consumer.evaluate_batch(forbidden, policy).slots[0].status !=
            supernova::Status::ok,
        "v1 cannot admit CPL");
  for (auto p : {supernova::ModelPointV2{cosmology::Model::flat_cpl_late_v1, .3,
                                         1, -1, 0},
                 {cosmology::Model::flat_lcdm_late_v1, .3, 0, -.9, 0},
                 {cosmology::Model::constant_q_flat_v1, .1, -1, -1, 0},
                 {cosmology::Model::flat_cpl_late_v1, .3, 0, -2.1, 0},
                 {static_cast<cosmology::Model>(99), .3, 0, -1, 0},
                 {cosmology::Model::flat_cpl_late_v1, .3, 0, -1,
                  std::numeric_limits<double>::infinity()}}) {
    auto s = consumer.evaluate_batch_v2(std::span(&p, 1), policy).slots[0];
    check(s.calculation.status != supernova::Status::ok,
          "invalid explicit model");
    check(s.source.model == p.model && s.source.w0 == p.w0 &&
              s.source.wa == p.wa,
          "attempt retained");
    check(s.calculation.shape_magnitudes.empty() &&
              s.calculation.profiled_residuals.empty(),
          "failure no vectors");
  }
  auto cap = policy;
  cap.maximum_models = 0;
  check(consumer.evaluate_batch_v2(points, cap).slots.empty(),
        "model cap early");
  auto work = policy;
  work.background.maximum_total_evaluations =
      result.slots[0].calculation.background_evaluations + 2;
  auto limited = consumer.evaluate_batch_v2(points, work);
  check(limited.slots[0].calculation.status == supernova::Status::ok &&
            limited.slots[1].calculation.status ==
                supernova::Status::work_limit,
        "mixed callback budget not reset");
  auto tight = policy;
  tight.maximum_forward_sensitivity = 1e-30;
  auto fail = consumer.evaluate_batch_v2(points, tight).slots[0].calculation;
  check(fail.status == supernova::Status::numerical_failure &&
            fail.numerical_status ==
                numerics::Status::conditioning_budget_exceeded,
        "wide cache conditioning cause");
  check(consumer.observations().source().uncertainty_matrix ==
            source.source().uncertainty_matrix,
        "raw asymmetry retained");
  std::printf("{\"suite\":\"supernova_cpl_owner\",\"checks\":%u,\"max_absolute_"
              "error\":%.17Lg,\"passed\":true}\n",
              checks, max_error);
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc == 1)
      ordinary();
    else if (argc == 4 && std::string(argv[3]) == "--verified-original-assets")
      actual(argv[1], argv[2]);
    else
      throw std::runtime_error(
          "usage: test_supernova_cpl [COV TABLE --verified-original-assets]; "
          "Rust SHA256 guard required");
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s checks=%u\n", e.what(), checks);
    return 1;
  }
}
