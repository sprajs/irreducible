// Owned deterministic magnitude-hypothesis tests. Exact 2x2 adjugate profile
// and constant-q=-1 analytic I=z; shared libm is disclosed. Score budget1e-10.
// Current Pantheon adaptation is not the historical 1820-object J06 replay.
#include "fixtures/ldlt_reference.hpp"
#include "fixtures/w01_grey_magnitude.hpp"
#include "fixtures/w01_historical.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
using namespace irred;
namespace {
int checks = 0;
long double worst = 0;
void check(bool v, const char *w) {
  ++checks;
  if (!v)
    throw std::runtime_error(w);
}
void near(long double a, long double b) {
  worst = std::max(worst, std::abs(a - b));
  check(std::abs(a - b) <= 1e-10, "rational profile");
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
long double reference_radial(long double z, const supernova::ModelPointV2 &p,
                             unsigned panels) {
  if (p.model == cosmology::Model::constant_q_flat_v1) {
    if (p.constant_q == 0)
      return std::log1p(z);
    return (1 - std::pow(1 + z, -(long double)p.constant_q)) / p.constant_q;
  }
  return radial(z, p.omega_m, p.w0, p.wa, panels);
}
// Geometry is supplied by native C++; this mode compares the older class's
// fixed-amplitude projection only, not an independent background provider.
void direct_native(const char *cov, const char *table) {
  auto p = supernova::Policy{};
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.maximum_forward_sensitivity = 1e-10;
  p.background.maximum_total_evaluations = 20000000;
  p.background.integration = {1e-12, 1e-12, 100000, 30};
  auto c = supernova::prepare(original(cov, table), p);
  check(c.status() == supernova::Status::ok, "direct preparation");
  std::array<supernova::ModelPointV2, 4> backgrounds{
      {{cosmology::Model::flat_lcdm_late_v1, .3, 0, -1, 0},
       {cosmology::Model::constant_q_flat_v1, 0, 0, -1, 0},
       {cosmology::Model::constant_q_flat_v1, 0, .5, -1, 0},
       {cosmology::Model::flat_cpl_late_v1, .35, 0, -.9, .2}}};
  std::vector<supernova::GreyMagnitudePoint> pts;
  for (auto point : backgrounds)
    for (double e : {-.2, 0., .2})
      pts.emplace_back(point, e);
  auto out = c.evaluate_grey_magnitude_batch(pts, {p, 256 * 1024 * 1024});
  check(out.slots.size() == 12, "direct rows");
  std::printf(
      "{\"suite\":\"grey_native_geometry_for_original_class\",\"asset_identity_"
      "provenance\":\"external Rust SHA guard; C++ does not hash\",\"rows\":[");
  for (size_t j = 0; j < out.slots.size(); ++j) {
    auto &r = out.slots[j];
    check(r.calculation.status == supernova::Status::ok, "direct finite");
    const auto &reference = irred::test_reference::grey_magnitude_reference[j];
    check(reference.model == (unsigned)r.source.background.model &&
              reference.omega_m == r.source.background.omega_m &&
              reference.constant_q == r.source.background.constant_q &&
              reference.w0 == r.source.background.w0 &&
              reference.wa == r.source.background.wa &&
              reference.epsilon_mag == r.source.epsilon_mag,
          "direct frozen identity");
    check(std::abs((long double)r.calculation.relative_profile_score -
                   reference.relative_score) <= 2e-7L,
          "native migrated original-class fixed-amplitude fact");
    check(std::abs((long double)r.calculation.relative_profile_score -
                   reference.original_basis_score) <= 2e-7L,
          "native migrated original-class zHD basis fact");
    if (j)
      std::printf(",");
    std::printf("{\"point\":%zu,\"epsilon_mag\":%.17g,\"native_score\":%.17g,"
                "\"geometry\":[",
                j, r.source.epsilon_mag, r.calculation.relative_profile_score);
    for (size_t k = 0; k < r.calculation.shape_magnitudes.size(); ++k)
      std::printf("%s%.17g", k ? "," : "", r.calculation.shape_magnitudes[k]);
    std::printf("],\"B\":[");
    for (size_t k = 0; k < r.magnitude_shifts.size(); ++k)
      std::printf("%s%.17g", k ? "," : "", r.magnitude_shifts[k]);
    std::printf("]}");
  }
  std::printf("]}\n");
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
  std::array<supernova::ModelPointV2, 4> backgrounds{
      {{cosmology::Model::flat_lcdm_late_v1, .3, 0, -1, 0},
       {cosmology::Model::constant_q_flat_v1, 0, 0, -1, 0},
       {cosmology::Model::constant_q_flat_v1, 0, .5, -1, 0},
       {cosmology::Model::flat_cpl_late_v1, .35, 0, -.9, .2}}};
  std::vector<supernova::GreyMagnitudePoint> grid;
  for (auto point : backgrounds)
    for (double epsilon : {-.2, 0., .2})
      grid.emplace_back(point, epsilon);
  supernova::GreyMagnitudePolicy grey_policy{policy, 256 * 1024 * 1024};
  std::array<supernova::GreyMagnitudeBatchResult, 3> batches;
  for (unsigned level = 0; level < 3; ++level) {
    auto eps = level == 0 ? 1e-8 : level == 1 ? 1e-10 : 1e-12;
    grey_policy.evaluation.background.integration = {eps, eps, 100000, 30};
    batches[level] = consumer.evaluate_grey_magnitude_batch(grid, grey_policy);
    check(batches[level].slots.size() == 12, "original model batch");
  }
  long double gamma = (2 * n + 2) * std::numeric_limits<long double>::epsilon();
  gamma /= 1 - gamma;
  long double eps64 = std::numeric_limits<double>::epsilon();
  bool all_pass = true;
  for (unsigned m = 0; m < 12; ++m) {
    const auto &p = grid[m].background;
    auto epsilon = grid[m].epsilon_mag;
    const auto &s = batches[2].slots[m].calculation;
    check(s.status == supernova::Status::ok, "original finite point");
    std::vector<long double> r32(n), r64(n), canonical(n);
    long double mu_error = 0;
    for (std::size_t i = 0; i < n; ++i) {
      auto z = in.zhd[indices[i]];
      auto pref = 1 + (long double)in.zhel[indices[i]];
      auto mu32 = 5 * std::log10(pref * reference_radial(z, p, 32)),
           mu64 = 5 * std::log10(pref * reference_radial(z, p, 64));
      const auto shift =
          (long double)epsilon * std::log1p((long double)z) / std::log(2.L);
      mu32 += shift;
      mu64 += shift;
      r32[i] = (long double)in.values[indices[i]] - mu32;
      r64[i] = (long double)in.values[indices[i]] - mu64;
      canonical[i] = s.base_residuals[i];
      mu_error = std::max(
          mu_error, std::abs((long double)s.shape_magnitudes[i] +
                             batches[2].slots[m].magnitude_shifts[i] - mu64));
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
    auto refinement = std::abs(
        (long double)batches[1].slots[m].calculation.relative_profile_score -
        s.relative_profile_score);
    auto score_error =
        std::abs((long double)s.relative_profile_score + b.q / 2);
    bool pass = reference_refinement <= 2e-7L && factor_error <= 3e-7L &&
                allowance <= 3e-7L && background_error <= 5e-7L &&
                refinement <= 5e-7L && mu_error <= 1e-8L &&
                score_error <= 1e-6L;
    if (epsilon == 0) {
      std::array<supernova::ModelPointV2, 1> baseline{p};
      auto old = consumer.evaluate_batch_v2(baseline, grey_policy.evaluation);
      check(s.relative_profile_score ==
                    old.slots[0].calculation.relative_profile_score &&
                s.base_residuals == old.slots[0].calculation.base_residuals,
            "original zero exact");
    }
    std::printf(
        "{\"point\":%u,\"model\":%u,\"constant_q\":%.17g,\"epsilon_mag\":%.17g,"
        "\"omega_m\":%.17g,\"w0\":%.17g,"
        "\"wa\":%.17g,\"relative_"
        "score\":%.17g,\"reference_refinement\":%.21Lg,\"factor_error\":%.21Lg,"
        "\"estimated_factor_allowance\":%.21Lg,\"background_error\":%.21Lg,"
        "\"production_refinement\":%.21Lg,\"magnitude_error\":%.21Lg,\"total_"
        "score_error\":%.21Lg,\"passed\":%s}\n",
        m, (unsigned)p.model, p.constant_q, epsilon, p.omega_m, p.w0, p.wa,
        s.relative_profile_score, reference_refinement, factor_error, allowance,
        background_error, refinement, mu_error, score_error,
        pass ? "true" : "false");
    const auto &frozen = irred::test_reference::grey_magnitude_reference[m];
    check(frozen.model == (unsigned)p.model && frozen.omega_m == p.omega_m &&
              frozen.constant_q == p.constant_q && frozen.w0 == p.w0 &&
              frozen.wa == p.wa && frozen.epsilon_mag == epsilon,
          "frozen grey fixedpoint identity");
    check(std::abs((long double)s.relative_profile_score -
                   frozen.relative_score) <= 2e-7L,
          "extracted original-class frozen relative score");
    check(std::abs((long double)s.relative_profile_score -
                   frozen.original_basis_score) <= 2e-7L,
          "extracted original-class frozen brightness basis");
    all_pass &= pass;
  }
  std::printf("{\"suite\":\"grey_magnitude_current_Pantheon_conditional_"
              "fixedpoints\",\"n\":%zu,"
              "\"reference_condition_estimate\":%.21Lg,\"passed\":true}\n",
              n, Cnorm * inverse_norm);
  check(all_pass, "all current-Pantheon conditional fixedpoint budgets");
}
supernova::Policy policy() {
  supernova::Policy p;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.maximum_forward_sensitivity = 1e-10;
  p.background.maximum_total_evaluations = 10000;
  return p;
}
std::array<cosmology::Query, 3> queries(bool same = false) {
  return {{{.5, .75, cosmology::Convention::released_zhd_zhel},
           {same ? .5 : 1., 1.2, cosmology::Convention::released_zhd_zhel},
           {.005, .005, cosmology::Convention::geometric_same_redshift}}};
}
supernova::Consumer consumer(bool same = false, bool third = false) {
  auto qs = queries(same);
  if (third)
    qs[2] = {.7, .9, cosmology::Convention::released_zhd_zhel};
  observations::Input in{};
  in.profile = observations::Profile::gaussian_fixture_v1;
  in.role = observations::Role::synthetic_control;
  in.unit = observations::Unit::magnitude;
  in.calibration = observations::Calibration::not_applicable;
  in.uncertainty = observations::Uncertainty::covariance;
  in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  in.component = observations::Component::total;
  in.ordering_provenance = "exact generated row order";
  in.table_sha256 = std::string(64, 'a');
  in.uncertainty_sha256 = std::string(64, 'b');
  in.measurement_ids = {"A", "B", "excluded"};
  in.event_ids = in.measurement_ids;
  in.uncertainty_axis_ids = in.measurement_ids;
  for (auto q : qs)
    in.values.push_back(5 * std::log10((1 + q.z_observer) * q.z_expansion) +
                        10 + .08 * std::log1p(q.z_expansion) / std::log(2.));
  in.uncertainty_matrix =
      third ? std::vector<double>{4, 1, 0, 1, 9, 0, 0, 0, 16}
            : std::vector<double>{4, 1, 8, 1, 9, 5, 7, 5, -100};
  in.missing = {0, 0, 0};
  in.source_selection = {1, 1, 1};
  in.quality = {0, 0, 0};
  auto prepared = observations::prepare(std::move(in), {3, 9, 4096});
  std::array<std::string, 3> ids{"A", "B", "excluded"};
  return supernova::prepare_synthetic(std::move(prepared), ids, qs, policy());
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc == 4 && std::string(argv[3]) == "--verified-original-assets") {
      if (std::getenv("IRRED_GREY_DIRECT_NATIVE"))
        direct_native(argv[1], argv[2]);
      else
        actual(argv[1], argv[2]);
      return 0;
    }
    if (argc != 1)
      throw std::runtime_error(
          "expected ordinary or SHA-guarded original mode");
    auto c = consumer();
    check(c.status() == supernova::Status::ok, "prepare");
    supernova::GreyMagnitudePolicy p{policy(), 1048576};
    supernova::ModelPointV2 model{cosmology::Model::constant_q_flat_v1, 0, -1,
                                  -1, 0};
    std::array<supernova::ModelPointV2, 1> base{model};
    auto old = c.evaluate_batch_v2(base, p.evaluation);
    for (double e : {-.5, -.2, -0., 0., .08, .2, .5}) {
      std::array<supernova::GreyMagnitudePoint, 1> pts{{{model, e}}};
      auto r = c.evaluate_grey_magnitude_batch(pts, p);
      check(r.status == supernova::Status::ok && r.slots.size() == 1, "batch");
      auto &s = r.slots[0];
      check(s.calculation.status == supernova::Status::ok, "finite");
      check(std::bit_cast<std::uint64_t>(s.source.epsilon_mag) ==
                std::bit_cast<std::uint64_t>(e),
            "epsilon source bits");
      auto f0 = std::log1p(.5L) / std::log(2.L), f1 = 1.L;
      near(s.calculation.quadratic,
           (.08L - e) * (.08L - e) * (f0 - f1) * (f0 - f1) / 11);
      near(s.magnitude_shifts[0], e * f0);
      near(s.magnitude_shifts[1], e);
      check(s.calculation.shape_magnitudes ==
                old.slots[0].calculation.shape_magnitudes,
            "geometry unchanged");
      if (e == 0) {
        check(s.calculation.relative_profile_score ==
                  old.slots[0].calculation.relative_profile_score,
              "zero score exact");
        check(s.calculation.base_residuals ==
                  old.slots[0].calculation.base_residuals,
              "zero residual exact");
        check(s.calculation.solve_diagnostics.adjusted_residuals ==
                  old.slots[0].calculation.solve_diagnostics.adjusted_residuals,
              "zero canonical residual exact");
        check(s.calculation.background_evaluations ==
                  old.slots[0].calculation.background_evaluations,
              "zero work exact");
      }
      check(s.magnitude_shifts[0] != e * std::log1p(.75) / std::log(2.) ||
                e == 0,
            "zHD not zHEL");
    }
    for (double e : {std::nextafter(.5, 1.), std::nextafter(-.5, -1.),
                     std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::infinity()}) {
      std::array<supernova::GreyMagnitudePoint, 1> pts{{{model, e}}};
      auto r = c.evaluate_grey_magnitude_batch(pts, p);
      auto &s = r.slots[0];
      check(s.calculation.status != supernova::Status::ok, "domain reject");
      check(s.magnitude_shifts.empty() &&
                s.calculation.shape_magnitudes.empty() &&
                s.calculation.base_residuals.empty(),
            "failed arrays absent");
    }
    auto constant = consumer(true);
    for (double e : {-.5, .08, .5}) {
      std::array<supernova::GreyMagnitudePoint, 1> pts{{{model, e}}};
      auto r = constant.evaluate_grey_magnitude_batch(pts, p);
      near(r.slots[0].calculation.quadratic, 0);
    }
    std::array<supernova::GreyMagnitudePoint, 1> pts{{{model, .2}}};
    auto no = p;
    no.maximum_native_bytes = 0;
    auto empty = c.evaluate_grey_magnitude_batch(pts, no);
    check(empty.status == supernova::Status::work_limit && empty.slots.empty(),
          "required bytecap");
    no = p;
    no.evaluation.maximum_models = 0;
    empty = c.evaluate_grey_magnitude_batch(pts, no);
    check(empty.status == supernova::Status::work_limit && empty.slots.empty(),
          "model cap");
    empty = c.evaluate_grey_magnitude_batch(
        std::span<const supernova::GreyMagnitudePoint>(
            static_cast<const supernova::GreyMagnitudePoint *>(nullptr),
            SIZE_MAX),
        p);
    check(empty.status == supernova::Status::work_limit && empty.slots.empty(),
          "impossible count before dereference");
    no = p;
    no.evaluation.background.maximum_total_evaluations = 0;
    auto r = c.evaluate_grey_magnitude_batch(pts, no);
    check(r.slots[0].calculation.status == supernova::Status::work_limit &&
              r.slots[0].magnitude_shifts.empty(),
          "work failure absent B");
    no = p;
    no.evaluation.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    empty = c.evaluate_grey_magnitude_batch(pts, no);
    check(empty.status == supernova::Status::incompatible_metadata &&
              empty.slots.empty(),
          "precision mismatch");
    auto three = consumer(false, true);
    std::vector<supernova::GreyMagnitudePoint> many;
    for (int j = 0; j < 7; ++j)
      many.emplace_back(model, .2);
    auto budget = p;
    constexpr size_t n = 3;
    budget.maximum_native_bytes =
        many.size() * sizeof(supernova::GreyMagnitudeSlot) +
        5 * many.size() * n * sizeof(double) +
        n * (sizeof(cosmology::Slot) + 16 * sizeof(double) +
             8 * sizeof(long double));
    auto bounded = three.evaluate_grey_magnitude_batch(many, budget);
    check(bounded.slots.size() == 7, "nonpower multirow exact admission");
    for (auto &row : bounded.slots) {
      check(row.calculation.status == supernova::Status::ok, "nonpower finite");
      check(row.calculation.shape_magnitudes.capacity() == n &&
                row.calculation.base_residuals.capacity() == n &&
                row.magnitude_shifts.capacity() == n,
            "requested vector payload no geometric growth");
    }
    --budget.maximum_native_bytes;
    auto rejected = three.evaluate_grey_magnitude_batch(many, budget);
    check(rejected.status == supernova::Status::work_limit &&
              rejected.slots.empty(),
          "byte boundary");
    for (double e : {std::numeric_limits<double>::denorm_min(),
                     std::numeric_limits<double>::min()}) {
      std::array<supernova::GreyMagnitudePoint, 1> tiny{{{model, e}}};
      auto result = c.evaluate_grey_magnitude_batch(tiny, p);
      check(result.slots[0].calculation.status ==
                    supernova::Status::numerical_failure &&
                result.slots[0].calculation.numerical_status ==
                    numerics::Status::outside_domain,
            "nonzero B underflow truthful");
      check(result.slots[0].magnitude_shifts.empty() &&
                result.slots[0].calculation.base_residuals.empty(),
            "range failed arrays absent");
    }
    std::printf("grey magnitude owner PASS %d maxerror %.18Lg\n", checks,
                worst);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "grey magnitude FAIL %s\n", e.what());
    return 1;
  }
}
