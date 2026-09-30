// Owned deterministic magnitude-hypothesis tests. Exact 2x2 adjugate profile
// and constant-q=-1 analytic I=z; shared libm is disclosed. Score budget1e-10.
// Current Pantheon adaptation is not the historical 1820-object J06 replay.
#include "fixtures/ldlt_reference.hpp"
#include "fixtures/w01_cpl.hpp"
#include "fixtures/w01_grey_magnitude.hpp"
#include "fixtures/w01_historical.hpp"
#include "fixtures/w01_piecewise.hpp"
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
void check(bool v, const char *w) {
  ++checks;
  if (!v)
    throw std::runtime_error(w);
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
// Independent split scale-factor GL8; multiplicative E, same canonical q/edges.
long double piecewise_radial(long double z, const std::array<double, 5> &q,
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
long double reference_radial(long double z, const cosmology::ExpansionSpec &p,
                             unsigned panels) {
  return std::visit(
      [&](const auto &m) -> long double {
        using T = std::decay_t<decltype(m)>;
        if constexpr (std::is_same_v<T, cosmology::ConstantQ>) {
          return m.q == 0 ? std::log1p(z)
                          : (1 - std::pow(1 + z, -(long double)m.q)) / m.q;
        } else if constexpr (std::is_same_v<T, cosmology::LCDM>)
          return radial(z, m.omega_m, -1, 0, panels);
        else if constexpr (std::is_same_v<T, cosmology::CPL>)
          return radial(z, m.omega_m, m.w0, m.wa, panels);
        else
          return piecewise_radial(z, m.q, panels);
      },
      p);
}
void actual(const char *cov, const char *table) {
  auto source = original(cov, table);
  check(source.status() == observations::Status::ok, "original typed source");
  supernova::PreparationPolicy preparation;
  preparation.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  preparation.maximum_selected_rows = 1701;
  preparation.maximum_matrix_elements = 1701 * 1701;
  preparation.maximum_forward_sensitivity = 1e-10;
  preparation.maximum_native_bytes = 512 * 1024 * 1024;
  auto shared =
      std::make_shared<const observations::Prepared>(std::move(source));
  auto consumer =
      supernova::prepare(supernova::pantheon_zhd_gt_001(shared), preparation);
  supernova::Policy policy;
  policy.arithmetic = preparation.arithmetic;
  policy.maximum_forward_sensitivity = 1e-10;
  policy.maximum_models = 12;
  policy.maximum_native_bytes = 256 * 1024 * 1024;
  policy.requested = 63;
  policy.background.maximum_queries = 1701;
  policy.background.maximum_callbacks = 20000000;
  policy.background.maximum_segment_visits = 50000;
  policy.background.maximum_native_bytes = 32 * 1024 * 1024;
  check(consumer.status() == supernova::Status::ok,
        "original wide preparation");
  auto indices = consumer.selected_source().source_indices;
  check(indices.size() == 1590, "original selection");
  auto n = indices.size();
  const auto &in = consumer.selected_source().source->source();
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
  std::array<cosmology::ExpansionSpec, 4> backgrounds{
      cosmology::LCDM{.3}, cosmology::ConstantQ{0}, cosmology::ConstantQ{.5},
      cosmology::CPL{.35, -.9, .2}};
  std::vector<supernova::ModelPoint> grid;
  for (auto point : backgrounds)
    for (double epsilon : {-.2, 0., .2})
      grid.emplace_back(point, cosmology::FlatFLRW{},
                        supernova::GreyLog1pMagnitude{epsilon});
  const bool all_named = std::getenv("IRRED_SN_ALL_NAMED") != nullptr;
  std::vector<size_t> batch_edges{0, grid.size()};
  if (all_named) {
    for (const auto &f : test_reference::w01_historical) {
      cosmology::ExpansionSpec e =
          f.lcdm ? cosmology::ExpansionSpec(cosmology::LCDM(f.parameter))
                 : cosmology::ExpansionSpec(cosmology::ConstantQ(f.parameter));
      grid.emplace_back(e, cosmology::FlatFLRW{},
                        supernova::NoMagnitudeEffect{});
    }
    batch_edges.push_back(grid.size());
    for (const auto &f : test_reference::w01_cpl_points)
      grid.emplace_back(cosmology::CPL(f.omega_m, f.w0, f.wa),
                        cosmology::FlatFLRW{}, supernova::NoMagnitudeEffect{});
    batch_edges.push_back(grid.size());
    for (const auto &q : test_reference::w01_piecewise_q)
      grid.emplace_back(cosmology::FixedFiveBinQ(q), cosmology::FlatFLRW{},
                        supernova::NoMagnitudeEffect{});
    batch_edges.push_back(grid.size());
  }
  std::array<supernova::BatchResult, 3> batches;
  for (unsigned level = 0; level < 3; ++level) {
    auto eps = level == 0 ? 1e-8 : level == 1 ? 1e-10 : 1e-12;
    policy.background.integration =
        numerics::IntegrationPolicy{eps, eps, 100000, 30};
    // Named comparison batches retain their frozen separate work policies;
    // one immutable source/factor/reference is reused across all of them.
    for (size_t k = 1; k < batch_edges.size(); ++k) {
      const auto begin = batch_edges[k - 1], count = batch_edges[k] - begin;
      auto batch = consumer.evaluate_batch(
          std::span<const supernova::ModelPoint>(grid.data() + begin, count),
          policy);
      check(batch.slots.size() == count, "original named model batch");
      const char *lineage = k == 1   ? "grey12_current1590_conditional_not1820"
                            : k == 2 ? "historical7"
                            : k == 3 ? "CPL4"
                                     : "fixedq6";
      std::printf("{\"batch\":\"%s\",\"level\":%u,\"models\":%zu,"
                  "\"callback_cap\":20000000,\"segment_cap\":50000,"
                  "\"callbacks\":%zu,\"segments\":%zu,"
                  "\"reference_allocation\":2e-7,\"factor_allocation\":3e-7,"
                  "\"background_allocation\":5e-7,\"score_allocation\":1e-6}\n",
                  lineage, level, count, batch.work.callbacks,
                  batch.work.segment_visits);
      for (auto &slot : batch.slots)
        batches[level].slots.push_back(std::move(slot));
    }
  }
  long double gamma = (2 * n + 2) * std::numeric_limits<long double>::epsilon();
  gamma /= 1 - gamma;
  long double eps64 = std::numeric_limits<double>::epsilon();
  bool all_pass = true;
  for (unsigned m = 0; m < grid.size(); ++m) {
    const auto &p = grid[m].expansion;
    const auto *grey =
        std::get_if<supernova::GreyLog1pMagnitude>(&grid[m].source_effect);
    const double epsilon = grey ? grey->epsilon_mag : 0;
    const bool analytic = std::holds_alternative<cosmology::FixedFiveBinQ>(p);
    const unsigned coarse_panels = analytic ? 64 : 32;
    const unsigned fine_panels = analytic ? 128 : 64;
    const auto &s = batches[2].slots[m];
    check(s.status == supernova::Status::ok, "original finite point");
    std::vector<long double> r32(n), r64(n), canonical(n);
    long double mu_error = 0;
    for (std::size_t i = 0; i < n; ++i) {
      auto z = in.zhd[indices[i]];
      auto pref = 1 + (long double)in.zhel[indices[i]];
      auto mu32 = 5 * std::log10(pref * reference_radial(z, p, coarse_panels)),
           mu64 = 5 * std::log10(pref * reference_radial(z, p, fine_panels));
      const auto shift =
          (long double)epsilon * std::log1p((long double)z) / std::log(2.L);
      mu32 += shift;
      mu64 += shift;
      r32[i] = (long double)in.values[indices[i]] - mu32;
      r64[i] = (long double)in.values[indices[i]] - mu64;
      canonical[i] = s.corrected_residuals[i];
      mu_error = std::max(
          mu_error, std::abs((long double)s.geometric_shape[i] +
                             batches[2].slots[m].magnitude_shifts[i] - mu64));
    }
    auto a = score(ref, r32, gram), b = score(ref, r64, gram),
         same = score(ref, canonical, gram);
    auto reference_refinement = std::abs(a.q - b.q) / 2,
         factor_error = std::abs((long double)s.score->quadratic - same.q) / 2,
         background_error = std::abs(same.q - b.q) / 2;
    const auto &d = *s.diagnostics;
    auto wx = consumer.profile_operator().cached_response_solution();
    long double wxnorm = 0, h = 0, hdot = 0, du1 = 0, duinf = 0;
    for (std::size_t i = 0; i < n; ++i) {
      wxnorm = std::max(wxnorm, std::abs((long double)wx[i]));
      h += (long double)wx[i] * s.profiled_residuals[i];
      hdot += std::abs((long double)wx[i] * s.profiled_residuals[i]);
      auto du = std::abs((long double)s.profiled_residuals[i] -
                         ((long double)s.corrected_residuals[i] -
                          s.score->offset_coefficient));
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
    auto da = (db + std::abs((long double)s.score->offset_coefficient) * dg) /
                  (G - dg) +
              eps64 / (1 - eps64) *
                  std::abs((long double)s.score->offset_coefficient);
    auto dh = d.adjusted_residual_l1 * ex + gamma * hdot;
    auto allowance =
        .5L * d.adjusted_residual_l1 * eu +
        gamma * d.adjusted_residual_l1 * d.adjusted_solution_norm_inf +
        (std::abs(h) + dh) * da + .5L * (G + dg) * da * da +
        du1 * d.adjusted_solution_norm_inf + .5L * du1 * inverse_norm * duinf +
        .5L * eps64 / (1 - eps64) * std::abs((long double)s.score->quadratic);
    auto refinement = std::abs(
        (long double)batches[1].slots[m].score->relative_profile_score -
        s.score->relative_profile_score);
    auto score_error =
        std::abs((long double)s.score->relative_profile_score + b.q / 2);
    bool pass = reference_refinement <= 2e-7L && factor_error <= 3e-7L &&
                allowance <= 3e-7L && background_error <= 5e-7L &&
                refinement <= 5e-7L && mu_error <= 1e-8L &&
                score_error <= 1e-6L;
    if (epsilon == 0) {
      std::array<supernova::ModelPoint, 1> baseline{{supernova::ModelPoint{
          p, cosmology::FlatFLRW{}, supernova::NoMagnitudeEffect{}}}};
      auto old = consumer.evaluate_batch(baseline, policy);
      check(s.score->relative_profile_score ==
                    old.slots[0].score->relative_profile_score &&
                s.corrected_residuals == old.slots[0].corrected_residuals,
            "original zero exact");
    }
    if (m < 12) {
      const auto &frozen = irred::test_reference::grey_magnitude_reference[m];
      std::printf(
          "{\"point\":%u,\"model\":%u,\"constant_q\":%.17g,\"epsilon_mag\":%."
          "17g,"
          "\"omega_m\":%.17g,\"w0\":%.17g,"
          "\"wa\":%.17g,\"relative_"
          "score\":%.17g,\"reference_refinement\":%.21Lg,\"factor_error\":%."
          "21Lg,"
          "\"estimated_factor_allowance\":%.21Lg,\"background_error\":%.21Lg,"
          "\"production_refinement\":%.21Lg,\"magnitude_error\":%.21Lg,\"total_"
          "score_error\":%.21Lg,\"passed\":%s}\n",
          m, frozen.model, frozen.constant_q, epsilon, frozen.omega_m,
          frozen.w0, frozen.wa, s.score->relative_profile_score,
          reference_refinement, factor_error, allowance, background_error,
          refinement, mu_error, score_error, pass ? "true" : "false");
      const bool identity = std::visit(
          [&](const auto &v) {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, cosmology::LCDM>)
              return frozen.model == 0 && frozen.omega_m == v.omega_m;
            else if constexpr (std::is_same_v<T, cosmology::ConstantQ>)
              return frozen.model == 1 && frozen.constant_q == v.q;
            else if constexpr (std::is_same_v<T, cosmology::CPL>)
              return frozen.model == 2 && frozen.omega_m == v.omega_m &&
                     frozen.w0 == v.w0 && frozen.wa == v.wa;
            else
              return false;
          },
          p);
      check(identity && frozen.epsilon_mag == epsilon,
            "frozen active model/effect identity");
      check(std::abs((long double)s.score->relative_profile_score -
                     frozen.relative_score) <= 2e-7L,
            "extracted original-class frozen relative score");
      check(std::abs((long double)s.score->relative_profile_score -
                     frozen.original_basis_score) <= 2e-7L,
            "extracted original-class frozen brightness basis");
    } else {
      const char *name = m < 19 ? "historical7" : m < 23 ? "CPL4" : "fixedq6";
      const double frozen_score =
          m < 19   ? test_reference::w01_historical[m - 12].stable_score
          : m < 23 ? test_reference::w01_cpl_points[m - 19].relative_score
                   : test_reference::w01_piecewise_scores[m - 23];
      check(std::abs((long double)s.score->relative_profile_score -
                     frozen_score) <= 2e-7L,
            "frozen named original score");
      if (m < 19)
        check(std::abs((long double)s.score->relative_profile_score -
                       test_reference::w01_historical[m - 12].full_score) <=
                  2e-7L,
              "historical full inverse named score");
      if (m >= 23) {
        check(
            std::abs((long double)s.score->relative_profile_score -
                     test_reference::w01_piecewise_historical_stable[m - 23]) <=
                2e-7L,
            "historical piecewise stable score");
        check(std::abs(
                  (long double)s.score->relative_profile_score -
                  test_reference::w01_piecewise_historical_full_inverse[m -
                                                                        23]) <=
                  2e-7L,
              "historical piecewise full inverse score");
      }
      std::printf(
          "{\"case\":\"%s\",\"point\":%u,\"relative_score\":%.17g,"
          "\"reference_refinement\":%.21Lg,\"factor_error\":%.21Lg,"
          "\"estimated_factor_allowance\":%.21Lg,\"background_error\":%.21Lg,"
          "\"production_refinement\":%.21Lg,\"analytic_provider\":%s,"
          "\"magnitude_error\":%.21Lg,\"total_score_error\":%.21Lg,"
          "\"reference_pass\":%s,\"factor_pass\":%s,\"background_pass\":%s,"
          "\"passed\":%s}\n",
          name, m, s.score->relative_profile_score, reference_refinement,
          factor_error, allowance, background_error, refinement,
          analytic ? "true" : "false", mu_error, score_error,
          reference_refinement <= 2e-7L ? "true" : "false",
          factor_error <= 3e-7L && allowance <= 3e-7L ? "true" : "false",
          background_error <= 5e-7L && refinement <= 5e-7L && mu_error <= 1e-8L
              ? "true"
              : "false",
          pass ? "true" : "false");
    }
    all_pass &= pass;
  }
  check(all_pass, "all current-Pantheon named fixedpoint budgets");
  std::printf("{\"suite\":\"%s\",\"models\":%zu,\"n\":%zu,"
              "\"reference_condition_estimate\":%.21Lg,\"passed\":true}\n",
              all_named ? "current1590_named29"
                        : "grey12_current1590_conditional_not1820",
              grid.size(), n, Cnorm * inverse_norm);
}
} // namespace
int main(int argc, char **argv) {
  if (argc == 4 && std::string_view(argv[3]) == "--verified-original-assets") {
    actual(argv[1], argv[2]);
    return 0;
  }
  throw std::runtime_error("optional exact-source guard required: COV TABLE "
                           "--verified-original-assets; C++ asserts external "
                           "guard provenance, does not hash files");
}
