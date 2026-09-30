// Optional native developer comparison, explicit original assets, not CI data.
// Frozen7 points/3 quadrature policies; score budgets: reference2e-7,
// factor+solve3e-7, background5e-7, total1e-6; magnitude screen1e-8.
// Same canonical binary64 selected C promoted into independent LDLT.
// Reference LCDM .3 fixed panels4096/8192/16384; others analytic limits.
// Shared libm and same statistics author disclosed: independent algorithm,
// not an independent author/interval proof. Never edits original files.
#include "fixtures/ldlt_reference.hpp"
#include "fixtures/w01_historical.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
using namespace irred;
namespace {
std::vector<std::string> words(const std::string &s) {
  std::istringstream in(s);
  std::vector<std::string> w;
  for (std::string x; in >> x;)
    w.push_back(x);
  return w;
}
long double integral(const supernova::ModelPoint &p, long double z,
                     unsigned panels) {
  if (p.model == cosmology::Model::constant_q_flat_v1) {
    auto power = -static_cast<long double>(p.constant_q);
    return power == 0 ? std::log1p(z)
                      : std::expm1(power * std::log1p(z)) / power;
  }
  if (p.omega_m == 0)
    return z;
  if (p.omega_m == 1)
    return 2 * (1 - 1 / std::sqrt(1 + z));
  auto f = [&](long double x) {
    auto u = 1 + x;
    return 1 / std::sqrt(p.omega_m * u * u * u + 1 - p.omega_m);
  };
  auto h = z / panels;
  long double sum = f(0) + f(z);
  for (unsigned i = 1; i < panels; ++i)
    sum += (i % 2 ? 4 : 2) * f(h * i);
  return h * sum / 3;
}
struct RefScore {
  long double q = 0, a = 0;
};
RefScore profile(const test_reference::LDLT &f, std::span<const long double> r,
                 std::span<const long double> wone, long double gram) {
  auto wr = test_reference::solve_longdouble(f, r);
  long double b = 0;
  for (auto x : wr)
    b += x;
  auto a = b / gram;
  std::vector<long double> u(r.size());
  for (std::size_t i = 0; i < r.size(); ++i)
    u[i] = r[i] - a;
  auto wu = test_reference::solve_longdouble(f, u);
  long double q = 0;
  for (std::size_t i = 0; i < r.size(); ++i)
    q += u[i] * wu[i];
  (void)wone;
  return {q, a};
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 3 && argc != 4)
      throw std::runtime_error(
          "usage: test_w01_scores COV TABLE; caller validates original hashes");
    const bool verified_historical =
        argc == 4 && std::string(argv[3]) == "--verified-original-assets";
    if (argc == 4 && !verified_historical)
      throw std::runtime_error("unknown verification flag");
    observations::Input in{};
    in.profile = observations::Profile::pantheon_plus_released_v1;
    in.role = observations::Role::released_fitted_summary;
    in.unit = observations::Unit::magnitude;
    in.calibration = observations::Calibration::released_corrected;
    in.uncertainty = observations::Uncertainty::covariance;
    in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
    in.component = observations::Component::total;
    in.ordering_provenance = "supplied release row ordering declaration; "
                             "independently selected symmetry audited";
    in.quality_dictionary = "absent decoded quality column; zero placeholders "
                            "are not all-clear flags";
    in.table_sha256 =
        "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8";
    in.uncertainty_sha256 =
        "abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc";
    std::ifstream table(argv[2]), cov(argv[1]);
    std::string line;
    std::getline(table, line);
    auto names = words(line);
    auto column = [&](const char *n) {
      auto p = std::find(names.begin(), names.end(), n);
      if (p == names.end())
        throw std::runtime_error("missing source column");
      return static_cast<std::size_t>(p - names.begin());
    };
    auto zm = column("zHD"), zh = column("zHEL"), zc = column("zCMB"),
         mag = column("m_b_corr"), cid = column("CID"),
         survey = column("IDSURVEY");
    while (std::getline(table, line)) {
      if (line.empty())
        continue;
      auto row = words(line);
      if (row.size() != names.size())
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
      throw std::runtime_error("source dimension");
    in.uncertainty_matrix.resize(n * n);
    for (auto &x : in.uncertainty_matrix)
      if (!(cov >> x) || !std::isfinite(x))
        throw std::runtime_error("matrix token");
    if (cov >> line)
      throw std::runtime_error("matrix trailing bytes");
    in.uncertainty_axis_ids = in.measurement_ids;
    in.missing.assign(n, 0);
    in.zhd_missing.assign(n, 0);
    in.zcmb_missing.assign(n, 0);
    in.zhel_missing.assign(n, 0);
    in.quality.assign(n, 0);
    in.source_selection.assign(n, 1);
    auto obs =
        observations::prepare(std::move(in), {1701, 1701 * 1701, 1048576});
    if (obs.status() != observations::Status::ok)
      throw std::runtime_error("typed source reject");
    auto selected = obs.select(observations::Selection::pantheon_zhd_gt_001);
    auto k = selected.source_indices.size();
    if (selected.status != observations::Status::ok || k != 1590)
      throw std::runtime_error("selection mismatch");
    std::vector<double> c(k * k);
    long double cnorm = 0;
    for (std::size_t i = 0; i < k; ++i) {
      long double rowsum = 0;
      for (std::size_t j = 0; j < k; ++j) {
        c[i * k + j] =
            obs.source().uncertainty_matrix[selected.source_indices[i] * n +
                                            selected.source_indices[j]];
        rowsum += std::abs(static_cast<long double>(c[i * k + j]));
      }
      cnorm = std::max(cnorm, rowsum);
    }
    auto ref = test_reference::factor(c, k);
    if (!ref.valid)
      throw std::runtime_error("reference factor reject");
    std::vector<long double> one(k, 1);
    auto wone = test_reference::solve_longdouble(ref, one);
    long double gref = 0;
    for (auto x : wone)
      gref += x;
    std::vector<long double> inverse_rows(k), basis(k);
    for (std::size_t j = 0; j < k; ++j) {
      basis[j] = 1;
      auto col = test_reference::solve_longdouble(ref, basis);
      basis[j] = 0;
      for (std::size_t i = 0; i < k; ++i)
        inverse_rows[i] += std::abs(col[i]);
    }
    auto current_inverse_norm_estimate =
        *std::max_element(inverse_rows.begin(), inverse_rows.end());
    std::printf("{\"reference_inverse_norm_estimate\":%.21Lg,\"reference_"
                "condition_estimate\":%.21Lg,\"meaning\":\"current same-input "
                "empirical LDLT basis solves\"}\n",
                current_inverse_norm_estimate,
                cnorm * current_inverse_norm_estimate);
    supernova::Policy policy;
    policy.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
    policy.maximum_forward_sensitivity = 1e-10;
    policy.background.maximum_total_evaluations = 20000000;
    auto consumer = supernova::prepare(std::move(obs), policy);
    std::printf("{\"preparation_status\":%u,\"arithmetic\":\"F02/"
                "longdouble-cpu/v1\",\"n\":%zu}\n",
                static_cast<unsigned>(consumer.status()), k);
    if (consumer.status() != supernova::Status::ok)
      return 1;
    std::array<supernova::ModelPoint, 7> grid{
        {{cosmology::Model::flat_lcdm_late_v1, 0, 0},
         {cosmology::Model::flat_lcdm_late_v1, .3, 0},
         {cosmology::Model::flat_lcdm_late_v1, 1, 0},
         {cosmology::Model::constant_q_flat_v1, 0, -1},
         {cosmology::Model::constant_q_flat_v1, 0, -.5},
         {cosmology::Model::constant_q_flat_v1, 0, 0},
         {cosmology::Model::constant_q_flat_v1, 0, .5}}};
    std::array<supernova::BatchResult, 3> batches;
    for (unsigned level = 0; level < 3; ++level) {
      auto eps = level == 0 ? 1e-8 : level == 1 ? 1e-10 : 1e-12;
      policy.background.integration = {eps, eps, 100000, 30};
      batches[level] = consumer.evaluate_batch(grid, policy);
      if (batches[level].slots.size() != grid.size())
        throw std::runtime_error("batch shape");
    }
    bool passed = true;
    for (std::size_t m = 0; m < grid.size(); ++m) {
      for (auto &b : batches)
        if (b.slots[m].status != supernova::Status::ok)
          throw std::runtime_error("candidate slot failed");
      auto &s = batches[2].slots[m];
      std::vector<long double> same(s.base_residuals.begin(),
                                    s.base_residuals.end());
      auto matched = profile(ref, same, wone, gref);
      std::array<RefScore, 3> references;
      long double magerror = 0;
      for (unsigned level = 0; level < 3; ++level) {
        std::vector<long double> r(k);
        for (std::size_t i = 0; i < k; ++i) {
          auto source = selected.source_indices[i];
          auto query = consumer.selected_queries()[i];
          auto mu =
              5 *
              std::log10((1 + static_cast<long double>(query.z_observer)) *
                         integral(grid[m], query.z_expansion, 4096u << level));
          r[i] = consumer.observations().source().values[source] - mu;
          if (level == 2)
            magerror = std::max(magerror, std::abs(mu - s.shape_magnitudes[i]));
        }
        references[level] = profile(ref, r, wone, gref);
      }
      auto factor_error = std::abs(-.5L * s.quadratic + .5L * matched.q);
      auto reference_refinement =
          .5L * std::abs(references[2].q - references[1].q);
      auto background_error = .5L * std::abs(matched.q - references[2].q);
      auto production_refinement =
          std::abs(s.relative_profile_score -
                   batches[1].slots[m].relative_profile_score);
      const auto &cache = consumer.profile_operator();
      auto wx = cache.cached_response_solution();
      long double wxinf = 0, sumwx = 0, h = 0, habs = 0;
      for (std::size_t i = 0; i < k; ++i) {
        wxinf = std::max(wxinf, std::abs(static_cast<long double>(wx[i])));
        sumwx += std::abs(static_cast<long double>(wx[i]));
        h += static_cast<long double>(wx[i]) * s.profiled_residuals[i];
        habs +=
            std::abs(static_cast<long double>(wx[i]) * s.profiled_residuals[i]);
      }
      auto gamma =
          (2 * k + 2) * std::numeric_limits<long double>::epsilon() /
          (1 - (2 * k + 2) * std::numeric_limits<long double>::epsilon());
      auto ex = cache.cached_response_forward_sensitivity() * wxinf;
      auto er = s.solve_diagnostics.coefficient_solve_forward_sensitivity *
                s.solve_diagnostics.solution_norm_inf;
      auto eu = s.solve_diagnostics.estimated_forward_sensitivity *
                s.solve_diagnostics.adjusted_solution_norm_inf;
      auto dg = k * ex + gamma * sumwx;
      auto db = k * er + gamma * k * s.solve_diagnostics.solution_norm_inf;
      auto gram = cache.gram();
      if (dg >= gram)
        throw std::runtime_error("Gram estimate failure");
      auto da = (db + std::abs(s.offset_coefficient) * dg) / (gram - dg) +
                std::numeric_limits<double>::epsilon() *
                    std::abs(s.offset_coefficient) /
                    (1 - std::numeric_limits<double>::epsilon());
      auto dh = s.solve_diagnostics.adjusted_residual_l1 * ex + gamma * habs;
      auto score_estimate =
          .5L * s.solve_diagnostics.adjusted_residual_l1 * eu +
          .5L * gamma * s.solve_diagnostics.adjusted_residual_l1 *
              s.solve_diagnostics.adjusted_solution_norm_inf +
          (std::abs(h) + dh) * da + .5L * (gram + dg) * da * da;
      long double delta_l1 = 0, delta_inf = 0;
      for (std::size_t i = 0; i < k; ++i) {
        auto delta =
            std::abs(static_cast<long double>(s.profiled_residuals[i]) -
                     (static_cast<long double>(s.base_residuals[i]) -
                      s.offset_coefficient));
        delta_l1 += delta;
        delta_inf = std::max(delta_inf, delta);
      }
      auto inverse_norm_estimate = current_inverse_norm_estimate;
      score_estimate +=
          .5L *
              (std::numeric_limits<double>::epsilon() /
               (1 - std::numeric_limits<double>::epsilon())) *
              std::abs(s.quadratic) +
          delta_l1 * s.solve_diagnostics.adjusted_solution_norm_inf +
          .5L * delta_l1 * inverse_norm_estimate * delta_inf;
      bool ok = factor_error <= 3e-7L && score_estimate <= 3e-7L &&
                reference_refinement <= 2e-7L && background_error <= 5e-7L &&
                production_refinement <= 5e-7L && magerror <= 1e-8L &&
                factor_error + reference_refinement + background_error <= 1e-6L;
      if (verified_historical) {
        const auto &historic = test_reference::w01_historical[m];
        if (historic.lcdm !=
                (grid[m].model == cosmology::Model::flat_lcdm_late_v1) ||
            historic.parameter !=
                (historic.lcdm ? grid[m].omega_m : grid[m].constant_q))
          throw std::runtime_error("historical fixture point mismatch");
        auto stable_error =
                 std::abs(s.relative_profile_score - historic.stable_score),
             full_error =
                 std::abs(s.relative_profile_score - historic.full_score),
             compressed_error =
                 std::abs(s.relative_profile_score - historic.compressed_score);
        ok &= stable_error <= 2e-7 && full_error <= 2e-7;
        std::printf("{\"historical_point\":%zu,\"stable_score_error\":%.17g,"
                    "\"full_projection_error\":%.17g,\"compressed_spline_"
                    "error\":%.17g,\"compressed_precise_reference_pass\":%s,"
                    "\"source_identity\":\"required external SHA2 wrapper\"}\n",
                    m, stable_error, full_error, compressed_error,
                    compressed_error <= 2e-7 ? "true" : "false");
      }
      passed &= ok;
      std::printf(
          "{\"point\":%zu,\"model\":%u,\"omega_m\":%.17g,\"q\":%.17g,"
          "\"relative_score\":%.17g,\"offset\":%.17g,\"factor_score_error\":%."
          "21Lg,\"estimated_factor_score_allowance\":%.21Lg,\"reference_"
          "refinement\":%.21Lg,\"background_score_error\":%.21Lg,\"coarse_"
          "medium_trend\":%.17g,\"production_"
          "refinement\":%.17g,\"max_magnitude_error\":%.21Lg,\"passed\":%s}\n",
          m, static_cast<unsigned>(grid[m].model), grid[m].omega_m,
          grid[m].constant_q, s.relative_profile_score, s.offset_coefficient,
          factor_error, score_estimate, reference_refinement, background_error,
          std::abs(batches[1].slots[m].relative_profile_score -
                   batches[0].slots[m].relative_profile_score),
          production_refinement, magerror, ok ? "true" : "false");
    }
    return passed ? 0 : 1;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
