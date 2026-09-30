// Optional developer benchmark: caller SHA2 wrapper verifies original assets.
// One serial process, frozen seven-point fine policy, three retained repeats.
// Timing is not qualification: accepted fixture/guard checks remain required.
// Factor-only mode reports factorization+condition together; opaque API offers
// no non-invasive split. No excluded measurements or automatic faster backend.
#include "fixtures/w01_historical.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
using namespace irred;
namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}
std::vector<std::string> words(const std::string &s) {
  std::istringstream in(s);
  std::vector<std::string> w;
  for (std::string x; in >> x;)
    w.push_back(x);
  return w;
}
observations::Prepared read(const char *cov_path, const char *table_path) {
  observations::Input in{};
  in.profile = observations::Profile::pantheon_plus_released_v1;
  in.role = observations::Role::released_fitted_summary;
  in.unit = observations::Unit::magnitude;
  in.calibration = observations::Calibration::released_corrected;
  in.uncertainty = observations::Uncertainty::covariance;
  in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  in.component = observations::Component::total;
  in.ordering_provenance = "supplied original release row order";
  in.quality_dictionary = "absent decoded quality; zeros are placeholders";
  in.table_sha256 = std::string(test_reference::w01_table_sha256);
  in.uncertainty_sha256 = std::string(test_reference::w01_cov_sha256);
  std::ifstream table(table_path), cov(cov_path);
  std::string line;
  std::getline(table, line);
  auto names = words(line);
  auto column = [&](const char *n) {
    auto p = std::find(names.begin(), names.end(), n);
    if (p == names.end())
      throw std::runtime_error("missing column");
    return static_cast<std::size_t>(p - names.begin());
  };
  auto zm = column("zHD"), zh = column("zHEL"), zc = column("zCMB"),
       mag = column("m_b_corr"), cid = column("CID"),
       survey = column("IDSURVEY");
  while (std::getline(table, line)) {
    if (line.empty())
      continue;
    auto r = words(line);
    if (r.size() != names.size())
      throw std::runtime_error("table shape");
    auto i = in.values.size();
    in.values.push_back(std::stod(r[mag]));
    in.zhd.push_back(std::stod(r[zm]));
    in.zhel.push_back(std::stod(r[zh]));
    in.zcmb.push_back(std::stod(r[zc]));
    in.measurement_ids.push_back(r[cid] + "/" + r[survey] +
                                 "/original-row:" + std::to_string(i));
    in.event_ids.push_back(r[cid]);
  }
  std::size_t n;
  if (!(cov >> n) || n != 1701 || in.values.size() != n)
    throw std::runtime_error("source shape");
  in.uncertainty_matrix.resize(n * n);
  for (auto &x : in.uncertainty_matrix)
    if (!(cov >> x) || !std::isfinite(x))
      throw std::runtime_error("matrix value");
  if (cov >> line)
    throw std::runtime_error("trailing matrix");
  in.uncertainty_axis_ids = in.measurement_ids;
  in.missing.assign(n, 0);
  in.zhd_missing.assign(n, 0);
  in.zcmb_missing.assign(n, 0);
  in.zhel_missing.assign(n, 0);
  in.quality.assign(n, 0);
  in.source_selection.assign(n, 1);
  auto out = observations::prepare(std::move(in), {1701, 1701 * 1701, 1048576});
  if (out.status() != observations::Status::ok)
    throw std::runtime_error("typed source failure");
  return out;
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 5 || std::string(argv[3]) != "--verified-original-assets")
      throw std::runtime_error(
          "usage: benchmark_w01 COV TABLE --verified-original-assets "
          "factor|consumer (SHA2 wrapper required)");
    auto start = Clock::now();
    auto source = read(argv[1], argv[2]);
    auto acquired = Clock::now();
    auto mode = std::string(argv[4]);
    if (mode == "factor") {
      auto selected =
          source.select(observations::Selection::pantheon_zhd_gt_001);
      auto k = selected.source_indices.size();
      if (k != 1590)
        throw std::runtime_error("selection");
      std::vector<double> c(k * k);
      for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j)
          c[i * k + j] =
              source.source()
                  .uncertainty_matrix[selected.source_indices[i] * 1701 +
                                      selected.source_indices[j]];
      auto ready = Clock::now();
      auto f = numerics::cholesky(c, k, k * k,
                                  numerics::Arithmetic::longdouble_cpu_v1);
      auto factored = Clock::now();
      std::vector<double> ones(k, 1);
      auto x = numerics::solve(f, ones, 1e-10);
      auto solved = Clock::now();
      if (f.status() != numerics::Status::ok ||
          x.status != numerics::Status::ok)
        throw std::runtime_error("factor quality failure");
      std::printf(
          "{\"mode\":\"factor\",\"acquisition_seconds\":%.9g,\"selection_"
          "seconds\":%.9g,\"factor_and_condition_seconds\":%.9g,\"response_"
          "solve_seconds\":%.9g,\"end_to_end_seconds\":%.9g,\"condition_"
          "estimate\":%.17g,\"returned_sensitivity\":%.17g,\"passed\":true}\n",
          seconds(start, acquired), seconds(acquired, ready),
          seconds(ready, factored), seconds(factored, solved),
          seconds(start, solved), f.condition_estimate_inf(),
          x.estimated_forward_sensitivity);
      return 0;
    }
    if (mode != "consumer")
      throw std::runtime_error("unknown phase");
    supernova::Policy policy;
    policy.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
    policy.maximum_forward_sensitivity = 1e-10;
    policy.background.integration = {1e-12, 1e-12, 100000, 30};
    policy.background.maximum_total_evaluations = 20000000;
    auto consumer = supernova::prepare(std::move(source), policy);
    auto prepared = Clock::now();
    if (consumer.status() != supernova::Status::ok)
      throw std::runtime_error("preparation quality failure");
    std::array<supernova::ModelPoint, 7> grid;
    for (std::size_t i = 0; i < grid.size(); ++i) {
      auto p = test_reference::w01_historical[i];
      grid[i] = {p.lcdm ? cosmology::Model::flat_lcdm_late_v1
                        : cosmology::Model::constant_q_flat_v1,
                 p.lcdm ? p.parameter : 0, p.lcdm ? 0 : p.parameter};
    }
    std::printf("{\"mode\":\"consumer\",\"acquisition_seconds\":%.9g,"
                "\"preparation_seconds\":%.9g,\"arithmetic\":\"F02/"
                "longdouble-cpu/v1\",\"guard\":1e-10,\"repeats\":3}\n",
                seconds(start, acquired), seconds(acquired, prepared));
    bool passed = true;
    for (unsigned repeat = 0; repeat < 3; ++repeat) {
      auto before = Clock::now();
      auto result = consumer.evaluate_batch(grid, policy);
      auto after = Clock::now();
      if (result.slots.size() != grid.size())
        throw std::runtime_error("batch shape");
      std::size_t callbacks = 0;
      for (std::size_t i = 0; i < grid.size(); ++i) {
        auto &s = result.slots[i];
        if (s.status != supernova::Status::ok ||
            s.numerical_status != numerics::Status::ok) {
          passed = false;
          callbacks += s.background_evaluations;
          std::printf(
              "{\"repeat\":%u,\"point\":%zu,\"status\":%u,\"passed\":false}\n",
              repeat, i, static_cast<unsigned>(s.status));
          continue;
        }
        auto err = std::abs(s.relative_profile_score -
                            test_reference::w01_historical[i].stable_score);
        bool ok = s.status == supernova::Status::ok &&
                  s.numerical_status == numerics::Status::ok && err <= 2e-7;
        passed &= ok;
        callbacks += s.background_evaluations;
        std::printf("{\"repeat\":%u,\"point\":%zu,\"relative_score\":%.17g,"
                    "\"reference_error\":%.17g,\"passed\":%s}\n",
                    repeat, i, s.relative_profile_score, err,
                    ok ? "true" : "false");
      }
      std::printf(
          "{\"repeat\":%u,\"retained_batch_seconds\":%.9g,\"callbacks\":%zu}\n",
          repeat, seconds(before, after), callbacks);
    }
    std::printf(
        "{\"mode\":\"consumer\",\"end_to_end_seconds\":%.9g,\"passed\":%s}\n",
        seconds(start, Clock::now()), passed ? "true" : "false");
    return passed ? 0 : 1;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
