// Optional developer reference; no external data dependency in ordinary CI.
// Same-input arithmetic comparison: parse binary64, select zHD>0.01, promote
// to long double LDLT. No averaging/repair. Source provenance is pinned in
// test_statistics.cpp; caller must check original asset hashes before run.
// Predeclared constructed-solve component screen1e-10, score reference2e-7.
// This r=Cv probe validates the reference, not a cosmological W01 target.
#include "fixtures/ldlt_reference.hpp"
#ifdef IRRED_PRODUCTION_PREPARATION_PROBE
#include "irred/statistics.hpp"
#endif
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error("usage: test_w01_reference COVARIANCE TABLE");
    std::ifstream cov(argv[1]), table(argv[2]);
    std::size_t n;
    if (!(cov >> n) || n != 1701)
      throw std::runtime_error("expected released1701 covariance");
    std::vector<double> full(n * n);
    for (auto &v : full)
      if (!(cov >> v) || !std::isfinite(v))
        throw std::runtime_error("invalid covariance");
    std::string line, token;
    if (cov >> token)
      throw std::runtime_error("trailing matrix");
    std::getline(table, line);
    std::istringstream header(line);
    std::vector<std::string> names;
    while (header >> token)
      names.push_back(token);
    auto pos = std::find(names.begin(), names.end(), "zHD");
    if (pos == names.end())
      throw std::runtime_error("missing zHD");
    std::size_t zi = pos - names.begin(), row = 0;
    std::vector<std::size_t> kept;
    while (std::getline(table, line)) {
      if (line.empty())
        continue;
      std::istringstream in(line);
      std::vector<std::string> fields;
      while (in >> token)
        fields.push_back(token);
      if (fields.size() != names.size() || row >= n)
        throw std::runtime_error("table shape");
      auto z = std::stod(fields[zi]);
      if (!std::isfinite(z))
        throw std::runtime_error("invalid redshift");
      if (z > .01)
        kept.push_back(row);
      ++row;
    }
    if (row != n || kept.size() != 1590)
      throw std::runtime_error("unexpected selection");
    auto k = kept.size();
    std::vector<double> c(k * k), r(k);
    std::vector<long double> v(k);
    for (std::size_t i = 0; i < k; ++i) {
      v[i] = static_cast<long double>(static_cast<int>(i % 9) - 4) / 8;
      for (std::size_t j = 0; j < k; ++j)
        c[i * k + j] = full[kept[i] * n + kept[j]];
    }
#ifdef IRRED_PRODUCTION_PREPARATION_PROBE
#ifdef IRRED_WIDE_PREPARATION_PROBE
    auto wide_factor = irred::numerics::cholesky(
        c, k, k * k, irred::numerics::Arithmetic::longdouble_cpu_v1);
    std::vector<double> wide_ones(k, 1);
    auto wide_response = irred::numerics::solve(wide_factor, wide_ones, 1e-10);
    std::printf("{\"suite\":\"W01_explicitwide_preparation_diagnostic\","
                "\"policy\":\"F02/longdouble-cpu/"
                "v1\",\"factor_status\":%u,\"solve_status\":%u,\"condition_"
                "estimate_inf\":%.17g,\"pre_cast_sensitivity\":%.17g,\"cast_"
                "inf\":%.17g,\"cast_relative\":%.17g,\"post_cast_backward\":%."
                "17g,\"returned_sensitivity\":%.17g,\"qualified\":false}\n",
                static_cast<unsigned>(wide_factor.status()),
                static_cast<unsigned>(wide_response.status),
                wide_factor.condition_estimate_inf(),
                wide_response.pre_cast_sensitivity_estimate,
                wide_response.output_rounding_error_inf,
                wide_response.output_rounding_error_relative,
                wide_response.post_cast_backward_residual,
                wide_response.estimated_forward_sensitivity);
    return 0;
#endif
    irred::statistics::Metadata metadata;
    for (std::size_t i = 0; i < k; ++i)
      metadata.ordered_ids.push_back("original-row:" + std::to_string(kept[i]));
    metadata.measure = "product d(magnitude)";
    metadata.ordering_provenance =
        "original source indices in strict ascending selected order";
    auto gaussian = irred::statistics::prepare_gaussian(
        c, irred::statistics::MatrixKind::covariance, std::move(metadata),
        k * k, 1e-10);
    std::vector<double> ones(k, 1);
    auto ordered_ids = gaussian.metadata().ordered_ids;
    auto profile =
        std::move(gaussian).prepare_offset_profile(ones, ordered_ids, 1e-10);
    std::printf("{\"suite\":\"W01_initial_production_preparation_probe\","
                "\"guard\":1e-10,\"profile_status\":%u,\"numerical_status\":%u,"
                "\"qualified\":false}\n",
                static_cast<unsigned>(profile.status()),
                static_cast<unsigned>(profile.numerical_status()));
    return 0;
#endif
    long double rounding = 0;
    for (std::size_t i = 0; i < k; ++i) {
      long double sum = 0;
      for (std::size_t j = 0; j < k; ++j)
        sum += static_cast<long double>(c[i * k + j]) * v[j];
      r[i] = static_cast<double>(sum);
      rounding = std::max(rounding, std::abs(sum - r[i]));
    }
    auto f = irred::test_reference::factor(c, k);
    if (!f.valid)
      throw std::runtime_error("selected reference not SPD");
    auto x = irred::test_reference::solve(f, r);
    long double err = 0, residual_inf = 0, matrix_inf = 0, rhs_inf = 0,
                solution_inf = 0, reconstruction_sample = 0;
    for (std::size_t i = 0; i < k; ++i)
      err = std::max(err, std::abs(x[i] - v[i]));
    for (std::size_t i = 0; i < k; ++i) {
      long double ax = 0, rowsum = 0;
      for (std::size_t j = 0; j < k; ++j) {
        ax += static_cast<long double>(c[i * k + j]) * x[j];
        rowsum += std::abs(static_cast<long double>(c[i * k + j]));
      }
      residual_inf = std::max(residual_inf, std::abs(ax - r[i]));
      matrix_inf = std::max(matrix_inf, rowsum);
      rhs_inf = std::max(rhs_inf, std::abs(static_cast<long double>(r[i])));
      solution_inf = std::max(solution_inf, std::abs(x[i]));
      // Every diagonal and one deterministic offdiagonal per row; not a full
      // reconstruction proof.
      for (auto j : {i, (i * 17 + 3) % k}) {
        long double value = 0;
        for (std::size_t t = 0; t <= std::min(i, j); ++t)
          value += f.lower[i * k + t] * f.diagonal[t] * f.lower[j * k + t];
        reconstruction_sample =
            std::max(reconstruction_sample, std::abs(value - c[i * k + j]));
      }
    }
    std::vector<long double> inverse_rows(k);
    std::vector<double> basis(k);
    for (std::size_t j = 0; j < k; ++j) {
      basis[j] = 1;
      auto col = irred::test_reference::solve(f, basis);
      basis[j] = 0;
      for (std::size_t i = 0; i < k; ++i)
        inverse_rows[i] += std::abs(col[i]);
    }
    auto inverse_inf =
        *std::max_element(inverse_rows.begin(), inverse_rows.end());
    auto condition = matrix_inf * inverse_inf;
    const auto rhs_rounding_solution_estimate = inverse_inf * rounding;
    std::printf("{\"reference_uncertainty\":\"empirical\",\"rhs_rounding_"
                "solution_estimate\":%.21Lg,\"long_double_digits\":%d}\n",
                rhs_rounding_solution_estimate,
                std::numeric_limits<long double>::digits);
    auto backward = residual_inf / (matrix_inf * solution_inf + rhs_inf);
    std::printf(
        "{\"suite\":\"W01_reference_diagnostics\",\"residual_inf\":%.21Lg,"
        "\"backward_residual\":%.21Lg,\"matrix_norm_inf\":%.21Lg,\"sampled_"
        "reconstruction_error\":%.21Lg,\"condition_estimate_inf\":%.21Lg}\n",
        residual_inf, backward, matrix_inf, reconstruction_sample, condition);
    std::printf("{\"suite\":\"W01_same_binary64_LDLT_constructed\",\"n\":%zu,"
                "\"max_solution_error\":%.21Lg,\"rhs_rounding_inf\":%.21Lg,"
                "\"logdet\":%.21Lg,\"passed\":%s}\n",
                k, err, rounding, f.logdet, err <= 1e-10L ? "true" : "false");
    return err <= 1e-10L ? 0 : 1;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
