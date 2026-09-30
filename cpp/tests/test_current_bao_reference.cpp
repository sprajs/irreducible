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
long double reference_integral(long double z, const cosmology::CPL &p,
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
                        const std::vector<bao::ModelPoint> &models,
                        const bao::DensityBatch &candidate,
                        const std::vector<double> &previous_scores) {
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
  std::vector<long double> prior_scores(models.size());
  std::vector<std::vector<long double>> prior_predictions(models.size());
  for (int panels : {16, 32, 64})
    for (std::size_t m = 0; m < models.size(); ++m) {
      std::vector<long double> predictions(n), residual(n);
      const auto &p = std::get<cosmology::CPL>(models[m].expansion);
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
      const auto score = -.5L * (q + f.logdet + normalization);
      if (!std::isfinite(score) || !std::isfinite(allowance))
        throw std::runtime_error("nonfinite independent reference");
      if (panels == 64) {
        const auto &c = candidate.slots.at(m);
        const auto &g = *c.result;
        const auto half_ulp = [](double value) {
          return .5L * std::max(std::abs((long double)std::nextafter(value,
                                                                     INFINITY) -
                                         value),
                                std::abs((long double)value -
                                         std::nextafter(value, -INFINITY)));
        };
        std::vector<long double> canonical(n);
        for (std::size_t i = 0; i < n; ++i)
          canonical[i] = c.residuals.at(i);
        const auto xc = test_reference::solve_longdouble(f, canonical);
        long double qc = 0, ac = 0, cl1 = 0, xinf = 0, canonical_residual = 0;
        long double dl1 = 0, dinf = 0, background = 0;
        const auto dx = inverse_bound * residual_bound;
        for (std::size_t i = 0; i < n; ++i) {
          qc += canonical[i] * xc[i];
          ac += std::abs(canonical[i] * xc[i]);
          cl1 += std::abs(canonical[i]);
          xinf = std::max(xinf, std::abs(xc[i]));
          long double cx = 0, abs_cx = 0;
          for (std::size_t j = 0; j < n; ++j) {
            auto t = (long double)input.covariance[i * n + j] * xc[j];
            cx += t;
            abs_cx += std::abs(t);
          }
          canonical_residual =
              std::max(canonical_residual,
                       std::abs(canonical[i] - cx) +
                           gamma * (abs_cx + std::abs(canonical[i])));
          const auto d =
              std::abs(canonical[i] - residual[i]) +
              std::abs(predictions[i] - prior_predictions[m][i]) +
              half_ulp(static_cast<double>(predictions[i])) +
              epsilon / (1 - epsilon) *
                  (std::abs(input.observed[i]) + std::abs(predictions[i]));
          dl1 += d;
          dinf = std::max(dinf, d);
          background += (std::abs(x[i]) + dx) * d;
        }
        background += .5L * inverse_bound * dl1 * dinf;
        const auto canonical_allowance =
            .5L * (cl1 * inverse_bound * canonical_residual + gamma * ac +
                   logdet_allowance + 64 * n * epsilon);
        const auto canonical_score = -.5L * (qc + f.logdet + normalization);
        const auto eta = (long double)g.estimated_forward_sensitivity;
        if (!(eta >= 0 && eta < 1) || !std::isfinite(eta))
          throw std::runtime_error("invalid candidate sensitivity");
        // Reported q/normalization casts are diagnostics, not density
        // dependencies.
        const auto dot_gamma =
            (2 * n + 2) * epsilon / (1 - (2 * n + 2) * epsilon);
        const auto factor =
            .5L * cl1 * xinf * eta / (1 - eta) + canonical_allowance +
            dot_gamma * ac / 2 + half_ulp(g.density.log_value) +
            .5L * std::abs((long double)g.log_determinant - f.logdet) +
            .5L * std::abs(n * 1.8378770664093454835606594728112352797228L -
                           normalization);
        const auto reference_refinement = std::abs(score - prior_scores[m]);
        const auto production_refinement =
            std::abs(g.density.log_value - previous_scores.at(m));
        const bool reference_pass = allowance + reference_refinement <= 2e-9L;
        const bool factor_pass =
            std::isfinite(factor) && factor <= 3e-9L &&
            std::abs((long double)g.density.log_value - canonical_score) +
                    canonical_allowance <=
                3e-9L;
        const bool background_pass =
            std::isfinite(background) && background <= 5e-9L;
        const bool refinement_pass = production_refinement <= 5e-9L;
        const bool total_pass =
            std::abs((long double)g.density.log_value - score) + allowance <=
            1e-8L;
        std::printf(
            "{\"point\":%zu,\"reference_allowance\":%.21Lg,\"reference_"
            "refinement\":%.21Lg,\"factor_allowance\":%.21Lg,\"background_"
            "allowance\":%.21Lg,\"production_refinement\":%.21g,\"reference_"
            "pass\":%s,\"factor_pass\":%s,\"background_pass\":%s,\"refinement_"
            "pass\":%s,\"total_pass\":%s}\n",
            m, allowance, reference_refinement, factor, background,
            production_refinement, reference_pass ? "true" : "false",
            factor_pass ? "true" : "false", background_pass ? "true" : "false",
            refinement_pass ? "true" : "false", total_pass ? "true" : "false");
        if (!(reference_pass && factor_pass && background_pass &&
              refinement_pass && total_pass))
          throw std::runtime_error(
              "original allocated independent gate failed");
      }
      prior_scores[m] = score;
      prior_predictions[m] = std::move(predictions);
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
  p.maximum_queries = 13;
  p.maximum_native_bytes = 1 << 20;
  p.background.maximum_queries = 13;
  p.background.maximum_callbacks = 2000000;
  p.background.maximum_segment_visits = 2000;
  p.background.maximum_native_bytes = 1 << 20;
  bao::PreparationPolicy preparation{
      13, 169, 10000, 1000000, 1e-10, numerics::Arithmetic::longdouble_cpu_v1};
  bao::DensityPolicy dp{
      p, 20, 1000000, 1e-10, numerics::Arithmetic::longdouble_cpu_v1, 7};
  auto prepared = bao::prepare_density(input, preparation);
  if (prepared.status() != statistics::DensityStatus::finite)
    throw std::runtime_error("original Gaussian preparation failed");
  std::vector<bao::ModelPoint> models;
  for (double om : {0., .3, 1.})
    models.push_back(
        {cosmology::CPL{om, -1, 0}, cosmology::FlatFLRW{}, bao::Ruler(10000)});
  models.push_back(
      {cosmology::CPL{.3, -.8, -.5}, cosmology::FlatFLRW{}, bao::Ruler(10000)});
  models.push_back({cosmology::CPL{.353, -.42, -1.75}, cosmology::FlatFLRW{},
                    bao::Ruler(10000)});
  for (double hrd : {5000., 15000.})
    for (double om : {0., .3, 1.})
      models.push_back(
          {cosmology::CPL{om, -1, 0}, cosmology::FlatFLRW{}, bao::Ruler(hrd)});
  bao::DensityBatch fine;
  std::vector<double> previous_scores;
  for (double tolerance : {1e-8, 1e-10, 1e-12, 1e-14}) {
    dp.observables.background.integration =
        numerics::IntegrationPolicy{tolerance, tolerance, 100000, 30};
    dp.observables.background.maximum_callbacks = 2000000;
    auto batch = prepared.evaluate(models, dp);
    if (batch.status != statistics::DensityStatus::finite ||
        batch.numerical_status != numerics::Status::ok ||
        batch.slots.size() != 11)
      throw std::runtime_error("original model count");
    for (std::size_t i = 0; i < batch.slots.size(); ++i) {
      const auto &slot = batch.slots[i];
      const auto available_ok = [](const bao::OutputState &state) {
        return state.availability == cosmology::Availability::available &&
               state.status == cosmology::Status::ok &&
               state.numerical_status == numerics::Status::ok;
      };
      if (!available_ok(slot.predictions_state) ||
          !available_ok(slot.residuals_state) ||
          !available_ok(slot.density_state) || !slot.result ||
          slot.numerical_status != numerics::Status::ok ||
          slot.result->density.numerical_status != numerics::Status::ok ||
          slot.result->density.status != statistics::DensityStatus::finite ||
          !std::isfinite(slot.result->density.log_value) ||
          !std::isfinite(slot.result->quadratic) ||
          !std::isfinite(slot.result->log_determinant) ||
          !std::isfinite(slot.result->normalization) ||
          slot.predictions.size() != 13 || slot.residuals.size() != 13)
        throw std::runtime_error("invalid original candidate payload");
      for (double value : slot.predictions)
        if (!std::isfinite(value))
          throw std::runtime_error("nonfinite prediction");
      for (double value : slot.residuals)
        if (!std::isfinite(value))
          throw std::runtime_error("nonfinite residual");
      std::printf(
          "BAO_ORIGINAL %.17g %zu %u %u %.17g %.17g %.17g %.17g %.17g %zu",
          tolerance, i, static_cast<unsigned>(slot.result->density.status),
          static_cast<unsigned>(slot.numerical_status),
          slot.result->density.log_value, slot.result->quadratic,
          slot.result->log_determinant, slot.result->normalization,
          slot.result->estimated_forward_sensitivity, slot.work.callbacks);
      for (double prediction : slot.predictions)
        std::printf(" %.17g", prediction);
      std::printf("\n");
      if (tolerance == 1e-14) {
        const auto &expected = test_reference::bao_points[i];
        if (!(std::abs(slot.result->density.log_value - expected.log_density) <=
              1e-8))
          throw std::runtime_error("frozen original normalized density budget");
        if (slot.predictions.size() != expected.prediction.size())
          throw std::runtime_error("original frozen prediction shape");
        for (std::size_t j = 0; j < expected.prediction.size(); ++j)
          if (!(std::abs(slot.predictions[j] - expected.prediction[j]) <=
                2e-12 + 2e-10 * std::abs(expected.prediction[j])))
            throw std::runtime_error("frozen original observable budget");
      }
      if (slot.result->density.status != statistics::DensityStatus::finite)
        throw std::runtime_error("original fixedpoint failure retained");
    }
    if (tolerance == 1e-12)
      for (const auto &slot : batch.slots)
        previous_scores.push_back(slot.result->density.log_value);
    if (tolerance == 1e-14)
      fine = std::move(batch);
  }
  original_reference(input, models, fine, previous_scores);
  std::printf("{\"suite\":\"current_BAO_original11\",\"models\":11,\"n\":13,"
              "\"passed\":true}\n");
  return 0;
}
int main(int argc, char **argv) {
  if (argc == 4 && std::string_view(argv[1]) == "--verified-original-assets")
    return original(argv[2], argv[3]);
  throw std::runtime_error(
      "optional external exactSHA guard required: --verified-original-assets "
      "MEAN COV; C++ does not hash inputs");
}
