// Frozen before implementation: coefficient/prediction 2e-11 abs+relative,
// H0 2e-10 abs+relative, quadratic 1e-9 absolute. Synthetic data only.
#include "irred/background.hpp"
#include "irred/calibration_ladder.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred::calibration;
using irred::statistics::DensityStatus;
namespace {
int checks = 0;
void check(bool x, const char *why) {
  ++checks;
  if (!x)
    throw std::runtime_error(why);
}
void near(double x, double y, const char *why) {
  check(std::abs(x - y) <= 2e-11 + 2e-11 * std::abs(y), why);
}
Model fixture() {
  Model m;
  m.ordered_host_ids = {"A", "B", "C"};
  m.magnitude_convention = "synthetic common photometric magnitude convention";
  m.metallicity_coordinate_identity = "synthetic log metal abundance dex";
  m.distance_shape_identity =
      "synthetic fixed FLRW Href70 supplied luminosity modulus";
  m.calibration_identity = "synthetic one delta measurement";
  m.dependence_identity = "known synthetic rank-one unrelated covariance";
  m.conditional_covariance_identity =
      "C=.01 I+.002 uuT conditional delta; u independent instrument mode";
  for (const auto *h : {"A", "B"})
    m.rows.push_back({RowKind::anchor_modulus, std::string("anchor/") + h, h,
                      "", 0, 0, 0, 0});
  for (int h = 0; h < 3; ++h)
    for (int k = 0; k < 3; ++k) {
      const double p[] = {.6, 1, 1.5}, z[] = {-.2, .3, -.1};
      const auto event = "cep/" + std::to_string(h) + "/" + std::to_string(k);
      m.rows.push_back({RowKind::cepheid, event, m.ordered_host_ids[h], event,
                        (k - 1) * .5, p[k], z[k], 0});
    }
  m.rows.push_back(
      {RowKind::calibrator_supernova, "cal/A", "A", "SN/A", 1, 0, 0, 0});
  m.rows.push_back(
      {RowKind::calibrator_supernova, "cal/C", "C", "SN/C", -.5, 0, 0, 0});
  m.rows.push_back(
      {RowKind::hubble_supernova, "hf/1", "", "SN/H1", .25, 0, 0, 36});
  m.rows.push_back(
      {RowKind::hubble_supernova, "hf/2", "", "SN/H2", 1.5, 0, 0, 38});
  m.rows.push_back({RowKind::calibration_measurement, "delta/measurement", "",
                    "", 1, 0, 0, 0});
  return m;
}
std::vector<double> truth() { return {31, 32, 33, -5, -3, .2, -19, .5, .03}; }
std::vector<std::string> rows(const Model &m) {
  std::vector<std::string> ids;
  for (const auto &r : m.rows)
    ids.push_back(r.row_id);
  return ids;
}
// Independent direct transcription, never obtains X or predictions from core.
std::vector<double> observations(const Model &m, const std::vector<double> &b) {
  std::vector<double> y;
  for (const auto &r : m.rows) {
    double mu = 0;
    if (!r.host_id.empty())
      mu = b[static_cast<std::size_t>(std::find(m.ordered_host_ids.begin(),
                                                m.ordered_host_ids.end(),
                                                r.host_id) -
                                      m.ordered_host_ids.begin())];
    switch (r.kind) {
    case RowKind::anchor_modulus:
      y.push_back(mu);
      break;
    case RowKind::cepheid:
      y.push_back(mu + b[3] + b[4] * (r.log10_period_days - 1) +
                  b[5] * (r.metallicity_dex - m.metallicity_reference_dex) +
                  r.calibration_response * b[8]);
      break;
    case RowKind::calibrator_supernova:
      y.push_back(mu + b[6] + r.calibration_response * b[8]);
      break;
    case RowKind::hubble_supernova:
      y.push_back(b[6] + r.reference_modulus_mag - b[7] +
                  r.calibration_response * b[8]);
      break;
    case RowKind::calibration_measurement:
      y.push_back(b[8]);
      break;
    }
  }
  return y;
}
irred::statistics::Gaussian gaussian(const Model &m, bool correlated = true) {
  const auto n = m.rows.size();
  std::vector<double> c(n * n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      auto row_mode = [&](std::size_t row) {
        unsigned value = 0;
        for (unsigned char ch : m.rows[row].row_id)
          value += ch;
        return static_cast<double>(static_cast<int>(value % 3) - 1);
      };
      const double u_i = row_mode(i), u_j = row_mode(j);
      c[i * n + j] = (i == j ? .01 : 0) + (correlated ? .002 * u_i * u_j : 0);
    }
  irred::statistics::Metadata md;
  md.ordered_ids = rows(m);
  md.measure = "product d(mag)";
  md.source_semantics = "synthetic controls";
  md.table_identity = "hand constructed ladder v1";
  md.uncertainty_identity = m.conditional_covariance_identity;
  md.ordering_provenance = "exact synthetic row order";
  md.calibration_provenance = m.calibration_identity;
  md.dependence_provenance = m.dependence_identity;
  return irred::statistics::prepare_gaussian(
      c, irred::statistics::MatrixKind::covariance, md, n * n, 1e-10);
}
void recovered(const Recovery &r, const std::vector<double> &b,
               double href = 70) {
  check(r.relative_fit.status == DensityStatus::finite,
        "relative recovery admitted");
  for (std::size_t j = 0; j < b.size(); ++j)
    near(r.relative_fit.coefficients[j], b[j], "injected coordinate recovery");
  const auto h = href * std::pow(10, b[7] / 5);
  check(std::abs(r.h0_km_s_Mpc - h) <= 2e-10 + 2e-10 * std::abs(h),
        "injected H0 projection");
  check(std::abs(r.relative_fit.quadratic) <= 1e-9,
        "noiseless quadratic budget");
}
} // namespace
int main() {
  try {
    auto m = fixture();
    auto b = truth();
    auto y = observations(m, b);
    auto ids = rows(m);
    auto g = gaussian(m);
    auto ladder = Ladder::prepare(std::move(g), m);
    check(ladder.status() == DensityStatus::finite &&
              g.status() == DensityStatus::invalid_input,
          "prepare consumes covariance once");
    check(ladder.rank() ==
              irred::statistics::DesignRank::full_within_conditioning_contract,
          "identified design");
    recovered(ladder.fit(y, ids), b);
    auto prediction = ladder.predict(b, parameter_ids(m));
    check(prediction.status == DensityStatus::finite, "forward admitted");
    for (std::size_t i = 0; i < y.size(); ++i)
      near(prediction.values_mag[i], y[i], "independent forward rows");
    auto again = ladder.fit(y, ids);
    recovered(again, b);
    check(ladder.design_metadata().shared_nuisance_ids ==
              std::vector<std::string>{"delta_shared_calibration"},
          "one shared coordinate");
    // Coordinate transformation preserves every forward value and physical H0.
    auto changed = m;
    changed.metallicity_reference_dex = .4;
    changed.h_reference_km_s_Mpc = 100;
    changed.distance_shape_identity = "same fixed shape at Href100";
    const auto shift = 5 * std::log10(100. / 70.);
    for (auto &r : changed.rows)
      if (r.kind == RowKind::hubble_supernova)
        r.reference_modulus_mag -= shift;
    auto bc = b;
    bc[3] += .4 * b[5];
    bc[7] -= shift;
    auto gc = gaussian(changed);
    auto lc = Ladder::prepare(std::move(gc), changed);
    recovered(lc.fit(y, rows(changed)), bc, 100);
    auto pc = lc.predict(bc, parameter_ids(changed));
    for (std::size_t i = 0; i < y.size(); ++i)
      near(pc.values_mag[i], y[i], "reference coordinate invariance");
    // Row order can change only along with observations and covariance axes.
    auto perm = m;
    std::reverse(perm.rows.begin(), perm.rows.end());
    auto gp = gaussian(perm);
    auto lp = Ladder::prepare(std::move(gp), perm);
    recovered(lp.fit(observations(perm, b), rows(perm)), b);
    auto hp = m;
    std::swap(hp.ordered_host_ids[0], hp.ordered_host_ids[2]);
    auto hb = b;
    std::swap(hb[0], hb[2]);
    auto gh = gaussian(hp);
    auto lh = Ladder::prepare(std::move(gh), hp);
    recovered(lh.fit(y, rows(hp)), hb);
    // Remove a Cepheid, refit remaining rows, then predict that unfit row.
    auto held = m;
    held.rows = {m.rows[5]};
    auto train = m;
    train.rows.erase(train.rows.begin() + 5);
    auto gt = gaussian(train);
    auto lt = Ladder::prepare(std::move(gt), train);
    const auto rt = lt.fit(observations(train, b), rows(train));
    recovered(rt, b);
    auto ph = predict(held, rt.relative_fit.coefficients, parameter_ids(held));
    check(ph.status == DensityStatus::finite && ph.values_mag.size() == 1,
          "independent heldout admitted");
    near(ph.values_mag[0], observations(held, b)[0],
         "remove refit heldout Cepheid");
    // Explicit existing background geometry produces one reference distance.
    using namespace irred::cosmology;
    auto background = prepare(ConstantQ{0}, FlatFLRW{});
    std::vector<Request> query;
    query.emplace_back(
        .1, static_cast<std::uint32_t>(Observable::flat_distances_volume),
        Observer{.1, Convention::geometric_same_redshift}, PhysicalScale{70});
    EvaluationPolicy bp;
    bp.integration =
        irred::numerics::IntegrationPolicy{1e-13, 1e-13, 10000, 30};
    bp.maximum_queries = 1;
    bp.maximum_callbacks = 10000;
    bp.maximum_segment_visits = 100;
    bp.maximum_native_bytes = 1 << 20;
    const auto geometry = background.evaluate(query, bp);
    check(geometry.slots.size() == 1 &&
              geometry.slots[0].physical.value.has_value(),
          "reference background distance");
    auto geometry_model = m;
    geometry_model.rows[13].reference_modulus_mag =
        5 * std::log10(geometry.slots[0].physical.value->luminosity_mpc) + 25;
    geometry_model.distance_shape_identity =
        "P01 flat constant-q0; geometric-same-redshift z=.1; Href70; "
        "luminosity Mpc";
    auto gb = gaussian(geometry_model);
    auto lb = Ladder::prepare(std::move(gb), geometry_model);
    recovered(lb.fit(observations(geometry_model, b), rows(geometry_model)), b);
    const auto analytic_mu =
        5 * std::log10((299792.458 / 70) * std::log1p(.1) * 1.1) + 25;
    near(geometry_model.rows[13].reference_modulus_mag, analytic_mu,
         "analytic constant-q reference shape");
    // Explicit absence of calibration observation and identical responses:
    // delta is M_Cep+M_SN and unidentified. No parameter is dropped or
    // regularized.
    auto deg = m;
    deg.rows.pop_back();
    for (auto &r : deg.rows)
      if (r.kind != RowKind::anchor_modulus)
        r.calibration_response = 1;
    auto gd = gaussian(deg);
    auto ld = Ladder::prepare(std::move(gd), deg);
    check(ld.status() != DensityStatus::finite &&
              gd.status() == DensityStatus::finite,
          "calibration degeneracy preserves source");
    auto unused = m;
    unused.ordered_host_ids.push_back("unused host");
    auto gu = gaussian(unused);
    auto lu = Ladder::prepare(std::move(gu), unused);
    check(lu.status() != DensityStatus::finite &&
              gu.status() == DensityStatus::finite,
          "unconnected host rank failure preserves source");
    auto near_rank = m;
    for (auto &r : near_rank.rows)
      if (r.kind == RowKind::cepheid)
        r.metallicity_dex = r.log10_period_days - 1;
    near_rank.rows[3].metallicity_dex += 1e-9;
    auto gn = gaussian(near_rank);
    auto ln = Ladder::prepare(std::move(gn), near_rank);
    check(ln.status() != DensityStatus::finite &&
              gn.status() == DensityStatus::finite,
          "near period metallicity rank failure preserves source");
    auto independent = m;
    independent.rows.pop_back();
    auto gi = gaussian(independent);
    auto li = Ladder::prepare(std::move(gi), independent);
    recovered(li.fit(observations(independent, b), rows(independent)), b);
    for (int failure = 0; failure < 11; ++failure) {
      auto bad = m;
      switch (failure) {
      case 0:
        bad.rows[1].row_id = bad.rows[0].row_id;
        break;
      case 1:
        bad.rows[3].event_id = bad.rows[2].event_id;
        break;
      case 2:
        bad.rows[2].host_id = "unknown";
        break;
      case 3:
        bad.rows[0].reference_modulus_mag = 1;
        break;
      case 4:
        bad.rows[2].metallicity_dex = NAN;
        break;
      case 5:
        bad.ordered_host_ids[1] = bad.ordered_host_ids[0];
        break;
      case 6:
        bad.conditional_covariance_identity.clear();
        break;
      case 7:
        bad.rows[0].kind = static_cast<RowKind>(999);
        break;
      case 8:
        bad.rows[0].calibration_response = 1;
        break;
      case 9:
        bad.magnitude_convention.clear();
        break;
      case 10:
        bad.metallicity_coordinate_identity.clear();
        break;
      }
      auto keep = gaussian(m);
      auto reject = Ladder::prepare(std::move(keep), bad);
      check(reject.status() != DensityStatus::finite &&
                keep.status() == DensityStatus::finite,
            "invalid model preserves covariance");
    }
    for (int mismatch = 0; mismatch < 4; ++mismatch) {
      auto mismatch_model = m;
      if (mismatch == 0)
        mismatch_model.calibration_identity = "conflicting calibration";
      if (mismatch == 1)
        mismatch_model.dependence_identity = "conflicting dependence";
      if (mismatch == 2)
        mismatch_model.conditional_covariance_identity =
            "conflicting covariance";
      auto original = gaussian(m);
      if (mismatch == 3) {
        auto metadata = original.metadata();
        metadata.source_semantics = "real observations";
        std::vector<double> covariance(original.covariance().begin(),
                                       original.covariance().end());
        original = irred::statistics::prepare_gaussian(
            covariance, irred::statistics::MatrixKind::covariance, metadata,
            covariance.size(), 1e-10);
      }
      auto incompatible = Ladder::prepare(std::move(original), mismatch_model);
      check(incompatible.status() == DensityStatus::incompatible_metadata &&
                original.status() == DensityStatus::finite,
            "conflicting source metadata preserves covariance");
    }
    auto wrong = ids;
    std::swap(wrong[0], wrong[1]);
    check(ladder.fit(y, wrong).relative_fit.status ==
              DensityStatus::incompatible_metadata,
          "ordered fit IDs");
    auto wrongp = parameter_ids(m);
    std::swap(wrongp[0], wrongp[1]);
    check(ladder.predict(b, wrongp).status ==
              DensityStatus::incompatible_metadata,
          "ordered coefficient IDs");
    auto invalid_y = y;
    invalid_y[0] = INFINITY;
    check(ladder.fit(invalid_y, ids).relative_fit.status ==
              DensityStatus::invalid_input,
          "nonfinite measurement");
    auto tiny = irred::statistics::DesignPolicy{};
    tiny.maximum_payload_bytes = 1;
    auto keep = gaussian(m);
    auto limited = Ladder::prepare(std::move(keep), m, tiny);
    check(limited.numerical_status() == irred::numerics::Status::work_limit &&
              keep.status() == DensityStatus::finite,
          "wrapper preparation budget");
    check(ladder.fit(y, ids, tiny).relative_fit.coefficients.empty(),
          "wrapper fit budget no outputs");
    check(ladder.predict(b, parameter_ids(m), tiny).values_mag.empty(),
          "wrapper prediction budget no outputs");
    for (const double eta : {2000., -2000.}) {
      auto enormous = b;
      enormous[7] = eta;
      const auto fail = ladder.fit(observations(m, enormous), ids);
      check(fail.relative_fit.status == DensityStatus::numerical_failure &&
                fail.relative_fit.coefficients.empty() && fail.h0_km_s_Mpc == 0,
            "H0 overflow underflow fail entire recovery");
    }
    auto moved = std::move(ladder);
    check(ladder.status() == DensityStatus::invalid_input,
          "moved source invalid");
    recovered(moved.fit(y, ids), b);
    auto *self = &moved;
    moved = std::move(*self);
    recovered(moved.fit(y, ids), b);
    check(ladder.model().rows.empty() &&
              ladder.model().ordered_host_ids.empty(),
          "successful moved source releases model");
    Ladder assigned;
    assigned = std::move(moved);
    check(moved.status() == DensityStatus::invalid_input &&
              moved.model().rows.empty(),
          "move assignment invalidates source and clears model");
    recovered(assigned.fit(y, ids), b);
    auto *assigned_self = &assigned;
    assigned = std::move(*assigned_self);
    recovered(assigned.fit(y, ids), b);
    Ladder failed_moved(std::move(limited));
    check(limited.status() == DensityStatus::invalid_input &&
              limited.numerical_status() ==
                  irred::numerics::Status::invalid_input &&
              failed_moved.numerical_status() ==
                  irred::numerics::Status::work_limit,
          "failed wrapper move clears source fallback status");
    Ladder failed_assigned;
    failed_assigned = std::move(failed_moved);
    check(failed_moved.status() == DensityStatus::invalid_input &&
              failed_assigned.numerical_status() ==
                  irred::numerics::Status::work_limit,
          "failed wrapper move assignment clears fallback status");
    // Saturated square control: eta=y_cal-y_anchorA-y_HF+mu_ref-.75*y_delta.
    // Independent diagonal .01 covariance gives Var(eta)=.01*(3+.75^2).
    auto square = m;
    square.ordered_host_ids = {"A", "B"};
    square.rows.clear();
    for (const auto index : {0, 1, 2, 3, 4, 11, 13, 15})
      square.rows.push_back(m.rows[index]);
    square.conditional_covariance_identity =
        "synthetic diagonal .01 I conditional delta";
    auto square_source = gaussian(square, false);
    auto square_ladder = Ladder::prepare(std::move(square_source), square);
    check(square_ladder.status() == DensityStatus::finite,
          "square eta sampling control");
    const std::vector<double> probabilities{.8413447460685429, .5,
                                            .15865525393145705};
    const auto law = square_ladder.h0_estimator_law(
        .5, "supplied synthetic generating eta=.5", probabilities);
    check(law.status == DensityStatus::finite,
          "explicit conditional H0 sampling law");
    check(std::abs(law.eta_estimator.variance - .035625) <=
              2e-12 + 2e-12 * .035625,
          "independent linear-noise eta variance");
    const auto a = std::log(10.) / 5;
    const auto nominal = 70 * std::pow(10., .1);
    const auto spread = a * std::sqrt(.035625);
    auto h_near = [&](double actual, double expected, const char *why) {
      check(std::abs(actual - expected) <= 2e-10 + 2e-10 * std::abs(expected),
            why);
    };
    h_near(law.nominal_h0_km_s_Mpc, nominal, "H0 nominal median projection");
    h_near(law.sampling_expectation_h0_km_s_Mpc,
           nominal * std::exp(a * a * .035625 / 2),
           "H0 nonlinear sampling expectation");
    check(law.sampling_expectation_h0_km_s_Mpc > law.nominal_h0_km_s_Mpc,
          "expectation distinct from nominal fit projection");
    h_near(law.h0_quantiles_km_s_Mpc[0], nominal * std::exp(spread),
           "near-one-sigma upper quantile");
    h_near(law.h0_quantiles_km_s_Mpc[1], nominal, "median quantile");
    h_near(law.h0_quantiles_km_s_Mpc[2], nominal * std::exp(-spread),
           "near-one-sigma lower quantile");
    check(law.probabilities == probabilities &&
              !law.eta_mean_identity.empty() && law.assumed_eta_mean == .5 &&
              law.eta_estimator.metadata.output_unit == "mag" &&
              law.maximum_relative_tail_probability_error <= 1e-12 &&
              law.maximum_quantile_log_error_estimate <= 1e-10 &&
              law.normal_cdf_evaluations <= 98 * probabilities.size(),
          "law exact probability order provenance and bounded diagnostics");
    const auto endpoints = square_ladder.h0_estimator_law(
        .5, "same synthetic mean", std::vector<double>{1e-12, 1 - 1e-12});
    check(endpoints.status == DensityStatus::finite &&
              endpoints.h0_quantiles_km_s_Mpc[0] < nominal &&
              endpoints.h0_quantiles_km_s_Mpc[1] > nominal,
          "bounded quantile endpoints");
    for (const double probability :
         {0., 1., std::numeric_limits<double>::quiet_NaN()})
      check(square_ladder
                    .h0_estimator_law(.5, "mean",
                                      std::vector<double>{probability})
                    .status == DensityStatus::invalid_input,
            "unsupported probability rejected");
    check(square_ladder.h0_estimator_law(.5, "", probabilities).status ==
              DensityStatus::invalid_input,
          "mean provenance required");
    check(square_ladder.h0_estimator_law(.5, "mean", probabilities, tiny)
              .h0_quantiles_km_s_Mpc.empty(),
          "sampling law byte quota no output");
    auto count_quota = irred::statistics::DesignPolicy{};
    count_quota.maximum_elements = probabilities.size() - 1;
    check(square_ladder.h0_estimator_law(.5, "mean", probabilities, count_quota)
                  .numerical_status == irred::numerics::Status::work_limit,
          "sampling law element quota");
    for (const double eta : {2000., -2000.}) {
      const auto failed_law =
          square_ladder.h0_estimator_law(eta, "extreme mean", probabilities);
      check(failed_law.status == DensityStatus::numerical_failure &&
                failed_law.h0_quantiles_km_s_Mpc.empty() &&
                failed_law.nominal_h0_km_s_Mpc == 0 &&
                failed_law.eta_estimator.metadata.weight_units.empty(),
            "unrepresentable law projection withholds all calculation payload");
    }
    check(ladder.h0_estimator_law(.5, "mean", probabilities).status ==
              DensityStatus::invalid_input,
          "moved sampling owner invalid");
    const auto round = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    check(assigned.predict(b, parameter_ids(m)).status ==
              DensityStatus::unsupported_domain,
          "rounding contract prediction");
    std::fesetround(round);
    std::printf("calibration ladder owner: %d checks passed; bounded synthetic "
                "relative controls\n",
                checks);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL: %s (%d checks)\n", e.what(), checks);
    return 1;
  }
}
