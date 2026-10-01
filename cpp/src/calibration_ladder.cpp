#include "irred/calibration_ladder.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <set>
#include <utility>

namespace irred::calibration {
namespace {
using statistics::DensityStatus;
bool arithmetic_supported() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384;
}
bool representable(double v) {
  return std::isfinite(v) && (v == 0 || std::isnormal(v));
}
bool add(std::size_t &a, std::size_t b) {
  if (b > SIZE_MAX - a)
    return false;
  a += b;
  return true;
}
bool product(std::size_t a, std::size_t b, std::size_t &out) {
  if (a && b > SIZE_MAX / a)
    return false;
  out = a * b;
  return true;
}
std::optional<std::size_t> payload(const Model &m) {
  std::size_t n = sizeof(Model), q;
  if (!product(m.rows.capacity(), sizeof(Row), q) || !add(n, q) ||
      !product(m.ordered_host_ids.capacity(), sizeof(std::string), q) ||
      !add(n, q))
    return {};
  auto string_charge = [&](const std::string &s) {
    return add(n, s.capacity() + 1);
  };
  for (const auto &s : m.ordered_host_ids)
    if (!string_charge(s))
      return {};
  for (const auto &r : m.rows)
    if (!string_charge(r.row_id) || !string_charge(r.host_id) ||
        !string_charge(r.event_id))
      return {};
  if (!string_charge(m.magnitude_convention) ||
      !string_charge(m.metallicity_coordinate_identity) ||
      !string_charge(m.distance_shape_identity) ||
      !string_charge(m.calibration_identity) ||
      !string_charge(m.dependence_identity) ||
      !string_charge(m.conditional_covariance_identity))
    return {};
  return n;
}
bool valid(const Model &m, bool recovery) {
  if (m.ordered_host_ids.size() < 2 || m.rows.empty() ||
      !representable(m.metallicity_reference_dex) ||
      !std::isnormal(m.h_reference_km_s_Mpc) || m.h_reference_km_s_Mpc <= 0 ||
      m.magnitude_convention.empty() ||
      m.metallicity_coordinate_identity.empty() ||
      m.distance_shape_identity.empty() || m.calibration_identity.empty() ||
      m.dependence_identity.empty() ||
      m.conditional_covariance_identity.empty())
    return false;
  std::set<std::string> hosts, rows, events, anchors;
  for (const auto &s : m.ordered_host_ids)
    if (s.empty() || !hosts.insert(s).second)
      return false;
  std::size_t cep = 0, cal = 0, hf = 0, measurements = 0;
  for (const auto &r : m.rows) {
    if (r.row_id.empty() || !rows.insert(r.row_id).second ||
        !representable(r.calibration_response) ||
        !representable(r.log10_period_days) ||
        !representable(r.metallicity_dex) ||
        !representable(r.reference_modulus_mag))
      return false;
    switch (r.kind) {
    case RowKind::anchor_modulus:
      if (!hosts.contains(r.host_id) || !r.event_id.empty() ||
          r.calibration_response != 0 || r.log10_period_days != 0 ||
          r.metallicity_dex != 0 || r.reference_modulus_mag != 0)
        return false;
      anchors.insert(r.host_id);
      break;
    case RowKind::cepheid:
      ++cep;
      if (!hosts.contains(r.host_id) || r.reference_modulus_mag != 0 ||
          r.event_id.empty() || !events.insert(r.event_id).second)
        return false;
      break;
    case RowKind::calibrator_supernova:
      ++cal;
      if (!hosts.contains(r.host_id) || r.log10_period_days != 0 ||
          r.metallicity_dex != 0 || r.reference_modulus_mag != 0 ||
          r.event_id.empty() || !events.insert(r.event_id).second)
        return false;
      break;
    case RowKind::hubble_supernova:
      ++hf;
      if (!r.host_id.empty() || r.log10_period_days != 0 ||
          r.metallicity_dex != 0 || r.event_id.empty() ||
          !events.insert(r.event_id).second)
        return false;
      break;
    case RowKind::calibration_measurement:
      if (++measurements > 1 || !r.host_id.empty() || !r.event_id.empty() ||
          r.calibration_response != 1 || r.log10_period_days != 0 ||
          r.metallicity_dex != 0 || r.reference_modulus_mag != 0)
        return false;
      break;
    default:
      return false;
    }
  }
  return !recovery || (anchors.size() >= 2 && cep && cal && hf);
}
std::size_t eta_index(const Model &m) { return m.ordered_host_ids.size() + 4; }
long double eta_log_h0_scale() { return std::log(10.L) / 5; }
long double log_h0_projection(const Model &m, double eta) {
  return std::log(static_cast<long double>(m.h_reference_km_s_Mpc)) +
         eta_log_h0_scale() * eta;
}
double row_offset(const Row &r) {
  return r.kind == RowKind::hubble_supernova ? r.reference_modulus_mag : 0;
}
// Single owner of each row's linear equation, used by prediction and
// preparation.
bool equation(const Model &m, const Row &r, std::span<double> x,
              double &offset) {
  std::fill(x.begin(), x.end(), 0);
  const auto h = m.ordered_host_ids.size();
  if (!r.host_id.empty()) {
    const auto it = std::find(m.ordered_host_ids.begin(),
                              m.ordered_host_ids.end(), r.host_id);
    if (it == m.ordered_host_ids.end())
      return false;
    x[static_cast<std::size_t>(it - m.ordered_host_ids.begin())] = 1;
  }
  offset = row_offset(r);
  switch (r.kind) {
  case RowKind::anchor_modulus:
    break;
  case RowKind::cepheid:
    x[h] = 1;
    x[h + 1] = r.log10_period_days - 1;
    x[h + 2] = r.metallicity_dex - m.metallicity_reference_dex;
    break;
  case RowKind::calibrator_supernova:
    x[h + 3] = 1;
    break;
  case RowKind::hubble_supernova:
    x[h + 3] = 1;
    x[h + 4] = -1;
    break;
  case RowKind::calibration_measurement:
    break;
  default:
    return false;
  }
  x[h + 5] = r.calibration_response;
  return std::all_of(x.begin(), x.end(), representable);
}
bool model_work(const Model &m, statistics::DesignPolicy p,
                std::size_t &charge) {
  if (m.ordered_host_ids.size() > SIZE_MAX - 6)
    return false;
  std::size_t cells, bytes;
  if (!product(m.rows.size(), m.ordered_host_ids.size() + 6, cells) ||
      cells > p.maximum_elements || !product(cells, sizeof(double), bytes))
    return false;
  const auto own = payload(m);
  if (!own)
    return false;
  charge = *own;
  // Local equation vector, output/residual vector and transient IDs/validation
  // sets: conservative string/object storage, excluding allocator bookkeeping.
  if (!add(charge, *own) || !add(charge, bytes) || !add(charge, bytes) ||
      charge > p.maximum_payload_bytes)
    return false;
  return true;
}
Recovery failed(DensityStatus s, numerics::Status ns) {
  Recovery v;
  v.relative_fit.status = s;
  v.relative_fit.numerical_status = ns;
  return v;
}
} // namespace
std::vector<std::string> parameter_ids(const Model &m) {
  std::vector<std::string> ids;
  ids.reserve(m.ordered_host_ids.size() + 6);
  for (const auto &s : m.ordered_host_ids)
    ids.push_back("mu_host/" + s);
  for (const auto *s : {"M_Cep", "b_log10_period", "gamma_metallicity", "M_SN",
                        "eta_5log10_H0_over_Href", "delta_shared_calibration"})
    ids.emplace_back(s);
  return ids;
}
Prediction predict(const Model &m, std::span<const double> b,
                   std::span<const std::string> ids,
                   statistics::DesignPolicy policy) {
  Prediction out;
  if (!arithmetic_supported()) {
    out.status = DensityStatus::unsupported_domain;
    return out;
  }
  std::size_t charge;
  if (!model_work(m, policy, charge)) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = numerics::Status::work_limit;
    return out;
  }
  if (!valid(m, false) || b.size() != m.ordered_host_ids.size() + 6 ||
      !std::all_of(b.begin(), b.end(), representable))
    return out;
  const auto expected = parameter_ids(m);
  if (ids.size() != expected.size() ||
      !std::equal(ids.begin(), ids.end(), expected.begin())) {
    out.status = DensityStatus::incompatible_metadata;
    return out;
  }
  std::vector<double> x(b.size()), values;
  values.reserve(m.rows.size());
  for (const auto &r : m.rows) {
    double offset;
    if (!equation(m, r, x, offset))
      return out;
    long double v = offset;
    for (std::size_t j = 0; j < b.size(); ++j)
      v += static_cast<long double>(x[j]) * b[j];
    const auto d = static_cast<double>(v);
    if (!representable(d) || (v != 0 && d == 0)) {
      out.status = DensityStatus::numerical_failure;
      out.numerical_status = numerics::Status::outside_domain;
      return out;
    }
    values.push_back(d);
  }
  out.values_mag = std::move(values);
  out.status = DensityStatus::finite;
  out.numerical_status = numerics::Status::ok;
  return out;
}
Ladder::Ladder(Ladder &&other) noexcept
    : preparation_status_(std::exchange(other.preparation_status_,
                                        DensityStatus::invalid_input)),
      preparation_numerical_status_(
          std::exchange(other.preparation_numerical_status_,
                        numerics::Status::invalid_input)),
      model_(std::move(other.model_)), profile_(std::move(other.profile_)) {
  other.model_ = Model{};
}
Ladder &Ladder::operator=(Ladder &&other) noexcept {
  if (this != &other) {
    preparation_status_ =
        std::exchange(other.preparation_status_, DensityStatus::invalid_input);
    preparation_numerical_status_ = std::exchange(
        other.preparation_numerical_status_, numerics::Status::invalid_input);
    model_ = std::move(other.model_);
    profile_ = std::move(other.profile_);
    other.model_ = Model{};
  }
  return *this;
}
Ladder Ladder::prepare(statistics::Gaussian &&g, Model m,
                       statistics::DesignPolicy policy) {
  Ladder out;
  std::size_t charge;
  // Invalid wrapper preparation returns an invalid owner without consuming g.
  if (!arithmetic_supported()) {
    out.preparation_status_ = DensityStatus::unsupported_domain;
    return out;
  }
  if (!model_work(m, policy, charge)) {
    out.preparation_status_ = DensityStatus::numerical_failure;
    out.preparation_numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  if (!valid(m, true))
    return out;
  const auto &gm = g.metadata();
  if (gm.source_semantics != "synthetic controls" ||
      gm.measure != "product d(mag)" ||
      gm.calibration_provenance != m.calibration_identity ||
      gm.dependence_provenance != m.dependence_identity ||
      gm.uncertainty_identity != m.conditional_covariance_identity ||
      gm.table_identity.empty() || gm.ordering_provenance.empty()) {
    out.preparation_status_ = DensityStatus::incompatible_metadata;
    return out;
  }

  auto ids = parameter_ids(m);
  statistics::DesignMetadata md;
  md.ordered_parameter_ids = ids;
  md.parameter_units.assign(ids.size(), "mag");
  md.parameter_units[m.ordered_host_ids.size() + 1] = "mag per log10(day)";
  md.parameter_units[m.ordered_host_ids.size() + 2] = "mag per dex";
  md.shared_nuisance_ids = {ids.back()};
  md.residual_unit = "mag";
  md.design_identity = "empirical-synthetic-absolute-ladder/supplied-shape/v1";
  md.dependence_identity = m.dependence_identity + "; conditional covariance=" +
                           m.conditional_covariance_identity;
  const auto bound =
      statistics::DesignProfile::preparation_payload_bound(g, ids.size(), md);
  if (!bound || !add(charge, *bound) || charge > policy.maximum_payload_bytes) {
    out.preparation_status_ = DensityStatus::numerical_failure;
    out.preparation_numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  std::vector<std::string> row_ids;
  row_ids.reserve(m.rows.size());
  std::vector<double> x(m.rows.size() * ids.size());
  for (std::size_t i = 0; i < m.rows.size(); ++i) {
    row_ids.push_back(m.rows[i].row_id);
    double offset;
    if (!equation(m, m.rows[i],
                  std::span<double>(x).subspan(i * ids.size(), ids.size()),
                  offset))
      return out;
  }
  out.profile_ = statistics::DesignProfile::prepare(std::move(g), x, row_ids,
                                                    std::move(md), policy);
  if (out.profile_.status() == DensityStatus::finite)
    out.model_ = std::move(m);
  return out;
}
Recovery Ladder::fit(std::span<const double> y,
                     std::span<const std::string> ids,
                     statistics::DesignPolicy policy) const {
  if (status() != DensityStatus::finite)
    return failed(status(), numerical_status());
  const auto bound = profile_.evaluation_payload_bound();
  std::size_t charge = 0, extra;
  if (!bound || !product(y.size(), sizeof(double), extra) ||
      !add(charge, *bound) || !add(charge, extra) ||
      charge > policy.maximum_payload_bytes)
    return failed(DensityStatus::numerical_failure,
                  numerics::Status::work_limit);
  if (y.size() != model_.rows.size() ||
      !std::all_of(y.begin(), y.end(), representable))
    return failed(DensityStatus::invalid_input,
                  numerics::Status::invalid_input);
  std::vector<double> residual(y.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    const auto offset = row_offset(model_.rows[i]);
    const long double r = static_cast<long double>(y[i]) - offset;
    residual[i] = static_cast<double>(r);
    if (!representable(residual[i]) || (r != 0 && residual[i] == 0))
      return failed(DensityStatus::numerical_failure,
                    numerics::Status::outside_domain);
  }
  Recovery out;
  out.relative_fit = profile_.evaluate(residual, ids, policy);
  if (out.relative_fit.status != DensityStatus::finite)
    return out;
  const long double exponent = log_h0_projection(
      model_, out.relative_fit.coefficients[eta_index(model_)]);
  const long double h0 = std::exp(exponent);
  out.h0_km_s_Mpc = static_cast<double>(h0);
  if (!std::isnormal(out.h0_km_s_Mpc) || !std::isfinite(h0))
    return failed(DensityStatus::numerical_failure,
                  numerics::Status::outside_domain);
  return out;
}
Prediction Ladder::predict(std::span<const double> b,
                           std::span<const std::string> ids,
                           statistics::DesignPolicy policy) const {
  if (status() != DensityStatus::finite) {
    Prediction p;
    p.status = status();
    p.numerical_status = numerical_status();
    return p;
  }
  return calibration::predict(model_, b, ids, policy);
}
H0EstimatorLaw Ladder::h0_estimator_law(double eta_mean,
                                        std::string mean_identity,
                                        std::span<const double> probabilities,
                                        statistics::DesignPolicy policy) const {
  H0EstimatorLaw out;
  auto fail = [](DensityStatus status, numerics::Status numerical_status) {
    H0EstimatorLaw result;
    result.status = status;
    result.numerical_status = numerical_status;
    return result;
  };
  if (!arithmetic_supported())
    return fail(DensityStatus::unsupported_domain,
                numerics::Status::outside_domain);
  if (status() != DensityStatus::finite)
    return fail(status(), numerical_status());
  if (!representable(eta_mean) || mean_identity.empty() ||
      probabilities.empty() ||
      !std::isfinite(policy.maximum_forward_sensitivity) ||
      policy.maximum_forward_sensitivity <= 0)
    return out;
  const auto &design = profile_.design_metadata();
  const auto p = design.ordered_parameter_ids.size();
  detail::PayloadAccounting bytes(sizeof(H0EstimatorLaw));
  bytes.string(mean_identity);
  bytes.add(probabilities.size(), 2 * sizeof(double));
  bytes.add(p, sizeof(double) + 2 * sizeof(long double) + sizeof(std::string));
  // Bound construction/result capacities of each declared weight unit, plus
  // fixed functional/unit identities. No model/covariance readback is needed.
  for (const auto &unit : design.parameter_units) {
    bytes.add(unit.size(), 2);
    bytes.add(64, 1);
  }
  bytes.add(256, 1);
  const auto bound = bytes.result();
  if (!bound || *bound > policy.maximum_payload_bytes ||
      p > policy.maximum_elements ||
      probabilities.size() > policy.maximum_elements ||
      probabilities.size() > std::numeric_limits<std::size_t>::max() / 98)
    return fail(DensityStatus::numerical_failure, numerics::Status::work_limit);
  if (std::any_of(probabilities.begin(), probabilities.end(),
                  [](double probability) {
                    return !std::isfinite(probability) || probability < 1e-12 ||
                           probability > 1 - 1e-12;
                  }))
    return out;
  statistics::LinearFunctionalMetadata contrast;
  contrast.functional_identity = "empirical-synthetic-ladder/eta-estimator/v1";
  contrast.output_unit = "mag";
  contrast.weight_units.reserve(p);
  for (const auto &unit : design.parameter_units)
    contrast.weight_units.push_back("mag / (" + unit + ")");
  std::vector<double> weights(p, 0);
  weights[eta_index(model_)] = 1;
  auto variance = profile_.estimator_variance(
      weights, design.ordered_parameter_ids, std::move(contrast), policy);
  if (variance.status != DensityStatus::finite)
    return fail(variance.status, variance.numerical_status);
  const long double a = eta_log_h0_scale();
  const long double log_nominal = log_h0_projection(model_, eta_mean);
  const long double spread =
      a * std::sqrt(static_cast<long double>(variance.variance));
  const long double log_expectation =
      log_nominal + a * a * variance.variance / 2;
  auto positive_output = [](long double value, double &rounded) {
    if (!std::isfinite(value) || !(value > 0) ||
        value > std::numeric_limits<double>::max())
      return false;
    rounded = static_cast<double>(value);
    return std::isnormal(rounded);
  };
  const auto nominal = std::exp(log_nominal),
             expectation = std::exp(log_expectation);
  if (!positive_output(nominal, out.nominal_h0_km_s_Mpc) ||
      !positive_output(expectation, out.sampling_expectation_h0_km_s_Mpc))
    return fail(DensityStatus::numerical_failure,
                numerics::Status::outside_domain);
  const auto epsilon = std::numeric_limits<long double>::epsilon();
  const auto expectation_error =
      a * a * variance.variance / 2 * variance.estimated_forward_sensitivity +
      8 * epsilon * (1 + std::abs(log_expectation)) +
      std::abs(expectation - out.sampling_expectation_h0_km_s_Mpc) /
          expectation;
  if (!std::isfinite(expectation_error) ||
      expectation_error > policy.maximum_forward_sensitivity)
    return fail(DensityStatus::numerical_failure,
                numerics::Status::conditioning_budget_exceeded);
  std::vector<double> quantiles;
  quantiles.reserve(probabilities.size());
  for (const double requested : probabilities) {
    const long double probability = requested;
    const long double tail = std::min(probability, 1 - probability);
    long double z = 0, achieved = .5L;
    if (requested != .5) {
      long double lower = 0, upper = 8;
      auto small_tail = [&](long double value) {
        ++out.normal_cdf_evaluations;
        return std::erfc(value / std::sqrt(2.L)) / 2;
      };
      if (small_tail(upper) > tail)
        return fail(DensityStatus::numerical_failure,
                    numerics::Status::outside_domain);
      for (unsigned iteration = 0; iteration < 96; ++iteration) {
        const auto middle = (lower + upper) / 2;
        if (middle == lower || middle == upper)
          break;
        if (small_tail(middle) > tail)
          lower = middle;
        else
          upper = middle;
      }
      z = (lower + upper) / 2;
      achieved = small_tail(z);
    }
    const auto tail_error = std::abs(achieved - tail) / tail;
    if (!std::isfinite(tail_error) || tail_error > 1e-12L)
      return fail(DensityStatus::numerical_failure,
                  numerics::Status::conditioning_budget_exceeded);
    out.maximum_relative_tail_probability_error =
        std::max(out.maximum_relative_tail_probability_error,
                 static_cast<double>(tail_error));
    const auto pdf = std::exp(-z * z / 2) / std::sqrt(2 * std::acos(-1.L));
    const auto quantile_normal_error = std::abs(achieved - tail) / pdf;
    if (requested < .5)
      z = -z;
    const auto log_quantile = log_nominal + spread * z;
    const auto quantile = std::exp(log_quantile);
    double rounded;
    if (!positive_output(quantile, rounded))
      return fail(DensityStatus::numerical_failure,
                  numerics::Status::outside_domain);
    const auto log_error =
        std::abs(spread * z) * variance.estimated_forward_sensitivity / 2 +
        std::abs(spread) * quantile_normal_error +
        8 * epsilon * (1 + std::abs(log_quantile)) +
        std::abs(quantile - rounded) / quantile;
    if (!std::isfinite(log_error) ||
        log_error > policy.maximum_forward_sensitivity)
      return fail(DensityStatus::numerical_failure,
                  numerics::Status::conditioning_budget_exceeded);
    out.maximum_quantile_log_error_estimate =
        std::max(out.maximum_quantile_log_error_estimate,
                 static_cast<double>(log_error));
    quantiles.push_back(rounded);
  }
  out.eta_estimator = std::move(variance);
  out.assumed_eta_mean = eta_mean;
  out.eta_mean_identity = std::move(mean_identity);
  out.probabilities.assign(probabilities.begin(), probabilities.end());
  out.h0_quantiles_km_s_Mpc = std::move(quantiles);
  out.status = DensityStatus::finite;
  out.numerical_status = numerics::Status::ok;
  return out;
}
} // namespace irred::calibration
