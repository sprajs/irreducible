#include "irred/calibration_ladder.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <set>

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
  const long double exponent =
      std::log(static_cast<long double>(model_.h_reference_km_s_Mpc)) +
      std::log(10.L) *
          out.relative_fit.coefficients[model_.ordered_host_ids.size() + 4] / 5;
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
} // namespace irred::calibration
