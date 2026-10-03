#include "irred/calibration_predictive.hpp"
#include "offset_translation.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
namespace irred::calibration {
namespace {
using D = statistics::DensityStatus;
using N = numerics::Status;
using detail::PayloadAccounting;
void model_bytes(PayloadAccounting &b, const Model &m) noexcept {
  b.vector(m.ordered_host_ids);
  for (const auto &s : m.ordered_host_ids)
    b.string(s);
  b.vector(m.rows);
  for (const auto &r : m.rows) {
    b.string(r.row_id);
    b.string(r.host_id);
    b.string(r.event_id);
  }
  for (const auto *s :
       {&m.magnitude_convention, &m.metallicity_coordinate_identity,
        &m.distance_shape_identity, &m.calibration_identity,
        &m.dependence_identity, &m.conditional_covariance_identity})
    b.string(*s);
}
// Conservative candidate linearization storage, including generated strings.
void linear_bound(PayloadAccounting &b, const Model &m) noexcept {
  b.add(linearization_preparation_payload_bound(m).value_or(SIZE_MAX), 1);
}
void linear_bytes(PayloadAccounting &b, const Linearization &l) noexcept {
  b.embedded(l.retained_payload_bound(), sizeof(Linearization));
}
bool source_matches(const statistics::Gaussian &g, const Model &m) {
  const auto &md = g.metadata();
  return md.source_semantics == "synthetic controls" &&
         md.measure == "product d(mag)" && !md.table_identity.empty() &&
         !md.ordering_provenance.empty() &&
         md.calibration_provenance == m.calibration_identity &&
         md.dependence_provenance == m.dependence_identity &&
         md.uncertainty_identity == m.conditional_covariance_identity;
}
bool same_models(const Model &a, const Model &b) {
  return a.ordered_host_ids == b.ordered_host_ids &&
         a.metallicity_reference_dex == b.metallicity_reference_dex &&
         a.h_reference_km_s_Mpc == b.h_reference_km_s_Mpc &&
         a.magnitude_convention == b.magnitude_convention &&
         a.metallicity_coordinate_identity ==
             b.metallicity_coordinate_identity &&
         a.distance_shape_identity == b.distance_shape_identity &&
         a.calibration_identity == b.calibration_identity &&
         a.dependence_identity == b.dependence_identity;
}
bool normal(double x) {
  return std::isfinite(x) && (x == 0 || std::isnormal(x));
}
using detail::center_exact;
template <class P> bool valid_policy(const P &p) {
  return p.maximum_elements && p.maximum_payload_bytes &&
         p.maximum_work_units && std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
template <class P> P tighter(P a, const P &b) {
  a.maximum_elements = std::min(a.maximum_elements, b.maximum_elements);
  a.maximum_payload_bytes =
      std::min(a.maximum_payload_bytes, b.maximum_payload_bytes);
  a.maximum_work_units = std::min(a.maximum_work_units, b.maximum_work_units);
  a.maximum_forward_sensitivity =
      std::min(a.maximum_forward_sensitivity, b.maximum_forward_sensitivity);
  return a;
}
std::optional<size_t> wrapper_work(const Model &m) noexcept {
  PayloadAccounting b(0);
  if (m.ordered_host_ids.size() > SIZE_MAX - 6)
    return {};
  const auto p = m.ordered_host_ids.size() + 6;
  if (p && m.rows.size() > SIZE_MAX / p)
    return {};
  b.add(m.rows.size() * p, 64);
  b.add(m.rows.size(), 64);
  return b.result();
}
template <class P> bool reserve_work(P &p, std::optional<size_t> work) {
  if (!work || *work >= p.maximum_work_units)
    return false;
  p.maximum_work_units -= *work;
  return true;
}
bool fits(std::optional<size_t> a, size_t limit) { return a && *a <= limit; }
statistics::GaussianResult density_fail(D d, N n) {
  statistics::GaussianResult r;
  r.density.status = d;
  r.density.numerical_status = n;
  return r;
}
statistics::PosteriorMean mean_fail(D d, N n) {
  statistics::PosteriorMean r;
  r.status = d;
  r.numerical_status = n;
  return r;
}
} // namespace
std::optional<size_t> LadderPosterior::preparation_payload_bound(
    const statistics::Gaussian &g, const Model &m,
    const statistics::ParameterPrior &p) noexcept {
  PayloadAccounting b(sizeof(LadderPosterior));
  model_bytes(b, m);
  linear_bound(b, m);
  b.add(statistics::GaussianPosterior::preparation_payload_bound(g, p).value_or(
            SIZE_MAX),
        1);
  // Linearization validation sets, temporary strings and equation construction.
  model_bytes(b, m);
  linear_bound(b, m);
  return b.result();
}
std::optional<size_t> LadderPosterior::retained_payload_bound() const noexcept {
  PayloadAccounting b(sizeof(LadderPosterior));
  model_bytes(b, model_);
  linear_bytes(b, linear_);
  if (posterior_)
    b.embedded(posterior_->retained_payload_bound(),
               sizeof(statistics::GaussianPosterior));
  return b.result();
}
std::optional<size_t>
LadderPosterior::evaluation_payload_bound() const noexcept {
  PayloadAccounting b(0);
  b.add(linear_.offsets_mag.size(), sizeof(double));
  if (posterior_)
    b.add(posterior_->evaluation_payload_bound().value_or(SIZE_MAX), 1);
  return b.result();
}
LadderPosterior::LadderPosterior(LadderPosterior &&o) noexcept {
  *this = std::move(o);
}
LadderPosterior &LadderPosterior::operator=(LadderPosterior &&o) noexcept {
  if (this != &o) {
    model_ = std::move(o.model_);
    linear_ = std::move(o.linear_);
    posterior_ = std::move(o.posterior_);
    policy_ = o.policy_;
    status_ = std::exchange(o.status_, D::invalid_input);
    numerical_status_ = std::exchange(o.numerical_status_, N::invalid_input);
    o.model_ = Model{};
    o.linear_ = Linearization{};
  }
  return *this;
}
LadderPosterior LadderPosterior::prepare(statistics::Gaussian &&g, Model m,
                                         statistics::ParameterPrior prior,
                                         statistics::PosteriorPolicy policy) {
  LadderPosterior out;
  auto fail = [&](D d, N n) {
    out = LadderPosterior{};
    out.status_ = d;
    out.numerical_status_ = n;
    return std::move(out);
  };
  if (!valid_policy(policy))
    return out;
  if (!fits(preparation_payload_bound(g, m, prior),
            policy.maximum_payload_bytes))
    return fail(D::numerical_failure, N::work_limit);
  auto parent_policy = policy;
  if (!reserve_work(parent_policy, wrapper_work(m)))
    return fail(D::numerical_failure, N::work_limit);
  try {
    auto l =
        linearize(m, {policy.maximum_elements, policy.maximum_payload_bytes,
                      policy.maximum_forward_sensitivity});
    if (l.status != D::finite)
      return fail(l.status, l.numerical_status);
    const auto &md = l.metadata;
    if (!source_matches(g, m) ||
        prior.ordered_parameter_ids != md.ordered_parameter_ids ||
        prior.parameter_units != md.parameter_units ||
        prior.shared_nuisance_ids != md.shared_nuisance_ids ||
        prior.residual_unit != md.residual_unit ||
        prior.design_identity != md.design_identity ||
        prior.dependence_identity != md.dependence_identity)
      return fail(D::incompatible_metadata, N::invalid_input);
    // All wrapper candidate allocations precede the parent's final source move.
    out.model_ = std::move(m);
    out.linear_ = std::move(l);
    out.policy_ = policy;
    out.posterior_ = statistics::GaussianPosterior::prepare(
        std::move(g), out.linear_.design, out.linear_.ordered_row_ids,
        std::move(prior), parent_policy);
    if (out.posterior_->status() != D::finite)
      return fail(out.posterior_->status(), out.posterior_->numerical_status());
    out.status_ = D::finite;
    out.numerical_status_ = N::ok;
  } catch (const std::bad_alloc &) {
    return fail(D::numerical_failure, N::work_limit);
  }
  return out;
}
statistics::PosteriorMean
LadderPosterior::condition(std::span<const double> y,
                           std::span<const std::string> rows,
                           statistics::PosteriorPolicy policy) const {
  if (status_ != D::finite)
    return mean_fail(status_, numerical_status_);
  if (!valid_policy(policy))
    return mean_fail(D::invalid_input, N::invalid_input);
  policy = tighter(policy, policy_);
  PayloadAccounting b(0);
  b.add(retained_payload_bound().value_or(SIZE_MAX), 1);
  b.add(evaluation_payload_bound().value_or(SIZE_MAX), 1);
  if (!fits(b.result(), policy.maximum_payload_bytes))
    return mean_fail(D::numerical_failure, N::work_limit);
  try {
    std::vector<double> r;
    auto ns = center_exact(y, linear_.offsets_mag, r);
    if (ns != N::ok)
      return mean_fail(D::numerical_failure, ns);
    if (!reserve_work(policy, linear_.offsets_mag.size()))
      return mean_fail(D::numerical_failure, N::work_limit);
    return posterior_->condition(r, rows, policy);
  } catch (const std::bad_alloc &) {
    return mean_fail(D::numerical_failure, N::work_limit);
  }
}
statistics::GaussianResult LadderPosterior::log_density(
    std::span<const double> y, std::span<const std::string> rows,
    std::span<const double> beta, std::span<const std::string> ids,
    statistics::PosteriorPolicy policy) const {
  if (status_ != D::finite)
    return density_fail(status_, numerical_status_);
  if (!valid_policy(policy))
    return density_fail(D::invalid_input, N::invalid_input);
  policy = tighter(policy, policy_);
  PayloadAccounting b(0);
  b.add(retained_payload_bound().value_or(SIZE_MAX), 1);
  b.add(evaluation_payload_bound().value_or(SIZE_MAX), 1);
  if (!fits(b.result(), policy.maximum_payload_bytes))
    return density_fail(D::numerical_failure, N::work_limit);
  try {
    std::vector<double> r;
    auto ns = center_exact(y, linear_.offsets_mag, r);
    if (ns != N::ok)
      return density_fail(D::numerical_failure, ns);
    if (!reserve_work(policy, linear_.offsets_mag.size()))
      return density_fail(D::numerical_failure, N::work_limit);
    return posterior_->log_density(r, rows, beta, ids, policy);
  } catch (const std::bad_alloc &) {
    return density_fail(D::numerical_failure, N::work_limit);
  }
}
H0PosteriorProjection
LadderPosterior::h0_projection(std::span<const double> y,
                               std::span<const std::string> rows,
                               statistics::PosteriorPolicy policy) const {
  H0PosteriorProjection out;
  if (status_ != D::finite) {
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!valid_policy(policy))
    return out;
  policy = tighter(policy, policy_);
  PayloadAccounting projection_work(32);
  const auto parameter_count = linear_.metadata.ordered_parameter_ids.size();
  if (parameter_count && parameter_count > SIZE_MAX / parameter_count)
    projection_work.add(SIZE_MAX, 2);
  else
    projection_work.add(parameter_count * parameter_count, 1);
  if (!reserve_work(policy, projection_work.result())) {
    out.status = D::numerical_failure;
    out.numerical_status = N::work_limit;
    return out;
  }
  const auto m = condition(y, rows, policy);
  out.status = m.status;
  out.numerical_status = m.numerical_status;
  if (m.status != D::finite)
    return out;
  const size_t i = model_.ordered_host_ids.size() + 4, p = m.value.size();
  const auto a = log_h0_scale();
  const double eta = m.value[i], v = posterior_->covariance()[i * p + i];
  long double vnorm = 0;
  const auto covariance = posterior_->covariance();
  for (size_t row = 0; row < p; ++row) {
    long double sum = 0;
    for (size_t col = 0; col < p; ++col)
      sum += std::abs(static_cast<long double>(covariance[row * p + col]));
    vnorm = std::max(vnorm, sum);
  }
  const long double variance_error =
      vnorm * posterior_->covariance_relative_error_estimate();
  const double median = project_h0(model_, eta);
  const long double expectation =
      static_cast<long double>(median) * std::exp(a * a * v / 2);
  const double e = static_cast<double>(expectation),
               sd = static_cast<double>(a *
                                        std::sqrt(static_cast<long double>(v)));
  if (!std::isnormal(median) || !std::isnormal(e) || !std::isnormal(sd) ||
      !std::isfinite(expectation)) {
    out = H0PosteriorProjection{};
    out.status = D::numerical_failure;
    out.numerical_status = N::outside_domain;
    return out;
  }
  const auto wide_eps = std::numeric_limits<long double>::epsilon();
  const long double median_error = a * m.absolute_error_estimates[i] +
                                   16 * wide_eps * (1 + std::abs(a * eta)) +
                                   2 * std::numeric_limits<double>::epsilon();
  const long double expectation_error = median_error +
                                        a * a * variance_error / 2 +
                                        16 * wide_eps * (1 + a * a * v / 2) +
                                        std::abs(expectation - e) / expectation;
  const long double sd_error = variance_error / (2 * v) + 16 * wide_eps +
                               2 * std::numeric_limits<double>::epsilon();
  const double limit = std::min(policy.maximum_forward_sensitivity,
                                policy_.maximum_forward_sensitivity);
  auto upward = [](long double x) {
    double d = static_cast<double>(x);
    if (static_cast<long double>(d) < x)
      d = std::nextafter(d, std::numeric_limits<double>::infinity());
    return d;
  };
  if (!std::isfinite(expectation_error) || !std::isfinite(sd_error) ||
      std::max({median_error, expectation_error, sd_error}) > limit ||
      !normal(upward(variance_error)) ||
      (variance_error > 0 && upward(variance_error) == 0)) {
    out = H0PosteriorProjection{};
    out.status = D::numerical_failure;
    out.numerical_status = N::conditioning_budget_exceeded;
    return out;
  }
  out.eta_mean_absolute_error_estimate = m.absolute_error_estimates[i];
  out.eta_variance_absolute_error_estimate = upward(variance_error);
  out.median_relative_error_estimate = upward(median_error);
  out.expectation_relative_error_estimate = upward(expectation_error);
  out.log_standard_deviation_relative_error_estimate = upward(sd_error);
  out.eta_mean = eta;
  out.eta_variance = v;
  out.median_h0_km_s_Mpc = median;
  out.expectation_h0_km_s_Mpc = e;
  out.log_standard_deviation = sd;
  return out;
}
std::optional<size_t> LadderPredictive::preparation_payload_bound(
    const LadderPosterior &p, const statistics::Gaussian &g, const Model &m,
    const statistics::PredictiveMetadata &md) noexcept {
  if (p.status() != D::finite)
    return {};
  PayloadAccounting b(sizeof(LadderPredictive));
  model_bytes(b, m);
  linear_bound(b, m);
  model_bytes(b, m);
  linear_bound(b, m);
  b.add(statistics::GaussianPredictive::preparation_payload_bound(p.posterior(),
                                                                  g, md)
            .value_or(SIZE_MAX),
        1);
  // Parent bound includes retained GaussianPosterior, not its wrapper storage.
  const auto all = p.retained_payload_bound(),
             parent = p.posterior().retained_payload_bound();
  if (!all || !parent || *all < *parent)
    b.add(SIZE_MAX, 2);
  else
    b.add(*all - *parent, 1);
  b.add(p.linearization().offsets_mag.size(), sizeof(double));
  b.add(m.rows.size(), 2 * sizeof(double));
  return b.result();
}
std::optional<size_t>
LadderPredictive::retained_payload_bound() const noexcept {
  PayloadAccounting b(sizeof(LadderPredictive));
  model_bytes(b, model_);
  linear_bytes(b, linear_);
  b.vector(mean_);
  b.vector(mean_errors_);
  b.embedded(predictive_.retained_payload_bound(),
             sizeof(statistics::GaussianPredictive));
  return b.result();
}
std::optional<size_t>
LadderPredictive::evaluation_payload_bound() const noexcept {
  PayloadAccounting b(0);
  b.add(linear_.offsets_mag.size(), sizeof(double));
  b.add(predictive_.evaluation_payload_bound().value_or(SIZE_MAX), 1);
  return b.result();
}
LadderPredictive::LadderPredictive(LadderPredictive &&o) noexcept {
  *this = std::move(o);
}
LadderPredictive &LadderPredictive::operator=(LadderPredictive &&o) noexcept {
  if (this != &o) {
    model_ = std::move(o.model_);
    linear_ = std::move(o.linear_);
    predictive_ = std::move(o.predictive_);
    mean_ = std::move(o.mean_);
    mean_errors_ = std::move(o.mean_errors_);
    policy_ = o.policy_;
    status_ = std::exchange(o.status_, D::invalid_input);
    numerical_status_ = std::exchange(o.numerical_status_, N::invalid_input);
    o.model_ = Model{};
    o.linear_ = Linearization{};
    o.mean_.clear();
    o.mean_errors_.clear();
  }
  return *this;
}
LadderPredictive LadderPredictive::prepare(
    const LadderPosterior &p, std::span<const double> training,
    std::span<const std::string> rows, const statistics::Gaussian &g, Model m,
    statistics::PredictiveMetadata md, statistics::PredictivePolicy policy) {
  return prepare_impl(p, training, rows, g, std::move(m), std::move(md), policy,
                      false);
}
LadderPredictive LadderPredictive::prepare_impl(
    const LadderPosterior &p, std::span<const double> training,
    std::span<const std::string> rows, const statistics::Gaussian &g, Model m,
    statistics::PredictiveMetadata md, statistics::PredictivePolicy policy,
    bool fixed) {
  LadderPredictive out;
  auto fail = [&](D d, N n) {
    out = LadderPredictive{};
    out.status_ = d;
    out.numerical_status_ = n;
    return std::move(out);
  };
  if (p.status() != D::finite)
    return fail(p.status(), p.numerical_status());
  if (!valid_policy(policy))
    return out;
  if (!fits(preparation_payload_bound(p, g, m, md),
            policy.maximum_payload_bytes))
    return fail(D::numerical_failure, N::work_limit);
  auto parent_policy = policy;
  PayloadAccounting work(0);
  work.add(wrapper_work(m).value_or(SIZE_MAX), 1);
  work.add(training.size(), 32);
  if (!reserve_work(parent_policy, work.result()))
    return fail(D::numerical_failure, N::work_limit);
  try {
    auto l =
        linearize(m, {policy.maximum_elements, policy.maximum_payload_bytes,
                      policy.maximum_forward_sensitivity});
    if (l.status != D::finite)
      return fail(l.status, l.numerical_status);
    if (!same_models(p.model(), m) || !source_matches(g, m) ||
        g.metadata().ordered_ids != l.ordered_row_ids ||
        md.training_event_ids != p.linearization().predictive_event_ids ||
        md.future_event_ids != l.predictive_event_ids ||
        md.ordered_parameter_ids != l.metadata.ordered_parameter_ids ||
        md.parameter_units != l.metadata.parameter_units ||
        md.response_identity != l.metadata.design_identity ||
        md.dependence_identity != l.metadata.dependence_identity ||
        md.conditional_noise_identity != m.conditional_covariance_identity)
      return fail(D::incompatible_metadata, N::invalid_input);
    std::vector<double> r;
    if (!fixed) {
      const auto ns = center_exact(training, p.linearization().offsets_mag, r);
      if (ns != N::ok)
        return fail(D::numerical_failure, ns);
    }
    out.predictive_ = statistics::GaussianPredictive::prepare_impl(
        p.posterior(), r, rows, g, l.design, md, parent_policy, fixed);
    if (out.predictive_.status() != D::finite)
      return fail(out.predictive_.status(), out.predictive_.numerical_status());
    if (!fixed) {
      const auto ns = detail::add_offsets(
          out.predictive_.mean(),
          out.predictive_.mean_absolute_error_estimates(), l.offsets_mag,
          policy.maximum_forward_sensitivity, out.mean_, out.mean_errors_);
      if (ns != N::ok)
        return fail(D::numerical_failure, ns);
    }
    out.model_ = std::move(m);
    out.linear_ = std::move(l);
    out.policy_ = policy;
    out.status_ = D::finite;
    out.numerical_status_ = N::ok;
  } catch (const std::bad_alloc &) {
    return fail(D::numerical_failure, N::work_limit);
  }
  return out;
}
statistics::GaussianResult
LadderPredictive::log_density(std::span<const double> y,
                              std::span<const std::string> rows,
                              statistics::PredictivePolicy policy) const {
  if (status_ != D::finite)
    return density_fail(status_, numerical_status_);
  if (!valid_policy(policy))
    return density_fail(D::invalid_input, N::invalid_input);
  policy = tighter(policy, policy_);
  PayloadAccounting b(0);
  b.add(retained_payload_bound().value_or(SIZE_MAX), 1);
  b.add(evaluation_payload_bound().value_or(SIZE_MAX), 1);
  if (!fits(b.result(), policy.maximum_payload_bytes))
    return density_fail(D::numerical_failure, N::work_limit);
  try {
    std::vector<double> r;
    auto ns = center_exact(y, linear_.offsets_mag, r);
    if (ns != N::ok)
      return density_fail(D::numerical_failure, ns);
    if (!reserve_work(policy, linear_.offsets_mag.size()))
      return density_fail(D::numerical_failure, N::work_limit);
    return predictive_.log_density(r, rows, policy);
  } catch (const std::bad_alloc &) {
    return density_fail(D::numerical_failure, N::work_limit);
  }
}

LadderPredictiveConditioning::LadderPredictiveConditioning(
    LadderPredictiveConditioning &&o) noexcept
    : posterior_(std::move(o.posterior_)), law_(std::move(o.law_)) {
  o.posterior_.reset();
}
LadderPredictiveConditioning &LadderPredictiveConditioning::operator=(
    LadderPredictiveConditioning &&o) noexcept {
  if (this != &o) {
    this->~LadderPredictiveConditioning();
    new (this) LadderPredictiveConditioning(std::move(o));
  }
  return *this;
}
LadderPredictiveConditioning LadderPredictiveConditioning::prepare(
    LadderPosterior &&p, const statistics::Gaussian &g, Model m,
    statistics::PredictiveMetadata md, statistics::PredictivePolicy policy) {
  LadderPredictiveConditioning out;
  if (p.status() != D::finite || g.status() != D::finite) {
    out.law_.status_ = p.status() != D::finite ? p.status() : g.status();
    out.law_.numerical_status_ =
        p.status() != D::finite ? p.numerical_status() : g.numerical_status();
    return out;
  }
  if (!fits(preparation_payload_bound(p, g, m, md),
            policy.maximum_payload_bytes)) {
    out.law_.numerical_status_ = N::work_limit;
    return out;
  }
  out.law_ = LadderPredictive::prepare_impl(
      p, {}, p.linearization().ordered_row_ids, g, std::move(m), std::move(md),
      policy, true);
  if (out.law_.status() == D::finite) {
    PayloadAccounting b(sizeof(out));
    b.embedded(p.retained_payload_bound(), sizeof(LadderPosterior));
    b.embedded(out.law_.retained_payload_bound(), sizeof(LadderPredictive));
    if (!fits(b.result(), policy.maximum_payload_bytes)) {
      out.law_ = LadderPredictive{};
      out.law_.status_ = D::numerical_failure;
      out.law_.numerical_status_ = N::work_limit;
    } else
      out.posterior_ = std::move(p);
  }
  return out;
}
std::optional<size_t> LadderPredictiveConditioning::preparation_payload_bound(
    const LadderPosterior &p, const statistics::Gaussian &g, const Model &m,
    const statistics::PredictiveMetadata &md) noexcept {
  PayloadAccounting b(sizeof(LadderPredictiveConditioning) -
                      sizeof(LadderPredictive));
  b.add(LadderPredictive::preparation_payload_bound(p, g, m, md)
            .value_or(SIZE_MAX),
        1);
  return b.result();
}
std::optional<size_t>
LadderPredictiveConditioning::retained_payload_bound() const noexcept {
  PayloadAccounting b(sizeof(*this));
  if (posterior_)
    b.embedded(posterior_->retained_payload_bound(), sizeof(LadderPosterior));
  b.embedded(law_.retained_payload_bound(), sizeof(LadderPredictive));
  return b.result();
}
std::optional<size_t> LadderPredictiveConditioning::batch_payload_bound(
    size_t count, statistics::PredictiveOutputs outputs) const noexcept {
  if (!posterior_)
    return {};
  const auto total = retained_payload_bound(),
             p = posterior_->posterior().retained_payload_bound(),
             l = law_.predictive_.retained_payload_bound();
  if (!total || !p || !l || *p > *total || *l > *total - *p)
    return {};
  return law_.predictive_.batch_payload_bound(posterior_->posterior(), count,
                                              outputs, *total - *p - *l);
}
statistics::PredictiveBatch LadderPredictiveConditioning::evaluate(
    std::span<const double> training, std::span<const std::string> rows,
    std::span<const double> future, std::span<const std::string> future_rows,
    size_t count, statistics::PredictiveOutputs outputs,
    statistics::PredictivePolicy policy) const {
  if (!posterior_) {
    statistics::PredictiveBatch out;
    out.status = status();
    out.numerical_status = numerical_status();
    return out;
  }
  const auto total = retained_payload_bound(),
             p = posterior_->posterior().retained_payload_bound(),
             l = law_.predictive_.retained_payload_bound();
  if (!total || !p || !l || *p > *total || *l > *total - *p) {
    statistics::PredictiveBatch out;
    out.numerical_status = N::work_limit;
    return out;
  }
  if (!valid_policy(policy)) {
    statistics::PredictiveBatch out;
    return out;
  }
  policy = tighter(policy, law_.policy_);
  return law_.predictive_.condition_batch(
      posterior_->posterior(), training, rows, future, future_rows, count,
      outputs, policy, posterior_->linearization().offsets_mag,
      law_.linear_.offsets_mag, *total - *p - *l);
}
} // namespace irred::calibration
