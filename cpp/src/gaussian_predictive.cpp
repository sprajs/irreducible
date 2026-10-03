#include "irred/gaussian_predictive.hpp"
#include "offset_translation.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
namespace irred::statistics {
namespace {
using W = long double;
using S = numerics::Status;
bool normal(double x) { return x == 0 || std::fpclassify(x) == FP_NORMAL; }
bool environment() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
bool valid(PredictivePolicy p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
bool same(std::span<const std::string> a, std::span<const std::string> b) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}
bool unique(std::span<const std::string> ids) {
  for (size_t i = 0; i < ids.size(); ++i)
    if (ids[i].empty() ||
        std::find(ids.begin(), ids.begin() + i, ids[i]) != ids.begin() + i)
      return false;
  return true;
}
bool disjoint(std::span<const std::string> a, std::span<const std::string> b) {
  for (const auto &id : a)
    if (std::find(b.begin(), b.end(), id) != b.end())
      return false;
  return true;
}
bool cast(W x, double &v, S &status) {
  if (!std::isfinite(x) || std::abs(x) > std::numeric_limits<double>::max()) {
    status = S::overflow;
    return false;
  }
  v = static_cast<double>(x);
  if (x != 0 && std::fpclassify(v) != FP_NORMAL) {
    status = S::outside_domain;
    return false;
  }
  return true;
}
bool dimensions(size_t n, size_t p, size_t k) {
  return n && p >= 2 && k && n <= SIZE_MAX / n && p <= SIZE_MAX / p &&
         k <= SIZE_MAX / k && n <= SIZE_MAX / p && k <= SIZE_MAX / p;
}
std::optional<size_t> work(size_t n, size_t p, size_t k) {
  if (!dimensions(n, p, k))
    return {};
  detail::PayloadAccounting b(0);
  b.add(n * n, 1);
  b.add(n * p, 1);
  b.add(p * p, 1);
  b.add(k * p, p);
  b.add(k * k, p);
  b.add(k * k, 8 * k);
  const auto count = b.result();
  if (!count || *count > SIZE_MAX / 32)
    return {};
  return *count * 32;
}
void charge(detail::PayloadAccounting &b, const Metadata &m) noexcept {
  b.strings(m.ordered_ids);
  for (const auto *s :
       {&m.arithmetic_id, &m.measure, &m.table_identity,
        &m.uncertainty_identity, &m.ordering_provenance,
        &m.calibration_provenance, &m.dependence_provenance,
        &m.source_semantics, &m.input_matrix_convention, &m.treatment})
    b.string(*s);
}
void charge(detail::PayloadAccounting &b,
            const PredictiveMetadata &m) noexcept {
  for (const auto *v :
       {&m.ordered_parameter_ids, &m.parameter_units, &m.response_units,
        &m.training_event_ids, &m.future_event_ids})
    b.strings(*v);
  for (const auto *s :
       {&m.future_unit, &m.future_covariance_unit, &m.future_measure,
        &m.response_identity, &m.conditioning_identity,
        &m.conditional_noise_identity, &m.dependence_identity})
    b.string(*s);
}
PredictivePolicy restrict(PredictivePolicy a, PredictivePolicy b) {
  return {
      std::min(a.maximum_elements, b.maximum_elements),
      std::min(a.maximum_payload_bytes, b.maximum_payload_bytes),
      std::min(a.maximum_work_units, b.maximum_work_units),
      std::min(a.maximum_forward_sensitivity, b.maximum_forward_sensitivity)};
}
PosteriorMean project_mean(const PosteriorMean &conditioned,
                           std::span<const double> a, size_t k,
                           PredictivePolicy policy) {
  PosteriorMean out;
  if (conditioned.status != DensityStatus::finite) {
    out.status = conditioned.status;
    out.numerical_status = conditioned.numerical_status;
    return out;
  }
  const auto p = conditioned.value.size();
  out.value.resize(k);
  out.absolute_error_estimates.resize(k);
  S ns = S::ok;
  for (size_t i = 0; i < k; ++i) {
    W mu = 0, absolute = 0, err = 0;
    for (size_t j = 0; j < p; ++j) {
      const W t = static_cast<W>(a[i * p + j]) * conditioned.value[j];
      mu += t;
      absolute += std::abs(t);
      err += std::abs(a[i * p + j]) * conditioned.absolute_error_estimates[j];
    }
    if (!cast(mu, out.value[i], ns)) {
      out.value.clear();
      out.absolute_error_estimates.clear();
      out.status = DensityStatus::numerical_failure;
      out.numerical_status = ns;
      return out;
    }
    err += (4 * p + 4) * std::numeric_limits<W>::epsilon() * absolute +
           std::abs(mu - out.value[i]);
    if (!cast(err, out.absolute_error_estimates[i], ns) ||
        err / (1 + std::abs(mu)) > policy.maximum_forward_sensitivity) {
      out.value.clear();
      out.absolute_error_estimates.clear();
      out.status = DensityStatus::numerical_failure;
      out.numerical_status = S::conditioning_budget_exceeded;
      return out;
    }
  }
  out.status = DensityStatus::finite;
  out.numerical_status = S::ok;
  return out;
}
void clear(Metadata &m) noexcept {
  m.ordered_ids.clear();
  for (auto *s :
       {&m.arithmetic_id, &m.measure, &m.table_identity,
        &m.uncertainty_identity, &m.ordering_provenance,
        &m.calibration_provenance, &m.dependence_provenance,
        &m.source_semantics, &m.input_matrix_convention, &m.treatment})
    s->clear();
}
} // namespace
void GaussianPredictive::invalidate() noexcept {
  predictive_.reset();
  metadata_ = {};
  clear(training_source_);
  clear(future_noise_);
  prior_identity_.clear();
  training_design_identity_.clear();
  parameter_measure_.clear();
  prior_dependence_identity_.clear();
  shared_nuisance_ids_.clear();
  training_values_.clear();
  response_.clear();
  mean_.clear();
  mean_error_.clear();
  covariance_error_.clear();
  inverse_norm_estimate_ = covariance_perturbation_ = 0;
  relative_covariance_error_ = 0;
  preparation_work_units_ = 0;
  status_ = DensityStatus::invalid_input;
  numerical_status_ = S::invalid_input;
}
GaussianPredictive::GaussianPredictive(GaussianPredictive &&o) noexcept
    : predictive_(std::move(o.predictive_)), metadata_(std::move(o.metadata_)),
      training_source_(std::move(o.training_source_)),
      future_noise_(std::move(o.future_noise_)),
      prior_identity_(std::move(o.prior_identity_)),
      training_design_identity_(std::move(o.training_design_identity_)),
      parameter_measure_(std::move(o.parameter_measure_)),
      prior_dependence_identity_(std::move(o.prior_dependence_identity_)),
      shared_nuisance_ids_(std::move(o.shared_nuisance_ids_)),
      training_values_(std::move(o.training_values_)),
      response_(std::move(o.response_)), mean_(std::move(o.mean_)),
      mean_error_(std::move(o.mean_error_)),
      covariance_error_(std::move(o.covariance_error_)), policy_(o.policy_),
      inverse_norm_estimate_(o.inverse_norm_estimate_),
      covariance_perturbation_(o.covariance_perturbation_),
      relative_covariance_error_(o.relative_covariance_error_),
      preparation_work_units_(o.preparation_work_units_), status_(o.status_),
      numerical_status_(o.numerical_status_) {
  o.invalidate();
}
GaussianPredictive &
GaussianPredictive::operator=(GaussianPredictive &&o) noexcept {
  if (this != &o) {
    this->~GaussianPredictive();
    new (this) GaussianPredictive(std::move(o));
  }
  return *this;
}
std::optional<size_t> GaussianPredictive::preparation_payload_bound(
    const GaussianPosterior &posterior, const Gaussian &noise,
    const PredictiveMetadata &m) noexcept {
  const auto n = posterior.source().metadata().ordered_ids.size(),
             p = posterior.prior().mean.size(),
             k = noise.metadata().ordered_ids.size();
  if (!work(n, p, k))
    return {};
  detail::PayloadAccounting b(sizeof(GaussianPredictive));
  // Already-retained immutable borrowed inputs are part of this simultaneous
  // envelope. Their factors remain in their original owners, without copies.
  const auto pr = posterior.retained_payload_bound(),
             nr = noise.retained_payload_bound(),
             pe = posterior.evaluation_payload_bound();
  const auto gp = gaussian_preparation_payload_bound(
      k, MatrixKind::covariance, noise.arithmetic(), noise.metadata());
  if (!pr || !nr || !pe || !gp)
    return {};
  b.add(*pr, 1);
  b.add(*nr, 1);
  b.add(*pe, 1);
  b.add(*gp, 1);
  for (unsigned i = 0; i < 4; ++i)
    charge(b, m);
  charge(b, posterior.source().metadata());
  charge(b, noise.metadata());
  const auto &prior = posterior.prior();
  b.strings(prior.shared_nuisance_ids);
  for (const auto *s : {&prior.prior_identity, &prior.design_identity,
                        &prior.parameter_measure, &prior.dependence_identity})
    b.string(*s);
  // Covers simultaneous wide products/errors, stored A/training/mean/error
  // arrays and candidate covariance. Constant allowance covers fixed labels.
  b.add(k * p, 8 * sizeof(W));
  b.add(k * k, 8 * sizeof(W));
  b.add(p * p, 8 * sizeof(W));
  b.add(n, 8 * sizeof(W));
  b.add(p, 8 * sizeof(W));
  b.add(k, 8 * sizeof(W));
  b.add(4096, 1);
  return b.result();
}
std::optional<size_t>
GaussianPredictive::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  if (predictive_)
    b.embedded(predictive_->retained_payload_bound(), sizeof(Gaussian));
  charge(b, metadata_);
  charge(b, training_source_);
  charge(b, future_noise_);
  for (const auto *s : {&prior_identity_, &training_design_identity_,
                        &parameter_measure_, &prior_dependence_identity_})
    b.string(*s);
  b.strings(shared_nuisance_ids_);
  for (const auto *v : {&training_values_, &response_, &mean_, &mean_error_,
                        &covariance_error_})
    b.vector(*v);
  return b.result();
}
std::optional<size_t>
GaussianPredictive::evaluation_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(GaussianResult));
  if (!predictive_)
    return {};
  const auto e = predictive_->evaluation_payload_bound(
      future_noise_.ordered_ids.size(), false);
  if (!e)
    return {};
  b.add(*e, 1);
  b.add(future_noise_.ordered_ids.size(), 8 * sizeof(W));
  return b.result();
}
GaussianPredictive GaussianPredictive::prepare(
    const GaussianPosterior &posterior, std::span<const double> training,
    std::span<const std::string> rows, const Gaussian &noise,
    std::span<const double> a, const PredictiveMetadata &m,
    PredictivePolicy policy) {
  return prepare_impl(posterior, training, rows, noise, a, m, policy, false);
}
GaussianPredictive GaussianPredictive::prepare_impl(
    const GaussianPosterior &posterior, std::span<const double> training,
    std::span<const std::string> rows, const Gaussian &noise,
    std::span<const double> a, const PredictiveMetadata &m,
    PredictivePolicy policy, bool fixed) {
  GaussianPredictive out;
  const auto fail = [&](DensityStatus ds, S ns) {
    out.invalidate();
    out.status_ = ds;
    out.numerical_status_ = ns;
    return std::move(out);
  };
  if (!environment())
    return fail(DensityStatus::unsupported_domain, S::outside_domain);
  if (!valid(policy))
    return out;
  if (posterior.status() != DensityStatus::finite)
    return fail(posterior.status(), posterior.numerical_status());
  if (noise.status() != DensityStatus::finite)
    return fail(noise.status(), noise.numerical_status());
  const auto n = posterior.source().metadata().ordered_ids.size(),
             p = posterior.prior().mean.size(),
             k = noise.metadata().ordered_ids.size();
  const auto units = work(n, p, k),
             peak = preparation_payload_bound(posterior, noise, m);
  if (!units || !peak || n * n > policy.maximum_elements ||
      n * p > policy.maximum_elements || p * p > policy.maximum_elements ||
      k * k > policy.maximum_elements || k * p > policy.maximum_elements ||
      *units > policy.maximum_work_units ||
      *peak > policy.maximum_payload_bytes)
    return fail(DensityStatus::invalid_input, S::work_limit);
  try {
    const auto &prior = posterior.prior();
    if (!same(rows, posterior.source().metadata().ordered_ids) ||
        !same(m.ordered_parameter_ids, prior.ordered_parameter_ids) ||
        !same(m.parameter_units, prior.parameter_units) ||
        m.training_event_ids.size() != n || m.future_event_ids.size() != k ||
        !unique(m.training_event_ids) || !unique(m.future_event_ids) ||
        !disjoint(rows, noise.metadata().ordered_ids) ||
        !disjoint(m.training_event_ids, m.future_event_ids))
      return fail(DensityStatus::incompatible_metadata, S::invalid_input);
    if ((!fixed && training.size() != n) || a.size() != k * p ||
        m.response_units.size() != p || !m.future_noise_independence_declared ||
        !m.noise_conditional_on_parameters_declared ||
        m.response_identity.empty() || m.conditioning_identity.empty() ||
        m.conditional_noise_identity.empty() || m.dependence_identity.empty() ||
        m.future_unit.empty() || m.future_unit != prior.residual_unit ||
        m.future_covariance_unit != m.future_unit + "^2" ||
        m.future_measure != "product d(" + m.future_unit + ")" ||
        m.future_measure != noise.metadata().measure ||
        posterior.source().metadata().source_semantics !=
            "synthetic controls" ||
        noise.metadata().source_semantics != "synthetic controls" ||
        noise.arithmetic() != posterior.source().arithmetic() ||
        !noise.priors().empty() || !noise.mean_shift().empty())
      return fail(DensityStatus::incompatible_metadata, S::invalid_input);
    for (size_t j = 0; j < p; ++j)
      if (m.response_units[j] != m.future_unit + "/" + m.parameter_units[j])
        return fail(DensityStatus::incompatible_metadata, S::invalid_input);
    for (auto values : {training, a})
      for (auto v : values)
        if (!normal(v))
          return fail(DensityStatus::invalid_input, S::nonfinite_input);
    PosteriorPolicy pp{policy.maximum_elements, policy.maximum_payload_bytes,
                       policy.maximum_work_units,
                       policy.maximum_forward_sensitivity};
    PosteriorMean projected;
    if (!fixed) {
      const auto conditioned = posterior.condition(training, rows, pp);
      projected = project_mean(conditioned, a, k, policy);
      if (projected.status != DensityStatus::finite)
        return fail(projected.status, projected.numerical_status);
    }
    const auto v = posterior.covariance();
    const auto eps = std::numeric_limits<W>::epsilon();
    W vnorm = 0;
    for (size_t i = 0; i < p; ++i) {
      W row = 0;
      for (size_t j = 0; j < p; ++j)
        row += std::abs(v[i * p + j]);
      vnorm = std::max(vnorm, row);
    }
    const W inherited = vnorm * posterior.covariance_relative_error_estimate();
    std::vector<W> b(k * p), be(k * p), error(k * k), row_norm(k);
    std::vector<double> covariance(k * k), covariance_error(k * k);
    S projection = S::ok;
    for (size_t i = 0; i < k; ++i) {
      for (size_t j = 0; j < p; ++j)
        row_norm[i] += std::abs(a[i * p + j]);
      for (size_t j = 0; j < p; ++j) {
        W absolute_product = 0;
        for (size_t t = 0; t < p; ++t) {
          const W product = static_cast<W>(a[i * p + t]) * v[t * p + j];
          b[i * p + j] += product;
          absolute_product += std::abs(product);
        }
        be[i * p + j] =
            row_norm[i] * inherited + (4 * p + 4) * eps * absolute_product;
      }
    }
    for (size_t i = 0; i < k; ++i)
      for (size_t j = i; j < k; ++j) {
        W latent = 0, absolute = 0, err = 0;
        for (size_t t = 0; t < p; ++t) {
          const W product = b[i * p + t] * a[j * p + t];
          latent += product;
          absolute += std::abs(product);
          err += be[i * p + t] * std::abs(a[j * p + t]);
        }
        const W noise_value = noise.covariance()[i * k + j],
                value = noise_value + latent;
        if (!cast(value, covariance[i * k + j], projection))
          return fail(DensityStatus::numerical_failure, projection);
        err += (4 * p + 4) * eps * absolute +
               4 * eps * (std::abs(noise_value) + std::abs(latent)) +
               std::abs(value - covariance[i * k + j]);
        if (!cast(err, covariance_error[i * k + j], projection))
          return fail(DensityStatus::numerical_failure,
                      S::conditioning_budget_exceeded);
        covariance[j * k + i] = covariance[i * k + j];
        covariance_error[j * k + i] = covariance_error[i * k + j];
        error[i * k + j] = error[j * k + i] =
            std::max(err, static_cast<W>(covariance_error[i * k + j]));
      }
    Metadata md;
    md.ordered_ids = noise.metadata().ordered_ids;
    md.measure = m.future_measure;
    md.table_identity = m.conditioning_identity;
    md.uncertainty_identity = "R*+A V A^T; proper-Gaussian joint prediction";
    md.ordering_provenance = noise.metadata().ordering_provenance;
    md.calibration_provenance = noise.metadata().calibration_provenance;
    md.dependence_provenance = m.dependence_identity;
    md.source_semantics =
        "derived conditional predictive Gaussian; synthetic controls";
    auto predictive = prepare_gaussian(covariance, MatrixKind::covariance,
                                       std::move(md), policy.maximum_elements,
                                       policy.maximum_forward_sensitivity,
                                       noise.arithmetic());
    if (predictive.status() != DensityStatus::finite)
      return fail(predictive.status(), predictive.numerical_status());
    W norm = 0, error_norm = 0;
    for (size_t i = 0; i < k; ++i) {
      W row = 0, e = 0;
      for (size_t j = 0; j < k; ++j) {
        row += std::abs(covariance[i * k + j]);
        e += error[i * k + j];
      }
      norm = std::max(norm, row);
      error_norm = std::max(error_norm, e);
    }
    const W inverse_norm =
                static_cast<W>(predictive.factor_.condition_estimate_inf()) /
                norm,
            eta = inverse_norm * error_norm;
    if (!(norm > 0) || !std::isfinite(inverse_norm) || !(inverse_norm > 0) ||
        !std::isfinite(eta) || eta >= .01L ||
        eta / (1 - eta) > policy.maximum_forward_sensitivity)
      return fail(DensityStatus::numerical_failure,
                  S::conditioning_budget_exceeded);
    double relative = 0;
    if (!cast(eta / (1 - eta), relative, projection))
      return fail(DensityStatus::numerical_failure,
                  S::conditioning_budget_exceeded);
    out.training_source_ = posterior.source().metadata();
    out.future_noise_ = noise.metadata();
    out.prior_identity_ = prior.prior_identity;
    out.training_design_identity_ = prior.design_identity;
    out.parameter_measure_ = prior.parameter_measure;
    out.prior_dependence_identity_ = prior.dependence_identity;
    out.shared_nuisance_ids_ = prior.shared_nuisance_ids;
    out.training_values_.assign(training.begin(), training.end());
    out.response_.assign(a.begin(), a.end());
    out.metadata_ = m;
    out.mean_ = std::move(projected.value);
    out.mean_error_ = std::move(projected.absolute_error_estimates);
    out.covariance_error_ = std::move(covariance_error);
    out.predictive_ = std::move(predictive);
    out.policy_ = policy;
    out.inverse_norm_estimate_ = inverse_norm;
    out.covariance_perturbation_ = eta;
    out.relative_covariance_error_ = relative;
    out.preparation_work_units_ = *units;
    const auto retained = out.retained_payload_bound();
    if (!retained || *retained > policy.maximum_payload_bytes)
      return fail(DensityStatus::invalid_input, S::work_limit);
    out.status_ = DensityStatus::finite;
    out.numerical_status_ = S::ok;
    return out;
  } catch (const std::bad_alloc &) {
    return fail(DensityStatus::invalid_input, S::work_limit);
  } catch (const std::length_error &) {
    return fail(DensityStatus::invalid_input, S::work_limit);
  }
}
PosteriorMean GaussianPredictive::project(const PosteriorMean &mean,
                                          PredictivePolicy policy) const {
  return project_mean(mean, response_, future_noise_.ordered_ids.size(),
                      policy);
}
GaussianResult
GaussianPredictive::log_density(std::span<const double> y,
                                std::span<const std::string> ids,
                                PredictivePolicy requested) const {
  if (status_ != DensityStatus::finite) {
    GaussianResult out;
    out.density.status = status_;
    out.density.numerical_status = numerical_status_;
    return out;
  }
  if (mean_.empty()) {
    GaussianResult out;
    out.density.status = DensityStatus::invalid_input;
    return out;
  }
  return log_density_at_mean(y, ids, mean_, mean_error_, requested);
}
GaussianResult GaussianPredictive::log_density_at_mean(
    std::span<const double> y, std::span<const std::string> ids,
    std::span<const double> mean, std::span<const double> mean_error,
    PredictivePolicy requested) const {
  GaussianResult out;
  const auto fail = [&](DensityStatus ds, S ns) {
    GaussianResult bad;
    bad.density.status = ds;
    bad.density.numerical_status = ns;
    return bad;
  };
  if (status_ != DensityStatus::finite)
    return fail(status_, numerical_status_);
  if (!environment())
    return fail(DensityStatus::unsupported_domain, S::outside_domain);
  if (!valid(requested))
    return out;
  const auto policy = restrict(requested, policy_);
  const auto k = mean.size();
  if (!same(ids, future_noise_.ordered_ids))
    return fail(DensityStatus::incompatible_metadata, S::invalid_input);
  if (y.size() != k)
    return out;
  const auto retained = retained_payload_bound(),
             scratch = evaluation_payload_bound();
  detail::PayloadAccounting peak(0);
  if (!retained || !scratch)
    return fail(DensityStatus::invalid_input, S::work_limit);
  peak.add(*retained, 1);
  peak.add(*scratch, 1);
  const auto bytes = peak.result();
  if (!bytes || *bytes > policy.maximum_payload_bytes ||
      k * k > policy.maximum_elements || k * k > policy.maximum_work_units / 32)
    return fail(DensityStatus::invalid_input, S::work_limit);
  try {
    std::vector<double> centered(k);
    W delta = 0, z_l1 = 0;
    S projection = S::ok;
    for (size_t i = 0; i < k; ++i) {
      if (!normal(y[i]))
        return fail(DensityStatus::invalid_input, S::nonfinite_input);
      const W wide = static_cast<W>(y[i]) - mean[i];
      if (!cast(wide, centered[i], projection))
        return fail(DensityStatus::numerical_failure, projection);
      delta = std::max(delta, static_cast<W>(mean_error[i]) +
                                  std::abs(wide - centered[i]) +
                                  4 * std::numeric_limits<W>::epsilon() *
                                      (std::abs(y[i]) + std::abs(mean[i])));
      z_l1 += std::abs(centered[i]);
    }
    out = predictive_->evaluate(centered, ids,
                                policy.maximum_forward_sensitivity);
    if (out.density.status != DensityStatus::finite)
      return out;
    const W eta = covariance_perturbation_,
            projection_error = .5L * inverse_norm_estimate_ *
                               (2 * z_l1 * delta + k * delta * delta),
            covariance_error =
                .5L * (out.quadratic * eta / (1 - eta) - k * std::log1p(-eta)),
            inherited = out.estimated_forward_sensitivity *
                        (1 + std::abs(static_cast<W>(out.density.log_value))),
            total = projection_error + covariance_error + inherited,
            scaled =
                total / (1 + std::abs(static_cast<W>(out.density.log_value)));
    if (!std::isfinite(scaled) || scaled > policy.maximum_forward_sensitivity)
      return fail(DensityStatus::numerical_failure,
                  S::conditioning_budget_exceeded);
    double reported = 0;
    if (!cast(scaled, reported, projection))
      return fail(DensityStatus::numerical_failure,
                  S::conditioning_budget_exceeded);
    out.estimated_forward_sensitivity = reported;
    return out;
  } catch (const std::bad_alloc &) {
    return fail(DensityStatus::invalid_input, S::work_limit);
  } catch (const std::length_error &) {
    return fail(DensityStatus::invalid_input, S::work_limit);
  }
}

std::optional<size_t>
GaussianPredictive::batch_payload_bound(const GaussianPosterior &p,
                                        size_t count, PredictiveOutputs outputs,
                                        size_t extra) const noexcept {
  if (status_ != DensityStatus::finite || p.status() != DensityStatus::finite ||
      !count || (!outputs.means && !outputs.densities))
    return {};
  const auto n = p.source().metadata().ordered_ids.size(),
             k = future_noise_.ordered_ids.size();
  detail::PayloadAccounting b(sizeof(PredictiveBatch));
  b.add(retained_payload_bound().value_or(SIZE_MAX), 1);
  b.add(p.retained_payload_bound().value_or(SIZE_MAX), 1);
  b.add(extra, 1);
  b.add(p.evaluation_payload_bound().value_or(SIZE_MAX), 1);
  b.add(evaluation_payload_bound().value_or(SIZE_MAX), 1);
  b.add(n + k, sizeof(double)); // one exact translation scratch, if used
  b.add(k, 4 * sizeof(double)); // projected and offset means/errors
  b.add(count, sizeof(PredictiveRow));
  if (outputs.means) {
    if (k && count > SIZE_MAX / k)
      return {};
    b.add(count * k, 2 * sizeof(double));
  }
  if (outputs.densities)
    b.add(count, sizeof(GaussianResult));
  return b.result();
}
PredictiveBatch GaussianPredictive::condition_batch(
    const GaussianPosterior &p, std::span<const double> training,
    std::span<const std::string> rows, std::span<const double> future,
    std::span<const std::string> future_rows, size_t count,
    PredictiveOutputs outputs, PredictivePolicy requested,
    std::span<const double> training_offsets,
    std::span<const double> future_offsets, size_t extra) const {
  const auto bad = [](DensityStatus ds, S ns) {
    PredictiveBatch b;
    b.status = ds;
    b.numerical_status = ns;
    return b;
  };
  if (status_ != DensityStatus::finite)
    return bad(status_, numerical_status_);
  if (!environment())
    return bad(DensityStatus::unsupported_domain, S::outside_domain);
  if (!valid(requested) || !count || (!outputs.means && !outputs.densities))
    return bad(DensityStatus::invalid_input, S::invalid_input);
  const auto policy = restrict(requested, policy_);
  const auto n = p.source().metadata().ordered_ids.size(),
             parameter_count = p.prior().mean.size(),
             k = future_noise_.ordered_ids.size();
  if (!same(rows, training_source_.ordered_ids) ||
      (outputs.densities && !same(future_rows, future_noise_.ordered_ids)))
    return bad(DensityStatus::incompatible_metadata, S::invalid_input);
  if (!dimensions(n, parameter_count, k) || count > SIZE_MAX / n ||
      count > SIZE_MAX / k || training.size() != count * n ||
      (outputs.densities && future.size() != count * k) ||
      (!training_offsets.empty() && training_offsets.size() != n) ||
      (!future_offsets.empty() && future_offsets.size() != k))
    return bad(DensityStatus::invalid_input, S::invalid_input);
  detail::PayloadAccounting work_count(0);
  work_count.add(n * n, 32);
  work_count.add(n * parameter_count, 32);
  work_count.add(parameter_count * parameter_count, 32);
  work_count.add(k * parameter_count, 32);
  if (outputs.densities)
    work_count.add(k * k, 32);
  work_count.add(training_offsets.size(), 1);
  work_count.add(future_offsets.size(), 1 + (outputs.densities ? 1 : 0));
  const auto one = work_count.result(),
             peak = batch_payload_bound(p, count, outputs, extra);
  if (!one || !peak || *one > SIZE_MAX / count ||
      *one * count > policy.maximum_work_units ||
      *peak > policy.maximum_payload_bytes ||
      count * n > policy.maximum_elements ||
      count * k > policy.maximum_elements || n * n > policy.maximum_elements ||
      n * parameter_count > policy.maximum_elements ||
      parameter_count * parameter_count > policy.maximum_elements ||
      k * k > policy.maximum_elements)
    return bad(DensityStatus::invalid_input, S::work_limit);
  try {
    PredictiveBatch out;
    out.rows.resize(count);
    if (outputs.means) {
      out.means.resize(count * k);
      out.mean_absolute_error_estimates.resize(count * k);
    }
    if (outputs.densities)
      out.densities.resize(count);
    out.work_units = *one * count;
    out.payload_bound = *peak;
    for (size_t i = 0; i < count; ++i) {
      auto &row = out.rows[i];
      try {
        std::vector<double> translated_training, translated_future;
        auto r = training.subspan(i * n, n);
        if (!training_offsets.empty()) {
          const auto ns =
              detail::center_exact(r, training_offsets, translated_training);
          if (ns != S::ok) {
            row = {DensityStatus::numerical_failure, ns};
            continue;
          }
          r = translated_training;
        }
        if (std::any_of(r.begin(), r.end(),
                        [](double v) { return !normal(v); })) {
          row = {DensityStatus::invalid_input, S::nonfinite_input};
          continue;
        }
        PosteriorPolicy pp{
            policy.maximum_elements, policy.maximum_payload_bytes,
            policy.maximum_work_units, policy.maximum_forward_sensitivity};
        auto mean = project(p.condition(r, rows, pp), policy);
        if (mean.status != DensityStatus::finite) {
          row = {mean.status, mean.numerical_status};
          continue;
        }
        // Match complete ladder preparation: original-coordinate mean/error
        // admission is required even when the mean is not requested for output.
        std::vector<double> offset_mean, offset_error;
        std::span<const double> stored_mean = mean.value,
                                stored_error = mean.absolute_error_estimates;
        if (!future_offsets.empty()) {
          const auto ns = detail::add_offsets(
              mean.value, mean.absolute_error_estimates, future_offsets,
              policy.maximum_forward_sensitivity, offset_mean, offset_error);
          if (ns != S::ok) {
            row = {DensityStatus::numerical_failure, ns};
            continue;
          }
          stored_mean = offset_mean;
          stored_error = offset_error;
        }
        GaussianResult density;
        if (outputs.densities) {
          auto y = future.subspan(i * k, k);
          if (!future_offsets.empty()) {
            const auto ns =
                detail::center_exact(y, future_offsets, translated_future);
            if (ns != S::ok) {
              row = {DensityStatus::numerical_failure, ns};
              continue;
            }
            y = translated_future;
          }
          density = log_density_at_mean(y, future_rows, mean.value,
                                        mean.absolute_error_estimates, policy);
          if (density.density.status != DensityStatus::finite) {
            row = {density.density.status, density.density.numerical_status};
            continue;
          }
        }
        row = {DensityStatus::finite, S::ok};
        if (outputs.means) {
          std::copy(stored_mean.begin(), stored_mean.end(),
                    out.means.begin() + i * k);
          std::copy(stored_error.begin(), stored_error.end(),
                    out.mean_absolute_error_estimates.begin() + i * k);
        }
        if (outputs.densities)
          out.densities[i] = density;
      } catch (const std::bad_alloc &) {
        row = {DensityStatus::invalid_input, S::work_limit};
      } catch (const std::length_error &) {
        row = {DensityStatus::invalid_input, S::work_limit};
      }
    }
    out.status = DensityStatus::finite;
    out.numerical_status = S::ok;
    return out;
  } catch (const std::bad_alloc &) {
    return bad(DensityStatus::invalid_input, S::work_limit);
  } catch (const std::length_error &) {
    return bad(DensityStatus::invalid_input, S::work_limit);
  }
}
GaussianPredictiveConditioning::GaussianPredictiveConditioning(
    GaussianPredictiveConditioning &&o) noexcept
    : posterior_(std::move(o.posterior_)), law_(std::move(o.law_)) {
  o.posterior_.reset();
}
GaussianPredictiveConditioning &GaussianPredictiveConditioning::operator=(
    GaussianPredictiveConditioning &&o) noexcept {
  if (this != &o) {
    this->~GaussianPredictiveConditioning();
    new (this) GaussianPredictiveConditioning(std::move(o));
  }
  return *this;
}
GaussianPredictiveConditioning GaussianPredictiveConditioning::prepare(
    GaussianPosterior &&p, const Gaussian &noise, std::span<const double> a,
    const PredictiveMetadata &md, PredictivePolicy policy) {
  GaussianPredictiveConditioning out;
  if (p.status() != DensityStatus::finite ||
      noise.status() != DensityStatus::finite) {
    out.law_.status_ =
        p.status() != DensityStatus::finite ? p.status() : noise.status();
    out.law_.numerical_status_ = p.status() != DensityStatus::finite
                                     ? p.numerical_status()
                                     : noise.numerical_status();
    return out;
  }
  const auto peak = preparation_payload_bound(p, noise, md);
  if (!peak || *peak > policy.maximum_payload_bytes) {
    out.law_.numerical_status_ = S::work_limit;
    return out;
  }
  out.law_ = GaussianPredictive::prepare_impl(
      p, {}, p.source().metadata().ordered_ids, noise, a, md, policy, true);
  if (out.law_.status() == DensityStatus::finite) {
    const auto pr = p.retained_payload_bound(),
               lr = out.law_.retained_payload_bound();
    detail::PayloadAccounting b(sizeof(out));
    b.embedded(pr, sizeof(GaussianPosterior));
    b.embedded(lr, sizeof(GaussianPredictive));
    if (!b.result() || *b.result() > policy.maximum_payload_bytes) {
      out.law_.invalidate();
      out.law_.numerical_status_ = S::work_limit;
    } else
      out.posterior_ = std::move(p);
  }
  return out;
}
std::optional<size_t> GaussianPredictiveConditioning::preparation_payload_bound(
    const GaussianPosterior &p, const Gaussian &noise,
    const PredictiveMetadata &md) noexcept {
  detail::PayloadAccounting b(sizeof(GaussianPredictiveConditioning) -
                              sizeof(GaussianPredictive));
  b.add(GaussianPredictive::preparation_payload_bound(p, noise, md)
            .value_or(SIZE_MAX),
        1);
  return b.result();
}
std::optional<size_t>
GaussianPredictiveConditioning::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  if (posterior_)
    b.embedded(posterior_->retained_payload_bound(), sizeof(GaussianPosterior));
  b.embedded(law_.retained_payload_bound(), sizeof(GaussianPredictive));
  return b.result();
}
std::optional<size_t> GaussianPredictiveConditioning::batch_payload_bound(
    size_t count, PredictiveOutputs outputs) const noexcept {
  if (!posterior_)
    return {};
  const auto overhead =
      sizeof(*this) - sizeof(GaussianPosterior) - sizeof(GaussianPredictive);
  return law_.batch_payload_bound(*posterior_, count, outputs, overhead);
}
PredictiveBatch GaussianPredictiveConditioning::evaluate(
    std::span<const double> training, std::span<const std::string> rows,
    std::span<const double> future, std::span<const std::string> future_rows,
    size_t count, PredictiveOutputs outputs, PredictivePolicy policy) const {
  if (!posterior_) {
    PredictiveBatch out;
    out.status = status();
    out.numerical_status = numerical_status();
    return out;
  }
  const auto overhead =
      sizeof(*this) - sizeof(GaussianPosterior) - sizeof(GaussianPredictive);
  return law_.condition_batch(*posterior_, training, rows, future, future_rows,
                              count, outputs, policy, {}, {}, overhead);
}
} // namespace irred::statistics
