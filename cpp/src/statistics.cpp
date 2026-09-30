#include "irred/statistics.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <unordered_set>
namespace irred::statistics {
namespace {
void metadata_payload(detail::PayloadAccounting &b,
                      const Metadata &m) noexcept {
  b.strings(m.ordered_ids);
  for (const auto *s :
       {&m.arithmetic_id, &m.measure, &m.table_identity,
        &m.uncertainty_identity, &m.ordering_provenance,
        &m.calibration_provenance, &m.dependence_provenance,
        &m.source_semantics, &m.input_matrix_convention, &m.treatment})
    b.string(*s);
}
void history_payload(detail::PayloadAccounting &b,
                     const std::vector<SelectionRecord> &h) noexcept {
  b.vector(h);
  for (const auto &r : h) {
    b.string(r.operation);
    b.strings(r.kept_row_ids);
    b.strings(r.complement_row_ids);
  }
}
} // namespace
std::optional<size_t> gaussian_preparation_payload_bound(
    size_t n, MatrixKind kind, numerics::Arithmetic arithmetic,
    const Metadata &metadata) noexcept {
  if ((arithmetic != numerics::Arithmetic::binary64_legacy_v1 &&
       arithmetic != numerics::Arithmetic::longdouble_cpu_v1) ||
      (kind != MatrixKind::covariance && kind != MatrixKind::precision) ||
      (n && n > SIZE_MAX / n))
    return {};
  detail::PayloadAccounting b(sizeof(Gaussian) + 2 * sizeof(numerics::Factorization));
  metadata_payload(b, metadata);
  // Supported libstdc++ validation-set envelope: string-bearing node plus
  // bucket storage <=128 bytes/ID, separately charging copied string storage.
  b.add(n, 128);
  for (const auto &id : metadata.ordered_ids)
    b.string(id);
  b.add(n * n, kind == MatrixKind::precision
                   ? 6 * sizeof(double) + 3 * sizeof(long double)
                   : 4 * sizeof(double) + 2 * sizeof(long double));
  b.add(n, 16 * sizeof(double) + 16 * sizeof(long double));
  b.add(1, 4096); // compiled identity/default metadata growth envelope
  return b.result();
}
std::optional<size_t> selected_gaussian_preparation_payload_bound(
    const observations::Prepared &source, size_t k,
    numerics::Arithmetic arithmetic) noexcept {
  const auto &s = source.source();
  const auto n = s.values.size();
  if (source.status() != observations::Status::ok || !k || k > n ||
      (arithmetic != numerics::Arithmetic::binary64_legacy_v1 &&
       arithmetic != numerics::Arithmetic::longdouble_cpu_v1) ||
      (s.uncertainty != observations::Uncertainty::covariance &&
       s.uncertainty != observations::Uncertainty::precision) ||
      (n && n > SIZE_MAX / n) || k > SIZE_MAX / k)
    return {};
  detail::PayloadAccounting b(2 * sizeof(Gaussian) + 3 * sizeof(numerics::Factorization));
  b.add(n, 4 * sizeof(std::string) + 128 + sizeof(uint8_t));
  b.add(k, sizeof(size_t));
  // Full ordered-ID metadata is initially copied even for selected covariance;
  // retained kept/complement history and validation copies coexist by phase.
  for (const auto &id : s.measurement_ids) {
    if (id.capacity() == SIZE_MAX)
      return {};
    b.add(id.capacity() + 1, 8);
  }
  for (const auto *text :
       {&s.table_sha256, &s.uncertainty_sha256, &s.ordering_provenance,
        &s.calibration_provenance, &s.dependence_provenance}) {
    if (text->capacity() == SIZE_MAX)
      return {};
    b.add(text->capacity() + 1, 3);
  }
  b.add(k * k, 6 * sizeof(double) + 3 * sizeof(long double));
  if (s.uncertainty == observations::Uncertainty::precision)
    // Full inverse/factor remains live during the selected marginal factor.
    b.add(n * n, 6 * sizeof(double) + 3 * sizeof(long double));
  b.add(n, 16 * sizeof(double) + 16 * sizeof(long double));
  b.add(k, 16 * sizeof(double) + 16 * sizeof(long double));
  b.add(1, 12288);
  return b.result();
}
std::optional<size_t>
Gaussian::proper_offset_payload_bound(std::string_view identity) const noexcept {
  const auto retained = retained_payload_bound();
  const auto n = factor_.size();
  const auto preparation = gaussian_preparation_payload_bound(
      n, MatrixKind::covariance, factor_.arithmetic(), metadata_);
  if (!retained || !preparation || (n && n > SIZE_MAX / n) ||
      identity.size() == SIZE_MAX)
    return {};
  detail::PayloadAccounting b(0);
  b.add(*retained, 3); // existing owner plus conservative history/prior growth
  b.add(*preparation, 1);
  b.add(n * n, sizeof(double)); // proper-offset local covariance update
  b.add(n, 8 * sizeof(double) + 2 * sizeof(std::string));
  for (const auto &id : metadata_.ordered_ids) {
    if (id.capacity() == SIZE_MAX)
      return {};
    b.add(id.capacity() + 1, 2);
  }
  b.add(std::max(identity.size(), std::string{}.capacity()) + 1, 6);
  b.add(1, 4096);
  return b.result();
}
std::optional<size_t> Gaussian::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  metadata_payload(b, metadata_);
  b.vector(covariance_);
  b.vector(mean_shift_);
  b.strings(latent_ids_);
  b.vector(priors_);
  for (const auto &r : priors_) {
    b.string(r.latent_identity);
    b.vector(r.response);
    b.strings(r.applied_row_ids);
  }
  history_payload(b, history_);
  b.embedded(factor_.retained_payload_bound(), sizeof(factor_));
  return b.result();
}
std::optional<size_t>
Gaussian::evaluation_payload_bound(size_t row_count, bool profile_mode) const noexcept {
  if (!row_count)
    return size_t{0};
  detail::PayloadAccounting b(sizeof(GaussianResult) + sizeof(ProfileResult) +
                              3 * sizeof(numerics::SolveResult));
  const auto work_size = factor_.arithmetic() == numerics::Arithmetic::longdouble_cpu_v1
                             ? sizeof(long double) : sizeof(double);
  // Density: centered + returned solve + triangular work. Profile additionally
  // retains response solve and residual solve across adjusted-residual solve.
  b.add(factor_.size(), (profile_mode ? 4 : 2) * sizeof(double) + work_size);
  return b.result();
}
std::optional<size_t> ProfileOperator::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  metadata_payload(b, metadata_);
  history_payload(b, history_);
  b.vector(response_);
  b.vector(wx_);
  b.embedded(factor_.retained_payload_bound(), sizeof(factor_));
  return b.result();
}
namespace {
constexpr long double log2pi = 1.8378770664093454835606594728112352797228L;
DensityResult finite(long double v) {
  DensityResult r;
  if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<double>::max()) {
    r.status = DensityStatus::numerical_failure;
    r.numerical_status = numerics::Status::overflow;
    return r;
  }
  r.status = DensityStatus::finite;
  r.numerical_status = numerics::Status::ok;
  r.log_value = static_cast<double>(v);
  return r;
}
bool valid_ids(const Metadata &m) {
  if (m.ordered_ids.empty() || m.measure.empty() ||
      m.ordering_provenance.empty())
    return false;
  std::unordered_set<std::string> s;
  for (const auto &x : m.ordered_ids)
    if (x.empty() || !s.insert(x).second)
      return false;
  return true;
}
// Reported nonzero coefficient/residual/quadratic terms require normal f64.
// A rounded zero is not an exact scientific zero.
bool normal_or_zero_output(long double source, double output) {
  return source == 0 ? output == 0
                     : output != 0 && std::fpclassify(output) == FP_NORMAL;
}
bool finite_span(std::span<const double> x) {
  return std::all_of(x.begin(), x.end(),
                     [](double v) { return std::isfinite(v); });
}
bool indices(std::span<const std::size_t> x, std::size_t n) {
  if (x.empty())
    return false;
  std::unordered_set<std::size_t> s;
  for (auto i : x)
    if (i >= n || !s.insert(i).second)
      return false;
  return true;
}
long double dot(std::span<const double> a, std::span<const double> b) {
  long double v = 0;
  for (std::size_t i = 0; i < a.size(); ++i)
    v += static_cast<long double>(a[i]) * b[i];
  return v;
}
std::vector<double> inverse(const numerics::Factorization &f, double budget,
                            numerics::Status *cause = nullptr) {
  const auto n = f.size();
  std::vector<double> out(n * n), e(n);
  for (std::size_t j = 0; j < n; ++j) {
    std::fill(e.begin(), e.end(), 0);
    e[j] = 1;
    auto s = numerics::solve(f, e, budget);
    if (s.status != numerics::Status::ok) {
      if (cause)
        *cause = s.status;
      return {};
    }
    for (std::size_t i = 0; i < n; ++i)
      out[i * n + j] = s.value[i];
  } // independent column roundoff can break exact symmetry; use one triangle,
    // not source repair.
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < i; ++j)
      out[j * n + i] = out[i * n + j];
  return out;
}
} // namespace
DensityResult normal_log_density(double x, double mean, double sigma) noexcept {
  if (!std::isfinite(x) || !std::isfinite(mean) || !std::isfinite(sigma) ||
      sigma <= 0)
    return {};
  const auto z = (static_cast<long double>(x) - mean) / sigma;
  return finite(-z * z / 2 - std::log(static_cast<long double>(sigma)) -
                log2pi / 2);
}
DensityResult poisson_log_mass(std::uint64_t k, double rate) noexcept {
  if (!std::isfinite(rate) || rate < 0)
    return {};
  if (rate == 0) {
    if (k == 0)
      return finite(0);
    DensityResult r;
    r.status = DensityStatus::outside_support;
    return r;
  }
  if (k > 99) {
    DensityResult r;
    r.status = DensityStatus::unsupported_domain;
    return r;
  }
  auto g = numerics::log_gamma_positive(static_cast<double>(k + 1));
  if (g.status != numerics::Status::ok) {
    DensityResult r;
    r.status = DensityStatus::numerical_failure;
    r.numerical_status = g.status;
    return r;
  }
  return finite(static_cast<long double>(k) *
                    std::log(static_cast<long double>(rate)) -
                rate - g.value);
}
DensityResult selected_standard_normal_positive(double x) noexcept {
  if (!std::isfinite(x))
    return {};
  if (x <= 0) {
    DensityResult r;
    r.status = DensityStatus::outside_support;
    return r;
  }
  auto r = normal_log_density(x, 0, 1);
  if (r.status == DensityStatus::finite)
    return finite(static_cast<long double>(r.log_value) +
                  std::numbers::ln2_v<long double>);
  return r;
}
Gaussian prepare_gaussian(std::span<const double> matrix, MatrixKind kind,
                          Metadata m, std::size_t cap, double budget,
                          numerics::Arithmetic arithmetic) {
  Gaussian g;
  g.metadata_ = std::move(m);
  if (arithmetic != numerics::Arithmetic::binary64_legacy_v1 &&
      arithmetic != numerics::Arithmetic::longdouble_cpu_v1) {
    g.factor_ = numerics::cholesky({}, 0, 0, arithmetic);
    g.metadata_.arithmetic_id = g.factor_.arithmetic_id();
    return g;
  }

  g.metadata_.arithmetic_id =
      arithmetic == numerics::Arithmetic::longdouble_cpu_v1
          ? "F02/longdouble-cpu/v1"
          : "F02/binary64-legacy/v1";
  if (g.metadata_.input_matrix_convention.empty())
    g.metadata_.input_matrix_convention =
        kind == MatrixKind::covariance ? "covariance" : "precision";
  const auto n = g.metadata_.ordered_ids.size();
  if (!valid_ids(g.metadata_)) {
    g.status_ = DensityStatus::incompatible_metadata;
    return g;
  }
  if (kind != MatrixKind::covariance && kind != MatrixKind::precision)
    return g;
  if (!std::isfinite(budget) || budget <= 0)
    return g;
  auto f = numerics::cholesky(matrix, n, cap, arithmetic);
  g.metadata_.arithmetic_id = f.arithmetic_id();
  if (f.status() != numerics::Status::ok) {
    g.status_ = DensityStatus::numerical_failure;
    g.numerical_status_ = f.status();
    return g;
  }
  if (kind == MatrixKind::precision) {
    g.covariance_ = inverse(f, budget, &g.numerical_status_);
    if (g.covariance_.empty()) {
      g.status_ = DensityStatus::numerical_failure;
      return g;
    }
    g.factor_ = numerics::cholesky(g.covariance_, n, cap, arithmetic);
  } else {
    g.covariance_.assign(matrix.begin(), matrix.end());
    g.factor_ = std::move(f);
  }
  g.numerical_status_ = g.factor_.status();
  g.status_ = g.factor_.status() == numerics::Status::ok
                  ? DensityStatus::finite
                  : DensityStatus::numerical_failure;
  return g;
}
GaussianResult Gaussian::evaluate(std::span<const double> r,
                                  std::span<const std::string> ids,
                                  double budget) const {
  GaussianResult out;
  out.density.status = status_;
  out.density.numerical_status = numerical_status_;
  if (status_ != DensityStatus::finite)
    return out;
  if (ids.size() != metadata_.ordered_ids.size() ||
      !std::equal(ids.begin(), ids.end(), metadata_.ordered_ids.begin())) {
    out.density.status = DensityStatus::incompatible_metadata;
    out.density.numerical_status = numerics::Status::invalid_input;
    return out;
  }
  if (r.size() != factor_.size()) {
    out.density.status = DensityStatus::invalid_input;
    out.density.numerical_status = numerics::Status::invalid_input;
    return out;
  }
  if (!std::all_of(r.begin(), r.end(), [](double v) { return std::isfinite(v); })) {
    out.density.status = DensityStatus::numerical_failure;
    out.density.numerical_status = numerics::Status::nonfinite_input;
    return out;
  }
  std::vector<double> centered(r.begin(), r.end());
  if (!mean_shift_.empty())
    for (std::size_t i = 0; i < centered.size(); ++i) {
      const auto v = static_cast<long double>(centered[i]) - mean_shift_[i];
      if (!std::isfinite(v) ||
          std::abs(v) > std::numeric_limits<double>::max()) {
        out.density.status = DensityStatus::numerical_failure;
        out.density.numerical_status = numerics::Status::overflow;
        return out;
      }
      centered[i] = static_cast<double>(v);
    }
  auto s = numerics::solve(factor_, centered, budget);
  if (s.status != numerics::Status::ok) {
    out.density.status = DensityStatus::numerical_failure;
    out.density.numerical_status = s.status;
    return out;
  }
  const auto q = dot(centered, s.value);
  if (q < 0 || !std::isfinite(q) || q > std::numeric_limits<double>::max()) {
    out.density.status = DensityStatus::numerical_failure;
    out.density.numerical_status = q < 0 ? numerics::Status::outside_domain
                                        : numerics::Status::overflow;
    return out;
  }
  const auto reported_q = static_cast<double>(q);
  if (!normal_or_zero_output(q, reported_q)) {
    out.density.status = DensityStatus::numerical_failure;
    out.density.numerical_status = numerics::Status::outside_domain;
    return out;
  }
  out.quadratic = reported_q;
  out.log_determinant = factor_.log_determinant();
  const auto norm = static_cast<long double>(r.size()) * log2pi;
  out.normalization = static_cast<double>(norm);
  out.density = finite(-(q + out.log_determinant + norm) / 2);
  out.backward_residual = s.backward_residual;
  out.estimated_forward_sensitivity = s.estimated_forward_sensitivity;
  return out;
}
Gaussian Gaussian::marginal(std::span<const std::size_t> keep, std::size_t cap,
                            double budget) const {
  if (status_ != DensityStatus::finite || !indices(keep, factor_.size()))
    return {};
  Metadata m = metadata_;
  m.treatment += "; marginal";
  m.ordered_ids.clear();
  for (auto i : keep)
    m.ordered_ids.push_back(metadata_.ordered_ids[i]);
  if (keep.size() > cap / keep.size())
    return {};
  std::vector<double> c;
  for (auto i : keep)
    for (auto j : keep)
      c.push_back(covariance_[i * factor_.size() + j]);
  auto out = prepare_gaussian(c, MatrixKind::covariance, std::move(m), cap,
                              budget, factor_.arithmetic());
  out.latent_ids_ = latent_ids_;
  out.priors_ = priors_;
  out.history_ = history_;
  SelectionRecord h;
  h.operation = "marginal";
  h.kept_row_ids = out.metadata_.ordered_ids;
  for (std::size_t i = 0; i < factor_.size(); ++i)
    if (std::find(keep.begin(), keep.end(), i) == keep.end())
      h.complement_row_ids.push_back(metadata_.ordered_ids[i]);
  out.history_.push_back(std::move(h));
  if (!mean_shift_.empty())
    for (auto i : keep)
      out.mean_shift_.push_back(mean_shift_[i]);
  return out;
}
Gaussian
Gaussian::conditional_zero_complement(std::span<const std::size_t> keep,
                                      std::size_t cap, double budget) const {
  if (status_ != DensityStatus::finite || !indices(keep, factor_.size()))
    return {};
  numerics::Status cause = numerics::Status::invalid_input;
  auto precision = inverse(factor_, budget, &cause);
  if (precision.empty()) {
    Gaussian failed;
    failed.status_ = DensityStatus::numerical_failure;
    failed.numerical_status_ = cause;
    return failed;
  }
  Metadata m = metadata_;
  m.treatment += "; conditional complement CENTERED residual fixed zero";
  m.ordered_ids.clear();
  for (auto i : keep)
    m.ordered_ids.push_back(metadata_.ordered_ids[i]);
  if (keep.size() > cap / keep.size())
    return {};
  std::vector<double> p;
  for (auto i : keep)
    for (auto j : keep)
      p.push_back(precision[i * factor_.size() + j]);
  auto out = prepare_gaussian(p, MatrixKind::precision, std::move(m), cap,
                              budget, factor_.arithmetic());
  out.latent_ids_ = latent_ids_;
  out.priors_ = priors_;
  out.history_ = history_;
  SelectionRecord h;
  h.operation = "conditional complement CENTERED residual fixed zero";
  h.kept_row_ids = out.metadata_.ordered_ids;
  for (std::size_t i = 0; i < factor_.size(); ++i)
    if (std::find(keep.begin(), keep.end(), i) == keep.end())
      h.complement_row_ids.push_back(metadata_.ordered_ids[i]);
  out.history_.push_back(std::move(h));
  if (!mean_shift_.empty())
    for (auto i : keep)
      out.mean_shift_.push_back(mean_shift_[i]);
  return out;
}
namespace {
ProfileResult evaluate_profile(const numerics::Factorization &factor,
                               const Metadata &metadata,
                               std::span<const double> r,
                               std::span<const double> x, long double gram,
                               std::span<const std::string> ids,
                               double budget) {
  ProfileResult out;
  if (ids.size() != metadata.ordered_ids.size() ||
      !std::equal(ids.begin(), ids.end(), metadata.ordered_ids.begin())) {
    out.status = DensityStatus::incompatible_metadata;
    return out;
  }
  if (r.size() != factor.size() || x.size() != r.size() || !finite_span(r) ||
      !finite_span(x) || !(gram > 0) || !std::isfinite(gram))
    return out;
  auto wr = numerics::solve(factor, r, budget);
  if (wr.status != numerics::Status::ok) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = wr.status;
    return out;
  }
  const auto a = dot(x, wr.value) / gram;
  if (!std::isfinite(a) || std::abs(a) > std::numeric_limits<double>::max()) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = numerics::Status::overflow;
    return out;
  }
  const auto reported_a = static_cast<double>(a);
  if (!normal_or_zero_output(a, reported_a)) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = numerics::Status::outside_domain;
    return out;
  }
  std::vector<double> adjusted(r.size());
  for (std::size_t i = 0; i < r.size(); ++i) {
    const auto v = static_cast<long double>(r[i]) - a * x[i];
    if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<double>::max()) {
      out.status = DensityStatus::numerical_failure;
      out.numerical_status = numerics::Status::overflow;
      return out;
    }
    adjusted[i] = static_cast<double>(v);
    if (!normal_or_zero_output(v, adjusted[i])) {
      out.status = DensityStatus::numerical_failure;
      out.numerical_status = numerics::Status::outside_domain;
      return out;
    }
  }
  auto wa = numerics::solve(factor, adjusted, budget);
  if (wa.status != numerics::Status::ok) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = wa.status;
    return out;
  }
  const auto q = dot(adjusted, wa.value);
  if (q < 0 || !std::isfinite(q) || q > std::numeric_limits<double>::max()) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = numerics::Status::overflow;
    return out;
  }
  const auto reported_q = static_cast<double>(q);
  if (!normal_or_zero_output(q, reported_q)) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = numerics::Status::outside_domain;
    return out;
  }
  long double residual_l1 = 0, adjusted_l1 = 0;
  for (std::size_t i = 0; i < r.size(); ++i) {
    residual_l1 += std::abs(static_cast<long double>(r[i]));
    adjusted_l1 += std::abs(static_cast<long double>(adjusted[i]));
    out.solution_norm_inf =
        std::max(out.solution_norm_inf, std::abs(wr.value[i]));
    out.adjusted_solution_norm_inf =
        std::max(out.adjusted_solution_norm_inf, std::abs(wa.value[i]));
  }
  if (residual_l1 > std::numeric_limits<double>::max() ||
      adjusted_l1 > std::numeric_limits<double>::max()) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = numerics::Status::overflow;
    return out;
  }
  out.status = DensityStatus::finite;
  out.numerical_status = numerics::Status::ok;
  out.coefficient = reported_a;
  out.quadratic = reported_q;
  out.backward_residual = wa.backward_residual;
  out.estimated_forward_sensitivity = wa.estimated_forward_sensitivity;
  out.coefficient_solve_backward_residual = wr.backward_residual;
  out.coefficient_solve_forward_sensitivity = wr.estimated_forward_sensitivity;
  out.residual_l1 = static_cast<double>(residual_l1);
  out.adjusted_residual_l1 = static_cast<double>(adjusted_l1);
  out.adjusted_residuals = std::move(adjusted);
  return out;
}
} // namespace
ProfileResult Gaussian::profile_offset(std::span<const double> r,
                                       std::span<const double> x,
                                       std::span<const std::string> ids,
                                       double budget) const {
  ProfileResult out;
  if (status_ != DensityStatus::finite) {
    out.status = status_;
    return out;
  }
  if (!priors_.empty() || ids.size() != metadata_.ordered_ids.size() ||
      !std::equal(ids.begin(), ids.end(), metadata_.ordered_ids.begin())) {
    out.status = DensityStatus::incompatible_metadata;
    return out;
  }
  if (x.size() != factor_.size() || !finite_span(x))
    return out;
  auto wx = numerics::solve(factor_, x, budget);
  if (wx.status != numerics::Status::ok) {
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = wx.status;
    return out;
  }
  return evaluate_profile(factor_, metadata_, r, x, dot(x, wx.value), ids,
                          budget);
}
ProfileOperator
Gaussian::prepare_offset_profile(std::span<const double> x,
                                 std::span<const std::string> ids,
                                 double budget) && {
  ProfileOperator out;
  if (status_ != DensityStatus::finite) {
    out.status_ = status_;
    out.numerical_status_ = numerical_status_;
    return out;
  }
  if (!priors_.empty() || ids.size() != metadata_.ordered_ids.size() ||
      !std::equal(ids.begin(), ids.end(), metadata_.ordered_ids.begin())) {
    out.status_ = DensityStatus::incompatible_metadata;
    return out;
  }
  if (x.size() != factor_.size() || !finite_span(x))
    return out;
  auto wx = numerics::solve(factor_, x, budget);
  if (wx.status != numerics::Status::ok) {
    out.status_ = DensityStatus::numerical_failure;
    out.numerical_status_ = wx.status;
    return out;
  }
  const auto gram = dot(x, wx.value);
  if (!(gram > 0) || !std::isfinite(gram))
    return out;
  out.response_.assign(x.begin(), x.end());
  out.wx_ = std::move(wx.value);
  out.gram_ = gram;
  out.backward_ = wx.backward_residual;
  out.sensitivity_ = wx.estimated_forward_sensitivity;
  status_ = DensityStatus::invalid_input; // consumed even if later metadata
                                          // allocation throws
  out.factor_ = std::move(factor_);
  out.metadata_ = std::move(metadata_);
  out.history_ = std::move(history_);
  out.metadata_.treatment +=
      "; retained fixed offset profile score (not density)";
  out.status_ = DensityStatus::finite;
  out.numerical_status_ = numerics::Status::ok;
  status_ = DensityStatus::invalid_input;
  covariance_.clear();
  return out;
}
ProfileResult ProfileOperator::evaluate(std::span<const double> r,
                                        std::span<const std::string> ids,
                                        double budget) const {
  if (status_ != DensityStatus::finite) {
    ProfileResult out;
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!std::isfinite(budget) || budget <= 0) {
    return {};
  }
  if (sensitivity_ > budget) {
    ProfileResult out;
    out.status = DensityStatus::numerical_failure;
    out.numerical_status = numerics::Status::conditioning_budget_exceeded;
    return out;
  }
  return evaluate_profile(factor_, metadata_, r, response_, gram_, ids, budget);
}
Gaussian Gaussian::proper_offset(std::span<const double> x,
                                 std::span<const std::string> ids, double mean,
                                 double variance, std::string identity,
                                 bool independent, std::size_t cap,
                                 double budget) const {
  if (ids.size() != metadata_.ordered_ids.size() ||
      !std::equal(ids.begin(), ids.end(), metadata_.ordered_ids.begin())) {
    Gaussian failed;
    failed.status_ = DensityStatus::incompatible_metadata;
    return failed;
  }
  if (status_ != DensityStatus::finite || x.size() != factor_.size() ||
      !finite_span(x) || !std::isfinite(mean) || !std::isfinite(variance) ||
      variance <= 0 || !independent || identity.empty() ||
      std::find(latent_ids_.begin(), latent_ids_.end(), identity) !=
          latent_ids_.end())
    return {};
  auto c = covariance_;
  for (std::size_t i = 0; i < x.size(); ++i)
    for (std::size_t j = 0; j <= i; ++j) {
      auto v = static_cast<long double>(c[i * x.size() + j]) +
               static_cast<long double>(variance) * x[i] * x[j];
      if (!std::isfinite(v) ||
          std::abs(v) > std::numeric_limits<double>::max()) {
        Gaussian failed;
        failed.status_ = DensityStatus::numerical_failure;
        return failed;
      }
      c[i * x.size() + j] = c[j * x.size() + i] = static_cast<double>(v);
    }
  auto m = metadata_;
  m.treatment += "; proper independent offset:" + identity;
  auto out = prepare_gaussian(c, MatrixKind::covariance, std::move(m), cap,
                              budget, factor_.arithmetic());
  out.latent_ids_ = latent_ids_;
  out.latent_ids_.push_back(identity);
  out.priors_ = priors_;
  out.priors_.push_back(PriorRecord{mean, variance, std::move(identity),
                                    std::vector<double>(x.begin(), x.end()),
                                    metadata_.ordered_ids, true});
  out.history_ = history_;
  out.mean_shift_ = mean_shift_;
  if (out.mean_shift_.empty())
    out.mean_shift_.assign(x.size(), 0);
  for (std::size_t i = 0; i < x.size(); ++i) {
    auto v = static_cast<long double>(out.mean_shift_[i]) +
             static_cast<long double>(mean) * x[i];
    if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<double>::max()) {
      Gaussian failed;
      failed.status_ = DensityStatus::numerical_failure;
      return failed;
    }
    out.mean_shift_[i] = static_cast<double>(v);
  }
  return out;
}
Gaussian
prepare_selected_observations(const observations::Prepared &p,
                              std::span<const std::size_t> selected_indices,
                              std::size_t cap, double budget,
                              numerics::Arithmetic arithmetic) {
  if (arithmetic != numerics::Arithmetic::binary64_legacy_v1 &&
      arithmetic != numerics::Arithmetic::longdouble_cpu_v1)
    return prepare_gaussian({}, MatrixKind::covariance, {}, cap, budget,
                            arithmetic);
  if (p.status() != observations::Status::ok)
    return {};
  const auto &s = p.source();
  const auto n = s.values.size();
  if (selected_indices.empty() || selected_indices.size() > n)
    return {};
  const auto k_limit = selected_indices.size();
  if (s.uncertainty == observations::Uncertainty::covariance &&
      k_limit > cap / k_limit)
    return {};
  if (s.uncertainty == observations::Uncertainty::precision && n > cap / n)
    return {};

  observations::SelectionResult selected;
  selected.status = observations::Status::ok;
  selected.mask.assign(n, 0);
  selected.source_indices.assign(selected_indices.begin(),
                                 selected_indices.end());
  std::size_t previous = 0;
  bool first = true;
  for (auto i : selected_indices) {
    if (i >= n || (!first && i <= previous) ||
        (!s.source_selection.empty() && !s.source_selection[i]) ||
        s.missing[i] || !std::isfinite(s.values[i]))
      return {};
    if (s.profile == observations::Profile::pantheon_plus_released_v1 &&
        (s.zhd_missing[i] || s.zhel_missing[i] || !std::isfinite(s.zhd[i]) ||
         !std::isfinite(s.zhel[i])))
      return {};
    selected.mask[i] = 1;
    previous = i;
    first = false;
  }
  if (s.uncertainty == observations::Uncertainty::none)
    return {};
  Metadata m{};
  m.ordered_ids = s.measurement_ids;
  m.measure = s.unit == observations::Unit::magnitude ? "product d(magnitude)"
                                                      : "product d(metre)";
  m.table_identity = s.table_sha256;
  m.uncertainty_identity = s.uncertainty_sha256;
  m.ordering_provenance = s.ordering_provenance;
  m.calibration_provenance = s.calibration_provenance;
  m.dependence_provenance = s.dependence_provenance;
  m.source_semantics =
      "profile=" + std::to_string(static_cast<unsigned>(s.profile)) +
      ";role=" + std::to_string(static_cast<unsigned>(s.role)) +
      ";unit=" + std::to_string(static_cast<unsigned>(s.unit)) +
      ";calibration=" + std::to_string(static_cast<unsigned>(s.calibration)) +
      ";component=" + std::to_string(static_cast<unsigned>(s.component)) +
      ";uncertainty_kind=" +
      std::to_string(static_cast<unsigned>(s.uncertainty)) +
      ";uncertainty_unit=" +
      std::to_string(static_cast<unsigned>(s.uncertainty_unit));
  if (s.uncertainty == observations::Uncertainty::covariance) {
    const auto k = selected.source_indices.size(), n = s.values.size();
    if (k > cap / k)
      return {};
    m.ordered_ids.clear();
    for (auto i : selected.source_indices)
      m.ordered_ids.push_back(s.measurement_ids[i]);
    m.matrix_validation_scope = MatrixValidationScope::selected_covariance_only;
    std::vector<double> block;
    block.reserve(k * k);
    for (auto i : selected.source_indices)
      for (auto j : selected.source_indices)
        block.push_back(s.uncertainty_matrix[i * n + j]);
    auto out = prepare_gaussian(block, MatrixKind::covariance, std::move(m),
                                cap, budget, arithmetic);
    SelectionRecord history;
    history.operation = "caller-declared selected principal covariance block; "
                        "full source probability validity not assessed";
    history.kept_row_ids = out.metadata_.ordered_ids;
    for (std::size_t i = 0; i < n; ++i)
      if (!selected.mask[i])
        history.complement_row_ids.push_back(s.measurement_ids[i]);
    out.history_.push_back(std::move(history));
    return out;
  }
  m.matrix_validation_scope =
      MatrixValidationScope::full_precision_then_marginal;
  auto full = prepare_gaussian(s.uncertainty_matrix, MatrixKind::precision,
                               std::move(m), cap, budget, arithmetic);
  if (full.status() != DensityStatus::finite)
    return full;
  auto out = full.marginal(selected.source_indices, cap, budget);
  if (!out.history_.empty())
    out.history_.back().operation = "caller-declared source-order selection "
                                    "after validated full precision marginal";
  return out;
}
Gaussian prepare_observations(const observations::Prepared &p,
                              observations::Selection choice, std::size_t cap,
                              double budget, numerics::Arithmetic arithmetic) {
  if (arithmetic != numerics::Arithmetic::binary64_legacy_v1 &&
      arithmetic != numerics::Arithmetic::longdouble_cpu_v1)
    return prepare_gaussian({}, MatrixKind::covariance, {}, cap, budget,
                            arithmetic);
  if (p.status() != observations::Status::ok)
    return {};
  auto selected = p.select(choice);
  if (selected.status != observations::Status::ok)
    return {};
  auto out = prepare_selected_observations(p, selected.source_indices, cap,
                                           budget, arithmetic);
  if (!out.history_.empty())
    out.history_.back().operation +=
        "; compiled selection=" + std::to_string(static_cast<unsigned>(choice));
  return out;
}
} // namespace irred::statistics
