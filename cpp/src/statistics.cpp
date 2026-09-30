#include "irred/statistics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <unordered_set>
namespace irred::statistics {
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
std::vector<double> inverse(const numerics::Factorization &f, double budget) {
  const auto n = f.size();
  std::vector<double> out(n * n), e(n);
  for (std::size_t j = 0; j < n; ++j) {
    std::fill(e.begin(), e.end(), 0);
    e[j] = 1;
    auto s = numerics::solve(f, e, budget);
    if (s.status != numerics::Status::ok)
      return {};
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
                          Metadata m, std::size_t cap, double budget) {
  Gaussian g;
  g.metadata_ = std::move(m);
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
  auto f = numerics::cholesky(matrix, n, cap);
  if (f.status() != numerics::Status::ok) {
    g.status_ = DensityStatus::numerical_failure;
    return g;
  }
  if (!std::isfinite(budget) || budget <= 0)
    return g;
  if (kind == MatrixKind::precision) {
    g.covariance_ = inverse(f, budget);
    if (g.covariance_.empty()) {
      g.status_ = DensityStatus::numerical_failure;
      return g;
    }
    g.factor_ = numerics::cholesky(g.covariance_, n, cap);
  } else {
    g.covariance_.assign(matrix.begin(), matrix.end());
    g.factor_ = std::move(f);
  }
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
  if (status_ != DensityStatus::finite)
    return out;
  if (ids.size() != metadata_.ordered_ids.size() ||
      !std::equal(ids.begin(), ids.end(), metadata_.ordered_ids.begin())) {
    out.density.status = DensityStatus::incompatible_metadata;
    return out;
  }
  std::vector<double> centered(r.begin(), r.end());
  if (centered.size() != factor_.size()) {
    out.density.status = DensityStatus::invalid_input;
    return out;
  }
  if (!mean_shift_.empty())
    for (std::size_t i = 0; i < centered.size(); ++i) {
      const auto v = static_cast<long double>(centered[i]) - mean_shift_[i];
      if (!std::isfinite(v) ||
          std::abs(v) > std::numeric_limits<double>::max()) {
        out.density.status = DensityStatus::numerical_failure;
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
    return out;
  }
  out.quadratic = static_cast<double>(q);
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
  auto out =
      prepare_gaussian(c, MatrixKind::covariance, std::move(m), cap, budget);
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
  auto precision = inverse(factor_, budget);
  if (precision.empty()) {
    Gaussian failed;
    failed.status_ = DensityStatus::numerical_failure;
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
  auto out =
      prepare_gaussian(p, MatrixKind::precision, std::move(m), cap, budget);
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
ProfileResult Gaussian::profile_offset(std::span<const double> r,
                                       std::span<const double> x,
                                       std::span<const std::string> ids,
                                       double budget) const {
  ProfileResult out;
  if (ids.size() != metadata_.ordered_ids.size() ||
      !std::equal(ids.begin(), ids.end(), metadata_.ordered_ids.begin())) {
    out.status = DensityStatus::incompatible_metadata;
    return out;
  }
  if (status_ != DensityStatus::finite || r.size() != factor_.size() ||
      x.size() != r.size() || !finite_span(r) || !finite_span(x))
    return out;
  std::vector<double> centered(r.begin(), r.end());
  if (!mean_shift_.empty())
    for (std::size_t i = 0; i < centered.size(); ++i) {
      auto v = static_cast<long double>(centered[i]) - mean_shift_[i];
      if (!std::isfinite(v) ||
          std::abs(v) > std::numeric_limits<double>::max()) {
        out.status = DensityStatus::numerical_failure;
        return out;
      }
      centered[i] = static_cast<double>(v);
    }
  auto wr = numerics::solve(factor_, centered, budget),
       wx = numerics::solve(factor_, x, budget);
  if (wr.status != numerics::Status::ok || wx.status != numerics::Status::ok) {
    out.status = DensityStatus::numerical_failure;
    return out;
  }
  const auto gram = dot(x, wx.value);
  if (!(gram > 0))
    return out;
  const auto a = dot(x, wr.value) / gram;
  if (!std::isfinite(a) || std::abs(a) > std::numeric_limits<double>::max()) {
    out.status = DensityStatus::numerical_failure;
    return out;
  }
  std::vector<double> adjusted(r.size());
  for (std::size_t i = 0; i < r.size(); ++i) {
    auto v = static_cast<long double>(centered[i]) - a * x[i];
    if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<double>::max()) {
      out.status = DensityStatus::numerical_failure;
      return out;
    }
    adjusted[i] = static_cast<double>(v);
  }
  auto wa = numerics::solve(factor_, adjusted, budget);
  if (wa.status != numerics::Status::ok) {
    out.status = DensityStatus::numerical_failure;
    return out;
  }
  auto q = dot(adjusted, wa.value);
  if (q < 0 || !std::isfinite(q) || q > std::numeric_limits<double>::max()) {
    out.status = DensityStatus::numerical_failure;
    return out;
  }
  out.status = DensityStatus::finite;
  out.coefficient = static_cast<double>(a);
  out.quadratic = static_cast<double>(q);
  return out;
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
  auto out =
      prepare_gaussian(c, MatrixKind::covariance, std::move(m), cap, budget);
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
Gaussian prepare_observations(const observations::Prepared &p,
                              observations::Selection choice, std::size_t cap,
                              double budget) {
  if (p.status() != observations::Status::ok)
    return {};
  auto selected = p.select(choice);
  if (selected.status != observations::Status::ok ||
      selected.source_indices.empty())
    return {};
  const auto &s = p.source();
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
                                cap, budget);
    SelectionRecord history;
    history.operation = "selected principal covariance block; full source probability "
                        "validity not assessed";
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
                               std::move(m), cap, budget);
  if (full.status() != DensityStatus::finite)
    return full;
  return full.marginal(selected.source_indices, cap, budget);
}
} // namespace irred::statistics
