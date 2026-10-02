// Original independent calculations: covariance conditioning through
// T=C+X S X^T with explicit2x2/3x3 cofactors, not precision Gram/whitening.
// Exact rational rank-one facts were derived before reading implementation:
// mu=(767/1342,-511/2684), V=(622/671,-529/1342;-529/1342,171/671).
// Mean/cov2e-12*(1+|ref|), normalized logdensity2e-11*(1+|ref|).
// Direct prior*likelihood product GL8 quadrature64/128 panels checks posterior
// normalization/moments with <=5% reference-refinement allocation. Prior
// mean+-16 rectangle: sigma<=sqrt2, omitted prior mass<3e-29; likelihood
// upper normalization and evidence>0.009 imply posterior tail<1e-26.
// Refinement/error diagnostics are empirical estimates, not rigorous bounds.
// No external implementation, owner fixture or production solver is an oracle.
#include "irred/gaussian_design.hpp"
#include "irred/gaussian_posterior.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred;
namespace {
using W = long double;
using M = std::vector<W>;
using DS = statistics::DensityStatus;
unsigned checks = 0;
W largest_mean_error = 0, largest_covariance_error = 0,
  largest_density_error = 0;
void need(bool ok, const char *message) {
  ++checks;
  if (!ok)
    throw std::runtime_error(message);
}
void near(W got, W ref, W allocation, const char *message) {
  need(std::abs(got - ref) <= allocation * (1 + std::abs(ref)), message);
}
M wide(std::span<const double> x) { return M(x.begin(), x.end()); }
M multiply(const M &a, const M &b, size_t n, size_t k, size_t p) {
  M c(n * p);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < p; ++j)
      for (size_t t = 0; t < k; ++t)
        c[i * p + j] += a[i * k + t] * b[t * p + j];
  return c;
}
M transpose(const M &a, size_t n, size_t p) {
  M b(n * p);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < p; ++j)
      b[j * n + i] = a[i * p + j];
  return b;
}
struct Inverse {
  M value;
  W determinant;
};
Inverse inverse(const M &a, size_t n) {
  need(n == 2 || n == 3, "reference cofactor scope");
  M b(n * n);
  W det;
  if (n == 2) {
    det = a[0] * a[3] - a[1] * a[2];
    b = {a[3], -a[1], -a[2], a[0]};
  } else {
    b = {a[4] * a[8] - a[5] * a[7], a[2] * a[7] - a[1] * a[8],
         a[1] * a[5] - a[2] * a[4], a[5] * a[6] - a[3] * a[8],
         a[0] * a[8] - a[2] * a[6], a[2] * a[3] - a[0] * a[5],
         a[3] * a[7] - a[4] * a[6], a[1] * a[6] - a[0] * a[7],
         a[0] * a[4] - a[1] * a[3]};
    det = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
  }
  need(det > 0, "reference positive determinant");
  for (auto &x : b)
    x /= det;
  return {b, det};
}
W lognormal(std::span<const W> x, std::span<const W> mean, const M &cov) {
  size_t p = x.size();
  auto inv = inverse(cov, p);
  W q = 0;
  for (size_t i = 0; i < p; ++i)
    for (size_t j = 0; j < p; ++j)
      q += (x[i] - mean[i]) * inv.value[i * p + j] * (x[j] - mean[j]);
  return -.5L * (q + std::log(inv.determinant) +
                 p * std::log(2 * std::numbers::pi_v<W>));
}
struct Case {
  size_t n, p;
  std::vector<double> C, S, X, m, r;
};
Case rational() {
  return {
      2,          2,          {1, .25, .25, 2}, {2, .5, .5, 1}, {1, 2, 2, 4},
      {.5, -.25}, {1.25, -.5}};
}
Case full() {
  return {3,
          2,
          {2, .25, -.125, .25, 1, .125, -.125, .125, 1.5},
          {1, .25, .25, 2},
          {1, 0, 1, 1, 0, 1},
          {.5, -.25},
          {1, .25, -.5}};
}
struct Reference {
  M mean, covariance, T;
};
Reference reference(const Case &c) {
  M S = wide(c.S), X = wide(c.X), Xt = transpose(X, c.n, c.p), m = wide(c.m),
    r = wide(c.r);
  M sx = multiply(S, Xt, c.p, c.p, c.n), T = multiply(X, sx, c.n, c.p, c.n);
  for (size_t j = 0; j < T.size(); ++j)
    T[j] += c.C[j];
  M K = multiply(sx, inverse(T, c.n).value, c.p, c.n, c.n);
  auto predicted = multiply(X, m, c.n, c.p, 1);
  for (size_t j = 0; j < c.n; ++j)
    r[j] -= predicted[j];
  auto update = multiply(K, r, c.p, c.n, 1);
  for (size_t j = 0; j < c.p; ++j)
    m[j] += update[j];
  auto subtract = multiply(multiply(K, X, c.p, c.n, c.p), S, c.p, c.p, c.p);
  for (size_t j = 0; j < S.size(); ++j)
    S[j] -= subtract[j];
  return {m, S, T};
}
std::vector<std::string> rows(size_t n) {
  std::vector<std::string> r;
  for (size_t j = 0; j < n; ++j)
    r.push_back("row" + std::to_string(j));
  return r;
}
statistics::ParameterPrior prior(const Case &c) {
  statistics::ParameterPrior p;
  for (size_t j = 0; j < c.p; ++j) {
    p.ordered_parameter_ids.push_back("beta" + std::to_string(j));
    p.parameter_units.push_back("declared-coordinate" + std::to_string(j));
  }
  p.shared_nuisance_ids = {p.ordered_parameter_ids.back()};
  p.mean = c.m;
  p.covariance = c.S;
  p.prior_identity = "original synthetic proper prior";
  p.design_identity = "explicit fixed linear design";
  p.residual_unit = "synthetic residual coordinate";
  p.parameter_measure =
      "product of declared parameter-coordinate Lebesgue measures";
  p.dependence_identity =
      "prior and supplied Gaussian observation noise independent";
  p.noise_independence_declared = true;
  return p;
}
statistics::Gaussian source(const Case &c) {
  statistics::Metadata md;
  md.ordered_ids = rows(c.n);
  md.measure = "product of observation residual coordinates";
  md.table_identity = "original synthetic residual source";
  md.uncertainty_identity = "explicit full SPD covariance";
  md.ordering_provenance = "fixed row order";
  md.calibration_provenance = "synthetic fixed calibration";
  md.dependence_provenance = "full supplied covariance";
  return statistics::prepare_gaussian(c.C, statistics::MatrixKind::covariance,
                                      std::move(md), 1024, 1e-10,
                                      numerics::Arithmetic::longdouble_cpu_v1);
}
statistics::GaussianPosterior owner(const Case &c) {
  auto g = source(c);
  return statistics::GaussianPosterior::prepare(std::move(g), c.X, rows(c.n),
                                                prior(c));
}
W check_case(const Case &c, std::span<const double> beta) {
  auto ref = reference(c);
  auto got = owner(c);
  need(got.status() == DS::finite, "proper posterior prepared");
  auto mean = got.condition(c.r, rows(c.n));
  need(mean.status == DS::finite && mean.value.size() == c.p, "mean accepted");
  for (size_t j = 0; j < c.p; ++j) {
    largest_mean_error =
        std::max(largest_mean_error, std::abs(W(mean.value[j]) - ref.mean[j]));
    near(mean.value[j], ref.mean[j], 2e-12,
         "independent covariance-route conditional mean");
    need(std::isfinite(mean.absolute_error_estimates[j]) &&
             mean.absolute_error_estimates[j] >= 0,
         "mean numerical diagnostic separate from variance");
  }
  for (size_t j = 0; j < c.p * c.p; ++j) {
    largest_covariance_error =
        std::max(largest_covariance_error,
                 std::abs(W(got.covariance()[j]) - ref.covariance[j]));
    near(got.covariance()[j], ref.covariance[j], 2e-12,
         "independent covariance-route V");
  }
  auto density =
      got.log_density(c.r, rows(c.n), beta, prior(c).ordered_parameter_ids);
  need(density.density.status == DS::finite,
       "normalized parameter density accepted");
  W logp = lognormal(wide(beta), ref.mean, ref.covariance);
  largest_density_error = std::max(
      largest_density_error, std::abs(W(density.density.log_value) - logp));
  near(density.density.log_value, logp, 2e-11,
       "independent cofactor parameter density");
  near(density.log_determinant,
       std::log(inverse(ref.covariance, c.p).determinant), 2e-11,
       "posterior determinant normalization");
  near(density.normalization, c.p * std::log(2 * std::numbers::pi_v<W>), 2e-11,
       "parameter Lebesgue normalization");
  return density.density.log_value;
}
void exact_and_axes() {
  auto c = rational();
  auto ref = reference(c);
  const W exactmean[]{767.L / 1342, -511.L / 2684},
      exactcov[]{622.L / 671, -529.L / 1342, -529.L / 1342, 171.L / 671};
  for (unsigned j = 0; j < 2; ++j)
    near(ref.mean[j], exactmean[j], 1e-18, "exact rational mean reference");
  for (unsigned j = 0; j < 4; ++j)
    near(ref.covariance[j], exactcov[j], 1e-18,
         "exact rational covariance reference");
  const double b[]{.75, -.125};
  check_case(c, b);
  c = full();
  W before = check_case(c, b);
  // Reorder every observed/design/covariance axis together; the parameter
  // measure is unchanged, so its normalized density is invariant.
  auto perm = c;
  const size_t order[]{2, 0, 1};
  for (size_t i = 0; i < 3; ++i) {
    perm.r[i] = c.r[order[i]];
    for (size_t j = 0; j < 2; ++j)
      perm.X[i * 2 + j] = c.X[order[i] * 2 + j];
    for (size_t j = 0; j < 3; ++j)
      perm.C[i * 3 + j] = c.C[order[i] * 3 + order[j]];
  }
  near(check_case(perm, b), before, 2e-11, "row permutation density");
  perm = c;
  std::swap(perm.m[0], perm.m[1]);
  perm.S = {c.S[3], c.S[2], c.S[1], c.S[0]};
  for (size_t i = 0; i < 3; ++i)
    std::swap(perm.X[i * 2], perm.X[i * 2 + 1]);
  const double swapped[]{b[1], b[0]};
  near(check_case(perm, swapped), before, 2e-11,
       "parameter permutation density");
  // Explicit new coordinates beta'=D beta. Unit labels do not perform this
  // conversion; X'=X D^-1, m'=D m, S'=D S D^T, dbeta'=|det D|dbeta.
  perm = c;
  const double scale[]{-2, 4};
  double changed[2];
  for (size_t i = 0; i < 2; ++i) {
    perm.m[i] *= scale[i];
    changed[i] = b[i] * scale[i];
    for (size_t j = 0; j < 2; ++j)
      perm.S[i * 2 + j] *= scale[i] * scale[j];
    for (size_t j = 0; j < 3; ++j)
      perm.X[j * 2 + i] /= scale[i];
  }
  near(check_case(perm, changed), before - std::log(8.L), 2e-11,
       "explicit unit Jacobian");
  // Observation-unit scaling has its own Jacobian in evidence; the parameter
  // posterior measure is unchanged.
  perm = c;
  for (auto &x : perm.C)
    x *= 4;
  for (auto &x : perm.X)
    x *= 2;
  for (auto &x : perm.r)
    x *= 2;
  near(check_case(perm, b), before, 2e-11,
       "observation units cancel from posterior");
  c = rational();
  std::fill(c.X.begin(), c.X.end(), 0);
  check_case(c, b);
  auto null = owner(c);
  auto mean = null.condition(c.r, rows(2));
  for (unsigned j = 0; j < 2; ++j)
    near(mean.value[j], c.m[j], 2e-12, "null design returns proper prior mean");
  // More parameters than rows: proper prior retains the unobserved directions.
  c = {2,
       3,
       {1, .25, .25, 2},
       {2, .25, .125, .25, 1, .25, .125, .25, 1.5},
       {1, 0, 1, 0, 1, 1},
       {.5, -.25, .125},
       {1.25, -.5}};
  const double b3[]{.75, -.125, .25};
  check_case(c, b3);
}
struct Node {
  W x, w;
};
std::vector<Node> nodes(W center, unsigned panels) {
  constexpr W x[]{
      .183434642495649804939476142360184L, .525532409916328985817739049189246L,
      .796666477413626739591553936475831L, .960289856497536231683560868569473L},
      w[]{.362683783378361982965150449277195L,
          .313706645877887287337962201986601L,
          .222381034453374470544355994426240L,
          .101228536290376259152531354309962L};
  std::vector<Node> v;
  W h = 16.L / panels;
  for (unsigned i = 0; i < panels; ++i)
    for (unsigned j = 0; j < 4; ++j) {
      W middle = center - 16 + (2 * i + 1) * h;
      v.push_back({middle - h * x[j], h * w[j]});
      v.push_back({middle + h * x[j], h * w[j]});
    }
  return v;
}
struct Integrated {
  W evidence, mean0, mean1, cov00, cov01, cov11, mass;
};
Integrated joint(const Case &c, unsigned panels) {
  auto ref = reference(c);
  auto ci = inverse(wide(c.C), 2), si = inverse(wide(c.S), 2);
  const W norm = 1 / (4 * std::numbers::pi_v<W> * std::numbers::pi_v<W> *
                      std::sqrt(ci.determinant * si.determinant));
  Integrated out{};
  for (auto a : nodes(c.m[0], panels))
    for (auto b : nodes(c.m[1], panels)) {
      W beta[]{a.x, b.x},
          e[]{W(c.r[0]) - c.X[0] * beta[0] - c.X[1] * beta[1],
              W(c.r[1]) - c.X[2] * beta[0] - c.X[3] * beta[1]},
          d[]{beta[0] - c.m[0], beta[1] - c.m[1]}, q = 0;
      for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j)
          q += d[i] * si.value[i * 2 + j] * d[j] +
               e[i] * ci.value[i * 2 + j] * e[j];
      W weight = a.w * b.w * norm * std::exp(-q / 2);
      out.evidence += weight;
      out.mean0 += weight * beta[0];
      out.mean1 += weight * beta[1];
      out.cov00 += weight * beta[0] * beta[0];
      out.cov01 += weight * beta[0] * beta[1];
      out.cov11 += weight * beta[1] * beta[1];
    }
  out.mean0 /= out.evidence;
  out.mean1 /= out.evidence;
  out.cov00 = out.cov00 / out.evidence - out.mean0 * out.mean0;
  out.cov01 = out.cov01 / out.evidence - out.mean0 * out.mean1;
  out.cov11 = out.cov11 / out.evidence - out.mean1 * out.mean1;
  out.mass = out.evidence /
             std::exp(lognormal(
                 wide(c.r), multiply(wide(c.X), wide(c.m), 2, 2, 1), ref.T));
  return out;
}
void integration() {
  auto c = rational();
  auto coarse = joint(c, 64), fine = joint(c, 128);
  need(fine.evidence > .009,
       "positive evidence bounds omitted tail amplification");
  // Exact rational observation-space determinant671/16 and quadratic1211/671
  // independently pin evidence; no production Gaussian is involved.
  const W exact_evidence = std::exp(-1211.L / 1342) /
                           (2 * std::numbers::pi_v<W> * std::sqrt(671.L / 16));
  near(fine.evidence, exact_evidence, .05L * 2e-11,
       "direct evidence against exact rational marginal");
  near(fine.mass, 1, 2e-11, "direct joint density normalizes proper posterior");
  near(fine.evidence, coarse.evidence, .05L * 2e-11,
       "direct evidence refinement");
  const auto ref = reference(c);
  const W values[]{fine.mean0, fine.mean1, fine.cov00, fine.cov01, fine.cov11},
      old[]{coarse.mean0, coarse.mean1, coarse.cov00, coarse.cov01,
            coarse.cov11},
      expected[]{ref.mean[0], ref.mean[1], ref.covariance[0], ref.covariance[1],
                 ref.covariance[3]};
  for (unsigned j = 0; j < 5; ++j) {
    near(values[j], old[j], .05L * 2e-12,
         "independent direct joint moment refinement");
    near(values[j], expected[j], 2e-12, "integrated conditional moment");
  }
  // At fixed parameter points use direct prior*likelihood/evidence, without a
  // posterior precision or Gaussian library result on the reference side.
  auto got = owner(c);
  const M zeros(2, 0);
  for (const auto &b :
       std::array<std::array<double, 2>, 3>{{{0, 0}, {.75, -.125}, {-1, 1}}}) {
    M predicted = multiply(wide(c.X), wide(b), 2, 2, 1);
    M e = wide(c.r);
    for (unsigned j = 0; j < 2; ++j)
      e[j] -= predicted[j];
    W expected = lognormal(wide(b), wide(c.m), wide(c.S)) +
                 lognormal(e, zeros, wide(c.C)) - std::log(fine.evidence);
    auto density =
        got.log_density(c.r, rows(2), b, prior(c).ordered_parameter_ids);
    need(density.density.status == DS::finite,
         "direct joint density point admitted");
    near(density.density.log_value, expected, 2e-11,
         "normalized Bayes density with direct joint integration");
  }
  std::cout << "independent quadrature evidence=" << double(fine.evidence)
            << " mass=" << double(fine.mass) << '\n';
}
void boundaries() {
  auto c = rational();
  auto got = owner(c);
  const double beta[]{.75, -.125};
  auto ids = rows(2);
  std::swap(ids[0], ids[1]);
  need(got.condition(c.r, ids).status == DS::incompatible_metadata,
       "exact row identity refusal");
  auto pids = prior(c).ordered_parameter_ids;
  std::swap(pids[0], pids[1]);
  need(got.log_density(c.r, rows(2), beta, pids).density.status ==
           DS::incompatible_metadata,
       "exact parameter identity refusal");
  for (unsigned mode = 0; mode < 5; ++mode) {
    auto g = source(c);
    auto p = prior(c);
    statistics::PosteriorPolicy policy;
    if (mode == 0)
      p.noise_independence_declared = false;
    if (mode == 1)
      p.covariance = {1, 1, 1, 1};
    if (mode == 2)
      policy.maximum_elements = 1;
    if (mode == 3)
      policy.maximum_forward_sensitivity = 1e-30;
    if (mode == 4)
      p.parameter_measure.clear();
    auto bad = statistics::GaussianPosterior::prepare(
        std::move(g), c.X, rows(2), std::move(p), policy);
    need(bad.status() != DS::finite && g.status() == DS::finite,
         "failed proper-prior admission preserves source");
  }
  // A relative profile cannot normalize a parameter law in unconstrained
  // directions; the same null design remains lawful under a proper prior.
  std::fill(c.X.begin(), c.X.end(), 0);
  auto g = source(c);
  statistics::DesignMetadata dm;
  dm.ordered_parameter_ids = prior(c).ordered_parameter_ids;
  dm.parameter_units = prior(c).parameter_units;
  dm.residual_unit = "synthetic";
  dm.design_identity = "null relative design";
  dm.dependence_identity = "supplied full covariance";
  auto relative =
      statistics::DesignProfile::prepare(std::move(g), c.X, rows(2), dm);
  need(relative.status() != DS::finite,
       "null design has no unique relative profile");
  need(owner(c).status() == DS::finite,
       "null design proper posterior boundary");
  c = rational();
  c.m = {1e16, -.5e16};
  auto cancel = owner(c);
  need(cancel.status() == DS::finite,
       "cancellation source structurally lawful");
  auto refused =
      cancel.log_density(c.r, rows(2), c.m, prior(c).ordered_parameter_ids);
  need(refused.density.status != DS::finite,
       "unresolved mean cancellation cannot produce accepted normalized "
       "density");
  auto rounding = std::fegetround();
  need(std::fesetround(FE_DOWNWARD) == 0, "set hostile rounding");
  auto unsupported = owner(rational());
  need(std::fesetround(rounding) == 0, "restore rounding");
  need(unsupported.status() != DS::finite,
       "unsupported rounding refuses posterior");
}
} // namespace
int main() {
  try {
    exact_and_axes();
    integration();
    boundaries();
    std::cout << "PASS " << checks
              << " independent Gaussian posterior controls; maximum absolute "
                 "mean/covariance/logdensity errors "
              << double(largest_mean_error) << ' '
              << double(largest_covariance_error) << ' '
              << double(largest_density_error) << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
