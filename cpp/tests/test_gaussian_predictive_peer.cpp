// Original covariance-route conditioning T=C+XSX^T, F=ASX^T,
// U=R+ASA^T, b=Am+F T^-1(r-Xm), W=U-F T^-1 F^T.
// Explicit cofactors/Fraction facts and prior*noise quadrature were frozen
// before reading production. No precision Gram or posterior V is an oracle.
// Common physical/statistical equations, long-double/libm inputs disclosed.
// Mean/cov2e-12*(1+|ref|), log/q/logdet2e-11*(1+|ref|), mass2e-11;
// quadrature refinement<=5% allocation. Bounds are empirical, not universal.
#include "irred/gaussian_predictive.hpp"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
namespace {
namespace s = irred::statistics;
using W = long double;
using M = std::vector<W>;
unsigned checks = 0;
W largest_error = 0;
void need(bool b, const char *message) {
  ++checks;
  if (!b)
    throw std::runtime_error(message);
}
void near(W a, W b, W allowance, const char *message) {
  W err = std::abs(a - b);
  largest_error = std::max(largest_error, err);
  need(err <= allowance * (1 + std::abs(b)), message);
}
M mul(const M &a, const M &b, size_t n, size_t k, size_t p) {
  M r(n * p);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < p; ++j)
      for (size_t z = 0; z < k; ++z)
        r[i * p + j] += a[i * k + z] * b[z * p + j];
  return r;
}
M tr(const M &a, size_t n, size_t p) {
  M r(n * p);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < p; ++j)
      r[j * n + i] = a[i * p + j];
  return r;
}
M add(M a, const M &b, W sign = 1) {
  for (size_t i = 0; i < a.size(); ++i)
    a[i] += sign * b[i];
  return a;
}
W det(const M &a) { return a[0] * a[3] - a[1] * a[2]; }
M inv(const M &a) {
  W d = det(a);
  need(d > 0, "independent cofactor positive determinant");
  return {a[3] / d, -a[1] / d, -a[2] / d, a[0] / d};
}
M wide(const std::vector<double> &a) { return {a.begin(), a.end()}; }
W lognormal(const M &y, const M &mean, const M &C) {
  M z = add(y, mean, -1), inverse = inv(C);
  W q = mul(tr(z, 2, 1), mul(inverse, z, 2, 2, 1), 1, 2, 1)[0];
  return -.5L *
         (q + std::log(det(C)) + 2 * std::log(2 * std::numbers::pi_v<W>));
}
struct Case {
  std::vector<double> C{1, .25, .25, 2}, S{2, .5, .5, 1}, X{1, 0, 1, 1},
      m{.5, -.25}, r{1.25, -.5}, A{1, 2, 1, -1}, R{.5, .125, .125, .75};
};
struct Ref {
  M b, Wcov, T, cross;
};
Ref reference(const Case &c) {
  M C = wide(c.C), S = wide(c.S), X = wide(c.X), A = wide(c.A), R = wide(c.R),
    m = wide(c.m), r = wide(c.r);
  M T = add(C, mul(mul(X, S, 2, 2, 2), tr(X, 2, 2), 2, 2, 2));
  M F = mul(mul(A, S, 2, 2, 2), tr(X, 2, 2), 2, 2, 2),
    gain = mul(F, inv(T), 2, 2, 2);
  M U = add(R, mul(mul(A, S, 2, 2, 2), tr(A, 2, 2), 2, 2, 2));
  return {add(mul(A, m, 2, 2, 1),
              mul(gain, add(r, mul(X, m, 2, 2, 1), -1), 2, 2, 1)),
          add(U, mul(gain, tr(F, 2, 2), 2, 2, 2), -1), T, F};
}
std::vector<std::string> trainids{"train-cal0", "train-cal1"},
    futureids{"heldout-cal2", "heldout-cal3"};
s::Metadata noise_metadata(bool future) {
  s::Metadata m;
  m.ordered_ids = future ? futureids : trainids;
  m.measure = "product d(mag)";
  m.table_identity = future ? "original synthetic future noise"
                            : "original synthetic training noise";
  m.uncertainty_identity = "explicit full SPD covariance";
  m.ordering_provenance = "declared axes";
  m.calibration_provenance = "synthetic fixed calibration";
  m.dependence_provenance = "full internally correlated noise";
  m.source_semantics = "synthetic controls";
  return m;
}
s::ParameterPrior prior(const Case &c) {
  s::ParameterPrior p;
  p.ordered_parameter_ids = {"offset", "slope"};
  p.parameter_units = {"mag", "mag"};
  p.shared_nuisance_ids = {"offset"};
  p.mean = c.m;
  p.covariance = c.S;
  p.prior_identity = "original synthetic proper prior";
  p.design_identity = "fixed synthetic linear training design";
  p.residual_unit = "mag";
  p.parameter_measure = "product d(mag)";
  p.dependence_identity = "prior independent of training noise";
  p.noise_independence_declared = true;
  return p;
}
s::PredictiveMetadata metadata() {
  s::PredictiveMetadata m;
  m.ordered_parameter_ids = {"offset", "slope"};
  m.parameter_units = {"mag", "mag"};
  m.response_units = {"mag/mag", "mag/mag"};
  m.training_event_ids = {"synthetic-train-event0", "synthetic-train-event1"};
  m.future_event_ids = {"synthetic-future-event2", "synthetic-future-event3"};
  m.future_unit = "mag";
  m.future_covariance_unit = "mag^2";
  m.future_measure = "product d(mag)";
  m.response_identity = "fixed future linear response";
  m.conditioning_identity = "fixed training vector";
  m.conditional_noise_identity = "future conditional independent full noise";
  m.dependence_identity =
      "future noise independent of training noise and prior";
  m.future_noise_independence_declared = true;
  m.noise_conditional_on_parameters_declared = true;
  return m;
}
s::GaussianPredictive native(const Case &c) {
  auto g = s::prepare_gaussian(c.C, s::MatrixKind::covariance,
                               noise_metadata(false), 1024, 1e-10,
                               irred::numerics::Arithmetic::longdouble_cpu_v1);
  auto p = s::GaussianPosterior::prepare(std::move(g), c.X, trainids, prior(c));
  auto noise = s::prepare_gaussian(
      c.R, s::MatrixKind::covariance, noise_metadata(true), 1024, 1e-10,
      irred::numerics::Arithmetic::longdouble_cpu_v1);
  auto out =
      s::GaussianPredictive::prepare(p, c.r, trainids, noise, c.A, metadata());
  need(p.status() == s::DensityStatus::finite &&
           noise.status() == s::DensityStatus::finite,
       "borrowed inputs remain usable");
  return out;
}
void check_case(const Case &c) {
  auto ref = reference(c);
  auto got = native(c);
  need(got.status() == s::DensityStatus::finite && got.mean().size() == 2 &&
           got.covariance().size() == 4,
       "predictive prepared");
  for (size_t i = 0; i < 2; ++i)
    near(got.mean()[i], ref.b[i], 2e-12L, "covariance-route predictive mean");
  for (size_t i = 0; i < 4; ++i)
    near(got.covariance()[i], ref.Wcov[i], 2e-12L,
         "covariance-route predictive covariance");
  for (const auto &y :
       std::array<std::vector<double>, 3>{{{.2, -.75}, {0, 0}, {2, -3}}}) {
    auto d = got.log_density(y, futureids);
    need(d.density.status == s::DensityStatus::finite,
         "joint future density admitted");
    near(d.density.log_value, lognormal(wide(y), ref.b, ref.Wcov), 2e-11L,
         "cofactor joint future density");
    M z = add(wide(y), ref.b, -1);
    W q = mul(tr(z, 2, 1), mul(inv(ref.Wcov), z, 2, 2, 1), 1, 2, 1)[0];
    near(d.quadratic, q, 2e-11L, "independent predictive quadratic");
    near(d.log_determinant, std::log(det(ref.Wcov)), 2e-11L,
         "independent predictive determinant");
  }
  auto perm = futureids;
  std::swap(perm[0], perm[1]);
  need(got.log_density(std::array{0., 0.}, perm).density.status !=
           s::DensityStatus::finite,
       "future order mismatch refuses");
  need(got.training_values().size() == 2 && got.response().size() == 4 &&
           got.metadata().future_event_ids == metadata().future_event_ids,
       "retained fixed training and response lineage");
}
struct Rule {
  std::array<W, 8> x, w;
  Rule() {
    for (unsigned i = 0; i < 8; ++i) {
      W z = std::cos(std::numbers::pi_v<W> * (i + .75L) / 8.5L), d = 0;
      for (unsigned j = 0; j < 64; ++j) {
        W a = 1, b = z;
        for (unsigned k = 2; k <= 8; ++k) {
          W next = ((2 * k - 1) * z * b - (k - 1) * a) / k;
          a = b;
          b = next;
        }
        d = 8 * (z * b - a) / (z * z - 1);
        W delta = b / d;
        z -= delta;
        if (std::abs(delta) < 4e-19L)
          break;
      }
      x[i] = z;
      w[i] = 2 / ((1 - z * z) * d * d);
    }
  }
};
struct Product {
  W evidence = 0, future = 0;
};
Product direct_product(const Case &c, unsigned panels, const M &y) {
  static const Rule rule;
  W width = 32.L / panels;
  Product out;
  M S = wide(c.S), C = wide(c.C), R = wide(c.R), X = wide(c.X), A = wide(c.A),
    r = wide(c.r), m = wide(c.m);
  auto iS = inv(S), iC = inv(C), iR = inv(R);
  W normalS = 1 / (2 * std::numbers::pi_v<W> * std::sqrt(det(S))),
    normalC = 1 / (2 * std::numbers::pi_v<W> * std::sqrt(det(C))),
    normalR = 1 / (2 * std::numbers::pi_v<W> * std::sqrt(det(R)));
  auto density = [](W x, W y, const M &inverse, W normal) {
    return normal * std::exp(-.5L * (x * x * inverse[0] +
                                     x * y * (inverse[1] + inverse[2]) +
                                     y * y * inverse[3]));
  };
  for (unsigned i = 0; i < panels; ++i)
    for (unsigned j = 0; j < 8; ++j) {
      W b0 = m[0] - 16 + (i + .5L) * width + width / 2 * rule.x[j];
      for (unsigned k = 0; k < panels; ++k)
        for (unsigned l = 0; l < 8; ++l) {
          W b1 = m[1] - 16 + (k + .5L) * width + width / 2 * rule.x[l],
            weight = width * width / 4 * rule.w[j] * rule.w[l];
          W d = weight * density(b0 - m[0], b1 - m[1], iS, normalS) *
                density(r[0] - X[0] * b0 - X[1] * b1,
                        r[1] - X[2] * b0 - X[3] * b1, iC, normalC);
          out.evidence += d;
          out.future += d * density(y[0] - A[0] * b0 - A[1] * b1,
                                    y[1] - A[2] * b0 - A[3] * b1, iR, normalR);
        }
    }
  return out;
}
void quadratures() {
  Case c;
  auto ref = reference(c);
  M y{.2L, -.75L};
  auto coarse = direct_product(c, 64, y), fine = direct_product(c, 128, y);
  need(fine.evidence > .001L,
       "evidence lower floor for declared quadrature tail bound");
  W a = std::log(coarse.future / coarse.evidence),
    b = std::log(fine.future / fine.evidence);
  near(a, b, .05L * 2e-11L, "prior-noise product refinement");
  near(b, lognormal(y, ref.b, ref.Wcov), 2e-11L,
       "direct prior-training-future noise integral");
  // Prior-mean+-16 box: marginal sigma<=sqrt2 gives prior tail<3e-29;
  // likelihood upper<=.115 and evidence>.001 yield omitted conditional
  // mass<4e-27.
  auto owner = native(c);
  Rule rule;
  W previous = 0;
  for (unsigned panels : {32u, 64u}) {
    W mass = 0, width = 48.L / panels;
    for (unsigned i = 0; i < panels; ++i)
      for (unsigned j = 0; j < 8; ++j)
        for (unsigned k = 0; k < panels; ++k)
          for (unsigned l = 0; l < 8; ++l) {
            std::array<double, 2> values{
                double(ref.b[0] - 24 + (i + .5L) * width +
                       width / 2 * rule.x[j]),
                double(ref.b[1] - 24 + (k + .5L) * width +
                       width / 2 * rule.x[l])};
            auto d = owner.log_density(values, futureids);
            need(d.density.status == s::DensityStatus::finite,
                 "retained native joint mass quadrature");
            mass += width * width / 4 * rule.w[j] * rule.w[l] *
                    std::exp(W(d.density.log_value));
          }
    need(std::abs(mass - 1) <= 2e-11L, "normalized joint predictive mass");
    if (panels == 64)
      need(std::abs(mass - previous) <= 1e-12L,
           "native future mass refinement");
    previous = mass;
    std::cout << "MASS panels=" << panels << " value=" << mass << '\n';
  }
  // Largest marginal variance<3.4: +/-24 omitted future rectangle mass<1e-37.
  std::cout << "PRODUCT evidence=" << fine.evidence << " log_density=" << b
            << " refinement=" << std::abs(a - b) << '\n';
}
} // namespace
int main() {
  try {
    std::cout << std::setprecision(21);
    Case c;
    auto r = reference(c);
    near(r.b[0], -129.L / 334, 1e-18L, "exact Fraction predictive mean0");
    near(r.b[1], 855.L / 668, 1e-18L, "exact Fraction predictive mean1");
    for (size_t i = 0; i < 4; ++i)
      near(r.Wcov[i],
           std::array<W, 4>{4444.L / 1336, -773.L / 1336, -773.L / 1336,
                            2618.L / 1336}[i],
           1e-18L, "exact Fraction predictive covariance");
    check_case(c);
    c.X = {1, 2, 2, 4};
    check_case(c);
    c.X = {0, 0, 0, 0};
    check_case(c);
    c = Case{};
    c.A = {0, 0, 0, 0};
    check_case(c);
    c.A = {1, 2, 1, 2};
    check_case(c);
    c = Case{};
    std::swap(c.A[0], c.A[2]);
    std::swap(c.A[1], c.A[3]);
    std::swap(c.R[0], c.R[3]);
    check_case(c);
    quadratures();
    std::cout << "PASS " << checks
              << " independent predictive controls; maximum absolute "
                 "comparison error="
              << largest_error << '\n';
  } catch (const std::exception &e) {
    std::cerr << "FAIL after " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
