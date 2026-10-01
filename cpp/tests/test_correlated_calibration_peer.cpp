// Independent synthetic normalized-density controls; no copied software/assets.
// Explicit 3x3 cofactors and 2D Simpson integration avoid production Cholesky.
// Proper latent beta=m+Lz with explicit 2x2 LDL/square-root transform.
// Cofactor q/logdet/logdensity budget2e-12*(1+|ref|); density integral1e-8.
// Independent refinement and Gaussian-tail omissions each <=5e-10 density.
#include "irred/correlated_calibration.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred::statistics;
using W = long double;
namespace {
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void near(W a, W b, W tolerance = 2e-12L) {
  need(std::abs(a - b) <= tolerance * (1 + std::abs(b)),
       "independent numerical allocation");
}
using M = std::array<W, 9>;
using V = std::array<W, 3>;
struct Density {
  W q, ld, logp;
};
Density cofactor(M c, V r) {
  W a = c[0], b = c[1], d = c[2], e = c[4], f = c[5], g = c[8];
  M adj{e * g - f * f, d * f - b * g, b * f - d * e,
        d * f - b * g, a * g - d * d, b * d - a * f,
        b * f - d * e, b * d - a * f, a * e - b * b};
  W det = a * adj[0] + b * adj[1] + d * adj[2];
  if (!(det > 0))
    throw std::runtime_error("cofactor determinant positive");
  W q = 0;
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      q += r[i] * adj[i * 3 + j] * r[j] / det;
  W ld = std::log(det);
  return {q, ld, -(q + ld + 3 * std::log(2 * std::numbers::pi_v<W>)) / 2};
}
struct Fixture {
  std::vector<double> c{2, .25, -.1, .25, 3, .2, -.1, .2, 4};
  std::vector<double> x{1, .2, -.3, 1, .5, -.4}, y{.4, -.2, 1.1};
  std::vector<std::string> ids{"a", "b", "c"};
  CalibrationPrior p;
};
Fixture fixture() {
  Fixture f;
  f.p.ordered_parameter_ids = {"cal1", "cal2"};
  f.p.parameter_units = {"u", "u"};
  f.p.mean = {.1, -.2};
  f.p.covariance = {.5, .125, .125, .75};
  f.p.prior_identity = "independent synthetic correlated proper prior";
  f.p.response_identity = "explicit row-major synthetic response";
  f.p.residual_unit = "u";
  f.p.calibration_identity = "single prior calibration";
  f.p.dependence_identity = "noise and latent explicitly independent";
  f.p.measure_identity = "product d(u)";
  f.p.noise_independence_declared = true;
  return f;
}
Metadata meta(const Fixture &f) {
  Metadata m;
  m.ordered_ids = f.ids;
  m.measure = f.p.measure_identity;
  m.table_identity = "synthetic data";
  m.uncertainty_identity = "synthetic base covariance conditional on beta";
  m.calibration_provenance = f.p.calibration_identity;
  m.dependence_provenance = f.p.dependence_identity;
  m.ordering_provenance = "explicit rows";
  m.source_semantics = "synthetic controls";
  return m;
}
Gaussian gaussian(const Fixture &f) {
  return prepare_gaussian(f.c, MatrixKind::covariance, meta(f), 10000, 1e-10);
}
Density expected(const Fixture &f) {
  M sigma;
  V centered;
  for (size_t i = 0; i < 3; ++i) {
    centered[i] = f.y[i];
    for (size_t a = 0; a < 2; ++a)
      centered[i] -= W(f.x[i * 2 + a]) * f.p.mean[a];
    for (size_t j = 0; j < 3; ++j) {
      sigma[i * 3 + j] = f.c[i * 3 + j];
      for (size_t a = 0; a < 2; ++a)
        for (size_t b = 0; b < 2; ++b)
          sigma[i * 3 + j] +=
              W(f.x[i * 2 + a]) * f.p.covariance[a * 2 + b] * f.x[j * 2 + b];
    }
  }
  return cofactor(sigma, centered);
}
W prior_integral(const Fixture &f, unsigned panels) {
  M c;
  std::copy(f.c.begin(), f.c.end(), c.begin());
  W l00 = std::sqrt(W(f.p.covariance[0])), l10 = f.p.covariance[2] / l00,
    l11 = std::sqrt(W(f.p.covariance[3]) - l10 * l10);
  const W h = 18.L / panels;
  W sum = 0;
  for (unsigned i = 0; i <= panels; ++i) {
    W z0 = -9 + i * h, w0 = i == 0 || i == panels ? 1 : i % 2 ? 4 : 2;
    for (unsigned j = 0; j <= panels; ++j) {
      W z1 = -9 + j * h, w1 = j == 0 || j == panels ? 1 : j % 2 ? 4 : 2;
      W beta0 = f.p.mean[0] + l00 * z0,
        beta1 = f.p.mean[1] + l10 * z0 + l11 * z1;
      V r;
      for (size_t k = 0; k < 3; ++k)
        r[k] = W(f.y[k]) - f.x[k * 2] * beta0 - f.x[k * 2 + 1] * beta1;
      sum += w0 * w1 * std::exp(cofactor(c, r).logp - (z0 * z0 + z1 * z1) / 2) /
             (2 * std::numbers::pi_v<W>);
    }
  }
  return sum * h * h / 9;
}
GaussianResult compare(const Fixture &f, bool integrate) {
  auto g = gaussian(f);
  auto op = CorrelatedCalibration::prepare(std::move(g), f.x, f.ids, f.p);
  need(op.status() == DensityStatus::finite,
       "proper correlated owner prepared");
  need(g.status() == DensityStatus::invalid_input,
       "successful source consumed once");
  auto got = op.evaluate(f.y, f.ids);
  need(got.density.status == DensityStatus::finite,
       "normalized observed density");
  auto ref = expected(f);
  near(got.quadratic, ref.q);
  near(got.log_determinant, ref.ld);
  near(got.density.log_value, ref.logp);
  near(got.normalization, 3 * std::log(2 * std::numbers::pi_v<W>));
  if (integrate) {
    W fine = prior_integral(f, 256), coarse = prior_integral(f, 128);
    need(std::abs(fine - coarse) <= 5e-10L, "independent 2D prior refinement");
    M c;
    std::copy(f.c.begin(), f.c.end(), c.begin());
    W max_conditional = std::exp(cofactor(c, {0, 0, 0}).logp);
    W tail = 2 * std::erfc(9 / std::sqrt(2.L)) * max_conditional;
    need(tail <= 5e-10L, "analytic union-bound prior tail allocation");
    need(std::abs(fine - std::exp(ref.logp)) <= 1e-8L,
         "normalized integral cofactor agreement");
    need(std::abs(fine - std::exp(W(got.density.log_value))) <= 1e-8L,
         "normalized integral product agreement");
    std::cout << std::setprecision(20) << "integration fine=" << fine
              << " coarse=" << coarse
              << " refinement=" << std::abs(fine - coarse)
              << " tail_bound=" << tail << " cofactor=" << std::exp(ref.logp)
              << "\n";
  }
  auto wrong = f.ids;
  std::swap(wrong[0], wrong[1]);
  need(op.evaluate(f.y, wrong).density.status ==
           DensityStatus::incompatible_metadata,
       "ordered rows exact");
  auto moved = std::move(op);
  need(op.status() == DensityStatus::invalid_input &&
           moved.status() == DensityStatus::finite,
       "move owner lifetime");
  auto *same_owner = &moved;
  moved = std::move(*same_owner);
  need(moved.evaluate(f.y, f.ids).density.status == DensityStatus::finite,
       "self move retains prior owner");
  return got;
}
} // namespace
int main() {
  try {
    need(std::numeric_limits<W>::digits >= 64,
         "wide cofactor/reference arithmetic");
    auto f = fixture();
    auto baseline = compare(f, true);
    auto perm = f;
    std::swap(perm.p.mean[0], perm.p.mean[1]);
    std::swap(perm.p.ordered_parameter_ids[0], perm.p.ordered_parameter_ids[1]);
    std::swap(perm.p.parameter_units[0], perm.p.parameter_units[1]);
    std::swap(perm.p.covariance[0], perm.p.covariance[3]);
    for (size_t i = 0; i < 3; ++i)
      std::swap(perm.x[i * 2], perm.x[i * 2 + 1]);
    near(compare(perm, false).density.log_value, baseline.density.log_value);
    perm = f;
    std::array<size_t, 3> order{2, 0, 1};
    for (size_t i = 0; i < 3; ++i) {
      perm.ids[i] = f.ids[order[i]];
      perm.y[i] = f.y[order[i]];
      for (size_t a = 0; a < 2; ++a)
        perm.x[i * 2 + a] = f.x[order[i] * 2 + a];
      for (size_t j = 0; j < 3; ++j)
        perm.c[i * 3 + j] = f.c[order[i] * 3 + order[j]];
    }
    near(compare(perm, false).density.log_value, baseline.density.log_value);
    auto scale = f;
    for (auto &v : scale.y)
      v *= 7;
    for (auto &v : scale.x)
      v *= 7;
    for (auto &v : scale.c)
      v *= 49;
    scale.p.residual_unit = "u/7";
    scale.p.measure_identity = "product d(u/7)";
    near(compare(scale, false).density.log_value,
         baseline.density.log_value - 3 * std::log(7.L));
    scale = f;
    for (size_t i = 0; i < 3; ++i)
      scale.x[i * 2] /= 8;
    scale.p.mean[0] *= 8;
    scale.p.covariance[0] *= 64;
    scale.p.covariance[1] *= 8;
    scale.p.covariance[2] *= 8;
    scale.p.parameter_units[0] = "u/8";
    near(compare(scale, false).density.log_value, baseline.density.log_value);
    auto null = f;
    std::fill(null.x.begin(), null.x.end(), 0);
    compare(null, true);
    auto rank = f;
    for (size_t i = 0; i < 3; ++i)
      rank.x[i * 2 + 1] = 2 * rank.x[i * 2];
    compare(rank, true);
    for (int mode = 0; mode < 6; ++mode) {
      auto bad = f;
      if (mode == 0)
        bad.p.covariance = {1, 1, 1, 1};
      if (mode == 1)
        bad.p.covariance = {1, 2, 2, 1};
      if (mode == 2)
        bad.p.noise_independence_declared = false;
      if (mode == 3)
        bad.p.ordered_parameter_ids[1] = bad.p.ordered_parameter_ids[0];
      if (mode == 4)
        bad.p.covariance[1] = .2;
      if (mode == 5)
        bad.x[0] = std::numeric_limits<double>::quiet_NaN();
      auto g = gaussian(bad);
      auto failed =
          CorrelatedCalibration::prepare(std::move(g), bad.x, bad.ids, bad.p);
      need(failed.status() != DensityStatus::finite &&
               g.status() == DensityStatus::finite,
           "hostile prior/response rejected source preserved");
    }
    auto g = gaussian(f);
    auto prior_bearing =
        g.proper_offset(std::array<double, 3>{1, 1, 1}, f.ids, 0, 1,
                        "already applied prior", true, 10000, 1e-10);
    auto reject = CorrelatedCalibration::prepare(std::move(prior_bearing), f.x,
                                                 f.ids, f.p);
    need(reject.status() != DensityStatus::finite &&
             prior_bearing.status() == DensityStatus::finite,
         "prior double counting rejected without consuming");
    CalibrationPolicy quota;
    quota.maximum_payload_bytes = 0;
    reject =
        CorrelatedCalibration::prepare(std::move(g), f.x, f.ids, f.p, quota);
    need(reject.status() != DensityStatus::finite &&
             g.status() == DensityStatus::finite,
         "preparation quota preserves source");
    auto op = CorrelatedCalibration::prepare(std::move(g), f.x, f.ids, f.p);
    need(op.evaluate(f.y, f.ids, quota).density.status != DensityStatus::finite,
         "evaluation quota withholds density");
    auto bady = f.y;
    bady[1] = std::numeric_limits<double>::infinity();
    need(op.evaluate(bady, f.ids).density.status != DensityStatus::finite,
         "nonfinite observation rejected");
    auto rounding = std::fegetround();
    std::fesetround(FE_UPWARD);
    auto failedround = op.evaluate(f.y, f.ids);
    std::fesetround(rounding);
    need(failedround.density.status == DensityStatus::unsupported_domain,
         "rounding gate prior density");
    std::cout << "PASS " << checks
              << " correlated calibration independent peer controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
