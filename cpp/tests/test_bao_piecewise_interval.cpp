// Optional test infrastructure only: GMP exact rational canonical C/r and
// MPFR directed endpoint arithmetic. No production dependency. Original mode
// requires external SHA-before/after guard; ordinary controls use no assets.
// Frozen widths2e-9, factor3e-9, background5e-9, full1e-8, sensitivity1e-10.
#include "fixtures/bao_reference.hpp"
#include "irred/bao.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <gmpxx.h>
#include <mpfr.h>
#include <sstream>
#include <stdexcept>
using namespace irred;
namespace {
mpfr_prec_t precision = 256;
mpq_class exact(double x) {
  mpq_class q(x);
  q.canonicalize();
  return q;
}
unsigned checks = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
struct F {
  mpfr_t v;
  F() {
    mpfr_init2(v, precision);
    mpfr_set_zero(v, 0);
  }
  F(const F &o) {
    mpfr_init2(v, mpfr_get_prec(o.v));
    mpfr_set(v, o.v, MPFR_RNDN);
  }
  F &operator=(const F &o) {
    if (this != &o) {
      mpfr_set_prec(v, mpfr_get_prec(o.v));
      mpfr_set(v, o.v, MPFR_RNDN);
    }
    return *this;
  }
  ~F() { mpfr_clear(v); }
};
struct I {
  F lo, hi;
  I() = default;
  explicit I(const mpq_class &q) {
    mpfr_set_q(lo.v, q.get_mpq_t(), MPFR_RNDD);
    mpfr_set_q(hi.v, q.get_mpq_t(), MPFR_RNDU);
  }
  explicit I(double x) : I(exact(x)) {}
  explicit I(long x) : I(mpq_class(x)) {}
};
I operator+(const I &a, const I &b) {
  I c;
  mpfr_add(c.lo.v, a.lo.v, b.lo.v, MPFR_RNDD);
  mpfr_add(c.hi.v, a.hi.v, b.hi.v, MPFR_RNDU);
  return c;
}
I operator-(const I &a, const I &b) {
  I c;
  mpfr_sub(c.lo.v, a.lo.v, b.hi.v, MPFR_RNDD);
  mpfr_sub(c.hi.v, a.hi.v, b.lo.v, MPFR_RNDU);
  return c;
}
I operator-(const I &a) {
  I c;
  mpfr_neg(c.lo.v, a.hi.v, MPFR_RNDD);
  mpfr_neg(c.hi.v, a.lo.v, MPFR_RNDU);
  return c;
}
I corners(const I &a, const I &b, bool divide) {
  if (divide)
    check(mpfr_sgn(b.lo.v) > 0 || mpfr_sgn(b.hi.v) < 0,
          "interval denominator excludes zero");
  I c;
  F low, high;
  const F *aa[] = {&a.lo, &a.hi};
  const F *bb[] = {&b.lo, &b.hi};
  bool first = true;
  for (auto x : aa)
    for (auto y : bb) {
      if (divide) {
        mpfr_div(low.v, x->v, y->v, MPFR_RNDD);
        mpfr_div(high.v, x->v, y->v, MPFR_RNDU);
      } else {
        mpfr_mul(low.v, x->v, y->v, MPFR_RNDD);
        mpfr_mul(high.v, x->v, y->v, MPFR_RNDU);
      }
      if (first || mpfr_cmp(low.v, c.lo.v) < 0)
        mpfr_set(c.lo.v, low.v, MPFR_RNDD);
      if (first || mpfr_cmp(high.v, c.hi.v) > 0)
        mpfr_set(c.hi.v, high.v, MPFR_RNDU);
      first = false;
    }
  return c;
}
I operator*(const I &a, const I &b) { return corners(a, b, false); }
I operator/(const I &a, const I &b) { return corners(a, b, true); }
I log(const I &a) {
  check(mpfr_sgn(a.lo.v) > 0, "positive log domain");
  I c;
  mpfr_log(c.lo.v, a.lo.v, MPFR_RNDD);
  mpfr_log(c.hi.v, a.hi.v, MPFR_RNDU);
  return c;
}
I exp(const I &a) {
  I c;
  mpfr_exp(c.lo.v, a.lo.v, MPFR_RNDD);
  mpfr_exp(c.hi.v, a.hi.v, MPFR_RNDU);
  return c;
}
I expm1(const I &a) {
  I c;
  mpfr_expm1(c.lo.v, a.lo.v, MPFR_RNDD);
  mpfr_expm1(c.hi.v, a.hi.v, MPFR_RNDU);
  return c;
}
I cbrt(const I &a) {
  check(mpfr_sgn(a.lo.v) >= 0, "nonnegative cube-root domain");
  I c;
  mpfr_cbrt(c.lo.v, a.lo.v, MPFR_RNDD);
  mpfr_cbrt(c.hi.v, a.hi.v, MPFR_RNDU);
  return c;
}
F width(const I &a) {
  F x;
  mpfr_sub(x.v, a.hi.v, a.lo.v, MPFR_RNDU);
  return x;
}
F absmax(const I &a) {
  F x, y;
  mpfr_abs(x.v, a.lo.v, MPFR_RNDU);
  mpfr_abs(y.v, a.hi.v, MPFR_RNDU);
  if (mpfr_cmp(y.v, x.v) > 0)
    x = y;
  return x;
}
bool within(const F &a, const char *budget) {
  mpq_class q(budget);
  return mpfr_cmp_q(a.v, q.get_mpq_t()) <= 0;
}
bool contains(const I &a, const I &b) {
  return mpfr_cmp(a.lo.v, b.lo.v) <= 0 && mpfr_cmp(a.hi.v, b.hi.v) >= 0;
}
std::string text(const F &a) {
  char s[256];
  mpfr_snprintf(s, sizeof(s), "%.110RUg", a.v);
  return s;
}
struct Matrix {
  size_t n;
  std::vector<mpq_class> inverse;
  mpq_class determinant = 1;
  explicit Matrix(const std::vector<double> &c)
      : n(static_cast<size_t>(std::sqrt(c.size()))) {
    check(n * n == c.size(), "square canonical matrix");
    std::vector<mpq_class> a(n * 2 * n);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j) {
        check(c[i * n + j] == c[j * n + i], "exact symmetric source");
        a[i * 2 * n + j] = exact(c[i * n + j]);
        a[i * 2 * n + n + j] = (i == j ? 1 : 0);
      }
    for (size_t k = 0; k < n; ++k) {
      mpq_class pivot = a[k * 2 * n + k];
      check(pivot > 0, "positive exact Schur pivot");
      determinant *= pivot;
      for (size_t j = 0; j < 2 * n; ++j)
        a[k * 2 * n + j] /= pivot;
      for (size_t i = 0; i < n; ++i)
        if (i != k) {
          mpq_class t = a[i * 2 * n + k];
          for (size_t j = 0; j < 2 * n; ++j)
            a[i * 2 * n + j] -= t * a[k * 2 * n + j];
        }
    }
    inverse.resize(n * n);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        inverse[i * n + j] = a[i * 2 * n + n + j];
    // Exact inverse identity is a second algebraic check independent of pivots.
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j) {
        mpq_class sum = 0;
        for (size_t k = 0; k < n; ++k)
          sum += exact(c[i * n + k]) * inverse[k * n + j];
        check(sum == (i == j ? 1 : 0), "exact C*inverse identity");
      }
  }
  mpq_class quadratic(const std::vector<double> &r) const {
    check(r.size() == n, "canonical RHS dimension");
    mpq_class q = 0;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        q += exact(r[i]) * inverse[i * n + j] * exact(r[j]);
    return q;
  }
  I quadratic(const std::vector<I> &r) const {
    check(r.size() == n, "interval RHS dimension");
    I q(0L);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        q = q + r[i] * I(inverse[i * n + j]) * r[j];
    return q;
  }
  I normalization() const {
    I pi;
    mpfr_const_pi(pi.lo.v, MPFR_RNDD);
    mpfr_const_pi(pi.hi.v, MPFR_RNDU);
    return I(static_cast<long>(n)) * log(I(2L) * pi);
  }
  I density(const I &q) const {
    return -(q + log(I(determinant)) + normalization()) / I(2L);
  }
};
struct Geometry {
  I E, radial;
};
Geometry geometry(double z, const std::array<double, 5> &q) {
  check(z >= 0 && z <= 2.5, "reference redshift");
  I E(1L), radial(0L);
  for (size_t k = 0; k < 5 && z > cosmology::piecewise_q_edges[k]; ++k) {
    const auto start = cosmology::piecewise_q_edges[k],
               end = std::min(z, cosmology::piecewise_q_edges[k + 1]);
    const I u = I(1L) + I(start), v = I(1L) + I(end), L = log(v / u),
            coefficient(q[k]);
    radial = radial +
             u / E * (q[k] == 0 ? L : -expm1(-coefficient * L) / coefficient);
    E = E * exp((I(1L) + coefficient) * L);
  }
  return {E, radial};
}
I predict(bao::Query query, const bao::PiecewiseModelPoint &p) {
  auto g = geometry(query.z, p.q);
  const I scale(mpq_class(299792458L) / 1000 / exact(p.ruler.h0_rd_km_s));
  auto dm = scale * g.radial, dh = scale / g.E;
  switch (query.observable) {
  case bao::Observable::transverse_over_ruler:
    return dm;
  case bao::Observable::hubble_over_ruler:
    return dh;
  case bao::Observable::volume_over_ruler:
    return cbrt(I(query.z) * dm * dm * dh);
  }
  throw std::runtime_error("unknown reference observable");
}
std::array<std::array<double, 5>, 8> grid{{{0, 0, 0, 0, 0},
                                           {-1, -1, -1, -1, -1},
                                           {.5, .5, .5, .5, .5},
                                           {-.4, -.4, -.2, .1, .3},
                                           {-1, 0, -1, 0, -1},
                                           {-3, 2, -3, 2, -3},
                                           {-3, -3, -3, -3, -3},
                                           {2, 2, 2, 2, 2}}};
void controls() {
  Matrix m({4, 1, 1, 9});
  check(m.determinant == 35 &&
            m.quadratic(std::vector<double>{2, -3}) == mpq_class(12, 5),
        "exact cofactor control");
  check(
      contains(m.quadratic(std::vector<I>{I(2L), I(-3L)}), I(mpq_class(12, 5))),
      "exact interval q control");
  for (double z : {0., .1, .3, .6, 1., 2.33, 2.5})
    for (double q : {-3., -1., 0., 2.}) {
      std::array<double, 5> a;
      a.fill(q);
      auto g = geometry(z, a);
      mpq_class u = 1 + exact(z);
      if (q == -3) {
        check(contains(g.E, I(1 / (u * u))),
              "allminus3 rational E containment");
        check(contains(g.radial, I((u * u * u - 1) / 3)),
              "allminus3 rational I containment");
      }
      if (q == 2) {
        check(contains(g.E, I(u * u * u)), "all2 rational E containment");
        check(contains(g.radial, I((1 - 1 / (u * u)) / 2)),
              "all2 rational I containment");
      }
      if (q == 0) {
        check(contains(g.radial, log(I(u))), "qzero removable log containment");
        check(contains(g.E, I(u)), "qzero rational expansion");
      }
      if (q == -1) {
        check(contains(g.E, I(1L)), "allminus1 rational expansion");
        check(contains(g.radial, I(exact(z))), "allminus1 rational radial");
      }
    }
  Matrix scalar({4});
  check(scalar.determinant == 4 &&
            scalar.quadratic(std::vector<double>{2}) == 1,
        "scalar variance/standardized quadratic control");
  I pi;
  mpfr_const_pi(pi.lo.v, MPFR_RNDD);
  mpfr_const_pi(pi.hi.v, MPFR_RNDU);
  auto normal = -I(mpq_class(1, 2)) - log(I(2L)) - log(I(2L) * pi) / I(2L);
  auto normalized = scalar.density(I(1L));
  check(mpfr_cmp(normal.lo.v, normalized.hi.v) <= 0 &&
            mpfr_cmp(normal.hi.v, normalized.lo.v) >= 0,
        "normalized scalar normal analytic identity");
  // Elementary sanity has independent rational series ancestry. Machin's
  // identity uses alternating atan tails; ln2 uses the positive atanh series.
  auto atan_bounds = [](long denominator) {
    mpq_class x(1, denominator), power = x, sum = 0;
    for (long k = 0; k < 48; ++k) {
      sum += (k % 2 ? -power : power) / (2 * k + 1);
      power *= x * x;
    }
    return std::pair<mpq_class, mpq_class>{sum, sum + power / 97};
  };
  auto a = atan_bounds(5), b = atan_bounds(239);
  mpq_class pi_lower = 16 * a.first - 4 * b.second,
            pi_upper = 16 * a.second - 4 * b.first;
  check(mpfr_cmp_q(pi.lo.v, pi_lower.get_mpq_t()) >= 0 &&
            mpfr_cmp_q(pi.hi.v, pi_upper.get_mpq_t()) <= 0,
        "pi enclosed by rational Machin series");
  mpq_class third(1, 3), power = third, ln2_lower = 0;
  for (long k = 0; k < 32; ++k) {
    ln2_lower += 2 * power / (2 * k + 1);
    power *= third * third;
  }
  mpq_class ln2_upper = ln2_lower + 2 * power / 65 / (1 - third * third);
  auto ln2 = log(I(2L));
  check(mpfr_cmp_q(ln2.lo.v, ln2_lower.get_mpq_t()) >= 0 &&
            mpfr_cmp_q(ln2.hi.v, ln2_upper.get_mpq_t()) <= 0,
        "log2 enclosed by rational positive series");
  const I neg(-3L), pos(2L);
  check(contains(neg * pos, I(-6L)), "negative product sign");
  check(contains(pos / neg, I(mpq_class(-2, 3))), "negative divisor sign");
  bool rejected = false;
  try {
    I crosses;
    mpfr_set_si(crosses.lo.v, -1, MPFR_RNDN);
    mpfr_set_si(crosses.hi.v, 1, MPFR_RNDN);
    (void)(pos / crosses);
  } catch (...) {
    rejected = true;
  }
  check(rejected, "zero-crossing denominator rejected");
}
bao::DensityInput read(const char *mean, const char *cov) {
  bao::DensityInput in;
  std::ifstream a(mean), b(cov);
  check(bool(a) && bool(b), "assets available");
  std::string line;
  while (std::getline(a, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream s(line);
    double z, value;
    std::string kind, extra;
    check(bool(s >> z >> value >> kind) && !(s >> extra),
          "three source fields");
    bao::Observable o;
    if (kind == "DM_over_rs")
      o = bao::Observable::transverse_over_ruler;
    else if (kind == "DH_over_rs")
      o = bao::Observable::hubble_over_ruler;
    else if (kind == "DV_over_rs")
      o = bao::Observable::volume_over_ruler;
    else
      throw std::runtime_error("source tag");
    in.queries.push_back({z, o});
    in.observed.push_back(value);
    in.ordered_ids.push_back("original13-row:" +
                             std::to_string(in.queries.size() - 1));
  }
  double x;
  while (b >> x)
    in.covariance.push_back(x);
  check(b.eof() && in.queries.size() == 13 && in.covariance.size() == 169,
        "source13 shape");
  check(in.queries[11].observable == bao::Observable::hubble_over_ruler &&
            in.queries[12].observable == bao::Observable::transverse_over_ruler,
        "last DH then DM order");
  in.role = bao::RowRole::released_fitted_distance_summary;
  in.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  in.table_identity = std::string(test_reference::bao_mean_sha256);
  in.covariance_identity = std::string(test_reference::bao_cov_sha256);
  in.ordering_provenance =
      "declared source order; external SHA guard assertion";
  in.calibration_provenance = "free ruler unknown early calibration";
  in.dependence_provenance = "full released covariance; other probes unknown";
  return in;
}
struct Saved {
  I canonical, ideal, background;
  std::vector<I> predictions;
};
std::vector<Saved>
references(const bao::DensityInput &in, const Matrix &m,
           const std::vector<bao::PiecewiseModelPoint> &models,
           const bao::PiecewiseDensityBatch &batch) {
  std::vector<Saved> out;
  for (size_t k = 0; k < models.size(); ++k) {
    std::vector<I> p, r;
    for (size_t j = 0; j < in.queries.size(); ++j) {
      p.push_back(predict(in.queries[j], models[k]));
      r.push_back(I(in.observed[j]) - p.back());
    }
    const I canonicalq(m.quadratic(batch.slots[k].residuals)),
        idealq = m.quadratic(r);
    out.push_back({m.density(canonicalq), m.density(idealq),
                   (canonicalq - idealq) / I(2L), p});
  }
  return out;
}
void original(const char *mean, const char *cov) {
  auto in = read(mean, cov);
  bao::DensityPolicy prep;
  prep.maximum_models = 24;
  prep.maximum_matrix_elements = 169;
  prep.maximum_string_bytes = 10000;
  prep.maximum_native_bytes = 10000000;
  prep.maximum_forward_sensitivity = 1e-10;
  prep.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  auto owner = bao::prepare_density(in, prep);
  check(owner.status() == statistics::DensityStatus::finite,
        "source factor qualified guard");
  std::vector<bao::PiecewiseModelPoint> models;
  for (auto q : grid)
    for (double r : {5000., 10000., 15000.})
      models.emplace_back(q, bao::Ruler(r));
  bao::PiecewiseDensityPolicy policy;
  policy.maximum_models = 24;
  policy.maximum_native_bytes = 10000000;
  policy.maximum_forward_sensitivity = 1e-10;
  policy.arithmetic = prep.arithmetic;
  policy.observables.background.maximum_segment_visits = 2000;
  auto batch = owner.evaluate_piecewise(models, policy);
  check(batch.slots.size() == 24, "all24 returned");
  for (const auto &s : batch.slots)
    check(s.result.density.status == statistics::DensityStatus::finite,
          "finite candidate");
  Matrix m(in.covariance);
  precision = 256;
  auto coarse = references(in, m, models, batch);
  precision = 384;
  auto fine = references(in, m, models, batch);
  unsigned failures = 0;
  for (size_t i = 0; i < 24; ++i) {
    auto &f = fine[i];
    check(contains(coarse[i].canonical, f.canonical) &&
              contains(coarse[i].ideal, f.ideal) &&
              contains(coarse[i].background, f.background),
          "precision nesting density");
    for (size_t j = 0; j < 13; ++j)
      check(contains(coarse[i].predictions[j], f.predictions[j]),
            "precision nesting prediction");
    auto factor =
             absmax(I(batch.slots[i].result.density.log_value) - f.canonical),
         bg = absmax(f.background),
         full = absmax(I(batch.slots[i].result.density.log_value) - f.ideal);
    bool ref = within(width(f.canonical), "1/500000000") &&
               within(width(f.ideal), "1/500000000") &&
               within(width(f.background), "1/500000000");
    F observable_fraction;
    bool observable_pass = true;
    for (size_t j = 0; j < 13; ++j) {
      F lower_abs, upper_abs;
      mpfr_abs(lower_abs.v, f.predictions[j].lo.v, MPFR_RNDD);
      mpfr_abs(upper_abs.v, f.predictions[j].hi.v, MPFR_RNDD);
      if (mpfr_cmp(upper_abs.v, lower_abs.v) < 0)
        lower_abs = upper_abs;
      if (mpfr_sgn(f.predictions[j].lo.v) <= 0 &&
          mpfr_sgn(f.predictions[j].hi.v) >= 0)
        mpfr_set_zero(lower_abs.v, 0);
      I minimum;
      minimum.lo = lower_abs;
      minimum.hi = lower_abs;
      auto allowed = I(mpq_class(1, 500000000000L)) +
                     I(mpq_class(1, 5000000000L)) * minimum;
      auto error = absmax(I(batch.slots[i].predictions[j]) - f.predictions[j]);
      if (mpfr_cmp(error.v, allowed.lo.v) > 0)
        observable_pass = false;
      F fraction;
      mpfr_div(fraction.v, error.v, allowed.lo.v, MPFR_RNDU);
      if (mpfr_cmp(fraction.v, observable_fraction.v) > 0)
        observable_fraction = fraction;
    }
    bool fp = within(factor, "3/1000000000"), bp = within(bg, "1/200000000"),
         ap = within(full, "1/100000000");
    if (!(ref && fp && bp && ap && observable_pass))
      ++failures;
    std::printf(
        "{\"point\":%zu,\"q\":[%.17g,%.17g,%.17g,%.17g,%.17g],\"h0rd\":%.17g,"
        "\"factor_upper\":\"%s\",\"background_upper\":\"%s\",\"full_upper\":\"%"
        "s\",\"reference_width\":\"%s\",\"reference_pass\":%s,\"factor_pass\":%"
        "s,\"background_pass\":%s,\"full_pass\":%s,\"observable_pass\":%s,"
        "\"observable_fraction\":\"%s\"}\n",
        i, models[i].q[0], models[i].q[1], models[i].q[2], models[i].q[3],
        models[i].q[4], models[i].ruler.h0_rd_km_s, text(factor).c_str(),
        text(bg).c_str(), text(full).c_str(), text(width(f.ideal)).c_str(),
        ref ? "true" : "false", fp ? "true" : "false", bp ? "true" : "false",
        ap ? "true" : "false", observable_pass ? "true" : "false",
        text(observable_fraction).c_str());
  }
  std::printf("{\"suite\":\"piecewise_BAO_directed_original24\",\"failed\":%u,"
              "\"segments\":%zu,\"gmp\":\"%s\",\"mpfr\":\"%s\"}\n",
              failures, batch.segment_visits, gmp_version, mpfr_get_version());
  check(failures == 0, "allocated directed gates");
}
} // namespace
int main(int argc, char **argv) {
  try {
    controls();
    if (argc == 4 && std::string(argv[1]) == "--verified-original-assets")
      original(argv[2], argv[3]);
    else
      check(argc == 1, "optional mode needs external SHA guard");
    std::printf("directed interval native PASS %u checks\n", checks);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "directed interval native: %s\n", e.what());
    return 1;
  }
}
