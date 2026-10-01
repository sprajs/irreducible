// Independently constructed synthetic PLZ ladder, no released assets/copied
// code. KKT [C X;X^T 0] elimination avoids production normal/Cholesky equation
// route. Frozen allocations: coefficients/predictions2e-11 abs+rel,H0
// 2e-10abs+rel,q1e-9. Compression retains qmin and conditional covariance; no
// prior normalization.
#include "irred/calibration_ladder.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <string_view>
#include <utility>
using namespace irred;
using namespace irred::calibration;
using W = long double;
namespace {
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void near(W a, W b, W abs = 2e-11L, W rel = 2e-11L) {
  need(std::abs(a - b) <= abs + rel * std::abs(b), "peer numerical allocation");
}
std::vector<W> solve(std::vector<W> a, std::vector<W> b) {
  size_t n = b.size();
  for (size_t k = 0; k < n; ++k) {
    size_t p = k;
    for (size_t i = k + 1; i < n; ++i)
      if (std::abs(a[i * n + k]) > std::abs(a[p * n + k]))
        p = i;
    need(a[p * n + k] != 0, "KKT nonsingular");
    for (size_t j = k; j < n; ++j)
      std::swap(a[p * n + j], a[k * n + j]);
    std::swap(b[p], b[k]);
    W v = a[k * n + k];
    for (size_t j = k; j < n; ++j)
      a[k * n + j] /= v;
    b[k] /= v;
    for (size_t i = 0; i < n; ++i)
      if (i != k) {
        W v2 = a[i * n + k];
        for (size_t j = k; j < n; ++j)
          a[i * n + j] -= v2 * a[k * n + j];
        b[i] -= v2 * b[k];
      }
  }
  return b;
}
struct Oracle {
  std::vector<W> beta;
  W q, var_eta;
};
Oracle oracle(const std::vector<double> &c, const std::vector<double> &x,
              const std::vector<double> &y, size_t p, bool variance = true) {
  size_t n = y.size(), d = n + p;
  std::vector<W> a(d * d), rhs(d);
  for (size_t i = 0; i < n; ++i) {
    rhs[i] = y[i];
    for (size_t j = 0; j < n; ++j)
      a[i * d + j] = c[i * n + j];
    for (size_t j = 0; j < p; ++j)
      a[i * d + n + j] = a[(n + j) * d + i] = x[i * p + j];
  }
  auto result = solve(a, rhs);
  Oracle out;
  out.beta.assign(result.begin() + n, result.end());
  out.q = 0;
  for (size_t i = 0; i < n; ++i) {
    W residual = y[i];
    for (size_t j = 0; j < p; ++j)
      residual -= x[i * p + j] * out.beta[j];
    out.q += residual * result[i];
  }
  out.var_eta = 0;
  if (variance) {
    std::fill(rhs.begin(), rhs.end(), 0);
    rhs[n + 7] = 1;
    out.var_eta = -solve(a, rhs)[n + 7];
  }
  return out;
}
struct Fixture {
  Model m;
  std::vector<double> x, c, y;
  std::vector<std::string> ids;
};
Fixture fixture() {
  Fixture f;
  f.m.ordered_host_ids = {"a", "b", "host"};
  f.m.h_reference_km_s_Mpc = 70;
  f.m.metallicity_reference_dex = 0;
  f.m.magnitude_convention =
      "synthetic common photometric magnitude convention";
  f.m.metallicity_coordinate_identity = "synthetic log metal abundance dex";
  f.m.distance_shape_identity = "synthetic supplied HF reference moduli";
  f.m.calibration_identity = "one observed zero point";
  f.m.dependence_identity = "synthetic common covariance";
  f.m.conditional_covariance_identity =
      "conditional on delta excludes marginalized zero point";
  auto add = [&](Row r, std::vector<double> v, double offset) {
    r.row_id = "r" + std::to_string(f.m.rows.size());
    f.ids.push_back(r.row_id);
    f.m.rows.push_back(r);
    f.x.insert(f.x.end(), v.begin(), v.end());
    const double truth[]{29, 30, 31, -5, -3, -.2, -19, .15, .03};
    W y = offset;
    for (size_t j = 0; j < 9; ++j)
      y += W(v[j]) * truth[j];
    f.y.push_back(double(y) + .001 * (int(f.y.size() % 5) - 2));
  };
  for (size_t h = 0; h < 3; ++h) {
    std::vector<double> v(9);
    v[h] = 1;
    Row r;
    r.kind = RowKind::anchor_modulus;
    r.host_id = f.m.ordered_host_ids[h];
    add(r, v, 0);
  }
  for (size_t h = 0; h < 3; ++h)
    for (int k = 0; k < 3; ++k) {
      std::vector<double> v(9);
      v[h] = 1;
      v[3] = 1;
      v[4] = k - .8;
      v[5] = (k == 1 ? .4 : -.3) + .1 * h;
      v[8] = .5 + .2 * h;
      Row r;
      r.kind = RowKind::cepheid;
      r.event_id = "cep" + std::to_string(h) + "-" + std::to_string(k);
      r.host_id = f.m.ordered_host_ids[h];
      r.log10_period_days = 1 + v[4];
      r.metallicity_dex = v[5];
      r.calibration_response = v[8];
      add(r, v, 0);
    }
  for (size_t h = 0; h < 3; ++h) {
    std::vector<double> v(9);
    v[h] = 1;
    v[6] = 1;
    v[8] = .3;
    Row r;
    r.kind = RowKind::calibrator_supernova;
    r.host_id = f.m.ordered_host_ids[h];
    r.event_id = "cal" + std::to_string(h);
    r.calibration_response = .3;
    add(r, v, 0);
  }
  for (int k = 0; k < 3; ++k) {
    std::vector<double> v(9);
    v[6] = 1;
    v[7] = -1;
    v[8] = -.2;
    Row r;
    r.kind = RowKind::hubble_supernova;
    r.event_id = "hf" + std::to_string(k);
    r.reference_modulus_mag = 35 + k;
    r.calibration_response = -.2;
    add(r, v, r.reference_modulus_mag);
  }
  std::vector<double> v(9);
  v[8] = 1;
  Row r;
  r.kind = RowKind::calibration_measurement;
  r.calibration_response = 1;
  add(r, v, 0);
  size_t n = f.y.size();
  f.c.resize(n * n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      f.c[i * n + j] = (i == j ? .04 : 0) + .002;
  return f;
}
statistics::Metadata meta(const Fixture &f) {
  statistics::Metadata m;
  m.ordered_ids = f.ids;
  m.measure = "product d(mag)";
  m.source_semantics = "synthetic controls";
  m.calibration_provenance = f.m.calibration_identity;
  m.dependence_provenance = f.m.dependence_identity;
  m.uncertainty_identity = f.m.conditional_covariance_identity;
  m.table_identity = "independent synthetic peer rows";
  m.ordering_provenance = "explicit synthetic rows";
  return m;
}
Oracle compare(Fixture f) {
  auto residual = f.y;
  for (size_t i = 0; i < residual.size(); ++i)
    if (f.m.rows[i].kind == RowKind::hubble_supernova)
      residual[i] -= f.m.rows[i].reference_modulus_mag;
  auto ref = oracle(f.c, f.x, residual, 9);
  need(ref.q > 1e-7L, "nonzero compression Schur constant required");
  auto g = statistics::prepare_gaussian(f.c, statistics::MatrixKind::covariance,
                                        meta(f), f.c.size(), 1e-10);
  auto ladder = Ladder::prepare(std::move(g), f.m);
  need(ladder.status() == statistics::DensityStatus::finite,
       "typed ladder prepared");
  auto got = ladder.fit(f.y, f.ids);
  need(got.relative_fit.status == statistics::DensityStatus::finite,
       "typed ladder fit");
  for (size_t j = 0; j < 9; ++j)
    near(got.relative_fit.coefficients[j], ref.beta[j]);
  near(got.relative_fit.quadratic, ref.q, 1e-9L, 0);
  near(got.relative_fit.relative_log_score, -ref.q / 2, 1e-9L, 0);
  near(got.h0_km_s_Mpc, 70 * std::pow(10.L, ref.beta[7] / 5), 2e-10L, 2e-10L);
  auto predicted =
      ladder.predict(got.relative_fit.coefficients, parameter_ids(f.m));
  need(predicted.status == statistics::DensityStatus::finite,
       "typed prediction");
  for (size_t i = 0; i < f.y.size(); ++i) {
    W expected = f.m.rows[i].kind == RowKind::hubble_supernova
                     ? f.m.rows[i].reference_modulus_mag
                     : 0;
    for (size_t j = 0; j < 9; ++j)
      expected += f.x[i * 9 + j] * ref.beta[j];
    near(predicted.values_mag[i], expected);
  }
  // Profile eta at displaced points, compare direct conditional KKT with
  // compressed quadratic qmin+(eta-etahat)^2/Var(eta), retaining qmin
  // explicitly.
  for (W delta : {-.03L, .02L}) {
    W eta = ref.beta[7] + delta;
    std::vector<double> xx, res = residual;
    for (size_t i = 0; i < res.size(); ++i) {
      res[i] -= double(f.x[i * 9 + 7] * eta);
      for (size_t j = 0; j < 9; ++j)
        if (j != 7)
          xx.push_back(f.x[i * 9 + j]);
    }
    auto fixed = oracle(f.c, xx, res, 8, false);
    near(fixed.q, ref.q + delta * delta / ref.var_eta, 1e-9L, 0);
  }
  return ref;
}
void estimator_law_controls(const Fixture &f){
 // Original 60/90-digit Decimal augmented KKT and exact linear-noise
 // contrast references were frozen before native evaluation. The supplied
 // fixed C=.04I+.00211^T is conditional on the one explicit delta column;
 // it is observation noise, not a second prior on that calibration.
 constexpr W variance=.05337175729789333492058992206202025997L;
 constexpr W nominal=75.00635136663244902683619522L;
 constexpr W expectation=75.43204820553021741161294790L;
 const std::array<double,7> probabilities{.5,.15865525393145705,.8413447460685429,
                                       1e-12,1-1e-12,.975,.025};
 // Original mpmath inverse-erf references at60/90digits, using exact binary64
 // input probabilities, not ideal complementary decimal endpoints. The two
 // rounded endpoint probabilities have distinct |z| and quantiles.
 constexpr std::array<W,7> quantiles{75.00635136663244902683620L,
  67.43624555932101417973107L,83.42624502109619850199308L,
  35.48723238641414936590057L,158.5346112775924196628594L,
  92.39692178296864305683107L,60.88896293049176993463111L};
 auto g=statistics::prepare_gaussian(f.c,statistics::MatrixKind::covariance,meta(f),f.c.size(),1e-10);
 auto ladder=Ladder::prepare(std::move(g),f.m);need(ladder.status()==statistics::DensityStatus::finite,"sampling-law ladder admitted");
 auto law=ladder.h0_estimator_law(.15,"explicit synthetic generating eta mean",probabilities);
 need(law.status==statistics::DensityStatus::finite&&law.eta_estimator.status==statistics::DensityStatus::finite,"conditional estimator sampling law admitted");
 near(law.eta_estimator.variance,variance,2e-12L,2e-12L);
 near(law.nominal_h0_km_s_Mpc,nominal,2e-10L,2e-10L);
 near(law.sampling_expectation_h0_km_s_Mpc,expectation,2e-10L,2e-10L);
 need(law.assumed_eta_mean==.15&&law.eta_mean_identity=="explicit synthetic generating eta mean","mean is a supplied identified assumption");
 need(law.h0_quantiles_km_s_Mpc.size()==probabilities.size()&&law.probabilities==std::vector<double>(probabilities.begin(),probabilities.end()),"quantile probability axis preserved");
 for(size_t i=0;i<probabilities.size();++i)near(law.h0_quantiles_km_s_Mpc[i],quantiles[i],2e-10L,2e-10L);
 need(std::string_view(law.method_id)=="retained-qr-eta-lognormal-quantiles/v1"&&std::string_view(law.sampling_law_id)=="fixed-design-gaussian-observation-noise/lognormal-H0/v1","H0 sampling law identity separate from posterior");
 need(law.sampling_expectation_h0_km_s_Mpc>law.nominal_h0_km_s_Mpc,"lognormal mean differs from median");
 need(law.maximum_relative_tail_probability_error<=1e-12&&law.maximum_quantile_log_error_estimate<=1e-10,"tail and log-quantile diagnostics satisfy declared gates");
 need(law.normal_cdf_evaluations<=98*probabilities.size(),"normal inversion bounded work");
 // Exact fixed-design noise contrasts, no random sampling: independent
 // diagonal and common Gaussian noise columns span supplied covariance C.
 const auto base=ladder.fit(f.y,f.ids);need(base.relative_fit.status==statistics::DensityStatus::finite,"linear response baseline");
 W summed=0;const size_t n=f.y.size();
 for(size_t k=0;k<=n;++k){auto y=f.y;
  if(k<n)y[k]+=std::sqrt(f.c[0]-f.c[1]);else for(auto &v:y)v+=std::sqrt(f.c[1]);
  const auto perturbed=ladder.fit(y,f.ids);need(perturbed.relative_fit.status==statistics::DensityStatus::finite,"deterministic noise response admitted");
  W delta=W(perturbed.relative_fit.coefficients[7])-base.relative_fit.coefficients[7];summed+=delta*delta;}
 near(summed,variance,2e-12L,2e-12L);
 for(const std::array<double,1> p:std::array<std::array<double,1>,3>{{{0},{1},{std::numeric_limits<double>::quiet_NaN()}}}){
  auto bad=ladder.h0_estimator_law(.15,"supplied mean",p);need(bad.status!=statistics::DensityStatus::finite&&bad.h0_quantiles_km_s_Mpc.empty(),"invalid probabilities withhold law payload");}
 need(ladder.h0_estimator_law(.15,"",probabilities).status!=statistics::DensityStatus::finite,"unidentified eta mean refused");
 need(ladder.h0_estimator_law(std::numeric_limits<double>::infinity(),"supplied mean",probabilities).status!=statistics::DensityStatus::finite,"nonfinite eta mean refused");
 statistics::DesignPolicy quota;quota.maximum_payload_bytes=0;
 need(ladder.h0_estimator_law(.15,"supplied mean",probabilities,quota).h0_quantiles_km_s_Mpc.empty(),"law byte quota withholds payload");
 // Preserve the previously observed moved-Ladder status bug: the repaired
 // consumer must transfer the model/profile and invalidate prep fallbacks.
 auto moved=std::move(ladder);need(ladder.status()!=statistics::DensityStatus::finite,"moved Ladder status invalid");
 need(ladder.fit(f.y,f.ids).relative_fit.status!=statistics::DensityStatus::finite,"moved Ladder fit refused");
 need(ladder.h0_estimator_law(.15,"supplied mean",probabilities).status!=statistics::DensityStatus::finite,"moved Ladder law refused");
 auto stable=moved.h0_estimator_law(.15,"supplied mean",probabilities);near(stable.eta_estimator.variance,variance,2e-12L,2e-12L);
 auto* self=&moved;moved=std::move(*self);stable=moved.h0_estimator_law(.15,"supplied mean",probabilities);
 need(stable.status==statistics::DensityStatus::finite,"Ladder self move preserves owner");near(stable.nominal_h0_km_s_Mpc,nominal,2e-10L,2e-10L);
 Ladder assigned;assigned=std::move(moved);need(moved.status()!=statistics::DensityStatus::finite,"move assignment invalidates Ladder source");
 auto retained=assigned.h0_estimator_law(.15,"supplied mean",probabilities);need(retained.status==statistics::DensityStatus::finite,"assigned Ladder retains law");near(retained.eta_estimator.variance,variance,2e-12L,2e-12L);
}
} // namespace
int main() {
  try {
    auto full = fixture();
    compare(full);
    estimator_law_controls(full);
    const size_t n = full.y.size();
    // Remove/refit a third anchor, retaining two distinct anchors; then a
    // separate host Cepheid row. Each prediction uses fresh held-out rows.
    for (size_t held : {size_t(2), size_t(11)}) {
      auto training = full;
      training.m.rows.erase(training.m.rows.begin() + held);
      training.ids.erase(training.ids.begin() + held);
      training.y.erase(training.y.begin() + held);
      training.x.erase(training.x.begin() + held * 9,
                       training.x.begin() + (held + 1) * 9);
      training.c.clear();
      for (size_t i = 0; i < n; ++i)
        if (i != held)
          for (size_t j = 0; j < n; ++j)
            if (j != held)
              training.c.push_back(full.c[i * n + j]);
      auto ref = compare(training);
      Model fresh = full.m;
      fresh.rows = {full.m.rows[held]};
      std::vector<double> beta(ref.beta.begin(), ref.beta.end());
      auto prediction = predict(fresh, beta, parameter_ids(fresh));
      need(prediction.status == statistics::DensityStatus::finite,
           "held out fresh row prediction");
      W expected = 0;
      for (size_t j = 0; j < 9; ++j)
        expected += full.x[held * 9 + j] * ref.beta[j];
      near(prediction.values_mag[0], expected);
      need(std::abs(expected - full.y[held]) > 1e-8,
           "held out is not forced training residual");
    }
    std::cout << "PASS " << checks
              << " calibration ladder independent peer controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
