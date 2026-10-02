#include "irred/calibration_predictive.hpp"
#include <array>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <numbers>
#include <stdexcept>
namespace {
std::atomic<bool> armed = false;
std::atomic<size_t> calls = 0, fail_on = 0;
std::atomic<long> live = 0;
} // namespace
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(size_t n) {
  if (armed && ++calls == fail_on)
    throw std::bad_alloc();
  if (void *p = std::malloc(n ? n : 1)) {
    ++live;
    return p;
  }
  throw std::bad_alloc();
}
NOINLINE void *operator new[](size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
NOINLINE void operator delete(void *p, size_t) noexcept {
  ::operator delete(p);
}
NOINLINE void operator delete[](void *p, size_t) noexcept {
  ::operator delete(p);
}
namespace {
using namespace irred::statistics;
using namespace irred::calibration;
using N = irred::numerics::Status;
unsigned checks = 0;
void need(bool v, const char *w) {
  ++checks;
  if (!v)
    throw std::runtime_error(w);
}
void near(long double a, long double b) {
  need(std::abs(a - b) <= 2e-12L * (1 + std::abs(b)),
       "independent scalar conjugate reference");
}
Model model(bool future = false) {
  Model m;
  m.ordered_host_ids = {"A", "B"};
  m.magnitude_convention = "common synthetic mag";
  m.metallicity_coordinate_identity = "synthetic dex";
  m.distance_shape_identity = "fixed supplied shape";
  m.calibration_identity = "single shared delta";
  m.dependence_identity = "noise independent of beta";
  m.conditional_covariance_identity = "explicit conditional noise";
  if (!future)
    m.rows = {{RowKind::anchor_modulus, "train-anchor", "A", "", 0, 0, 0, 0}};
  else
    m.rows = {
        {RowKind::cepheid, "future-cep", "A", "cep-new", 1, 1, 0, 0},
        {RowKind::hubble_supernova, "future-hf", "", "sn-new", 1, 0, 0, 40}};
  return m;
}
Gaussian noise(const Model &m) {
  Metadata md;
  for (const auto &r : m.rows)
    md.ordered_ids.push_back(r.row_id);
  md.measure = "product d(mag)";
  md.source_semantics = "synthetic controls";
  md.table_identity = "original owner control";
  md.ordering_provenance = "declared exact rows";
  md.calibration_provenance = m.calibration_identity;
  md.dependence_provenance = m.dependence_identity;
  md.uncertainty_identity = m.conditional_covariance_identity;
  std::vector<double> c = m.rows.size() == 1
                              ? std::vector<double>{1}
                              : std::vector<double>{1, .25, .25, 1};
  return prepare_gaussian(c, MatrixKind::covariance, md, 1000, 1e-10,
                          irred::numerics::Arithmetic::longdouble_cpu_v1);
}
ParameterPrior prior(const Model &m) {
  auto l = linearize(m);
  ParameterPrior p;
  p.ordered_parameter_ids = l.metadata.ordered_parameter_ids;
  p.parameter_units = l.metadata.parameter_units;
  p.shared_nuisance_ids = l.metadata.shared_nuisance_ids;
  p.residual_unit = "mag";
  p.design_identity = l.metadata.design_identity;
  p.dependence_identity = l.metadata.dependence_identity;
  p.prior_identity = "independent scalar normal controls";
  p.parameter_measure = "product of declared parameter coordinates";
  p.noise_independence_declared = true;
  p.mean = {31, 32, -5, -3, .25, -19, .5, .125};
  p.covariance.assign(64, 0);
  for (size_t i = 0; i < 8; ++i)
    p.covariance[i * 8 + i] = 2;
  return p;
}
PredictiveMetadata metadata(const LadderPosterior &p, const Model &f) {
  auto l = linearize(f);
  PredictiveMetadata md;
  md.ordered_parameter_ids = l.metadata.ordered_parameter_ids;
  md.parameter_units = l.metadata.parameter_units;
  for (const auto &u : md.parameter_units)
    md.response_units.push_back("mag/" + u);
  md.training_event_ids = p.linearization().predictive_event_ids;
  md.future_event_ids = l.predictive_event_ids;
  md.future_unit = "mag";
  md.future_covariance_unit = "mag^2";
  md.future_measure = "product d(mag)";
  md.response_identity = l.metadata.design_identity;
  md.dependence_identity = l.metadata.dependence_identity;
  md.conditioning_identity = "explicit original training values";
  md.conditional_noise_identity = f.conditional_covariance_identity;
  md.future_noise_independence_declared = true;
  md.noise_conditional_on_parameters_declared = true;
  return md;
}
void run() {
  const auto m = model(), f = model(true);
  auto g = noise(m);
  auto pr = prior(m);
  const auto lb = linearization_preparation_payload_bound(m);
  need(lb.has_value(), "standalone linearization preparation bound");
  DesignPolicy lp;
  lp.maximum_payload_bytes = *lb - 1;
  calls = 0;
  fail_on = SIZE_MAX;
  armed = true;
  auto linear_refused = linearize(m, lp);
  armed = false;
  need(linear_refused.status != DensityStatus::finite && calls == 0,
       "linearization quota precedes candidate allocation");
  lp.maximum_payload_bytes = *lb;
  need(linearize(m, lp).status == DensityStatus::finite,
       "linearization exact preparation boundary");
  auto input_m = m;
  auto input_prior = pr;
  const auto bound =
      LadderPosterior::preparation_payload_bound(g, input_m, input_prior);
  need(bound.has_value(), "posterior bound");
  PosteriorPolicy pp;
  pp.maximum_payload_bytes = *bound - 1;
  auto refused = LadderPosterior::prepare(std::move(g), std::move(input_m),
                                          std::move(input_prior), pp);
  need(refused.status() != DensityStatus::finite &&
           g.status() == DensityStatus::finite,
       "source preserved below bound");
  input_m = m;
  input_prior = pr;
  pp.maximum_payload_bytes =
      *LadderPosterior::preparation_payload_bound(g, input_m, input_prior);
  auto p = LadderPosterior::prepare(std::move(g), std::move(input_m),
                                    std::move(input_prior), pp);
  need(p.status() == DensityStatus::finite &&
           g.status() != DensityStatus::finite,
       "proper rank deficient training admitted");
  auto y = std::array<double, 1>{32};
  auto rows = p.linearization().ordered_row_ids;
  auto mu = p.condition(y, rows);
  need(mu.status == DensityStatus::finite, "condition");
  near(mu.value[0], 95.L / 3);
  near(p.posterior().covariance()[0], 2.L / 3);
  near(p.posterior().covariance()[9], 2);
  auto h = p.h0_projection(y, rows);
  need(h.status == DensityStatus::finite, "H0 projection");
  near(h.median_h0_km_s_Mpc, 70 * std::exp(std::log(10.L) / 10));
  near(h.expectation_h0_km_s_Mpc,
       70 * std::exp(std::log(10.L) / 10 + std::pow(std::log(10.L) / 5, 2)));
  auto ng = noise(f);
  auto md = metadata(p, f);
  auto input_f = f;
  auto input_md = md;
  auto pb =
      LadderPredictive::preparation_payload_bound(p, ng, input_f, input_md);
  need(pb.has_value(), "predictive bound");
  PredictivePolicy qp;
  qp.maximum_payload_bytes = *pb - 1;
  auto qbad = LadderPredictive::prepare(p, y, rows, ng, std::move(input_f),
                                        std::move(input_md), qp);
  need(qbad.status() != DensityStatus::finite, "predictive below bound");
  input_f = f;
  input_md = md;
  qp.maximum_payload_bytes =
      *LadderPredictive::preparation_payload_bound(p, ng, input_f, input_md);
  auto q = LadderPredictive::prepare(p, y, rows, ng, std::move(input_f),
                                     std::move(input_md), qp);
  need(q.status() == DensityStatus::finite, "prepare joint future");
  near(q.mean()[0], 643.L / 24);
  near(q.mean()[1], 165.L / 8);
  near(q.covariance()[0], 17.L / 3);
  near(q.covariance()[1], 9.L / 4);
  near(q.covariance()[3], 7);
  const std::array<double, 2> fy{27, 21};
  auto fr = q.linearization().ordered_row_ids;
  need(q.log_density(fy, fr).density.status == DensityStatus::finite,
       "whole joint density");
  PredictivePolicy ep;
  auto rb = q.retained_payload_bound(), eb = q.evaluation_payload_bound();
  need(rb && eb, "evaluation bounds");
  ep.maximum_payload_bytes = *rb + *eb - 1;
  need(q.log_density(fy, fr, ep).density.status != DensityStatus::finite,
       "evaluation byte boundary");
  ep.maximum_payload_bytes = *rb + *eb;
  need(q.log_density(fy, fr, ep).density.status == DensityStatus::finite,
       "evaluation exact byte boundary");
  auto rounded = fy;
  rounded[1] =
      std::nextafter(40., 0.); // subtraction by 40 remains exact: Sterbenz.
  need(q.log_density(rounded, fr).density.status == DensityStatus::finite,
       "exact near-offset subtraction");
  rounded[1] = .1;
  need(q.log_density(rounded, fr).density.status != DensityStatus::finite,
       "rounded subtraction refused");
  auto badmd = md;
  badmd.future_event_ids[0] = md.training_event_ids[0];
  need(LadderPredictive::prepare(p, y, rows, ng, f, badmd).status() !=
           DensityStatus::finite,
       "lineage mismatch");
  auto overlap = f;
  overlap.rows[0].row_id = m.rows[0].row_id;
  need(LadderPredictive::prepare(p, y, rows, ng, overlap, metadata(p, overlap))
               .status() != DensityStatus::finite,
       "overlap refused");
  auto qmove = std::move(q);
  need(q.status() != DensityStatus::finite && q.mean().empty(),
       "moved future invalid");
  auto *qself = &qmove;
  qmove = std::move(*qself);
  need(qmove.log_density(fy, fr).density.status == DensityStatus::finite,
       "future self move");
  auto pmove = std::move(p);
  need(p.status() != DensityStatus::finite, "moved posterior invalid");
  auto *pself = &pmove;
  pmove = std::move(*pself);
  need(pmove.condition(y, rows).status == DensityStatus::finite,
       "posterior self move");
  auto anisotropic = pr;
  anisotropic.covariance[6 * 8 + 6] = .0001;
  auto decimal_source = noise(m);
  auto decimal_refused =
      LadderPosterior::prepare(std::move(decimal_source), m, anisotropic);
  need(decimal_refused.status() != DensityStatus::finite &&
           decimal_source.status() == DensityStatus::finite,
       "preserved anisotropic decimal prerequisite refusal");
  anisotropic.covariance[6 * 8 + 6] = 1. / 4096;
  auto ag = noise(m);
  auto ap = LadderPosterior::prepare(std::move(ag), m, anisotropic);
  need(ap.status() != DensityStatus::finite &&
           ag.status() == DensityStatus::finite,
       "preserved small dyadic prerequisite refusal");
  anisotropic.covariance[6 * 8 + 6] = 1. / 16;
  auto admitted_source = noise(m);
  ap = LadderPosterior::prepare(std::move(admitted_source), m, anisotropic);
  need(ap.status() == DensityStatus::finite,
       "resolved anisotropic proper prior admitted");
  auto ah = ap.h0_projection(y, rows);
  if (ah.status == DensityStatus::finite) {
    long double norm = 0;
    auto cv = ap.posterior().covariance();
    for (size_t i = 0; i < 8; ++i) {
      long double row = 0;
      for (size_t j = 0; j < 8; ++j)
        row += std::abs(cv[i * 8 + j]);
      norm = std::max(norm, row);
    }
    need(ah.eta_variance_absolute_error_estimate >=
             norm * ap.posterior().covariance_relative_error_estimate(),
         "H0 absolute variance error uses full covariance norm");
  } else
    need(ah.numerical_status == N::conditioning_budget_exceeded,
         "anisotropic H0 conservative refusal");
  need(!LadderPredictive::preparation_payload_bound(LadderPosterior{}, ng, f,
                                                    md),
       "invalid optional owner bound safe");
  // Enumerate every actual preparation allocation with immutable arguments
  // built before arming, preserving the source on each failed move preparation.
  size_t allocated = 0;
  {
    auto src = noise(m);
    auto input = m;
    auto pin = pr;
    calls = 0;
    fail_on = SIZE_MAX;
    armed = true;
    auto ok = LadderPosterior::prepare(std::move(src), std::move(input),
                                       std::move(pin));
    armed = false;
    allocated = calls;
    need(ok.status() == DensityStatus::finite,
         "measured posterior allocation witness");
  }
  need(allocated > 0, "posterior allocation sites measured");
  for (size_t i = 1; i <= allocated; ++i) {
    auto src = noise(m);
    auto input = m;
    auto pin = pr;
    calls = 0;
    fail_on = i;
    armed = true;
    auto bad = LadderPosterior::prepare(std::move(src), std::move(input),
                                        std::move(pin));
    armed = false;
    need(bad.status() != DensityStatus::finite &&
             src.status() == DensityStatus::finite,
         "allocation failure preserves source");
  }
  allocated = 0;
  {
    auto input = f;
    auto min = md;
    calls = 0;
    fail_on = SIZE_MAX;
    armed = true;
    auto ok = LadderPredictive::prepare(pmove, y, rows, ng, std::move(input),
                                        std::move(min));
    armed = false;
    allocated = calls;
    need(ok.status() == DensityStatus::finite,
         "measured future allocation witness");
  }
  for (size_t i = 1; i <= allocated; ++i) {
    auto input = f;
    auto min = md;
    calls = 0;
    fail_on = i;
    armed = true;
    auto bad = LadderPredictive::prepare(pmove, y, rows, ng, std::move(input),
                                         std::move(min));
    armed = false;
    need(bad.status() != DensityStatus::finite &&
             pmove.status() == DensityStatus::finite &&
             ng.status() == DensityStatus::finite,
         "future allocation preserves owners");
  }
}
void repeated_conditioning() {
  auto m = model(), f = model(true);
  auto g = noise(m);
  auto p = LadderPosterior::prepare(std::move(g), m, prior(m));
  auto ng = noise(f);
  auto md = metadata(p, f);
  auto moved_f = f;
  auto moved_md = md;
  const auto bound = LadderPredictiveConditioning::preparation_payload_bound(
      p, ng, moved_f, moved_md);
  PredictivePolicy limit;
  limit.maximum_payload_bytes = *bound - 1;
  calls = 0;
  fail_on = SIZE_MAX;
  armed = true;
  auto refused = LadderPredictiveConditioning::prepare(
      std::move(p), ng, std::move(moved_f), std::move(moved_md), limit);
  armed = false;
  need(refused.status() != DensityStatus::finite && calls == 0 &&
           p.status() == DensityStatus::finite,
       "repeated ladder peak pre-admission preserves source and allocates "
       "nothing");
  auto retained =
      LadderPredictiveConditioning::prepare(std::move(p), ng, f, md);
  need(retained.status() == DensityStatus::finite &&
           p.status() != DensityStatus::finite,
       "repeated ladder consumes posterior once");
  ng = Gaussian{};
  const auto rows = retained.posterior().linearization().ordered_row_ids;
  const auto future_rows = retained.linearization().ordered_row_ids;
  const std::vector<double> y{32, 30, 31},
      fy{24, 22, 26, 24, 24, std::nextafter(16., 0.)};
  auto batch = retained.evaluate(y, rows, fy, future_rows, 3);
  need(batch.status == DensityStatus::finite &&
           batch.rows[0].status == DensityStatus::finite &&
           batch.rows[1].status == DensityStatus::finite &&
           batch.rows[2].status == DensityStatus::numerical_failure &&
           batch.rows[2].numerical_status == N::conditioning_budget_exceeded,
       "original future rounded offset refusal retained without clipping");
  near(batch.means[0], (31.L + 2 * 32) / 3 - 7.875L);
  near(batch.means[1], 21.625L);
  near(retained.covariance()[0], 23.L / 3);
  near(retained.covariance()[1], 2.25L);
  near(retained.covariance()[3], 7.L);
  auto original_ng = noise(f);
  for (size_t i = 0; i < 3; ++i) {
    auto complete = LadderPredictive::prepare(retained.posterior(),
                                              std::span(y).subspan(i, 1), rows,
                                              original_ng, f, md);
    auto d = complete.log_density(std::span(fy).subspan(i * 2, 2), future_rows);
    need(batch.rows[i].status == d.density.status &&
             batch.rows[i].numerical_status == d.density.numerical_status,
         "complete/retained ladder original refusal equality");
    if (batch.rows[i].status == DensityStatus::finite) {
      for (size_t j = 0; j < 2; ++j)
        need(batch.means[i * 2 + j] == complete.mean()[j] &&
                 batch.mean_absolute_error_estimates[i * 2 + j] ==
                     complete.mean_absolute_error_estimates()[j],
             "exact shared offset-mean arithmetic");
      need(batch.densities[i].density.log_value == d.density.log_value &&
               batch.densities[i].estimated_forward_sensitivity ==
                   d.estimated_forward_sensitivity,
           "exact same normalized ladder density and diagnostic");
    }
  }
  auto means = retained.evaluate(y, rows, {}, {}, 3, {true, false});
  const size_t means_work = 3 * (32 * (1 + 8 + 64 + 16) + 1 + 2);
  need(means.work_units == means_work &&
           batch.work_units == means_work + 3 * (32 * 4 + 2) &&
           means.rows[2].status == DensityStatus::finite &&
           means.densities.empty(),
       "means-only charges training and mean translations; densities charge "
       "future translation separately");
  limit = {};
  limit.maximum_work_units = means_work - 1;
  need(retained.evaluate(y, rows, {}, {}, 3, {true, false}, limit).rows.empty(),
       "means translation work cannot bypass cumulative quota");
  limit.maximum_work_units = means_work;
  need(retained.evaluate(y, rows, {}, {}, 3, {true, false}, limit).status ==
           DensityStatus::finite,
       "exact means-only cumulative work boundary");
  limit = {};
  limit.maximum_forward_sensitivity = std::numeric_limits<double>::infinity();
  need(retained.evaluate(y, rows, fy, future_rows, 3, {}, limit).rows.empty(),
       "infinite requested policy refused before restriction");
  const auto peak = retained.batch_payload_bound(3);
  limit = {};
  limit.maximum_payload_bytes = *peak - 1;
  armed = true;
  calls = 0;
  fail_on = SIZE_MAX;
  auto quota = retained.evaluate(y, rows, fy, future_rows, 3, {}, limit);
  armed = false;
  need(quota.rows.empty() && quota.numerical_status == N::work_limit &&
           calls == 0,
       "ladder batch payload pre-admission");
  limit.maximum_payload_bytes = *peak;
  need(retained.evaluate(y, rows, fy, future_rows, 3, {}, limit).status ==
           DensityStatus::finite,
       "exact ladder batch byte boundary");
  auto moved = std::move(retained);
  need(retained.status() != DensityStatus::finite,
       "moved repeated ladder owner invalid");
  moved = std::move(moved);
  need(moved.evaluate(y, rows, fy, future_rows, 3).rows[0].status ==
           DensityStatus::finite,
       "self move retains ladder owner after original noise destruction");
  armed = true;
  calls = 0;
  fail_on = SIZE_MAX;
  auto measured = moved.evaluate(std::span(y).first(2), rows,
                                 std::span(fy).first(4), future_rows, 2);
  armed = false;
  const auto sites = calls.load();
  for (size_t point = 1; point <= sites; ++point) {
    const auto before = live.load();
    {
      armed = true;
      calls = 0;
      fail_on = point;
      auto rejected = moved.evaluate(std::span(y).first(2), rows,
                                     std::span(fy).first(4), future_rows, 2);
      armed = false;
      bool failed = rejected.status != DensityStatus::finite;
      for (auto &row : rejected.rows)
        failed |= row.status != DensityStatus::finite &&
                  row.numerical_status == N::work_limit;
      need(failed, "every measured ladder batch allocation failure explicit");
    }
    need(live.load() == before, "failed ladder batches leak no storage");
  }
}
} // namespace
int main() {
  try {
    run();
    repeated_conditioning();
    std::cout << "calibration predictive owner: " << checks << " checks\n";
    return 0;
  } catch (const std::exception &e) {
    armed = false;
    std::cerr << e.what() << "\n";
    return 1;
  }
}
