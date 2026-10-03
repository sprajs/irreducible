// Whole consumer admission/lifetime controls; native parity is interface
// evidence.
#include "gaussian_predictive_abi_fixture.hpp"
#include "irred/gaussian_predictive.hpp"
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>
using namespace irred;
using namespace predictive_fixture;
#ifdef IRRED_PREDICTIVE_WORK_OBSERVATION
namespace observed {
inline bool active = false;
inline size_t factors = 0, whitenings = 0, solves = 0;
} // namespace observed
extern "C" irred::numerics::Factorization
__real__ZN5irred8numerics8choleskyESt4spanIKdLm18446744073709551615EEmmNS0_10ArithmeticE(
    std::span<const double>, size_t, size_t, irred::numerics::Arithmetic);
extern "C" irred::numerics::Factorization
__wrap__ZN5irred8numerics8choleskyESt4spanIKdLm18446744073709551615EEmmNS0_10ArithmeticE(
    std::span<const double> a, size_t n, size_t cap,
    irred::numerics::Arithmetic arithmetic) {
  if (observed::active)
    ++observed::factors;
  return __real__ZN5irred8numerics8choleskyESt4spanIKdLm18446744073709551615EEmmNS0_10ArithmeticE(
      a, n, cap, arithmetic);
}
#endif
// Track requested heap bytes, excluding this header/allocator overhead/RSS.
// Atomic observable state and noinline operators survive Release optimization.
namespace allocation {
struct alignas(std::max_align_t) Header {
  size_t bytes;
};
std::atomic<size_t> live = 0, calls = 0, peak = 0, baseline = 0, fail_on = 0;
std::atomic<bool> armed = false;
void start(size_t fail = 0) {
  baseline = live.load();
  calls = 0;
  peak = 0;
  fail_on = fail;
  armed = true;
}
} // namespace allocation
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(size_t n) {
  if (allocation::armed && ++allocation::calls == allocation::fail_on)
    throw std::bad_alloc();
  if (n > SIZE_MAX - sizeof(allocation::Header))
    throw std::bad_alloc();
  auto *h = static_cast<allocation::Header *>(
      std::malloc(sizeof(allocation::Header) + (n ? n : 1)));
  if (!h)
    throw std::bad_alloc();
  h->bytes = n;
  const auto current = allocation::live.fetch_add(n) + n;
  if (allocation::armed && current >= allocation::baseline &&
      current - allocation::baseline > allocation::peak)
    allocation::peak = current - allocation::baseline;
  return h + 1;
}
NOINLINE void *operator new[](size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept {
  if (p) {
    auto *h = static_cast<allocation::Header *>(p) - 1;
    allocation::live -= h->bytes;
    std::free(h);
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
void need(bool x, const char *s) {
  if (!x)
    throw std::runtime_error(s);
}
irred_gaussian_predictive_view view(const Result &r) {
  irred_gaussian_predictive_view v{};
  need(irred_gaussian_predictive_result_view(r.p, &v) == IRRED_OK,
       "owned view");
  return v;
}
void near(double a, long double b, long double allocation = 2e-12L) {
  need(std::isfinite(a) && std::abs(a - b) <= allocation * (1 + std::abs(b)),
       "frozen named allocation");
}
void run(const irred_gaussian_predictive_batch &b,
         const irred_gaussian_posterior_policy &q, Result &r) {
  need(irred_gaussian_predictive_evaluate(&b, &q, &r.p) == IRRED_OK,
       "one coarse ABI call");
}
void withheld(const irred_gaussian_predictive_view &v) {
  need(v.status != IRRED_GAUSSIAN_STATUS_FINITE && v.case_count == 0 && !v.rows,
       "global refusal has no science rows");
}
} // namespace
int main() {
  try {
    Fixture f;
    auto b = f.batch();
    const auto q = policy();
    Result original;
#ifdef IRRED_PREDICTIVE_WORK_OBSERVATION
    observed::active = true;
    observed::factors = 0;
#endif
    run(b, q, original);
#ifdef IRRED_PREDICTIVE_WORK_OBSERVATION
    observed::active = false;
    need(observed::factors == 6, "one C/S/P/V/R/W factorization each");
#endif
    const auto v = view(original);
    need(v.status == IRRED_GAUSSIAN_STATUS_FINITE && v.numerical_status == 0 &&
             v.case_count == 2 && v.phase == 5,
         "original coarse success");
    need(v.training_noise_attempted == 1 && v.training_noise_prepared == 1 &&
             v.posterior_attempted == 1 && v.posterior_prepared == 1 &&
             v.future_noise_attempted == 1 && v.future_noise_prepared == 1 &&
             v.predictive_attempted == 1 && v.predictive_prepared == 1 &&
             v.batch_called == 1 && v.batch_admission_completed == 1,
         "complete causal phase receipt");
    need(v.declared_training_noise_work == 256 &&
             v.declared_posterior_work == 2560 &&
             v.declared_future_noise_work == 256 &&
             v.declared_predictive_work == 2944 &&
             v.declared_batch_work == 1280 && v.declared_total_work == 7296 &&
             v.pooled_numeric_elements == 50,
         "fixed whole-request reservations");
    near(v.rows[0].mean.data[0], -258.L / 668);
    near(v.rows[0].mean.data[1], 855.L / 668);
    near(v.rows[1].mean.data[0], 0);
    near(v.rows[1].mean.data[1], .75L);
    const long double det = 11036863.L / 1784896;
    near(v.rows[1].log_determinant, std::log(det), 2e-11L);
    near(v.rows[1].quadratic, 0, 2e-11L);
    near(v.rows[1].log_density,
         -.5L * (std::log(det) +
                 2 * std::log(2 * std::numbers::pi_v<long double>)),
         2e-11L);
    const auto scalar_product =
        -.5L * (std::log(4444.L * 2618 / (1336.L * 1336)) +
                2 * std::log(2 * std::numbers::pi_v<long double>));
    need(v.rows[1].log_density - scalar_product > 1e-3L,
         "whole future law retains covariance witness");
    // Standalone repeated-conditioning API; same owner ancestry, exact output
    // parity.
    using namespace statistics;
    Metadata md;
    md.ordered_ids = {"train-cal0", "train-cal1"};
    md.source_semantics = "synthetic controls";
    md.measure = "product d(mag)";
    md.table_identity = "synthetic training conditional noise";
    md.calibration_provenance = "offset and slope represented once in beta";
    md.uncertainty_identity = "conditional Gaussian noise";
    md.dependence_provenance = "full training covariance; no cross-case law";
    md.ordering_provenance = "train-cal0,train-cal1 exact order";
    auto source =
        prepare_gaussian(f.C, MatrixKind::covariance, md, 1000000, 1e-10,
                         numerics::Arithmetic::longdouble_cpu_v1);
    ParameterPrior prior;
    prior.ordered_parameter_ids = {"offset", "slope"};
    prior.parameter_units = {"mag", "mag"};
    prior.shared_nuisance_ids = {"offset"};
    prior.mean = {.5, -.25};
    prior.covariance = {2, .5, .5, 1};
    prior.prior_identity = "original proper prior";
    prior.design_identity = "synthetic fixed X";
    prior.residual_unit = "mag";
    prior.parameter_measure = "d(offset-mag) d(slope-mag)";
    prior.dependence_identity = "prior independent of both noises";
    prior.noise_independence_declared = true;
    auto post = GaussianPosterior::prepare(std::move(source), f.X,
                                           md.ordered_ids, std::move(prior));
    md.ordered_ids = {"heldout-cal2", "heldout-cal3"};
    md.table_identity = "synthetic future conditional noise";
    md.ordering_provenance = "heldout-cal2,heldout-cal3 exact order";
    md.dependence_provenance =
        "full future noise independent of prior and training noise";
    auto noise =
        prepare_gaussian(f.R, MatrixKind::covariance, md, 1000000, 1e-10,
                         numerics::Arithmetic::longdouble_cpu_v1);
    PredictiveMetadata pm;
    pm.ordered_parameter_ids = {"offset", "slope"};
    pm.parameter_units = {"mag", "mag"};
    pm.response_units = {"mag/mag", "mag/mag"};
    pm.training_event_ids = {"synthetic-train-0", "synthetic-train-1"};
    pm.future_event_ids = {"synthetic-future-2", "synthetic-future-3"};
    pm.future_unit = "mag";
    pm.future_covariance_unit = "mag^2";
    pm.future_measure = "product d(mag)";
    pm.response_identity = "synthetic fixed A";
    pm.conditioning_identity = "fixed two-case family";
    pm.conditional_noise_identity = "synthetic future conditional noise";
    pm.dependence_identity =
        "shared beta and full independent future noise; no cross-case law";
    pm.future_noise_independence_declared = true;
    pm.noise_conditional_on_parameters_declared = true;
    auto owner = GaussianPredictiveConditioning::prepare(std::move(post), noise,
                                                         f.A, pm);
    need(owner.status() == DensityStatus::finite, "standalone retained owner");
    const std::array<std::string, 2> training_ids{"train-cal0", "train-cal1"};
    auto batch =
        owner.evaluate(f.training, training_ids, f.future, md.ordered_ids, 2);
    need(batch.status == DensityStatus::finite && batch.rows.size() == 2,
         "standalone borrowed coarse batch");
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        need(v.rows[i].mean.data[j] == batch.means[i * 2 + j],
             "native mean parity");
        need(v.rows[i].absolute_error_estimates.data[j] ==
                 batch.mean_absolute_error_estimates[i * 2 + j],
             "native error parity");
      }
      need(v.rows[i].log_density == batch.densities[i].density.log_value &&
               v.rows[i].quadratic == batch.densities[i].quadratic &&
               v.rows[i].log_determinant == batch.densities[i].log_determinant,
           "native joint density parity");
    }
    for (unsigned mask = 1; mask <= 3; ++mask) {
      auto changed = b;
      changed.requested_outputs = mask;
      if (mask == 1) {
        changed.future_vectors = values(std::array<double, 0>{});
        changed.future_vector_row_ids = {};
        changed.future_vector_event_ids = {};
      }
      Result r;
      run(changed, q, r);
      auto x = view(r);
      need(x.status == IRRED_GAUSSIAN_STATUS_FINITE && x.case_count == 2,
           "requested masks");
      for (size_t i = 0; i < 2; ++i)
        need(x.rows[i].mean.length == ((mask & 1) ? 2 : 0) &&
                 x.rows[i].absolute_error_estimates.length ==
                     ((mask & 1) ? 2 : 0) &&
                 x.rows[i].joint_density_available == ((mask & 2) ? 1 : 0),
             "unused outputs not retained");
    }
    const std::array<double, 4> zero{0, 0, 0, 0}, bad{1, 2, 2, 1};
    for (unsigned which = 0; which < 2; ++which) {
      auto changed = b;
      if (which == 0)
        changed.future_response = values(zero);
      else
        changed.training.design = values(zero);
      Result r;
      run(changed, q, r);
      const auto x = view(r);
      need(x.status == IRRED_GAUSSIAN_STATUS_FINITE,
           "null response/design remain proper");
      near(x.rows[1].mean.data[0], 0);
      near(x.rows[1].mean.data[1], which == 0 ? 0 : .75L);
      if (which == 0)
        near(x.rows[0].log_determinant, std::log(.5L * .75L - .125L * .125L),
             2e-11L);
    }
    for (unsigned stage = 0; stage < 3; ++stage) {
      auto changed = b;
      if (stage == 0)
        changed.training.noise_covariance = values(bad);
      if (stage == 1)
        changed.training.prior_covariance = values(bad);
      if (stage == 2)
        changed.future_noise_covariance = values(bad);
      Result r;
      run(changed, q, r);
      auto x = view(r);
      withheld(x);
      need(x.numerical_status == IRRED_NUMERICAL_STATUS_NOT_POSITIVE_DEFINITE &&
               x.phase == stage + 1 && !x.batch_called,
           "failed SPD phase causal receipt");
      need(x.training_noise_prepared == (stage > 0) &&
               x.posterior_prepared == (stage > 1) &&
               x.future_noise_prepared == 0,
           "failed stage not marked completed");
    }
    for (unsigned variant = 0; variant < 5; ++variant) {
      auto changed = b;
      if (variant == 0)
        changed.future_noise_independence_declared = 0;
      if (variant == 1)
        changed.noise_conditional_on_parameters_declared = 0;
      if (variant == 2)
        changed.future_covariance_unit = text("wrong^2");
      if (variant == 3)
        changed.future_noise_event_ids = changed.future_vector_event_ids =
            ids(f.events);
      const std::array<irred_bytes, 2> reversed{f.parameters[1],
                                                f.parameters[0]};
      if (variant == 4)
        changed.future_parameter_ids = ids(reversed);
      Result r;
      run(changed, q, r);
      const auto x = view(r);
      withheld(x);
      need(x.status == IRRED_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA &&
               !x.training_noise_attempted,
           "lineage/domain admission before native equations");
    }
    for (unsigned variant = 0; variant < 6; ++variant) {
      auto changed = b;
      if (variant == 0)
        changed.reserved = 1;
      if (variant == 1)
        changed.requested_outputs = 0;
      if (variant == 2)
        changed.future_response.byte_length -= 1;
      if (variant == 3)
        changed.future_response.data = nullptr;
      if (variant == 4)
        changed.future_response.reserved = 1;
      const std::array<uint8_t, 2> invalid{0xc0, 0x80};
      if (variant == 5)
        changed.future_noise_identity = {invalid.data(), invalid.size()};
      Result r;
      need(irred_gaussian_predictive_evaluate(&changed, &q, &r.p) ==
                   IRRED_INVALID_ARGUMENT &&
               !r.p,
           "invalid ABI wire never executed");
    }
    {
      auto changed = b;
      changed.training.case_count = 0;
      changed.training.case_ids = {};
      changed.training.conditioning_vectors = values(std::array<double, 0>{});
      changed.future_vectors = values(std::array<double, 0>{});
      Result r;
      run(changed, q, r);
      auto x = view(r);
      withheld(x);
      need(x.phase == 0 && !x.training_noise_attempted && !x.batch_called,
           "empty batch is causal admission refusal");
    }
    for (unsigned cap = 0; cap < 5; ++cap) {
      auto limited = q;
      if (cap == 0)
        limited.maximum_work_units = v.declared_total_work - 1;
      if (cap == 1)
        limited.maximum_elements = v.pooled_numeric_elements - 1;
      if (cap == 2)
        limited.maximum_cases = 1;
      if (cap == 3)
        limited.maximum_string_bytes = 0;
      if (cap == 4)
        limited.maximum_native_bytes = v.minimum_result_bytes;
      Result r;
      run(b, limited, r);
      auto x = view(r);
      withheld(x);
      need(x.numerical_status == IRRED_NUMERICAL_STATUS_WORK_LIMIT &&
               !x.training_noise_attempted,
           "original quota before scientific allocation");
    }
    for (uint64_t cap : {uint64_t(0), v.minimum_result_bytes - 1}) {
      auto limited = q;
      limited.maximum_native_bytes = cap;
      Result r;
      allocation::start();
      const auto code = irred_gaussian_predictive_evaluate(&b, &limited, &r.p);
      allocation::armed = false;
      need(code == IRRED_GAUSSIAN_PREDICTIVE_QUOTA_REFUSED && !r.p &&
               allocation::calls == 0 && allocation::peak == 0,
           "minimal refusal allocation free");
    }
    {
      auto changed = b;
      const std::array<double, 4> dynamic{
          1.25, -.5, std::numeric_limits<double>::denorm_min(), 0};
      changed.training.conditioning_vectors = values(dynamic);
      Result r;
      run(changed, q, r);
      const auto x = view(r);
      need(x.status == IRRED_GAUSSIAN_STATUS_FINITE &&
               x.rows[0].status == IRRED_GAUSSIAN_STATUS_FINITE &&
               x.rows[1].status != IRRED_GAUSSIAN_STATUS_FINITE &&
               !x.rows[1].mean.data &&
               !x.rows[1].absolute_error_estimates.data &&
               x.rows[1].joint_density_available == 0,
           "retained nondropped per-case refusal");
    }
    // Stored policy must reserve setup AND every case, not only setup work.
    {
      std::array<double, 64> training{}, future{};
      std::array<irred_bytes, 32> cases{};
      std::array<std::string, 32> names{};
      for (size_t i = 0; i < 32; ++i) {
        training[2 * i] = .5;
        training[2 * i + 1] = .25;
        future[2 * i] = 0;
        future[2 * i + 1] = .75;
        names[i] = "case-" + std::to_string(i);
        cases[i] = text(names[i]);
      }
      auto changed = b;
      changed.training.case_count = 32;
      changed.training.case_ids = ids(cases);
      changed.training.conditioning_vectors = values(training);
      changed.future_vectors = values(future);
      Result r;
#ifdef IRRED_PREDICTIVE_WORK_OBSERVATION
      observed::active = true;
      observed::factors = 0;
#endif
      run(changed, q, r);
#ifdef IRRED_PREDICTIVE_WORK_OBSERVATION
      observed::active = false;
      need(observed::factors == 6,
           "larger coarse batch does not refactor per case");
#endif
      const auto x = view(r);
      need(x.status == IRRED_GAUSSIAN_STATUS_FINITE && x.case_count == 32 &&
               x.declared_batch_work > x.declared_predictive_work,
           "retained policy admits larger legal batch");
      for (size_t i = 0; i < 32; ++i)
        need(x.rows[i].status == IRRED_GAUSSIAN_STATUS_FINITE,
             "every larger-batch case retained");
    }
    {
      const std::string maximum_text(256, 'x');
      auto changed = b;
      changed.future_noise_identity = changed.conditioning_identity =
          text(maximum_text);
      Result r;
      run(changed, q, r);
      need(view(r).status == IRRED_GAUSSIAN_STATUS_FINITE,
           "maximum text constructor envelope");
    }
    {
      auto changed = b;
      const std::array<double, 4> rank_deficient{1, 1, 2, 2};
      changed.training.design = values(rank_deficient);
      Result r;
      run(changed, q, r);
      need(view(r).status == IRRED_GAUSSIAN_STATUS_FINITE,
           "rank deficient design remains proper");
    }
    {
      const auto original_rounding = std::fegetround();
      need(original_rounding != -1, "rounding readable");
      need(std::fesetround(FE_DOWNWARD) == 0, "adversarial rounding available");
      Result r;
      run(b, q, r);
      const auto x = view(r);
      const auto restored = std::fesetround(original_rounding);
      need(restored == 0, "rounding restored");
      withheld(x);
      need(!x.batch_called, "unsupported arithmetic never proceeds to batch");
    }
    size_t total_calls = 0, measured_peak = 0, measured_retained = 0;
    const auto live_before = allocation::live.load();
    {
      Result r;
      allocation::start();
      const auto code = irred_gaussian_predictive_evaluate(&b, &q, &r.p);
      allocation::armed = false;
      total_calls = allocation::calls;
      measured_peak = allocation::peak;
      measured_retained = allocation::live - live_before;
      need(code == IRRED_OK && total_calls > 1 && total_calls <= 4096,
           "bounded measured allocation baseline");
      const auto x = view(r);
      need(x.status == IRRED_GAUSSIAN_STATUS_FINITE &&
               measured_peak <= x.peak_payload_bytes &&
               measured_retained <= x.retained_payload_bytes,
           "actual native payload within charged lifetime bounds");
    }
    need(allocation::live == live_before,
         "owner destruction releases original success");
    for (size_t point = 1; point <= total_calls; ++point) {
      {
        Result r;
        allocation::start(point);
        const auto code = irred_gaussian_predictive_evaluate(&b, &q, &r.p);
        allocation::armed = false;
        need(code == IRRED_OK || code == IRRED_ALLOCATION_FAILURE,
             "injected failure contained in ABI");
        if (code == IRRED_ALLOCATION_FAILURE)
          need(!r.p, "first owner failure transport null");
        else {
          const auto x = view(r);
          bool failed = x.status != IRRED_GAUSSIAN_STATUS_FINITE;
          if (failed)
            withheld(x);
          else
            for (size_t i = 0; i < x.case_count; ++i)
              if (x.rows[i].status != IRRED_GAUSSIAN_STATUS_FINITE) {
                failed = true;
                need(!x.rows[i].mean.data &&
                         !x.rows[i].absolute_error_estimates.data &&
                         !x.rows[i].joint_density_available,
                     "failed row all requested outputs absent");
              }
          need(failed, "injected allocation never silently accepted");
          need(allocation::peak <= x.peak_payload_bytes &&
                   allocation::live - live_before <= x.retained_payload_bytes,
               "every failed prefix payload charged");
        }
      }
      need(allocation::live == live_before,
           "every failure prefix releases all owners");
    }
    {
      auto changed = b;
      auto limited = q;
      limited.maximum_native_bytes = v.peak_payload_bytes - 1;
      Result r;
      allocation::start();
      const auto code =
          irred_gaussian_predictive_evaluate(&changed, &limited, &r.p);
      allocation::armed = false;
      need(code == IRRED_OK, "owned reduced-byte quota");
      const auto x = view(r);
      withheld(x);
      need(allocation::peak <= limited.maximum_native_bytes &&
               x.peak_payload_bytes <= limited.maximum_native_bytes,
           "no allocation prefix exceeds original smaller byte cap");
    }
    std::cout << "{\"control\":\"gaussian-predictive-abi\",\"declared_work\":"
              << v.declared_total_work
              << ",\"measured_allocation_calls\":" << total_calls
              << ",\"measured_peak_payload\":" << measured_peak
              << ",\"measured_retained_payload\":" << measured_retained
              << ",\"failure_prefixes\":" << total_calls << "}\n";
    return 0;
  } catch (const std::exception &e) {
    allocation::armed = false;
    std::cerr << e.what() << '\n';
    return 1;
  }
}
