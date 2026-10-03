// Named original rational/cofactor fixture, then standalone/ABI parity.
// Interface parity shares native ancestry; it is not independent physics.
#include "irred/abi.h"
#include "irred/gaussian_posterior.hpp"
#include <array>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
using namespace irred;
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
void need(bool x, const char *m) {
  if (!x)
    throw std::runtime_error(m);
}
irred_bytes text(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
template <size_t N> irred_strings ids(const std::array<irred_bytes, N> &a) {
  return {a.data(), N, sizeof(a)};
}
template <size_t N> irred_f64_buffer values(const std::array<double, N> &a) {
  return {sizeof(irred_f64_buffer),
          IRRED_ABI_VERSION,
          2,
          0,
          a.data(),
          N,
          sizeof(a)};
}
struct Result {
  irred_gaussian_posterior_result *p = nullptr;
  ~Result() { irred_gaussian_posterior_result_destroy(p); }
  irred_gaussian_posterior_view view() {
    irred_gaussian_posterior_view v{};
    need(irred_gaussian_posterior_result_view(p, &v) == IRRED_OK, "view");
    return v;
  }
};
void near(double a, long double b) {
  need(std::abs(a - b) <= 2e-12L * (1 + std::abs(b)),
       "original frozen moment allocation");
}
} // namespace
int main() {
  try {
    const std::array<double, 4> C{2, .25, .25, 1.5}, X{1, .5, -.25, 1},
        S{1, .375, .375, .5}, R{1.25, -.5, .375, -.375};
    const std::array<double, 2> m{.5, -.25};
    const std::array<irred_bytes, 2> row{text("r0"), text("r1")},
        event{text("e0"), text("e1")}, parameter{text("zero"), text("colour")},
        unit{text("mag"), text("mag")},
        column{text("mag/mag"), text("mag/mag")},
        cases{text("cofactor"), text("mean-response")};
    const std::array<irred_bytes, 1> shared{text("zero")};
    irred_gaussian_posterior_batch b{};
    b.struct_size = sizeof(b);
    b.abi_version = IRRED_ABI_VERSION;
    b.source_semantics = 1;
    b.noise_independence_declared = 1;
    b.case_count = 2;
    b.noise_covariance = values(C);
    b.design = values(X);
    b.prior_mean = values(m);
    b.prior_covariance = values(S);
    b.conditioning_vectors = values(R);
    b.noise_row_ids = b.design_row_ids = b.conditioning_row_ids = ids(row);
    b.noise_event_ids = b.conditioning_event_ids = ids(event);
    b.prior_parameter_ids = b.design_parameter_ids = ids(parameter);
    b.prior_parameter_units = b.design_parameter_units = ids(unit);
    b.column_units = ids(column);
    b.shared_nuisance_ids = ids(shared);
    b.case_ids = ids(cases);
    b.residual_unit = text("mag");
    b.noise_identity = text("synthetic covariance");
    b.calibration_identity = text("synthetic fixed calibration");
    b.noise_dependence_identity = text("full row covariance");
    b.ordering_provenance = text("explicit");
    b.design_identity = text("fixed X");
    b.prior_identity = text("proper synthetic prior");
    b.parameter_measure = text("d(zero) d(colour)");
    b.prior_dependence_identity = text("independent original prior and noise");
    irred_gaussian_posterior_policy q{
        sizeof(q), IRRED_ABI_VERSION, 1,         0,         1000000,
        65536,     1048576,           268435456, 100000000, 1e-10};
    Result result;
    need(irred_gaussian_posterior_evaluate(&b, &q, &result.p) == IRRED_OK,
         "one coarse call");
    auto v = result.view();
    need(v.status == IRRED_GAUSSIAN_STATUS_FINITE && v.case_count == 2 &&
             v.covariance.length == 4,
         "moments");
    need(v.declared_work_units == 3584 && v.pooled_numeric_elements == 30,
         "cumulative admission arithmetic");
    const long double cov[]{15160.L / 25511, 4466.L / 25511, 4466.L / 25511,
                            8592.L / 25511};
    for (size_t i = 0; i < 4; ++i)
      near(v.covariance.data[i], cov[i]);
    near(v.rows[0].mean.data[0], 20604.L / 25511);
    near(v.rows[0].mean.data[1], -7125.L / 51022);
    near(v.rows[1].mean.data[0], .5L);
    near(v.rows[1].mean.data[1], -.25L);
    need(v.peak_payload_bytes <= q.maximum_native_bytes &&
             v.retained_payload_bytes <= v.peak_payload_bytes,
         "charged phase bounds");
    need(v.noise_preparation_attempted == 1 && v.noise_prepared == 1 &&
             v.posterior_preparation_attempted == 1 &&
             v.posterior_prepared == 1 && v.conditioning_cases_attempted == 2,
         "native attempted and completed phases");
    statistics::Metadata md;
    md.ordered_ids = {"r0", "r1"};
    md.measure = "product d(mag)";
    md.ordering_provenance = "explicit";
    auto g = statistics::prepare_gaussian(
        C, statistics::MatrixKind::covariance, md, 1000000, 1e-10,
        numerics::Arithmetic::longdouble_cpu_v1);
    need(g.status() == statistics::DensityStatus::finite,
         "standalone source accepted original order");
    statistics::ParameterPrior prior{{"zero", "colour"},
                                     {"mag", "mag"},
                                     {"zero"},
                                     {.5, -.25},
                                     {1, .375, .375, .5},
                                     "proper synthetic prior",
                                     "fixed X",
                                     "mag",
                                     "d(zero) d(colour)",
                                     "independent original prior and noise",
                                     true};
    auto owner = statistics::GaussianPosterior::prepare(
        std::move(g), X, md.ordered_ids, std::move(prior));
    need(owner.status() == statistics::DensityStatus::finite,
         "standalone proper posterior finite");
    for (size_t k = 0; k < 2; ++k) {
      auto mean = owner.condition({R.data() + 2 * k, 2}, md.ordered_ids);
      need(mean.status == statistics::DensityStatus::finite &&
               mean.value.size() == 2 &&
               mean.absolute_error_estimates.size() == 2,
           "standalone finite mean and diagnostics before indexing");
      need(v.rows[k].status == static_cast<uint32_t>(mean.status),
           "native status parity");
      for (size_t j = 0; j < 2; ++j) {
        need(v.rows[k].mean.data[j] == mean.value[j],
             "native exact mean parity");
        need(v.rows[k].absolute_error_estimates.data[j] ==
                 mean.absolute_error_estimates[j],
             "native exact error parity");
      }
    }
    need(v.covariance_relative_error_estimate ==
             owner.covariance_relative_error_estimate(),
         "native covariance diagnostic parity");
    for (auto cap : {uint64_t(0), v.minimum_result_bytes - 1}) {
      auto limited = q;
      limited.maximum_native_bytes = cap;
      Result r;
      allocation::start();
      const auto transport =
          irred_gaussian_posterior_evaluate(&b, &limited, &r.p);
      allocation::armed = false;
      need(transport == IRRED_GAUSSIAN_POSTERIOR_QUOTA_REFUSED && !r.p,
           "no allocation envelope admitted");
      need(allocation::calls == 0 && allocation::peak == 0,
           "below minimum constructor really performs zero heap allocations");
    }
    {
      auto limited = q;
      limited.maximum_native_bytes = v.minimum_result_bytes;
      Result r;
      allocation::start();
      const auto transport =
          irred_gaussian_posterior_evaluate(&b, &limited, &r.p);
      allocation::armed = false;
      need(transport == IRRED_OK, "exact failure owner quota");
      auto x = r.view();
      need(x.numerical_status == IRRED_NUMERICAL_STATUS_WORK_LIMIT &&
               x.covariance.length == 0 && x.case_count == 0 &&
               x.retained_payload_bytes == x.minimum_result_bytes,
           "charged minimal refusal");
      need(
          allocation::calls == 1 &&
              allocation::peak == x.minimum_result_bytes &&
              allocation::live - allocation::baseline == x.minimum_result_bytes,
          "disengaged minimal owner has no hidden Gaussian string allocations");
    }
    for (unsigned field = 0; field < 4; ++field) {
      auto limited = q;
      if (field == 0)
        limited.maximum_work_units = 3583;
      if (field == 1)
        limited.maximum_elements = 29;
      if (field == 2)
        limited.maximum_cases = 1;
      if (field == 3)
        limited.maximum_string_bytes = 0;
      Result r;
      need(irred_gaussian_posterior_evaluate(&b, &limited, &r.p) == IRRED_OK,
           "owned quota refusal");
      auto x = r.view();
      need(x.numerical_status == IRRED_NUMERICAL_STATUS_WORK_LIMIT &&
               x.covariance.length == 0 && x.case_count == 0,
           "global moments withheld");
    }
    {
      auto changed = b;
      changed.noise_independence_declared = 0;
      Result r;
      need(irred_gaussian_posterior_evaluate(&changed, &q, &r.p) == IRRED_OK,
           "scientific declaration refusal");
      need(r.view().status == IRRED_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA,
           "independence not inferred");
    }
    {
      auto changed = b;
      const std::array<irred_bytes, 2> reverse{row[1], row[0]};
      changed.conditioning_row_ids = ids(reverse);
      Result r;
      need(irred_gaussian_posterior_evaluate(&changed, &q, &r.p) == IRRED_OK &&
               r.view().status == IRRED_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA,
           "exact original row order");
    }
    {
      auto changed = b;
      const std::array<double, 4> badR{
          1.25, -.5, std::numeric_limits<double>::denorm_min(), 0};
      changed.conditioning_vectors = values(badR);
      Result r;
      need(irred_gaussian_posterior_evaluate(&changed, &q, &r.p) == IRRED_OK,
           "ordered per-case refusal");
      auto x = r.view();
      need(x.status == IRRED_GAUSSIAN_STATUS_FINITE &&
               x.rows[0].status == IRRED_GAUSSIAN_STATUS_FINITE &&
               x.rows[1].status != IRRED_GAUSSIAN_STATUS_FINITE &&
               x.rows[1].mean.length == 0 &&
               x.rows[1].absolute_error_estimates.length == 0,
           "failed numeric views withheld and neighbor retained");
    }
    {
      auto changed = b;
      const std::array<double, 4> zero{0, 0, 0, 0};
      changed.design = values(zero);
      Result r;
      need(irred_gaussian_posterior_evaluate(&changed, &q, &r.p) == IRRED_OK &&
               r.view().status == IRRED_GAUSSIAN_STATUS_FINITE,
           "null design proper law");
      auto x = r.view();
      for (size_t i = 0; i < 4; ++i)
        near(x.covariance.data[i], S[i]);
    }
    {
      auto changed = b;
      const std::array<double, 4> bad{1, 2, 2, 1};
      changed.prior_covariance = values(bad);
      Result r;
      need(irred_gaussian_posterior_evaluate(&changed, &q, &r.p) == IRRED_OK &&
               r.view().numerical_status ==
                   IRRED_NUMERICAL_STATUS_NOT_POSITIVE_DEFINITE,
           "no prior jitter");
      auto x = r.view();
      need(x.noise_prepared == 1 && x.posterior_preparation_attempted == 1 &&
               x.posterior_prepared == 0 && x.conditioning_cases_attempted == 0,
           "failed posterior attempt retained as attempted, not completed");
    }
    // Entire native call: result constructor, metadata, preparation, all
    // conditioning and retained owners. Bookkeeping headers are excluded.
    size_t total_calls = 0;
    const size_t live_before = allocation::live;
    {
      Result measured;
      allocation::start();
      const auto transport =
          irred_gaussian_posterior_evaluate(&b, &q, &measured.p);
      allocation::armed = false;
      total_calls = allocation::calls;
      need(transport == IRRED_OK, "whole-phase allocation run");
      const auto measured_view = measured.view();
      need(measured_view.status == IRRED_GAUSSIAN_STATUS_FINITE &&
               total_calls > 1,
           "whole-phase allocation witness completes actual consumer");
      need(allocation::peak <= measured_view.peak_payload_bytes &&
               allocation::live - live_before <=
                   measured_view.retained_payload_bytes,
           "measured requested peak/retained heap fits declared payload scope");
    }
    need(allocation::live == live_before,
         "whole owner destruction releases heap");
    for (size_t point = 1; point <= total_calls; ++point) {
      {
        Result rejected;
        allocation::start(point);
        const auto transport =
            irred_gaussian_posterior_evaluate(&b, &q, &rejected.p);
        allocation::armed = false;
        need(transport == IRRED_ALLOCATION_FAILURE || transport == IRRED_OK,
             "injected allocation failure stays inside ABI");
        if (transport == IRRED_ALLOCATION_FAILURE)
          need(!rejected.p, "allocation transport has no owner");
        else {
          const auto x = rejected.view();
          bool refused = x.status != IRRED_GAUSSIAN_STATUS_FINITE;
          if (refused)
            need(x.case_count == 0 && x.covariance.length == 0,
                 "global allocation refusal withholds all moments");
          else
            for (size_t i = 0; i < x.case_count; ++i)
              if (x.rows[i].status != IRRED_GAUSSIAN_STATUS_FINITE) {
                refused = true;
                need(x.rows[i].mean.length == 0 &&
                         x.rows[i].absolute_error_estimates.length == 0,
                     "failed conditioning allocation has no numeric row "
                     "payload");
              }
          need(refused, "injected failure is not ignored");
          need(allocation::peak <= x.peak_payload_bytes &&
                   allocation::live - live_before <= x.retained_payload_bytes,
               "refused attempt actual retained capacity still charged");
        }
      }
      need(allocation::live == live_before,
           "every failure site releases complete owner");
    }
    {
      auto changed = b;
      const std::array<double, 4> diagonal_C{1, 0, 0, 1},
          diagonal_S{100, 0, 0, 100}, diagonal_X{.25, 0, 0, .25},
          overflow_R{1.25, -.5, std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};
      changed.noise_covariance = values(diagonal_C);
      changed.prior_covariance = values(diagonal_S);
      changed.design = values(diagonal_X);
      changed.conditioning_vectors = values(overflow_R);
      auto late_source = statistics::prepare_gaussian(
          diagonal_C, statistics::MatrixKind::covariance, md, 1000000, 1e-10,
          numerics::Arithmetic::longdouble_cpu_v1);
      need(late_source.status() == statistics::DensityStatus::finite,
           "late-capacity source finite");
      statistics::ParameterPrior late_prior{
          {"zero", "colour"},
          {"mag", "mag"},
          {"zero"},
          {.5, -.25},
          {100, 0, 0, 100},
          "proper synthetic prior",
          "fixed X",
          "mag",
          "d(zero) d(colour)",
          "independent original prior and noise",
          true};
      auto late_owner = statistics::GaussianPosterior::prepare(
          std::move(late_source), diagonal_X, md.ordered_ids,
          std::move(late_prior));
      need(late_owner.status() == statistics::DensityStatus::finite,
           "late-capacity posterior finite");
      const auto failed_mean =
          late_owner.condition({overflow_R.data() + 2, 2}, md.ordered_ids);
      need(failed_mean.status != statistics::DensityStatus::finite &&
               failed_mean.value.empty() &&
               failed_mean.absolute_error_estimates.capacity() >= 2,
           "original late cast failure retains allocated error capacity "
           "witness");
      Result r;
      allocation::start();
      const auto transport =
          irred_gaussian_posterior_evaluate(&changed, &q, &r.p);
      allocation::armed = false;
      need(transport == IRRED_OK, "late mean cast refusal transported");
      const auto x = r.view();
      need(x.status == IRRED_GAUSSIAN_STATUS_FINITE && x.case_count == 2 &&
               x.rows[0].status == IRRED_GAUSSIAN_STATUS_FINITE &&
               x.rows[1].status != IRRED_GAUSSIAN_STATUS_FINITE &&
               x.rows[1].mean.length == 0 &&
               x.rows[1].absolute_error_estimates.length == 0 &&
               x.conditioning_cases_attempted == 2,
           "last late failure preserves original neighbor/order and withholds "
           "failed numerics");
      need(allocation::peak <= x.peak_payload_bytes &&
               allocation::live - allocation::baseline <=
                   x.retained_payload_bytes,
           "late failed mean retained vectors included in measured bytes");
      auto limited = q;
      limited.maximum_native_bytes = x.peak_payload_bytes - 1;
      Result quota;
      need(irred_gaussian_posterior_evaluate(&changed, &limited, &quota.p) ==
                   IRRED_OK &&
               quota.view().status != IRRED_GAUSSIAN_STATUS_FINITE &&
               quota.view().case_count == 0 &&
               quota.view().covariance.length == 0,
           "unchanged late-failure inputs under low global cap withhold every "
           "moment");
    }
    {
      auto changed = b;
      changed.abi_version = 0;
      Result r;
      need(irred_gaussian_posterior_evaluate(&changed, &q, &r.p) ==
                   IRRED_ABI_MISMATCH &&
               !r.p,
           "ABI mismatch");
      changed = b;
      changed.design.byte_length -= 1;
      need(irred_gaussian_posterior_evaluate(&changed, &q, &r.p) ==
                   IRRED_INVALID_INPUT &&
               !r.p,
           "byte shape hostile");
    }
    for (unsigned mutation = 0; mutation < 7; ++mutation) {
      auto changed = b;
      auto policy = q;
      if (mutation == 0)
        policy.reserved = 1;
      if (mutation == 1)
        changed.design.reserved = 1;
      if (mutation == 2)
        changed.design.length = UINT64_MAX;
      if (mutation == 3)
        changed.noise_row_ids.length = UINT64_MAX;
      if (mutation == 4)
        changed.case_count = UINT64_MAX;
      if (mutation == 5)
        changed.noise_row_ids.byte_length -= 1;
      if (mutation == 6)
        changed.design.data = reinterpret_cast<const double *>(uintptr_t(1));
      Result r;
      allocation::start();
      const auto transport =
          irred_gaussian_posterior_evaluate(&changed, &policy, &r.p);
      allocation::armed = false;
      need(transport == IRRED_INVALID_INPUT && !r.p && allocation::calls == 0,
           "hostile count/overflow/reserved/alignment rejected before "
           "allocation or dereference");
    }
    {
      Result r;
      need(irred_gaussian_posterior_evaluate(
               reinterpret_cast<const irred_gaussian_posterior_batch *>(
                   uintptr_t(1)),
               &q, &r.p) == IRRED_INVALID_INPUT &&
               !r.p,
           "misaligned top-level descriptor never dereferenced");
      need(irred_gaussian_posterior_evaluate(
               &b, &q,
               reinterpret_cast<irred_gaussian_posterior_result **>(
                   uintptr_t(1))) == IRRED_INVALID_INPUT,
           "misaligned output never dereferenced");
    }
    std::fesetround(FE_DOWNWARD);
    {
      Result r;
      need(irred_gaussian_posterior_evaluate(&b, &q, &r.p) == IRRED_OK &&
               r.view().status != IRRED_GAUSSIAN_STATUS_FINITE,
           "arithmetic refusal");
      need(r.view().noise_preparation_attempted == 1 &&
               r.view().noise_prepared == 0 &&
               r.view().posterior_preparation_attempted == 0,
           "unsupported environment preserves attempt without claiming factor "
           "complete");
    }
    std::fesetround(FE_TONEAREST);
    std::cout << "proper Gaussian posterior coarse ABI controls passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::fesetround(FE_TONEAREST);
    std::cerr << e.what() << '\n';
    return 1;
  }
}
