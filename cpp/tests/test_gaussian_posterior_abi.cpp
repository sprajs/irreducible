// Named original rational/cofactor fixture, then standalone/ABI parity.
// Interface parity shares native ancestry; it is not independent physics.
#include "irred/abi.h"
#include "irred/gaussian_posterior.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
using namespace irred;
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
    statistics::Metadata md;
    md.ordered_ids = {"r0", "r1"};
    md.measure = "product d(mag)";
    auto g = statistics::prepare_gaussian(
        C, statistics::MatrixKind::covariance, md, 1000000, 1e-10,
        numerics::Arithmetic::longdouble_cpu_v1);
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
    for (size_t k = 0; k < 2; ++k) {
      auto mean = owner.condition({R.data() + 2 * k, 2}, md.ordered_ids);
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
      need(irred_gaussian_posterior_evaluate(&b, &limited, &r.p) ==
                   IRRED_GAUSSIAN_POSTERIOR_QUOTA_REFUSED &&
               !r.p,
           "no allocation envelope admitted");
    }
    {
      auto limited = q;
      limited.maximum_native_bytes = v.minimum_result_bytes;
      Result r;
      need(irred_gaussian_posterior_evaluate(&b, &limited, &r.p) == IRRED_OK,
           "exact failure owner quota");
      auto x = r.view();
      need(x.numerical_status == IRRED_NUMERICAL_STATUS_WORK_LIMIT &&
               x.covariance.length == 0 && x.case_count == 0 &&
               x.retained_payload_bytes == x.minimum_result_bytes,
           "charged minimal refusal");
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
    std::fesetround(FE_DOWNWARD);
    {
      Result r;
      need(irred_gaussian_posterior_evaluate(&b, &q, &r.p) == IRRED_OK &&
               r.view().status != IRRED_GAUSSIAN_STATUS_FINITE,
           "arithmetic refusal");
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
