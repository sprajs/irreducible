// Original synthetic fixture; transport parity shares native ancestry.
#pragma once
#include "irred/abi.h"
#include <array>
#include <string_view>
namespace predictive_fixture {
inline irred_bytes text(std::string_view s) {
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
          N ? a.data() : nullptr,
          N,
          N * sizeof(double)};
}
struct Fixture {
  const std::array<double, 4> C{1, .25, .25, 2}, X{1, 0, 1, 1}, S{2, .5, .5, 1},
      R{.5, .125, .125, .75}, A{1, 2, 1, -1}, training{1.25, -.5, .5, .25},
      future{.25, -.5, 0, .75};
  const std::array<double, 2> m{.5, -.25};
  const std::array<irred_bytes, 2> rows{text("train-cal0"), text("train-cal1")},
      events{text("synthetic-train-0"), text("synthetic-train-1")},
      future_rows{text("heldout-cal2"), text("heldout-cal3")},
      future_events{text("synthetic-future-2"), text("synthetic-future-3")},
      parameters{text("offset"), text("slope")},
      units{text("mag"), text("mag")},
      columns{text("mag/mag"), text("mag/mag")},
      cases{text("original-rational-control"),
            text("prior-mean-response-control")};
  const std::array<irred_bytes, 1> shared{text("offset")};
  irred_gaussian_predictive_batch batch() const {
    irred_gaussian_predictive_batch b{};
    b.struct_size = sizeof(b);
    b.abi_version = IRRED_ABI_VERSION;
    b.requested_outputs = 3;
    b.future_noise_independence_declared = 1;
    b.noise_conditional_on_parameters_declared = 1;
    auto &t = b.training;
    t.struct_size = sizeof(t);
    t.abi_version = IRRED_ABI_VERSION;
    t.source_semantics = 1;
    t.noise_independence_declared = 1;
    t.case_count = 2;
    t.noise_covariance = values(C);
    t.design = values(X);
    t.prior_mean = values(m);
    t.prior_covariance = values(S);
    t.conditioning_vectors = values(training);
    t.noise_row_ids = t.design_row_ids = t.conditioning_row_ids = ids(rows);
    t.noise_event_ids = t.conditioning_event_ids = ids(events);
    t.prior_parameter_ids = t.design_parameter_ids = ids(parameters);
    t.prior_parameter_units = t.design_parameter_units = ids(units);
    t.column_units = ids(columns);
    t.shared_nuisance_ids = ids(shared);
    t.case_ids = ids(cases);
    t.residual_unit = text("mag");
    t.noise_identity = text("synthetic training conditional noise");
    t.calibration_identity = text("offset and slope represented once in beta");
    t.noise_dependence_identity =
        text("full training covariance; no cross-case law");
    t.ordering_provenance = text("train-cal0,train-cal1 exact order");
    t.design_identity = text("synthetic fixed X");
    t.prior_identity = text("original proper prior");
    t.parameter_measure = text("d(offset-mag) d(slope-mag)");
    t.prior_dependence_identity = text("prior independent of both noises");
    b.future_noise_covariance = values(R);
    b.future_response = values(A);
    b.future_vectors = values(future);
    b.future_noise_row_ids = b.future_response_row_ids =
        b.future_vector_row_ids = ids(future_rows);
    b.future_noise_event_ids = b.future_vector_event_ids = ids(future_events);
    b.future_parameter_ids = ids(parameters);
    b.future_parameter_units = ids(units);
    b.future_column_units = ids(columns);
    b.future_unit = text("mag");
    b.future_noise_identity = text("synthetic future conditional noise");
    b.future_calibration_identity = t.calibration_identity;
    b.future_noise_dependence_identity =
        text("full future noise independent of prior and training noise");
    b.future_ordering_provenance =
        text("heldout-cal2,heldout-cal3 exact order");
    b.future_response_identity = text("synthetic fixed A");
    b.future_covariance_unit = text("mag^2");
    b.future_measure = text("product d(mag)");
    b.conditioning_identity = text("fixed two-case family");
    b.prediction_dependence_identity = text(
        "shared beta and full independent future noise; no cross-case law");
    return b;
  }
};
inline irred_gaussian_posterior_policy policy() {
  return {sizeof(irred_gaussian_posterior_policy),
          IRRED_ABI_VERSION,
          1,
          0,
          1000000,
          65536,
          1048576,
          268435456,
          100000000,
          1e-10};
}
struct Result {
  irred_gaussian_predictive_result *p = nullptr;
  Result() = default;
  Result(const Result &) = delete;
  Result &operator=(const Result &) = delete;
  ~Result() { irred_gaussian_predictive_result_destroy(p); }
};
} // namespace predictive_fixture
