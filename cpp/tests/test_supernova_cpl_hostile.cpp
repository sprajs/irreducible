// Independent author from SN v2 implementation. A-variable fixed-panel Simpson
// reference differs from production adaptive-z and owner's redshift GL8 route;
// conservation equations/libm shared. Synthetic diagonal weighted profile is
// derived analytically, never an absolute Gaussian density.
#include "irred/supernova.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <type_traits>
using namespace irred;
namespace {
int checks = 0;
long double max_error = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
void near(long double a, long double b) {
  max_error = std::max(max_error, std::abs(a - b));
  check(std::abs(a - b) < 1e-10, "synthetic assembled score/offset budget");
}
long double integral(double z, double om, double w0, double wa, int panels) {
  const auto lower = 1 / (1 + (long double)z), h = (1 - lower) / panels;
  long double sum = 0;
  for (int i = 0; i <= panels; ++i) {
    auto a = lower + i * h;
    auto E2 = om / (a * a * a) +
              (1 - om) * std::pow(a, -3 * (1 + (long double)w0 + wa)) *
                  std::exp(-3 * (long double)wa * (1 - a));
    sum += (i == 0 || i == panels ? 1
            : i % 2               ? 4
                                  : 2) /
           (a * a * std::sqrt(E2));
  }
  return sum * h / 3;
}
} // namespace
int main() {
  try {
    static_assert(!std::is_default_constructible_v<supernova::ModelPointV2>);
    observations::Input input;
    input.profile = observations::Profile::gaussian_fixture_v1;
    input.role = observations::Role::synthetic_control;
    input.unit = observations::Unit::magnitude;
    input.calibration = observations::Calibration::not_applicable;
    input.uncertainty = observations::Uncertainty::covariance;
    input.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
    input.component = observations::Component::total;
    input.table_sha256 = std::string(64, 'c');
    input.uncertainty_sha256 = std::string(64, 'd');
    input.ordering_provenance = "independent synthetic diagonal row order";
    input.measurement_ids = {"first", "second", "third", "excluded"};
    input.event_ids = input.measurement_ids;
    input.uncertainty_axis_ids = input.measurement_ids;
    input.values = {12, 15, 18, 99};
    input.uncertainty_matrix = {4, 0, 0,  8, 0, 9, 0, 7,
                                0, 0, 16, 6, 1, 2, 3, -1};
    input.missing = {0, 0, 0, 0};
    input.quality = {0, 0, 0, 0};
    input.source_selection = {1, 1, 1, 1};
    auto source = observations::prepare(input, {4, 16, 4096});
    std::array<cosmology::Query, 4> queries{
        {{.1, .11, cosmology::Convention::released_zhd_zhel},
         {.7, .71, cosmology::Convention::released_zhd_zhel},
         {2, 2.2, cosmology::Convention::released_zhd_zhel},
         {.01, .01, cosmology::Convention::released_zhd_zhel}}};
    supernova::Policy policy;
    policy.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
    policy.maximum_forward_sensitivity = 1e-10;
    policy.maximum_matrix_elements = 16;
    policy.background.integration = {1e-14, 1e-13, 100000, 30};
    auto consumer = supernova::prepare_synthetic(source, input.measurement_ids,
                                                 queries, policy);
    check(consumer.status() == supernova::Status::ok, "synthetic source valid");
    check(consumer.selected_source_indices().size() == 3,
          "strict threshold excludes fourth");
    std::array<supernova::ModelPointV2, 4> models{
        {{cosmology::Model::flat_cpl_late_v1, 0, 0, -2. / 3, 0},
         {cosmology::Model::flat_cpl_late_v1, .3, 0, -.9, .4},
         {cosmology::Model::flat_cpl_late_v1, .3, 0, -1.1, -.4},
         {cosmology::Model::flat_cpl_late_v1, 1, 0, -2, 2}}};
    auto batch = consumer.evaluate_batch_v2(models, policy);
    check(batch.status == supernova::Status::ok && batch.slots.size() == 4,
          "CPL batch owns slots");
    for (std::size_t m = 0; m < models.size(); ++m) {
      const auto &slot = batch.slots[m];
      check(slot.calculation.status == supernova::Status::ok, "CPL finite");
      check(slot.source.w0 == models[m].w0 && slot.source.wa == models[m].wa,
            "explicit source survives");
      long double residual[3], weighted = 0, gram = 0;
      const long double variances[]{4, 9, 16};
      for (std::size_t i = 0; i < 3; ++i) {
        auto fine = integral(queries[i].z_expansion, models[m].omega_m,
                             models[m].w0, models[m].wa, 4096);
        auto coarse = integral(queries[i].z_expansion, models[m].omega_m,
                               models[m].w0, models[m].wa, 2048);
        check(std::abs(fine - coarse) < 2e-13,
              "independent reference refinement");
        const auto prediction =
            5 * std::log10((1 + queries[i].z_observer) * fine);
        near(slot.calculation.shape_magnitudes[i], prediction);
        residual[i] = input.values[i] - prediction;
        weighted += residual[i] / variances[i];
        gram += 1 / variances[i];
      }
      const auto coefficient = weighted / gram;
      long double q = 0;
      for (std::size_t i = 0; i < 3; ++i) {
        q += (residual[i] - coefficient) * (residual[i] - coefficient) /
             variances[i];
        near(slot.calculation.profiled_residuals[i], residual[i] - coefficient);
      }
      near(slot.calculation.offset_coefficient, coefficient);
      near(slot.calculation.quadratic, q);
      near(slot.calculation.relative_profile_score, -q / 2);
      check(slot.calculation.profile_status ==
                statistics::DensityStatus::finite,
            "profile finite, no densityclaim");
    }
    std::array<supernova::ModelPoint, 1> old{
        {{cosmology::Model::constant_q_flat_v1, 0, -.5}}};
    auto constant = consumer.evaluate_batch(old, policy);
    near(batch.slots[0].calculation.relative_profile_score,
         constant.slots[0].relative_profile_score);
    for (double om : {0., .3, 1.}) {
      std::array<supernova::ModelPointV2, 1> lambda{
          {{cosmology::Model::flat_cpl_late_v1, om, 0, -1, 0}}};
      std::array<supernova::ModelPoint, 1> lcdm{
          {{cosmology::Model::flat_lcdm_late_v1, om, 0}}};
      auto a = consumer.evaluate_batch_v2(lambda, policy).slots[0].calculation;
      auto b = consumer.evaluate_batch(lcdm, policy).slots[0];
      check(a.shape_magnitudes == b.shape_magnitudes, "Lambda bitwise shapes");
      check(a.profiled_residuals == b.profiled_residuals,
            "Lambda canonical residual parity");
      near(a.quadratic, b.quadratic);
    }
    for (auto bad : {supernova::ModelPointV2{
                         cosmology::Model::flat_lcdm_late_v1, .3, 0, -.8, 0},
                     {cosmology::Model::flat_cpl_late_v1, .3, 1, -1, 0},
                     {cosmology::Model::constant_q_flat_v1, .1, -.5, -1, 0},
                     {cosmology::Model::flat_cpl_late_v1, .3, 0, -1, 3},
                     {static_cast<cosmology::Model>(88), .3, 0, -1, 0}}) {
      auto result =
          consumer.evaluate_batch_v2(std::span(&bad, 1), policy).slots[0];
      check(result.calculation.status != supernova::Status::ok,
            "inactive/domain/enum rejects");
      check(result.source.model == bad.model && result.source.w0 == bad.w0 &&
                result.source.wa == bad.wa,
            "failed attempt retained");
      check(result.calculation.shape_magnitudes.empty() &&
                result.calculation.profiled_residuals.empty(),
            "failed finite vectors omitted");
    }
    std::array<supernova::ModelPoint, 1> alias{
        {{cosmology::Model::flat_cpl_late_v1, .3, 0}}};
    check(consumer.evaluate_batch(alias, policy).slots[0].status !=
              supernova::Status::ok,
          "v1 CPL alias prohibited");
    auto cap = policy;
    cap.maximum_models = 0;
    check(consumer.evaluate_batch_v2(models, cap).slots.empty(), "model cap");
    auto work = policy;
    work.background.maximum_total_evaluations =
        batch.slots[0].calculation.background_evaluations;
    auto limited = consumer.evaluate_batch_v2(models, work);
    check(limited.slots[0].calculation.status == supernova::Status::ok,
          "first globalbudget slot finite");
    check(limited.slots[1].calculation.status ==
                  supernova::Status::work_limit &&
              limited.slots[1].calculation.shape_magnitudes.empty(),
          "global callback budget not reset");
    auto tight = policy;
    tight.maximum_forward_sensitivity = 1e-30;
    auto fail = consumer.evaluate_batch_v2(models, tight).slots[0].calculation;
    check(fail.numerical_status ==
                  numerics::Status::conditioning_budget_exceeded &&
              fail.profiled_residuals.empty(),
          "precise cache failcause no payload");
    check(consumer.observations().source().uncertainty_matrix ==
              input.uncertainty_matrix,
          "full raw source unchanged");
    std::printf("SN CPL peer %d PASS max_absolute_error %.17Lg\n", checks,
                max_error);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL %s after%d\n", e.what(), checks);
    return 1;
  }
}
