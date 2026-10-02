// Exact EdS and ownership/domain controls. Independent non-EdS scientific
// references belong to the separately authored peer test.
#include "irred/growth_amplitude.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
namespace {
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
FixedSigma8Source amplitude(double s = 1, double ref = 1) {
  return {
      s,
      ref,
      Sigma8Convention::linear_pressureless_total_matter_top_hat_8_over_h_mpc,
      AmplitudeTreatment::fixed_supplied,
      "synthetic-amplitude",
      "exact supplied control"};
}
} // namespace
int main() {
  try {
    const double a[]{1e-8, .25, .5, 1, 1};
    auto eds = prepare_gr_growth(prepare(LCDM(1), FlatFLRW{}));
    auto result = evaluate_growth_amplitude(eds, amplitude(), a, 3);
    need(result.status == irred::numerics::Status::ok &&
             result.rows.size() == 5,
         "EdS coarse batch admitted");
    for (size_t j = 0; j < 5; ++j) {
      need(result.rows[j].sigma8.value == a[j] &&
               result.rows[j].f_sigma8.value == a[j],
           "EdS exact normalized amplitude");
      need(result.rows[j].sigma8.status == irred::numerics::Status::ok &&
               result.rows[j].f_sigma8.status == irred::numerics::Status::ok,
           "EdS numerical acceptance");
    }
    need(result.callbacks == 0, "EdS exact branch has no callbacks");
    need(result.rows[3].sigma8.absolute_error_estimate == 0 && result.reference,
         "same-reference sigma8 cancellation");
    auto growth = prepare_gr_growth(prepare(LCDM(.315), FlatFLRW{}));
    auto ordinary = evaluate_growth_amplitude(growth, amplitude(.811), a, 3);
    need(ordinary.reference && ordinary.reference->d.value &&
             ordinary.reference->f.value,
         "retained reference D and f diagnostics");
    size_t callbacks = ordinary.reference->callbacks;
    for (const auto &r : ordinary.rows)
      callbacks += r.callbacks;
    need(callbacks == ordinary.callbacks && ordinary.rows[3].callbacks == 0 &&
             ordinary.rows[4].callbacks == 0,
         "reference evaluated once and joint work charged");
    auto doubled = evaluate_growth_amplitude(growth, amplitude(1.622), a, 2);
    for (size_t j = 0; j < 5; ++j) {
      need(ordinary.rows[j].f_sigma8.value && doubled.rows[j].f_sigma8.value &&
               *doubled.rows[j].f_sigma8.value ==
                   2 * *ordinary.rows[j].f_sigma8.value,
           "power-of-two amplitude scaling");
      need(!doubled.rows[j].sigma8.value &&
               doubled.rows[j].sigma8.availability ==
                   Availability::not_requested,
           "omitted sigma8 stays absent");
    }
    // H0 is a dimensional background projection input. At fixed fractions and
    // supplied 8/h-Mpc amplitude identity it does not enter this growth law.
    auto background = prepare(LCDM(.315), FlatFLRW{});
    const Request h50[]{
        Request(1, unsigned(Observable::expansion), {}, PhysicalScale(50))};
    const Request h100[]{
        Request(1, unsigned(Observable::expansion), {}, PhysicalScale(100))};
    const EvaluationPolicy ep{{}, 1, 1024, 1024, 1024 * 1024};
    auto lower_h = background.evaluate(h50, ep),
         upper_h = background.evaluate(h100, ep);
    need(lower_h.slots[0].expansion.value && upper_h.slots[0].expansion.value &&
             lower_h.slots[0].expansion.value->h_km_s_mpc.value &&
             upper_h.slots[0].expansion.value->h_km_s_mpc.value &&
             *upper_h.slots[0].expansion.value->h_km_s_mpc.value ==
                 2 * *lower_h.slots[0].expansion.value->h_km_s_mpc.value,
         "dimensional H changes with supplied H0");
    auto fixed_fraction = evaluate_growth_amplitude(
        prepare_gr_growth(background), amplitude(.811), a, 2);
    for (size_t j = 0; j < 5; ++j)
      need(fixed_fraction.rows[j].f_sigma8.value ==
               ordinary.rows[j].f_sigma8.value,
           "conditional fixed-fraction supplied-amplitude H0 invariance");
    GrowthAmplitudePolicy p;
    p.growth.maximum_total_callbacks = 0;
    auto zero = evaluate_growth_amplitude(growth, amplitude(0), a, 3, p);
    for (const auto &r : zero.rows)
      need(r.sigma8.value == 0 && r.f_sigma8.value == 0 &&
               r.f_sigma8.absolute_error_estimate == 0,
           "declared exact zero amplitude");
    need(zero.callbacks == 0 && !zero.reference,
         "zero requires no reference quadrature");
    const double bad[]{.5,   -1,   std::numeric_limits<double>::quiet_NaN(),
                       1e-9, 1.01, 1};
    auto mixed_zero =
        evaluate_growth_amplitude(growth, amplitude(0), bad, 2, p);
    need(mixed_zero.rows[0].f_sigma8.value == 0 &&
             !mixed_zero.rows[1].f_sigma8.value &&
             !mixed_zero.rows[2].f_sigma8.value &&
             !mixed_zero.rows[3].f_sigma8.value &&
             !mixed_zero.rows[4].f_sigma8.value &&
             mixed_zero.rows[5].f_sigma8.value == 0,
         "zero still validates every query");
    auto limited = evaluate_growth_amplitude(growth, amplitude(), a, 3, p);
    need(limited.callbacks == 0 && !limited.rows[0].f_sigma8.value &&
             limited.rows[3].sigma8.value == 1 &&
             !limited.rows[3].f_sigma8.value,
         "quota failure preserves independent exact endpoint sigma8");
    p = {};
    p.growth.maximum_total_callbacks = 17;
    auto capped = evaluate_growth_amplitude(growth, amplitude(), a, 3, p);
    need(capped.callbacks <= 17 && !capped.rows[0].f_sigma8.value,
         "failed reference callbacks obey original total cap");
    p = {};
    const auto bound = growth_amplitude_payload_bound(5);
    need(bound.has_value(), "representable payload");
    p.maximum_native_bytes = *bound - 1;
    need(evaluate_growth_amplitude(growth, amplitude(), a, 3, p).rows.empty(),
         "payload refused before returned allocation");
    p.maximum_native_bytes = *bound;
    need(evaluate_growth_amplitude(growth, amplitude(), a, 3, p).rows.size() ==
             5,
         "exact payload boundary admits");
    p = {};
    p.growth.maximum_points = 4;
    need(evaluate_growth_amplitude(growth, amplitude(0), a, 3, p).rows.empty(),
         "zero still obeys point quota");
    p = {};
    p.maximum_string_bytes = 0;
    need(evaluate_growth_amplitude(growth, amplitude(), a, 3, p).rows.empty(),
         "source string quota");
    p = {};
    p.absolute_tolerance = 0;
    p.relative_tolerance = 1e-30;
    auto tight = evaluate_growth_amplitude(growth, amplitude(), a, 3, p);
    need(!tight.rows[0].f_sigma8.value && tight.rows[3].sigma8.value == 1,
         "composed budget failure does not erase exact endpoint");
    for (auto s :
         {amplitude(-1), amplitude(std::numeric_limits<double>::denorm_min()),
          amplitude(1, 0), amplitude(1, 1.01)})
      need(evaluate_growth_amplitude(growth, s, a, 3).rows.empty(),
           "source domain refusal");
    auto source = amplitude();
    source.amplitude_identity.clear();
    need(evaluate_growth_amplitude(growth, source, a, 3).rows.empty(),
         "source identity refusal");
    source = amplitude();
    source.convention = Sigma8Convention::unknown;
    need(evaluate_growth_amplitude(growth, source, a, 3).rows.empty(),
         "amplitude convention refusal");
    source = amplitude();
    source.treatment = AmplitudeTreatment::unknown;
    need(evaluate_growth_amplitude(growth, source, a, 3).rows.empty(),
         "uncertainty declaration refusal");
    for (auto spec :
         {ExpansionSpec(CPL(.3, -1, 0)), ExpansionSpec(LCDM(1e-7))}) {
      auto incompatible = prepare_gr_growth(prepare(spec, FlatFLRW{}));
      need(evaluate_growth_amplitude(incompatible, amplitude(0), a, 3)
               .rows.empty(),
           "zero cannot repair unsupported physical closure");
    }
    auto tiny = evaluate_growth_amplitude(
        eds, amplitude(std::numeric_limits<double>::min()), std::span(a, 1), 2);
    need(!tiny.rows[0].f_sigma8.value, "positive underflow never becomes zero");
    auto huge = evaluate_growth_amplitude(
        eds, amplitude(std::numeric_limits<double>::max(), 1e-8),
        std::span(a + 3, 1), 2);
    need(!huge.rows[0].f_sigma8.value, "positive overflow refused");
    need(!growth_amplitude_payload_bound(SIZE_MAX),
         "payload overflow explicit");
    need(evaluate_growth_amplitude(growth, amplitude(), a, 0).rows.empty() &&
             evaluate_growth_amplitude(growth, amplitude(), a, 4).rows.empty(),
         "unknown output masks");
    const auto rounding = std::fegetround();
    need(std::fesetround(FE_DOWNWARD) == 0, "set unsupported rounding");
    auto wrong = evaluate_growth_amplitude(growth, amplitude(0), a, 3);
    need(std::fesetround(rounding) == 0, "restore rounding");
    need(wrong.rows.empty(), "zero validates arithmetic environment");
    std::cout << "PASS " << checks << " growth amplitude owner controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
