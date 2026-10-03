#include "irred/effective_fluid.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void require(bool x, const char *why) { ++checks; if (!x) throw std::runtime_error(why); }
EffectiveFluidRequest fixture(double xc = .5) {
  return {{70, .5, .125, .0625, .0625, {}}, xc, .5, .5,
          "synthetic-dyadic-emitted-reference", "fixed-synthetic-half-scale"};
}
}
int main() {
  ThermalObservablePolicy p;
  auto r = fixture();
  auto owner = prepare_effective_fluid(r, p);
  require(owner.status() == S::ok, "synthetic emitted reference prepares");
  const auto transition = owner.fluid_state(.5L);
  require(transition.status == S::ok && transition.density == .5L &&
          transition.equation_of_state == -.25L && transition.log_scale_derivative == -1.125L,
          "exact transition density/equation/conservation witness");
  const auto early = owner.fluid_state(0);
  require(early.status == S::ok && early.density == 1 && early.equation_of_state == -1 &&
          early.log_scale_derivative == 0, "exact early density/conservation limit");
  for (long double a : {.125L, .25L, .5L, .75L, 1.L}) {
    const auto s = owner.fluid_state(a);
    require(s.status == S::ok &&
            std::abs(s.log_scale_derivative + 3 * (1 + s.equation_of_state) * s.density) <=
            s.derivative_estimate + 3 * s.density * s.equation_of_state_estimate,
            "continuity law including actual arithmetic estimates");
  }
  require(owner.omega_lambda() && *owner.omega_lambda() < .25L,
          "finite today fluid compensated from actual wide Lambda");
  auto zero = prepare_effective_fluid(fixture(0), p);
  require(zero.status() == S::ok, "zero amplitude owner prepares");
  for (long double a : {0.L, .125L, .5L, 1.L}) {
    const auto ref = zero.reference_background().scaled_expansion(a, p.thermal),
               got = zero.scaled_expansion(a, p.thermal);
    require(ref.status == got.status && ref.a4_e2 == got.a4_e2 &&
            ref.error_estimate == got.error_estimate && ref.callbacks == got.callbacks,
            "zero amplitude preserves exact reference coordinate/error/work");
    if (a == 0 || a == 1) {
      const auto q = owner.scaled_expansion(a, p.thermal);
      require(q.status == ref.status && q.a4_e2 == ref.a4_e2 &&
              q.error_estimate == ref.error_estimate && q.callbacks == ref.callbacks,
              "nonzero fluid origin/today exact scaled reference identity");
    }
  }
  const double factors[]{1, .25, .5, .25, 0, -1, std::numeric_limits<double>::infinity()};
  const auto b = owner.evaluate(factors, thermal_e | thermal_h | effective_fluid_ruler, p);
  require(b.status == S::ok && b.rows.size() == 7 && b.ruler && b.ruler->status == S::ok,
          "ordered coarse batch with one ruler");
  require(b.rows[0].e.value == 1 && b.rows[0].e.error_estimate == 0 &&
          b.rows[0].h_km_s_mpc.value == 70 && b.rows[0].h_km_s_mpc.error_estimate == 0,
          "today exact E/H normalization");
  require(b.rows[1].e.value == b.rows[3].e.value && b.rows[1].scale_factor == factors[1] &&
          b.rows[4].e.status == S::outside_domain && b.rows[5].e.status == S::outside_domain &&
          b.rows[6].e.status == S::nonfinite_input, "duplicates/order and every invalid row retained");
  require(b.callbacks == b.outer_callbacks + b.momentum_callbacks && b.outer_callbacks > 0 &&
          b.momentum_callbacks == 0 && b.callbacks <= p.maximum_total_callbacks,
          "combined work categories exact for massless owner");
  const auto baseline = zero.evaluate(factors, thermal_e | effective_fluid_ruler, p);
  require(baseline.ruler && baseline.ruler->value && *b.ruler->value < *baseline.ruler->value,
          "matched finite endpoint ruler decreases for positive fluid");
  // The dyadic endpoint is the EXACT same wide value in the existing z route.
  // This checks literal helper parity; a rounded general z->a is not identical.
  const ThermalPhysicalModel physical{70, .02, .1, 2.7255, 1.7e-5, {}};
  const auto mapping = map_thermal_physical_model(physical, p.thermal);
  require(mapping.status == S::ok && mapping.model, "literal parity source map");
  const auto previous = prepare_thermal_observables(
      {physical, 1, "exact-z-one-half-scale", "literal-public-parity-source"}, p);
  const auto matched = prepare_effective_fluid(
      {*mapping.model, 0, .5, .5, "literal-public-parity-source", "exact-z-one-half-scale"}, p);
  const auto old_ruler = previous.evaluate({}, thermal_ruler_mask, p);
  const auto exact_ruler = matched.evaluate({}, effective_fluid_ruler, p);
  require(previous.status() == S::ok && matched.status() == S::ok && old_ruler.ruler &&
          exact_ruler.ruler && old_ruler.ruler->status == S::ok && exact_ruler.ruler->status == S::ok &&
          old_ruler.ruler->value == exact_ruler.ruler->value &&
          old_ruler.ruler->error_estimate == exact_ruler.ruler->error_estimate &&
          old_ruler.callbacks == exact_ruler.callbacks && old_ruler.outer_callbacks == exact_ruler.outer_callbacks &&
          old_ruler.momentum_callbacks == exact_ruler.momentum_callbacks,
          "zero fluid literal previous ruler value/diagnostic/work parity");
  for (unsigned i = 0; i < 4; ++i)
    require(b.rows[i].e.value && baseline.rows[i].e.value &&
            *b.rows[i].e.value >= *baseline.rows[i].e.value, "matched E cannot decrease");
  const auto omitted = owner.evaluate({factors, 2}, thermal_e, p);
  require(!omitted.ruler && !omitted.rows[0].h_km_s_mpc.value && omitted.outer_callbacks == 0,
          "omitted outputs perform no ruler/dimensional work");
  const auto only = owner.evaluate(factors, effective_fluid_ruler, p),
             empty = owner.evaluate({}, effective_fluid_ruler, p);
  require(only.ruler->value == empty.ruler->value && only.callbacks == empty.callbacks &&
          !only.rows[0].e.value && only.rows[0].callbacks == 0,
          "one fixed ruler independent of ignored rows");
  for (unsigned field = 0; field < 7; ++field) {
    auto bad = r;
    if (field == 0) bad.density_at_transition = -1;
    if (field == 1) bad.transition_scale = 0;
    if (field == 2) bad.transition_scale = 1;
    if (field == 3) bad.ruler_endpoint = 0;
    if (field == 4) bad.ruler_endpoint = std::nextafter(1., 2.);
    if (field == 5) bad.source_origin.clear();
    if (field == 6) bad.endpoint_origin.clear();
    require(prepare_effective_fluid(bad, p).status() == S::outside_domain,
            "source/domain refusal has no fallback");
  }
  auto bad = r; bad.density_at_transition = 100;
  const auto closure = prepare_effective_fluid(bad, p);
  require(closure.status() == S::outside_domain && closure.request().density_at_transition == 100 &&
          closure.reference_background().status() == S::ok, "negative closure retains prepared source prefix");
  bad = r; bad.density_at_transition = std::numeric_limits<double>::quiet_NaN();
  require(prepare_effective_fluid(bad, p).status() == S::nonfinite_input, "nonfinite source refusal");
  auto q = p; q.maximum_native_bytes = 0;
  require(prepare_effective_fluid(r, q).status() == S::work_limit, "payload refusal precedes source copy");
  const auto bytes = effective_fluid_payload_bound(7, 0, r.source_origin.size() + r.endpoint_origin.size() + 2);
  require(bytes.has_value(), "finite whole simultaneous payload law");
  q = p; q.maximum_native_bytes = *bytes - 1;
  require(owner.evaluate(factors, thermal_e, q).rows.empty(), "evaluation payload admission before rows");
  q = p; q.maximum_points = 0;
  require(owner.evaluate(factors, thermal_e, q).status == S::work_limit, "point cap");
  q = p; q.maximum_total_callbacks = 0;
  const auto limited = owner.evaluate({factors, 1}, thermal_e | effective_fluid_ruler, q);
  require(limited.ruler->status == S::work_limit && limited.rows[0].e.value == 1 && limited.callbacks == 0,
          "failed ruler retains independent today scalar without extra work");
  q = p; q.relative_tolerance = 1e-30;
  require(owner.evaluate({factors + 1, 1}, thermal_e, q).rows[0].e.status == S::conditioning_budget_exceeded,
          "unattainable scalar budget is not relaxed");
  q = p; q.thermal.momentum_method = ThermalMomentumMethod::nested_clenshaw_curtis;
  require(owner.evaluate({}, effective_fluid_ruler, q).status == S::invalid_input,
          "retained method selector checked even empty batch");
  require(owner.evaluate({}, 0, p).status == S::invalid_input &&
          owner.evaluate({}, 8, p).status == S::invalid_input, "strict output masks");
  require(owner.scaled_expansion(2, p.thermal).status == S::outside_domain,
          "future scale is not silently continued");
  require(owner.scaled_expansion(1e-2000L, p.thermal).status != S::ok,
          "unresolved positive fluid powers are refused");
  static_assert(std::is_nothrow_move_constructible_v<EffectiveFluidBackground>);
  auto copy = owner; auto moved = std::move(owner);
  require(owner.status() == S::invalid_input && owner.request().source_origin.empty() &&
          owner.reference_background().status() == S::invalid_input && owner.evaluate(factors, thermal_e).rows.empty(),
          "move invalidates complete old owner");
  EffectiveFluidBackground assigned; assigned = std::move(moved);
  auto *alias = &assigned; assigned = std::move(*alias);
  require(assigned.evaluate({factors, 2}, thermal_e).rows[1].e.value ==
          copy.evaluate({factors, 2}, thermal_e).rows[1].e.value, "copy/move/self-move preserve retained source");
  const int rounding = std::fegetround(); std::fesetround(FE_UPWARD);
  const auto unsupported = assigned.evaluate(factors, thermal_e, p);
  std::fesetround(rounding);
  require(unsupported.status == S::invalid_input && unsupported.rows.empty(), "nonnearest profile refusal");
  std::cout << "effective_fluid owner checks=" << checks << " callbacks=" << b.callbacks << '\n';
}
