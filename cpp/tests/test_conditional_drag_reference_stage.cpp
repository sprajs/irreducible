#include "conditional_drag_reference_stage.hpp"
#include <cstdio>
#include <type_traits>
namespace ref = conditional_drag_reference_certificate;
namespace {
unsigned checks = 0, failures = 0;
bool need(bool answer, const char *label) {
  ++checks;
  if (!answer) { ++failures; std::fprintf(stderr, "FAILED %s\n", label); }
  return answer;
}
constexpr ref::Profile assumptions{true, true, true, false, false, false};
bool setup(ref::Owner &owner, ref::Arithmetic &arithmetic, ref::StageSupport &support) {
  if (!need(!support.prepared() && owner.served()[unsigned(ref::Counter::destinations)] == 0,
            "automatic brace construction leaves stage unavailable before admission")) return false;
  if (!need(owner.admit_arithmetic(), "stage conditional arithmetic environment") ||
      !need(arithmetic.prepare_directions(), "stage once-owned direction source")) return false;
  // Fixture source literals are fixed declared inputs, not a numerical oracle.
  const ref::StageVector center{0.5L, 0.5L, 1, 0.5L, 0.5L, 1};
  if (!need(ref::prepare_stage_support(arithmetic, support, center),
            "single existing-candidate support box acquired")) return false;
  return need(support.prepared() && support.weight(0).lower == 0.03125L &&
              support.weight(2).lower == 0.0625L && support.box(0).lower > 0 &&
              support.box(0).upper < 1,
              "fraction and temperature weights preserve the declared one-sixteenth law");
}
bool linear_inputs(ref::Arithmetic &a, ref::StageMatrix &M, ref::StageMatrix &J,
                   std::array<ref::Interval, 6> &r, bool zero_M, bool uncertain_J,
                   bool zero_residual) {
  ref::SourceNumbers numbers(a);
  for (unsigned i = 0; i != 6; ++i) {
    if (!(zero_residual || i != 0 ? numbers.integer(r[i], 0) : numbers.rational(r[i], 1, 1024)))
      return false;
    for (unsigned j = 0; j != 6; ++j) {
      if (!numbers.integer(M[6 * i + j], !zero_M && i == j ? 1 : 0)) return false;
      if (uncertain_J && i == j) {
        const ref::W lo = 0.5L, hi = 1.5L;
        if (!a.enclose(J[6 * i + j], lo, hi)) return false;
      } else if (!numbers.integer(J[6 * i + j], i == j ? 1 : 0)) return false;
    }
  }
  return true;
}
void incoming_once_and_tight_radius() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  ref::StageSupport support{};
  if (!setup(owner, arithmetic, support)) return;
  ref::StageMatrix M, J;
  std::array<ref::Interval, 6> residual;
  if (!need(linear_inputs(arithmetic, M, J, residual, false, false, false),
            "exact identity Jacobian and retained preconditioner fixture")) return;
  const std::array<ref::W, 3> incoming{0.001953125L, 0, 0}; // 2/1024.
  ref::StageTransport result;
  if (!need(ref::certify_stage_transport(arithmetic, result, support, M, J, residual, incoming) &&
            result.valid, "actual production contraction kernel admits fixed linear system")) return;
  // Independent exact algebra: J=M=I, residual0=1/1024, E repeats
  // the incoming coordinate. Maximum stage displacement is3/1024.
  // A doubled incoming contribution would be5/1024 and fail this dyadic limit.
  need(result.radius[0].upper >= 0.0029296875L && result.radius[0].upper < 0.00390625L &&
       result.radius[0].upper < support.weight(0).lower &&
       result.contraction < 0.125L && result.scaled_displacement < 0.125L,
       "incoming radius contributes exactly once and tight radius stays below support weight");
}
void uniform_J_and_genuine_zero() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  ref::StageSupport support{};
  if (!setup(owner, arithmetic, support)) return;
  ref::StageMatrix M, J;
  std::array<ref::Interval, 6> residual;
  if (!need(linear_inputs(arithmetic, M, J, residual, false, true, false),
            "whole uniform Jacobian support supplied")) return;
  const std::array<ref::W, 3> zero{0, 0, 0};
  ref::StageTransport result;
  if (!need(ref::certify_stage_transport(arithmetic, result, support, M, J, residual, zero) &&
            result.valid, "uniform derivative defect is admitted")) return;
  // J ranges from1/2 to3/2, so its exact contraction modulus is1/2;
  // the admitted bound must cover the worst inverse response2/1024.
  need(result.contraction >= 0.5L && result.contraction < 0.625L &&
       result.radius[0].upper >= 0.001953125L && result.radius[0].upper < 0.0029296875L,
       "uniform Jacobian uncertainty changes the actual response bound");

  ref::Owner zero_owner(assumptions);
  ref::Arithmetic za(zero_owner);
  ref::StageSupport zs{};
  if (!setup(zero_owner, za, zs)) return;
  ref::StageMatrix zM, zJ;
  std::array<ref::Interval, 6> zr;
  if (!need(linear_inputs(za, zM, zJ, zr, false, false, true), "exact zero residual fixture")) return;
  ref::StageTransport genuine_zero;
  if (!need(ref::certify_stage_transport(za, genuine_zero, zs, zM, zJ, zr, zero) &&
            genuine_zero.valid, "zero residual and zero incoming retain exact-zero witness")) return;
  need(genuine_zero.scaled_displacement == 0 && genuine_zero.radius[0].lower == 0 &&
       genuine_zero.radius[0].upper == 0 && genuine_zero.radius[5].upper == 0,
       "positive support weights are not incorrectly reported as a zero-source error");
}
void refusals_and_one_box() {
  ref::Owner owner(assumptions);
  ref::Arithmetic arithmetic(owner);
  ref::StageSupport support{};
  if (!setup(owner, arithmetic, support)) return;
  ref::StageMatrix M, J;
  std::array<ref::Interval, 6> residual;
  if (!need(linear_inputs(arithmetic, M, J, residual, true, false, false),
            "zero preconditioner adversarial fixture")) return;
  const std::array<ref::W, 3> incoming{0, 0, 0};
  ref::StageTransport refused;
  need(!ref::certify_stage_transport(arithmetic, refused, support, M, J, residual, incoming) &&
       !refused.valid && owner.failure().cause == ref::Cause::unresolved_positive,
       "q at or above one refuses instead of enlarging support or Newton budgets");
  const auto first = owner.failure();
  need(!ref::certify_stage_transport(arithmetic, refused, support, M, J, residual, incoming) &&
       owner.failure().cause == first.cause && owner.failure().served_prefix == first.served_prefix,
       "later attempted inclusion retains the original causal failure");

  ref::Owner one_owner(assumptions);
  ref::Arithmetic one_arithmetic(one_owner);
  ref::StageSupport one_support{};
  if (!setup(one_owner, one_arithmetic, one_support)) return;
  const ref::StageVector changed_center{0.25L, 0.25L, 1, 0.25L, 0.25L, 1};
  need(!ref::prepare_stage_support(one_arithmetic, one_support, changed_center) &&
       one_owner.failure().cause == ref::Cause::invalid_input && one_support.prepared() &&
       one_support.center(0).lower == 0.5L,
       "prepared support cannot be replaced by a second candidate box");

  ref::Owner alias_owner(assumptions);
  ref::Arithmetic aa(alias_owner);
  ref::StageSupport as{};
  if (!setup(alias_owner, aa, as)) return;
  ref::StageMatrix aM, aJ;
  std::array<ref::Interval, 6> ar;
  if (!need(linear_inputs(aa, aM, aJ, ar, false, false, true), "output-alias fixture")) return;
  ref::StageTransport acquired;
  ref::SourceNumbers numbers(aa);
  if (!need(numbers.integer(acquired.radius[0], 1), "acquired overlapping residual field")) return;
  need(!ref::certify_stage_transport(aa, acquired, as, aM, aJ, acquired.radius, incoming) &&
       alias_owner.failure().cause == ref::Cause::invalid_ownership &&
       acquired.radius[0].valid && acquired.radius[0].lower == 1,
       "overlapping output refuses before altering the acquired input");
}
}
static_assert(!std::is_copy_constructible_v<ref::StageSupport>);
static_assert(!std::is_move_constructible_v<ref::StageTransport>);
int main() {
  incoming_once_and_tight_radius(); uniform_J_and_genuine_zero(); refusals_and_one_box();
  std::printf("conditional_drag_reference_stage checks=%u failures=%u "
              "true_J_caller_qualified=false complete_arithmetic_qualified=false whole_payload_qualified=false\n",
              checks, failures);
  return failures ? 1 : 0;
}
