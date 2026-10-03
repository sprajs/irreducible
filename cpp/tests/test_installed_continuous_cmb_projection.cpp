// Compile separately against a fresh installed public header/archive prefix.
#include "irred/continuous_cmb_projection.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <type_traits>
#include <utility>

namespace {
using namespace irred::projection;
using S = irred::numerics::Status;
void need(bool ok, const char *why) {
  if (!ok) { std::cerr << "FAIL installed " << why << '\n'; std::exit(1); }
}
long double j0(long double x) { return x == 0 ? 1 : std::sin(x) / x; }
ContinuousCmbSource source() {
  ContinuousCmbSource s;
  s.k_mpc_inverse = {.5, 1}; s.eta_mpc = {1, 3}; s.observer_eta_mpc = 12;
  s.t0 = {1, 1, 1, 1}; s.t1 = s.t2 = s.polarization = {0, 0, 0, 0};
  s.producer_id = "installed-synthetic-constant-monopole";
  s.signed_mode_id = "positive-unit-synthetic-mode";
  s.normalization_id = "dimensionless-unit-mode";
  return s;
}
} // namespace
int main() {
  static_assert(!std::is_copy_constructible_v<ContinuousCmbProjection>);
  auto input = source();
  const auto *original_buffer = input.t0.data();
  auto prepared = prepare_continuous_cmb_projection(std::move(input));
  need(prepared.status() == S::ok && prepared.source()->t0.data() == original_buffer,
       "one acquired original source buffer");
  const auto *identity = prepared.source();
  auto moved = std::move(prepared);
  need(prepared.status() != S::ok && !prepared.source() && moved.source() == identity,
       "move-only retained source identity");
  const std::array<unsigned, 1> ell{1};
  ContinuousCmbPolicy policy;
  policy.absolute_tolerance = 2e-12; policy.relative_tolerance = 2e-10;
  const auto result = project_continuous_cmb(moved, ell, continuous_temperature, policy);
  need(result.status == S::ok && result.rows.size() == 2 &&
       result.source_owner.get() == identity, "public archive supplies finite transfer");
  for (std::size_t i = 0; i < result.rows.size(); ++i) {
    const auto &r = result.rows[i]; const long double k = r.k_mpc_inverse;
    // d[j0(k*(eta0-eta))]/deta=k*j1, so every finite boundary is retained.
    const long double exact = (j0(9*k) - j0(11*k)) / k;
    need(r.status == S::ok && r.k_index == i && r.multipole_index == 0 &&
         r.ell == 1 && r.temperature && !r.e_mode &&
         std::abs(static_cast<long double>(*r.temperature) - exact) < 2e-10L &&
         !r.radial_error_estimate && !r.source_grid_error_estimate,
         "independent trig endpoint/mask/order/absent external errors");
  }
  auto limited = policy; limited.maximum_kernel_evaluations = 0;
  const auto refusal = project_continuous_cmb(moved, ell, continuous_temperature, limited);
  need(refusal.status == S::work_limit && refusal.source_owner.get() == identity &&
       refusal.rows.size() == 2 && !refusal.rows[0].temperature &&
       !refusal.rows[1].temperature, "retained work refusal without dependent output");
  const std::array<unsigned, 1> outside{65};
  const auto domain = project_continuous_cmb(moved, outside, continuous_temperature, policy);
  need(domain.status == S::outside_domain && domain.source_owner.get() == identity,
       "public ell domain refusal");
  moved = {};
  need(result.source_owner.get() == identity && result.source_owner->t0[0] == 1 &&
       refusal.source_owner.get() == identity, "success and refusal outlive prepared owner");
  std::cout << "PASS installed continuous finite-support scalar transfer consumer\n";
}
