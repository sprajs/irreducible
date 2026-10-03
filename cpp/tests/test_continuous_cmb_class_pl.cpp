#include "fixtures/continuous_cmb_class_pl_604.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

namespace {
using namespace irred::projection;
namespace f = irred::test_fixtures::continuous_cmb_class_pl;
using S = irred::numerics::Status;
void need(bool condition, const char *message) {
  if (!condition) { std::cerr << "FAIL " << message << '\n'; std::exit(1); }
}
// Local nearest-binary64 acceptance assembly only. One outward adjacent value
// covers each basic rounded operation; the original MPFR reference interval is
// represented by exact portable bit endpoints. This is no production theorem.
double up(double value) {
  need(std::isfinite(value) && value >= 0, "finite positive error assembly");
  return std::nextafter(value, std::numeric_limits<double>::infinity());
}
double down(double value) {
  need(std::isfinite(value) && value > 0, "positive unchanged allocation");
  return std::nextafter(value, 0.);
}
void quality(double value, const f::ReferenceInterval &reference,
             double time, double arithmetic, const ContinuousCmbPolicy &policy) {
  need(std::isfinite(value) && std::isfinite(time) && std::isfinite(arithmetic) &&
       time >= 0 && arithmetic >= 0 && reference.lower() < reference.upper(),
       "complete finite reference and native diagnostics");
  const double distance = std::max(up(std::abs(value-reference.lower())),
                                   up(std::abs(value-reference.upper())));
  const double error = up(up(distance+time)+arithmetic);
  const double relative = policy.relative_tolerance*std::abs(value);
  need(relative >= std::numeric_limits<double>::min(), "normal allocation product");
  const double allocation = down(policy.absolute_tolerance+down(relative));
  need(error <= allocation, "same PL geometry combined independent reference budget");
}
void actual_original_table() {
  const ContinuousCmbPolicy policy;
  need(policy.absolute_tolerance == 1e-9 && policy.relative_tolerance == 1e-5 &&
       policy.maximum_kernel_evaluations == 2000000,
       "unchanged default native numerical policy");
  auto source = f::make_source();
  const auto *original_eta = source.eta_mpc.data();
  auto owner = prepare_continuous_cmb_projection(std::move(source));
  need(owner.status() == S::ok && owner.source()->eta_mpc.data() == original_eta,
       "one original source acquisition moves its buffers");
  const auto *identity = owner.source();
  const std::array<unsigned,2> ell{2,19};
  const auto result = project_continuous_cmb(owner,ell,
      continuous_temperature|continuous_e_mode,policy);
  need(result.status == S::ok && result.rows.size() == 2 &&
       result.source_owner.get() == identity &&
       result.kernel_evaluations > 0 && result.kernel_evaluations <= 2000000,
       "original full-table evaluation and source ownership");
  for (std::size_t index=0; index<result.rows.size(); ++index) {
    const auto &row=result.rows[index];
    need(row.status == S::ok && row.ell == ell[index] &&
         row.multipole_index == index && row.k_index == 0 &&
         row.k_mpc_inverse == f::k_mpc_inverse &&
         row.completed_source_cells == 603 && row.temperature && row.e_mode &&
         !row.radial_error_estimate && !row.source_grid_error_estimate,
         "exact source order, finite support, mask and absent physical errors");
    need(f::reference[2*index].ell==row.ell && !f::reference[2*index].e_mode &&
         f::reference[2*index+1].ell==row.ell && f::reference[2*index+1].e_mode,
         "explicit reference ell/output order");
    quality(*row.temperature,f::reference[2*index],row.temperature_quadrature_estimate,
            row.temperature_arithmetic_estimate,policy);
    quality(*row.e_mode,f::reference[2*index+1],row.e_quadrature_estimate,
            row.e_arithmetic_estimate,policy);
  }
  // Required source failures retain the actual owner and attempted nodes.
  auto capped=policy;capped.maximum_kernel_evaluations=2;
  const auto refused=project_continuous_cmb(owner,ell,continuous_temperature,capped);
  need(refused.status == S::work_limit && refused.kernel_evaluations == 2 &&
       refused.source_owner.get() == identity && !refused.rows.empty(),
       "same full-table failed work is retained");
  for (const auto &row:refused.rows)
    need(!row.temperature && !row.e_mode, "refused dependent values withheld");
  owner=ContinuousCmbProjection{};
  need(result.source_owner.get() == identity &&
       result.source_owner->eta_mpc.size() == f::source_rows.size(),
       "completed result retains original source after owner release");
}
ContinuousCmbSource control_source(bool endpoint, double amplitude) {
  ContinuousCmbSource s;
  s.k_mpc_inverse={endpoint?1.:.5};
  s.eta_mpc=endpoint?std::vector<double>{0,1e-12}:std::vector<double>{1,3};
  s.observer_eta_mpc=endpoint?1e-12:12;
  s.t0={0,0};s.t1=endpoint?std::vector<double>{0,0}:std::vector<double>{amplitude,amplitude};
  s.t2=endpoint?std::vector<double>{1e12,1e12}:std::vector<double>{0,0};
  s.polarization=s.t2;
  s.producer_id=endpoint?"synthetic-tiny-observer-endpoint":"synthetic-signed-constant-Doppler";
  s.signed_mode_id=endpoint?"positive-regular-limit":"declared-signed-amplitude";
  s.normalization_id=endpoint?"explicit-1e12-amplitude":"dimensionless";
  return s;
}
void signed_finite_endpoints() {
  ContinuousCmbPolicy policy;policy.maximum_kernel_evaluations=1024;
  const std::array<unsigned,1> ell{2};
  for (unsigned index=0;index<2;++index) {
    const double amplitude=index?-3.:1.;
    auto owner=prepare_continuous_cmb_projection(control_source(false,amplitude));
    need(owner.status()==S::ok,"signed finite endpoint acquisition");
    const auto result=project_continuous_cmb(owner,ell,continuous_temperature,policy);
    need(result.status==S::ok && result.rows.size()==1 &&
         result.rows[0].temperature && !result.rows[0].e_mode &&
         result.rows[0].completed_source_cells==1 && result.kernel_evaluations<=1024,
         "signed Doppler finite endpoint result");
    const auto &row=result.rows[0];
    quality(*row.temperature,f::synthetic_reference[index],
            row.temperature_quadrature_estimate,row.temperature_arithmetic_estimate,policy);
    need(index?*row.temperature>0:*row.temperature<0,
         "actual outward-phase Doppler sign");
  }
  auto owner=prepare_continuous_cmb_projection(control_source(true,0));
  need(owner.status()==S::ok,"observer endpoint acquisition");
  const auto result=project_continuous_cmb(owner,ell,continuous_temperature|continuous_e_mode,policy);
  need(result.status==S::ok && result.rows.size()==1 && result.rows[0].temperature &&
       result.rows[0].e_mode && result.rows[0].completed_source_cells==1 &&
       result.kernel_evaluations<=1024,"regular observer endpoint result");
  const auto &row=result.rows[0];
  quality(*row.temperature,f::synthetic_reference[2],row.temperature_quadrature_estimate,
          row.temperature_arithmetic_estimate,policy);
  quality(*row.e_mode,f::synthetic_reference[3],row.e_quadrature_estimate,
          row.e_arithmetic_estimate,policy);
  need(*row.temperature>0 && *row.e_mode>0,
       "positive ell2 quadrupole/E regular one-fifth boundary limit");
}
} // namespace
int main() {
  need(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53 &&
       std::fegetround()==FE_TONEAREST,"portable binary64 reference acceptance profile");
  actual_original_table();signed_finite_endpoints();
  std::cout << "continuous_cmb_class_pl_contract original604/MPFRbounds/sign/endpoints/ownership passed\n";
}
