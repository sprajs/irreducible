#include "irred/hydrogen_helium_history.hpp"
#include "../src/thermal_ruler.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
using namespace irred::cosmology;
using S = irred::numerics::Status;
using W = long double;
namespace {
unsigned checks=0;
void need(bool b,const char *message){++checks;if(!b)throw std::runtime_error(message);}
void foundation_controls() {
  const auto bg = prepare_thermal_background({70, 3./16, 5./16, .25, .25, {}});
  const auto r = detail::thermal_baryon_loading(bg);
  need(r.status == S::ok && r.ratio_today == 1, "analytic loading R0=1");
  need(r.arithmetic_estimate == 32 * std::numeric_limits<W>::epsilon(), "frozen exact ratio arithmetic allowance");
  for (W u : {1.L, 301.L, 1101.L, 2701.L}) {
    const auto q = detail::thermal_baryon_loading_at_inverse_scale(r, u);
    const W original = 3 * W(bg.source().omega_b) / (4 * W(bg.source().omega_gamma) * u);
    need(q.status == S::ok && q.value == original, "pure-H original grouped center preserved");
    need(q.arithmetic_estimate <= 128 * std::numeric_limits<double>::epsilon() * q.value,
         "helper operation estimate inside inherited pure-H floor");
  }
  const auto small_bg = prepare_thermal_background({70, .25, 0, 1e-300, .5, {}});
  const auto small = detail::thermal_baryon_loading(small_bg);
  need(small.status == S::ok && small.ratio_today > 0, "TF02 positive loading admitted");
  need(detail::thermal_baryon_loading_at_inverse_scale(small, std::numeric_limits<W>::max()/2).status != S::ok,
       "TF02 positive inverse-scale underflow refuses");
  const auto zero = detail::thermal_baryon_loading(prepare_thermal_background({70,.25,0,0,.5,{}}));
  const auto zero_at = detail::thermal_baryon_loading_at_inverse_scale(zero, std::numeric_limits<W>::max()/2);
  need(zero.status == S::ok && zero_at.status == S::ok && zero_at.value == 0 && zero_at.arithmetic_estimate == 0,
       "exact zero-baryon witness retained");
  need(detail::thermal_baryon_sound(r,std::numeric_limits<W>::max()/2).status == S::outside_domain,
       "TF03 sound scale outside declared a<=1 refuses");
  const auto sound=detail::thermal_baryon_sound(r,1);
  need(sound.status==S::ok && sound.sound_speed_over_c>0 && sound.sound_speed_squared_over_c_squared>0 &&
       sound.sound_speed_estimate_over_c>0 && sound.sound_speed_squared_estimate_over_c_squared>0,
       "positive sound outputs and diagnostics in admitted domain");
  auto tiny = r; tiny.arithmetic_estimate = std::ldexp(1.L,-100);
  const W a = std::ldexp(1.L,-10);
  const auto xi = detail::thermal_ruler_internal::loading_xi(tiny,a);
  need(xi && *xi > 0 && std::isnormal(*xi), "tiny positive rationalized xi survives");
  need(std::abs(*xi/(tiny.arithmetic_estimate*a/(2*(1+a)))-1) < 1e-16L,
       "tiny xi independent first-order analytic limit");
  tiny.arithmetic_estimate=0;
  need(detail::thermal_ruler_internal::loading_xi(tiny,a)==0, "exact zero loading error xi witness");
  need(detail::thermal_ruler_internal::loading_xi(r,0)==0, "exact zero scale xi witness");
}

void retained_controls(){
 const ThermalPhysicalModel source{67.4,.02237,.12,2.7255,0,{}};
 const auto mapped=map_thermal_physical_model(source);
 need(mapped.status==S::ok && mapped.model && mapped.scalar_witnesses,"one complete map owns witnesses");
 const double emitted[]{mapped.model->omega_gamma,mapped.model->omega_b,mapped.model->omega_cdm,mapped.model->omega_massless_nonphoton};
 for(unsigned i=0;i<4;++i){const auto &w=(*mapped.scalar_witnesses)[i];
  need(w.emitted_value==emitted[i] && w.measured_absolute_cast_loss==std::abs(w.wide_value-W(emitted[i])),"actual ordered map cast loss");
  need(w.wide_operation_estimate==(i==0?128.L:32.L)*std::numeric_limits<W>::epsilon()*std::abs(w.wide_value),"exact admitted map arithmetic law");
 }
 need((*mapped.scalar_witnesses)[3].wide_value==0 && (*mapped.scalar_witnesses)[3].wide_operation_estimate==0 &&
  (*mapped.scalar_witnesses)[3].measured_absolute_cast_loss==0,"exact sourcezero map witness");
 auto bad=source;bad.species={{0,-1,2}};const auto failed=map_thermal_physical_model(bad);
 need(failed.status!=S::ok && !failed.model && !failed.scalar_witnesses,"failed map cannot publish partial witnesses");
 HydrogenHeliumHistoryPolicy p;p.maximum_total_work=1;
 auto h=prepare_hydrogen_helium_history({source,.19,.015,"source-only map retention control",2700,300},p);
 need(h.status()==S::work_limit && h.thermal_mapping_witnesses(),"later refused history retains original map");
 const auto before=*h.thermal_mapping_witnesses();auto copy=h;auto moved=std::move(h);
 need(!h.thermal_mapping_witnesses() && h.status()==S::invalid_input,"history move clears source map metadata");
 for(unsigned i=0;i<4;++i)need(copy.thermal_mapping_witnesses()->at(i).emitted_value==before[i].emitted_value &&
  moved.thermal_mapping_witnesses()->at(i).wide_value==before[i].wide_value,"history copy/move preserves map bits and wide witnesses");
 std::size_t spent=0,outer=SIZE_MAX,momentum=0;
 detail::ThermalRulerWorkBudget malformed(spent,SIZE_MAX,outer,momentum,SIZE_MAX,SIZE_MAX);
 need(malformed.remaining()==0 && !malformed.charge_outer() && spent==0 && outer==SIZE_MAX,"private bounded precondition rejects malformed counters");
 spent=2;outer=1;momentum=1;detail::ThermalRulerWorkBudget ledger(spent,4,outer,momentum,2,4);
 need(ledger.charge_outer() && ledger.charge_momentum(1) && spent==4 && outer==2 && momentum==2,"checked aggregate work retains previous attempts");
 need(!ledger.charge_outer() && spent==4,"exhaustion does not wrap or erase costs");
 const auto bg=prepare_thermal_background({70,3./16,5./16,.25,.25,{}});
 const auto load=detail::thermal_baryon_loading(bg);
 spent=outer=momentum=0;detail::ThermalRulerWorkBudget zero_budget(spent,200000,outer,momentum,200000,200000);
 detail::ThermalRulerAllowance zero_allowance;zero_allowance.absolute_tolerance_mpc=1e-8;
 need(detail::integrate_thermal_ruler(bg,load,.5L,zero_allowance,zero_budget).status==S::work_limit && spent==0,
  "named zero outer allowance refuses before callback");
 // Rational diagnostic law: Eq + (|I|+Eq)*(rho+eta)/(1-eta),
 // then the loading cross term. These call the actual integral's owner.
 detail::ThermalRulerIntegral composed;
 const auto composition=detail::thermal_ruler_internal::complete_diagnostics(composed,2,.25L,true,.5L,.25L,.5L,0,false);
 const W op=2*64*std::numeric_limits<double>::epsilon();
 need(composition==S::ok && std::abs(composed.error_estimate_mpc-(4.75L+1.5L*op))<1e-17L,
  "actual production response includes cast denominator, quadrature transfer and loading cross");
 need(composed.background_estimate_mpc==1.5L && composed.outer_quadrature_estimate_mpc==.25L,
  "production category decomposition independent rational values");
 detail::ThermalRulerIntegral underflow;
 const W tiny_xi=16*std::numeric_limits<W>::min();
 need(detail::thermal_ruler_internal::complete_diagnostics(underflow,1e-10L,0,false,0,0,tiny_xi,0,false)==S::conditioning_budget_exceeded,
  "actual production positive loading contribution underflow refuses");
 need(detail::thermal_ruler_internal::complete_diagnostics(underflow,1e-10L,0,false,0,0,0,0,false)==S::ok && underflow.loading_estimate_mpc==0,
  "genuine zero loading correction admitted separately");
 need(detail::thermal_ruler_internal::complete_diagnostics(underflow,1,0,true,0,0,0,0,false)==S::conditioning_budget_exceeded,
  "positive quadrature witness cannot be replaced by represented zero");
 for(W a:{.5L,.8L}){
  spent=outer=momentum=0;detail::ThermalRulerWorkBudget b(spent,200000,outer,momentum,200000,200000);
  detail::ThermalRulerAllowance allowance;allowance.absolute_tolerance_mpc=1e-8;allowance.relative_tolerance=2e-10;allowance.maximum_outer_callbacks=200000;
  const auto q=detail::integrate_thermal_ruler(bg,load,a,allowance,b);
  const W expected=(299792.458L/70)*std::sqrt(2.L/3)*std::log1p(a);
  need(q.status==S::ok && q.error_estimate_mpc>0 && std::abs(q.value_mpc-expected)<=q.error_estimate_mpc+1e-9L,"shared ruler independent closed expression");
  need(spent==outer+momentum && spent<=200000,"ruler allattempt work categories close");
 }
}
}
int main(){try{foundation_controls();retained_controls();}catch(const std::exception&e){std::cerr<<"FAILED: "<<e.what()<<" after "<<checks<<" checks\n";return 1;}std::cout<<checks<<" shared thermal foundation checks passed\n";}
