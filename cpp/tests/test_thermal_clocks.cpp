#include "irred/thermal_clocks.hpp"
#include "irred/quantities.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace {
using namespace irred::cosmology;
using S=irred::numerics::Status;
using W=long double;
constexpr unsigned all=(1u<<thermal_clock_output_count)-1;
unsigned checks=0;
void need(bool ok,const char *message) { ++checks;if(!ok)throw std::runtime_error(message); }
void near(const ThermalBackgroundValue &out,W reference,bool clock,const char *message) {
  need(out.status==S::ok&&out.value.has_value(),message);
  const W budget=1e-8L+2e-10L*std::abs(reference);
  need(std::abs(W(*out.value)-reference)<=budget,message);
  need(out.error_estimate<=budget&&out.error_estimate>=0,"reported output admission");
  (void)clock;
}
constexpr W seconds_gyr=365.25L*86400*1000000000;
W age_scale() { return irred::megaparsec_in_metres_wide()/70000/seconds_gyr; }
W horizon_scale() { return W(irred::speed_of_light_m_per_s)/70000; }
void analytic() {
  const double a[]{0,.001,.1,.5,1,std::nextafter(1.,0.)};
  for(auto model:{ThermalFlatModel{70,1,0,0,0,{}},ThermalFlatModel{70,0,0,1,0,{}}}) {
    const bool radiation=model.omega_gamma>0;
    auto owner=prepare_thermal_clocks(prepare_thermal_background(model));
    need(owner.status()==S::ok,"analytic clock prepare");
    auto result=owner.evaluate(a,all);need(result.status==S::ok&&result.rows.size()==6,"analytic ordered rows");
    for(auto &row:result.rows) {
      W x=row.scale_factor;
      const W t=radiation?x*x/2:2*std::pow(x,1.5L)/3;
      const W eta=radiation?x:2*std::sqrt(x);
      // Direct near-one expressions avoid a cancellation-limited reference.
      const W lookback=radiation?(1-x)*(1+x)/2:
          -2*std::expm1(1.5L*std::log(x))/3;
      near(row.outputs[0],age_scale()*t,true,"analytic proper age");
      near(row.outputs[1],age_scale()*lookback,true,"analytic direct lookback");
      near(row.outputs[2],horizon_scale()*eta,false,"analytic comoving horizon");
      near(row.outputs[3],horizon_scale()*x*eta,false,"analytic proper horizon");
      need(row.callbacks==row.outer_callbacks+row.momentum_callbacks&&row.momentum_callbacks==0,"analytic work receipt");
      if(x>0&&radiation)need(row.age_tail&&row.conformal_tail,"full radiation support retained");
      if(x==0)need(row.outputs[0].value==0&&row.outputs[2].value==0,"finite big-bang endpoint");
    }
  }
  auto lcdm=prepare_thermal_clocks(prepare_thermal_background({70,0,0,.3,0,{}}));
  const auto rows=lcdm.evaluate(a,thermal_clock_mask(ThermalClockOutput::age_gyr));
  for(auto &row:rows.rows) {
    const W t=2/(3*std::sqrt(.7L))*std::asinh(std::sqrt(.7L/.3L)*std::pow(W(row.scale_factor),1.5L));
    near(row.outputs[0],age_scale()*t,true,"independent flat matter-Lambda asinh age");
    need(!row.outputs[1].value&&!row.conformal_tail,"unrequested integrals omitted");
  }
  auto vacuum=prepare_thermal_clocks(prepare_thermal_background({70,0,0,0,0,{}}));
  auto result=vacuum.evaluate(a,all);
  for(auto &row:result.rows) {
    need(row.outputs[0].status==S::outside_domain&&!row.outputs[0].value&&
         row.outputs[2].status==S::outside_domain&&!row.outputs[2].value,"pure Lambda divergence retained");
    if(row.scale_factor>0)near(row.outputs[1],-age_scale()*std::log(W(row.scale_factor)),true,"finite pure Lambda interval");
    else need(!row.outputs[1].value,"pure Lambda zero endpoint direct lookback divergent");
  }
}
void controls() {
  auto owner=prepare_thermal_clocks(prepare_thermal_background({70,.0001,0,.05,.25,{{.06,.0002,2}}}));
  need(owner.status()==S::ok,"FD clock preparation");
  const double factors[]{.001,.1,1};
  auto result=owner.evaluate(factors,all);
  need(result.status==S::ok,"FD batch admission");
  for(const auto &row:result.rows) {
    for(unsigned i=0;i<row.outputs.size();++i) {
      const auto &value=row.outputs[i];
      if(value.status!=S::ok||!value.value)std::cerr<<"FD refusal a="<<row.scale_factor
          <<" output="<<i<<" status="<<static_cast<unsigned>(value.status)
          <<" callbacks="<<row.callbacks<<" outer="<<row.outer_callbacks<<" momentum="<<row.momentum_callbacks<<'\n';
      need(value.status==S::ok&&value.value,"FD clock output admission");
    }
    need(row.age_tail&&row.conformal_tail&&row.age_tail->midpoint>0&&row.conformal_tail->midpoint>0,"positive early support retained");
    need(row.callbacks==row.outer_callbacks+row.momentum_callbacks&&row.callbacks>0,"all FD and tail work counted");
    near(row.outputs[3],W(row.scale_factor)**row.outputs[2].value,false,"proper-comoving horizon identity");
  }
  // These independently refined finite differences have their own allocations.
  auto derivative=[&](double a,W step,unsigned output) {
    std::array<double,4> x{double(a*std::exp(-2*step)),double(a*std::exp(-step)),
                           double(a*std::exp(step)),double(a*std::exp(2*step))};
    const auto rows=owner.evaluate(x,1u<<output);
    for(auto &row:rows.rows)need(row.outputs[output].value.has_value(),"derivative supporting rows");
    return (W(*rows.rows[0].outputs[output].value)-8*W(*rows.rows[1].outputs[output].value)+
            8*W(*rows.rows[2].outputs[output].value)-W(*rows.rows[3].outputs[output].value))/(12*step);
  };
  for(double a:{.01,.5})for(unsigned output:{0u,2u}) {
    const W coarse=derivative(a,.002L,output),fine=derivative(a,.001L,output);
    const auto d=owner.background().scaled_expansion(a);
    need(d.status==S::ok,"derivative background dependency");
    const W reference=(output==0?age_scale()*a*a:horizon_scale()*a)/std::sqrt(d.a4_e2);
    const W budget=(output==0?2e-6L:2e-4L)+2e-7L*std::abs(reference);
    need(std::abs(coarse-fine)<=.05L*budget,"derivative independent step refinement");
    need(std::abs(fine-reference)<=budget,"clock log-a derivative identity");
  }
}
void adversarial() {
  auto owner=prepare_thermal_clocks(prepare_thermal_background({70,.0001,0,.05,.25,{{.06,.0002,2}}}));
  const double factors[]{.1,std::numeric_limits<double>::quiet_NaN(),-1,1.1,1};
  auto result=owner.evaluate(factors,all);
  need(result.rows.size()==5&&result.rows[1].outputs[0].status==S::nonfinite_input&&
       result.rows[2].outputs[0].status==S::outside_domain&&result.rows[3].outputs[0].status==S::outside_domain&&
       result.rows[4].outputs[0].status==S::ok,"failed rows retained without dropping");
  auto policy=ThermalClockPolicy{};policy.maximum_total_callbacks=1;
  auto failed=owner.evaluate(std::span(factors,1),all,policy);
  need(!failed.rows[0].outputs[0].value&&failed.callbacks<=1,"original work ceiling");
  policy={};policy.maximum_native_bytes=1;
  need(owner.evaluate(factors,all,policy).status==S::work_limit,"native bytes admission");
  policy={};policy.maximum_points=1;
  need(owner.evaluate(factors,all,policy).status==S::work_limit,"ordered points admission");
  policy={};policy.maximum_tail_refinements=0;
  failed=owner.evaluate(std::span(factors,1),all,policy);
  need(!failed.rows[0].outputs[0].value&&failed.rows[0].outputs[1].value,"early-tail refusal preserves independent direct lookback");
  policy={};policy.absolute_tolerance_gyr=1e-30;policy.relative_tolerance=0;
  failed=owner.evaluate(std::span(factors,1),thermal_clock_mask(ThermalClockOutput::lookback_gyr),policy);
  need(!failed.rows[0].outputs[1].value,"unattainable budget refusal");
  policy={};policy.thermal.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis;
  need(owner.evaluate(factors,all,policy).status==S::invalid_input,"method identity");
  need(owner.evaluate(factors,0).status==S::invalid_input&&owner.evaluate(factors,1u<<31).status==S::invalid_input,"output mask admission");
  auto copy=owner;auto moved=std::move(owner);
  need(owner.status()==S::invalid_input&&copy.status()==S::ok&&moved.status()==S::ok,"owned copy and move invalidation");
  moved=std::move(moved);need(moved.status()==S::ok,"self move preserved");
  std::fesetround(FE_UPWARD);need(copy.evaluate(factors,all).status==S::invalid_input,"rounding admission");std::fesetround(FE_TONEAREST);
  need(!thermal_clocks_payload_bound(std::numeric_limits<std::size_t>::max(),16),"payload overflow");
  auto raw=prepare_thermal_background({70,.0001,0,.05,.25,{}});
  auto uncertain=prepare_thermal_clocks(raw,{}, {0,1.});
  failed=uncertain.evaluate(std::span(factors,1),all);
  need(!failed.rows[0].outputs[0].value,"unresolved positive denominator refused");
}
}
int main() { try {analytic();controls();adversarial();std::cout<<"PASS thermal clocks "<<checks<<" checks\n";}
  catch(const std::exception &e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;} }
