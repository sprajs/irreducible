// Original high-precision original-a integral and 3x3 cofactor facts.
// Physical growing-mode equation shared; no production quadrature/factor oracle.
#include "irred/growth_amplitude.hpp"
#include "irred/rsd.hpp"
#include "growth_rsd_peer_facts.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
namespace c=irred::cosmology;namespace s=irred::statistics;namespace f=growth_rsd_peer_facts;
using W=long double;unsigned checks=0;W maximum_relative=0;
void need(bool b,const char*w){++checks;if(!b)throw std::runtime_error(w);}
void near(W x,W y,W allocation,const char*w){maximum_relative=std::max(maximum_relative,std::abs(x/y-1));need(std::abs(x-y)<=allocation*std::abs(y),w);}
c::FixedSigma8Source amplitude(double sigma=.811,double reference=1){return {sigma,reference,c::Sigma8Convention::linear_pressureless_total_matter_top_hat_8_over_h_mpc,c::AmplitudeTreatment::fixed_supplied,"original synthetic supplied sigma8 normalization","fixed top-hat8/h Mpc linear pressureless RMS; no spectrum prediction"};}
std::array<W,2> ode(double matter,double scale,unsigned n){
 W a0=std::min(1e-7L,W(scale)/1000),x0=std::log(a0),step=(std::log(W(scale))-x0)/n,q0=(1-W(matter))/matter*a0*a0*a0;
 std::array<W,2> y{1-2*q0/11,1-8*q0/11};
 auto rhs=[&](W x,std::array<W,2> g){W a=std::exp(x),om=1/(1+(1-W(matter))/matter*a*a*a);return std::array<W,2>{g[1]-g[0],-(3-1.5L*om)*g[1]+1.5L*om*g[0]};};
 auto add=[](auto y,auto k,W h){return std::array<W,2>{y[0]+h*k[0],y[1]+h*k[1]};};
 for(unsigned i=0;i<n;++i){W x=x0+i*step;auto k1=rhs(x,y),k2=rhs(x+step/2,add(y,k1,step/2)),k3=rhs(x+step/2,add(y,k2,step/2)),k4=rhs(x+step,add(y,k3,step));for(unsigned j=0;j<2;++j)y[j]+=step*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;}
 return {W(scale)*y[0],y[1]/y[0]};
}
void growth(){for(std::size_t i=0;i<f::input.size();++i){auto in=f::input[i];auto op=c::prepare_gr_growth(c::prepare(c::LCDM(in[0]),c::FlatFLRW{}));const std::array<double,1> a{in[1]};auto r=c::evaluate_growth_amplitude(op,amplitude(in[3],in[2]),a,3);
 need(r.status==irred::numerics::Status::ok&&r.rows.size()==1&&r.rows[0].sigma8.value&&r.rows[0].f_sigma8.value,"growth amplitude admitted");
 near(*r.rows[0].sigma8.value,f::values[i][2],1e-8L,"independent sigma8 integral");near(*r.rows[0].f_sigma8.value,f::values[i][3],1e-8L,"independent fsigma8 integral");
 auto fine=ode(in[0],in[1],16384),coarse=ode(in[0],in[1],8192),reference=ode(in[0],in[2],16384);W expected=in[3]*fine[0]/reference[0]*fine[1];
 near(*r.rows[0].f_sigma8.value,expected,1e-8L,"original GR ODE amplitude comparison");for(unsigned j=0;j<2;++j)near(fine[j],coarse[j],.05L*1e-8L,"ODE reference refinement allocation");
 if(in[1]==in[2])need(*r.rows[0].sigma8.value==in[3],"exact shared-reference cancellation");
 }}
void density(){irred::rsd::DensityInput in;auto&d=in.source;d.scale_factors={.5,2./3,1};d.observed={.4,.46,.43};d.ordered_ids={"RSD-a","RSD-b","RSD-c"};d.event_ids={"synthetic-a","synthetic-b","synthetic-c"};d.covariance_axis_ids=d.ordered_ids;d.role=irred::rsd::RowRole::synthetic_control;d.covariance_unit=irred::rsd::CovarianceUnit::dimensionless_f_sigma8_squared;d.amplitude_convention=amplitude().convention;d.table_identity="original synthetic RSD vector";d.covariance_identity="original full correlated3x3";d.ordering_provenance="three exact ordered synthetic axes";d.calibration_provenance="fixed synthetic amplitude";d.dependence_provenance="single full joint covariance";in.covariance={.01,.002,-.001,.002,.015,.0005,-.001,.0005,.02};
 auto op=irred::rsd::prepare_density(std::move(in),{3,9,4096,8*1024*1024,1e-10});need(op.status()==s::DensityStatus::finite,"retained correlated RSD source");
 irred::rsd::ModelPoint point{c::LCDM(.315),c::FlatFLRW{},amplitude()};irred::rsd::DensityPolicy p;p.maximum_models=1;p.maximum_queries=3;p.maximum_string_bytes=4096;p.maximum_native_bytes=16*1024*1024;p.maximum_total_callbacks=1000000;p.maximum_forward_sensitivity=1e-10;p.requested=7;
 auto batch=op.evaluate(std::span(&point,1),p);need(batch.status==s::DensityStatus::finite&&batch.slots.size()==1&&batch.slots[0].result,"synthetic normalized density accepted");auto&r=*batch.slots[0].result;need(r.density.status==s::DensityStatus::finite,"density separate projection gate");
 for(std::size_t i=0;i<3;++i){need(batch.slots[0].predictions[i].value.has_value(),"ordered predictions retained");near(*batch.slots[0].predictions[i].value,f::rsd_predictions[i],1e-8L,"independent RSD prediction");}
 const std::array<double,4> values{r.quadratic,r.log_determinant,r.normalization,r.density.log_value};for(std::size_t i=0;i<4;++i)need(std::abs(W(values[i])-f::rsd_density[i])<=1e-8L,"independent cofactor density component");
 p.maximum_projection_log_density_error=1e-30;auto fail=op.evaluate(std::span(&point,1),p);need(!fail.slots[0].result&&fail.slots[0].predictions[0].value,"density refusal preserves prediction");
}
}
int main(){try{growth();density();std::cout<<"PASS "<<checks<<" independent growth/RSD controls; max relative="<<double(maximum_relative)<<'\n';}catch(const std::exception&e){std::cerr<<"FAIL "<<checks<<": "<<e.what()<<'\n';return 1;}}
