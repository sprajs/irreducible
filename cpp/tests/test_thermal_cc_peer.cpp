// Independent original tanh-sinh60/90-digit scalar facts and GL16/32
// fixed Gauss-Legendre composite quadrature. Shared physical equations/constants
// and libm ancestry; direct production is differential evidence only.
#include "irred/thermal_neutrino.hpp"
#include "hydrogen_relic_reference.hpp"
#include <array>
#include <iostream>
namespace {
using namespace irred::cosmology;using S=irred::numerics::Status;using W=long double;
unsigned checks=0;void need(bool b,const char* m){++checks;if(!b)throw std::runtime_error(m);}
void near(W a,W b,W budget,const char* m){need(std::abs(a-b)<=budget,m);}
ThermalPolicy policy(){ThermalPolicy p;p.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis;return p;}
struct Anchor{double y;W rho,P;};
constexpr std::array<Anchor,7> facts{{
{0,5.682196976983475505459019L,1.894065658994491835153006L},
{.001,5.682197388216540604845209L,1.894065521917083721577896L},
{.1,5.686292609702345663989564L,1.892709606544609801033724L},
{1,6.045645184985879384949935L,1.791793950630836860129931L},
{10,19.12519061523234664911632L,.6895865848064602173201004L},
{100,180.4251007630648469226090L,.07765096968166704910351685L},
{10000,18030.85471393754948386860L,.0007776956972463193579237046L}}};
void moments(){
 for(auto f:facts){auto r=evaluate_thermal_moments(f.y,policy());W s=std::hypot(1.L,W(f.y));need(r.status==S::ok&&r.rho_moment&&r.pressure_moment,"independent CC moments admitted");near(*r.rho_moment/s,f.rho/s,1e-10L+2e-10L*f.rho/s,"tanh-sinh rho");near(*r.pressure_moment*s,f.P*s,1e-10L+2e-10L*f.P*s,"tanh-sinh pressure");need(r.rho_error_estimate/s<=1e-12L+2e-12L*f.rho/s&&r.pressure_error_estimate*s<=1e-12L+2e-12L*f.P*s,"unchanged normalized diagnostic");}
 hydrogen_relic_reference::MomentumRule<16> coarse;hydrogen_relic_reference::MomentumRule<32> fine;
 W previous=0;
 for(int j=-10;j<=10;++j){double y=std::pow(10.,j/2.);W s=std::hypot(1.L,W(y)),a=coarse.moment(y)/s,b=fine.moment(y)/s,allocation=1e-10L+2e-10L*std::abs(b);near(a,b,allocation*.05L,"independent GL refinement <=5%");auto r=evaluate_thermal_moments(y,policy());need(r.status==S::ok&&r.rho_moment&&r.pressure_moment,"logarithmic CC grid admitted");near(*r.rho_moment/s,b,allocation,"independent GL density");need(*r.rho_moment>=previous&&*r.pressure_moment>0&&*r.pressure_moment<=*r.rho_moment/3*(1+1e-14),"monotonic and causal moments");previous=*r.rho_moment;}
 auto p=policy();p.maximum_callbacks_per_evaluation=1;auto f=evaluate_thermal_moments(1,p);need(f.status==S::work_limit&&f.callbacks<=1&&!f.rho_moment&&!f.pressure_moment,"odd paired budget truthful refusal");
}
void background(){
 ThermalFlatModel m{67.4,5.443781491563183e-5,3.751e-5,.02237/(.674*.674),.12/(.674*.674),{{.06,.00016818966050500852,2}}};
 auto p=policy();auto op=prepare_thermal_background(m,p);need(op.status()==S::ok,"CC retained background");std::array<double,7> z{0,.1,1,10,100,1000,1e6},a{};for(size_t i=0;i<7;++i)a[i]=1/(1+z[i]);
 // CLASSv3.3.0 commit0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18
 // original adapter, explicit matched FD amplitude and legacy SI conversion;
 // frozen trapezoid4096/8192+ODErefinement<=5.972e-12rel (<5%allocation).
 constexpr std::array<W,7> E{1,1.0508317697470424L,1.7902809950671763L,20.519876107077128L,578.34844510805067L,20477.937796412822L,10236606364.874552L};
 auto r=op.evaluate(a,thermal_e|thermal_h,p);need(r.status==S::ok&&r.rows.size()==7,"CC background batch");for(size_t i=0;i<7;++i){need(r.rows[i].e.value&&r.rows[i].h_km_s_mpc.value,"requested background outputs");near(*r.rows[i].e.value,E[i],2e-10L*E[i],"matched independent CLASS E");near(*r.rows[i].h_km_s_mpc.value,67.4L*E[i],2e-10L*67.4L*E[i],"matched independent CLASS H");}
 need(op.evaluate(a,thermal_e).status!=S::ok,"unreported mixed method refused");
}
}
int main(){try{moments();background();std::cout<<"PASS "<<checks<<" CC peer controls\n";}catch(const std::exception&e){std::cerr<<"FAIL "<<checks<<": "<<e.what()<<'\n';return 1;}}
