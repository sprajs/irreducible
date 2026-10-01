// Original synthetic peer controls; no copied CLASS/CAMB implementation.
// Independent references frozen before native evaluation in primary evidence
// prerequisite-validation-20261001/peer: original mpmath 1.3.0 tanh-sinh,
// 60/90 decimal digits, tails 128/160, normalized refinement <2.38e-48.
// Product normalized moments: 1e-10 + 2e-10*abs(ref); estimator budgets
// 1e-12 + 2e-12*abs(normalized ref). Physical state and E/H: 2e-10 relative.
// Analytic limits and continuity allocations are stated at their controls.
// Zero chemical potential, g counts states (g=2 includes nu+antinu),
// T0 and m in eV, natural density/pressure in eV^4. No Planck substitution.
#include "irred/thermal_neutrino.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
namespace {
using namespace irred::cosmology;
using S=irred::numerics::Status;
using W=long double;
unsigned checks=0;
void check(bool ok,const char* name){++checks;if(!ok)throw std::runtime_error(name);}
void moment_near(double got,W ref,const char* name){check(std::abs(W(got)-ref)<=1e-10L+2e-10L*std::abs(ref),name);}
void relative_near(double got,W ref,const char* name){check(std::abs(W(got)-ref)<=2e-10L*std::abs(ref),name);}
struct Anchor{double y;W rho,pressure;};
constexpr std::array<Anchor,7> anchors{{
 {0,5.682196976983475505459019L,1.894065658994491835153006L},
 {.001,5.682197388216540604845209L,1.894065521917083721577896L},
 {.1,5.686292609702345663989564L,1.892709606544609801033724L},
 {1,6.045645184985879384949935L,1.791793950630836860129931L},
 {10,19.12519061523234664911632L,.6895865848064602173201004L},
 {100,180.4251007630648469226090L,.07765096968166704910351685L},
 {10000,18030.85471393754948386860L,.0007776956972463193579237046L}
}};
void moments(){
 W previous_rho=0,previous_pressure=2;
 for(const auto& a:anchors){auto got=evaluate_thermal_moments(a.y);
  check(got.status==S::ok&&got.rho_moment&&got.pressure_moment,"FD moments admitted");
  W scale=std::hypot(1.L,W(a.y));
  moment_near(*got.rho_moment/scale,a.rho/scale,"independent normalized FD rho");
  moment_near(*got.pressure_moment*scale,a.pressure*scale,"independent normalized FD pressure");
  check(got.rho_error_estimate>0&&got.rho_error_estimate/scale<=1e-12L+2e-12L*a.rho/scale,"rho normalized diagnostic allocation");
  check(got.pressure_error_estimate>0&&got.pressure_error_estimate*scale<=1e-12L+2e-12L*a.pressure*scale,"pressure normalized diagnostic allocation");
  check(*got.rho_moment>=previous_rho&&*got.pressure_moment<=previous_pressure,"mass monotonic density and pressure");
  check(*got.pressure_moment>0&&*got.pressure_moment<=*got.rho_moment/3*(1+1e-14),"0<w<=one third");
  previous_rho=*got.rho_moment;previous_pressure=*got.pressure_moment;
 }
 auto zero=evaluate_thermal_moments(0);const W pi=std::numbers::pi_v<W>,rho=7*pi*pi*pi*pi/120;
 moment_near(*zero.rho_moment,rho,"analytic massless rho");moment_near(*zero.pressure_moment,rho/3,"analytic massless pressure");
 // Asymptotic leading limits at y=10000 have O(y^-2) fractional correction;
 // frozen 1e-6 allocation is for this truncation, not a numerical budget.
 constexpr W zeta3=1.202056903159594285399738L,zeta5=1.036927755143369926331365L;
 auto heavy=evaluate_thermal_moments(10000);W nr=1.5L*zeta3*10000,np=7.5L*zeta5/10000;
 check(std::abs(W(*heavy.rho_moment)/nr-1)<1e-6L,"analytic nonrelativistic rho leading term");
 check(std::abs(W(*heavy.pressure_moment)/np-1)<1e-6L,"analytic nonrelativistic pressure leading term");
 auto tiny=evaluate_thermal_moments(1e-6);check(tiny.status==S::ok,"small mass continuity admitted");
 check(std::abs(W(*tiny.rho_moment)-rho)<1e-10L,"massless continuity");
 check(evaluate_thermal_moments(-1).status!=S::ok,"negative y refused");
 check(evaluate_thermal_moments(std::numeric_limits<double>::infinity()).status!=S::ok,"nonfinite y refused");
 ThermalPolicy limited;limited.maximum_callbacks_per_evaluation=4;
 auto failed=evaluate_thermal_moments(1,limited);check(failed.status==S::work_limit&&!failed.rho_moment&&!failed.pressure_moment,"callback refusal no partial values");
}
void species(){
 ThermalSpecies s{1,1,2};auto got=evaluate_thermal_species(s,1);
 check(got.status==S::ok&&got.rho_ev4&&got.pressure_ev4,"explicit species admitted");
 const W pi=std::numbers::pi_v<W>,factor=1/(pi*pi);
 relative_near(*got.rho_ev4,anchors[3].rho*factor,"state g2 rho convention");
 relative_near(*got.pressure_ev4,anchors[3].pressure*factor,"state g2 pressure convention");
 auto one=evaluate_thermal_species({1,1,1},1);relative_near(*one.rho_ev4,W(*got.rho_ev4)/2,"one state not paired twice");
 auto scaled=evaluate_thermal_species({2,2,2},1);relative_near(*scaled.rho_ev4,16*W(*got.rho_ev4),"temperature fourth power same y");
 relative_near(*scaled.pressure_ev4,16*W(*got.pressure_ev4),"pressure temperature fourth power");
 // Collisionless conservation d rho/d ln(a)=-3(rho+P), central log step
 // 1e-4. Frozen 1e-7 fractional allocation includes O(step^2) truncation
 // and subtraction of independently integrated density outputs.
 constexpr double a=.5,h=1e-4;
 auto center=evaluate_thermal_species(s,a),plus=evaluate_thermal_species(s,a*std::exp(h)),minus=evaluate_thermal_species(s,a*std::exp(-h));
 check(center.status==S::ok&&plus.status==S::ok&&minus.status==S::ok,"continuity species admitted");
 W derivative=(W(*plus.rho_ev4)-W(*minus.rho_ev4))/(2*h),expected=-3*(W(*center.rho_ev4)+W(*center.pressure_ev4));
 check(std::abs(derivative-expected)<=1e-7L*std::abs(expected),"independent continuity equation");
 for(const ThermalSpecies invalid:std::array<ThermalSpecies,3>{{{-1,1,2},{1,0,2},{1,1,0}}}){
  auto bad=evaluate_thermal_species(invalid,1);check(bad.status!=S::ok&&!bad.rho_ev4&&!bad.pressure_ev4,"invalid species no values");}
 check(evaluate_thermal_species(s,0).status!=S::ok&&evaluate_thermal_species(s,1.01).status!=S::ok,"scale factor domain");
}
void background(){
 ThermalFlatModel model{70,5e-5,2e-5,.05,.25,{{.06,.000168,2}}};
 auto retained=prepare_thermal_background(model);check(retained.status()==S::ok,"explicit synthetic background admitted");
 relative_near(*retained.omega_species_today(),.001310257272267406042397192L,"independent SI critical density omega nu");
 relative_near(*retained.omega_lambda(),.6986197427277325939576028L,"flat closure excludes nu from cb");
 std::array<double,4> a{1,.1,.01,.001};
 constexpr std::array<W,4> E{1,17.39871883575367278847808L,555.7022615469616192741175L,19536.23533993629128441921L};
 constexpr std::array<W,4> rho{5.197777235343074523098948e-14L,5.223721787446074351696885e-11L,7.157659011299627146846238e-8L,.0004627291507663699674916468L};
 constexpr std::array<W,4> pressure{1.757328858352840311867935e-18L,1.737096785403971994310422e-13L,1.065865484856115094392202e-8L,.0001515647298014865362802269L};
 auto result=retained.evaluate(a,thermal_e|thermal_h);check(result.status==S::ok&&result.rows.size()==a.size(),"E/H batch admitted");
 for(size_t i=0;i<a.size();++i){const auto& row=result.rows[i];check(row.e.status==S::ok&&row.h_km_s_mpc.status==S::ok&&row.e.value&&row.h_km_s_mpc.value,"requested per output statuses");
  relative_near(*row.e.value,E[i],"independent highprecision E");relative_near(*row.h_km_s_mpc.value,70*E[i],"independent highprecision H");
  auto state=evaluate_thermal_species(model.species[0],a[i]);check(state.status==S::ok,"background species state");
  relative_near(*state.rho_ev4,rho[i],"independent physical rho eV4");relative_near(*state.pressure_ev4,pressure[i],"independent physical pressure eV4");
 }
 model.species[0].mass_ev=0;auto lighter=prepare_thermal_background(model);
 auto le=lighter.evaluate(a,thermal_e);check(le.status==S::ok,"massless supplied species background");
 check(std::abs(*le.rows[0].e.value-1)<1e-14,"flat normalization independent species mass");
 // All supplied components zero except dust: exact matter-only E=a^-3/2.
 auto dust=prepare_thermal_background({70,0,0,.05,.95,{}});
 auto dustrow=dust.evaluate(a,thermal_e);check(dustrow.status==S::ok,"species-free matter limit");
 for(size_t i=0;i<a.size();++i)relative_near(*dustrow.rows[i].e.value,std::pow(W(a[i]),-1.5L),"analytic matter E limit");
 auto impossible=prepare_thermal_background({70,0,0,1,1,{}});check(impossible.status()!=S::ok,"negative lambda refused no compensation");
 ThermalPolicy payload;payload.maximum_native_bytes=0;
 check(prepare_thermal_background(model,payload).status()!=S::ok,"payload admission respected");
}
void matched_class(){
 // Independent original adapter calls CLASS v3.3.0 background_init and
 // background_at_z, pinned commit 0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18.
 // Explicit zero-chemical-potential FD family: CLASS deg=1 includes both
 // states, whereas this API uses g=2. All input fractions below are explicit.
 // Fixed legacy CLASS h/eV/G/Mpc differ from this API. At fixed supplied
 // T_eV=Tcmb*T_ncdm*kB_CLASS/eV_CLASS, amplitude is matched with
 // deg_CLASS=(Gnew/Gold)*(hold/hnew)^3*(eVnew/eVold)^4*(Mpcnew/Mpcold)^2
 // =1.000002824932919, independently derived from pinned primary source.
 // CLASS manual trapezoid qmax80 bins4096/8192, tightened background ODE,
 // frozen E reference refinement <=5.972e-12 relative, less than 5% of the
 // fixed 2e-10 relative comparison. Earlier failed adaptive-mesh arbitrary-y
 // reference preserved. These are synthetic matched E values, not a Planck
 // posterior, automatic Neff mapping, perturbations or distance qualification.
 ThermalFlatModel m{67.4,5.443781491563183e-5,3.751e-5,
                    .02237/(.674*.674),.12/(.674*.674),
                    {{.06,.00016818966050500852,2}}};
 auto op=prepare_thermal_background(m);check(op.status()==S::ok,"matched CLASS model admitted");
 relative_near(*op.omega_species_today(),.0014180872495669578L,"matched CLASS explicit nu fraction");
 std::array<double,7> z{0,.1,1,10,100,1000,1e6},a{};
 for(size_t i=0;i<a.size();++i)a[i]=1/(1+z[i]);
 constexpr std::array<W,7> ref{1,1.0508317697470424L,1.7902809950671763L,
   20.519876107077128L,578.34844510805067L,20477.937796412822L,10236606364.874552L};
 auto got=op.evaluate(a,thermal_e|thermal_h);check(got.status==S::ok&&got.rows.size()==ref.size(),"matched CLASS E/H batch");
 for(size_t i=0;i<ref.size();++i){const auto& row=got.rows[i];check(row.e.status==S::ok&&row.e.value&&row.h_km_s_mpc.status==S::ok&&row.h_km_s_mpc.value,"matched per-output admission");
  relative_near(*row.e.value,ref[i],"external CLASS matched E");
  relative_near(*row.h_km_s_mpc.value,67.4L*ref[i],"external CLASS matched H");}
}
}
int main(){try{moments();species();background();matched_class();
 const auto old=std::fegetround();check(std::fesetround(FE_DOWNWARD)==0,"set hostile rounding");
 auto bad=irred::cosmology::evaluate_thermal_moments(1);check(std::fesetround(old)==0,"restore rounding");
 check(bad.status!=S::ok&&!bad.rho_moment&&!bad.pressure_moment,"unsupported arithmetic no accepted moment");
 std::cout<<"PASS "<<checks<<" thermal neutrino peer controls\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
