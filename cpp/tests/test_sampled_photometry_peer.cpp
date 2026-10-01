// Synthetic controls authored before implementation outputs. No copied external code.
// Independent observed-frequency Simpson reference: Fnu=s*Lnu(s nu)/(4 pi D^2),
// Lnu=lambda_rest^2 Llambda/c. Photon weighting is 1/(h nu), T enters once.
// Split every rest/observer interpolation knot. Reference refinement allocation
// 2e-13 relative; product 2e-12 relative+1e-300 with separate positivity.
// Shared exact SI definitions and system pi are ancestry, not independent physics.
#include "irred/sampled_photometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>
namespace {
using W=long double;
constexpr W c=299792458.L,h=6.62607015e-34L;
unsigned checks=0;
void check(bool b,const char* label){++checks;if(!b)throw std::runtime_error(label);}
W interp(std::span<const double> x,std::span<const double> y,W at){
 if(at<x.front()||at>x.back())return 0;
 for(size_t i=1;i<x.size();++i)if(at<=x[i])return y[i-1]+(at-x[i-1])*(W(y[i])-y[i-1])/(W(x[i])-x[i-1]);
 return y.back();
}
std::array<W,3> oracle(const irred::photometry::SampledInput& in,unsigned panels){
 const W s=1+W(in.redshift),d=in.luminosity_distance_metre;
 W lo=std::max(s*in.spectrum.wavelength_metre.front(),W(in.passband.wavelength_metre.front()));
 W hi=std::min(s*in.spectrum.wavelength_metre.back(),W(in.passband.wavelength_metre.back()));
 if(hi<=lo)return {};
 std::vector<W> knots{lo,hi};
 for(double x:in.spectrum.wavelength_metre)if(s*x>lo&&s*x<hi)knots.push_back(s*x);
 for(double x:in.passband.wavelength_metre)if(x>lo&&x<hi)knots.push_back(x);
 std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
 std::array<W,3> out{};
 for(size_t k=1;k<knots.size();++k){
  const W a=c/knots[k],b=c/knots[k-1],step=(b-a)/panels;
  std::array<W,3> sum{};
  for(unsigned j=0;j<=panels;++j){
   W nu=a+j*step,lr=c/(s*nu),observed=c/nu;
   // Clamp roundoff at support edges in this reference only; named inputs have broad support.
   lr=std::clamp(lr,W(in.spectrum.wavelength_metre.front()),W(in.spectrum.wavelength_metre.back()));
   observed=std::clamp(observed,W(in.passband.wavelength_metre.front()),W(in.passband.wavelength_metre.back()));
   const W lnu=lr*lr*interp(in.spectrum.wavelength_metre,in.spectrum.luminosity_watt_per_metre,lr)/c;
   const W fnu=s*lnu/(4*std::numbers::pi_v<W>*d*d),t=interp(in.passband.wavelength_metre,in.passband.optical_transmission,observed);
   const W weight=(j==0||j==panels)?1:(j%2?4:2);
   sum[0]+=weight*fnu;sum[1]+=weight*fnu*t;sum[2]+=weight*fnu*t/(h*nu);
  }
  for(unsigned i=0;i<3;++i)out[i]+=step*sum[i]/3;
 }
 out[1]*=W(in.collecting_area_square_metre)*in.observer_exposure_second;
 out[2]*=W(in.collecting_area_square_metre)*in.observer_exposure_second;
 return out;
}
void near(const irred::photometry::Outcome& out,W expected,const char* label){
 using namespace irred::photometry;
 check(out.availability==Availability::available&&out.value.has_value(),label);
 if(expected>0)check(*out.value>0,label);
 check(std::abs(W(*out.value)-expected)<=2e-12L*std::abs(expected)+1e-300L,label);
}
void compare(const irred::photometry::SampledInput& in){
 auto a=oracle(in,8192),b=oracle(in,16384);
 for(unsigned i=0;i<3;++i){if(std::abs(a[i]-b[i])>2e-13L*std::abs(b[i]))std::cerr<<"reference component "<<i<<" relative refinement "<<double(std::abs(a[i]-b[i])/std::abs(b[i]))<<" z "<<in.redshift<<"\n";check(std::abs(a[i]-b[i])<=2e-13L*std::abs(b[i]),"frequency refinement");}
 auto out=irred::photometry::evaluate_sampled(in,{7,64,64});
 check(out.admission_status==irred::numerics::Status::ok,"sample admission");
 near(out.flux_watt_per_square_metre,b[0],"frequency incident");near(out.energy_joule,b[1],"frequency energy");near(out.expected_photons,b[2],"frequency photons");
}
}
int main(){try{
 using namespace irred::photometry;
 constexpr double m=1e-6;
 std::array<double,2> x{m,3*m},l{2,6},t{.25,.75};
 SampledInput in{{x,l},{x,t},2,0,3,5};
 auto out=evaluate_sampled(in,{7,64,64});
 const W norm=2/(4*std::numbers::pi_v<W>*4),scale=W(m);
 near(out.flux_watt_per_square_metre,4*norm*scale,"linear analytic incident");
 near(out.energy_joule,15*13.L/6*norm*scale,"linear analytic energy");
 near(out.expected_photons,15*5*norm*scale*scale/(h*c),"linear analytic photons");
 compare(in);
 std::array<double,2> ox{2*m,6*m};in.passband.wavelength_metre=ox;in.redshift=1;
 auto red=evaluate_sampled(in,{7,64,64});
 near(red.flux_watt_per_square_metre,*out.flux_watt_per_square_metre.value,"bolometric redshift incident");
 near(red.energy_joule,*out.energy_joule.value,"observer time redshift energy");
 near(red.expected_photons,2*W(*out.expected_photons.value),"redshift photon energy");compare(in);
 std::array<double,4> sx{m,1.5*m,2.25*m,4*m},sl{1,7,.5,3};
 std::array<double,5> px{1.2*m,1.75*m,2*m,2.8*m,3.7*m},pt{0,.8,.1,1,.3};
 in={{sx,sl},{px,pt},1e12,0,2,17};compare(in);
 in.redshift=.7;compare(in);
 auto old=evaluate_sampled(in,{7,64,64});pt.fill(0);
 auto zero=evaluate_sampled(in,{7,64,64});near(zero.flux_watt_per_square_metre,*old.flux_watt_per_square_metre.value,"T independent incident");near(zero.energy_joule,0,"zero T energy");near(zero.expected_photons,0,"zero T photons");
 sl[1]=-1;check(evaluate_sampled(in,{7,64,64}).admission_status!=irred::numerics::Status::ok,"zero T cannot mask invalid source");sl[1]=7;
 px[2]=px[1];check(evaluate_sampled(in,{7,64,64}).admission_status!=irred::numerics::Status::ok,"duplicate wavelengths");px[2]=2*m;
 pt[2]=std::numeric_limits<double>::quiet_NaN();check(evaluate_sampled(in,{1,64,64}).admission_status!=irred::numerics::Status::ok,"unrequested transmission still admission");pt.fill(1);
 // Binary-exact thin support, no redshift/product edge uncertainty.
 std::array<double,2> thin{1,std::nextafter(1.,2.)},ones{1,1};
 in={{thin,ones},{thin,ones},1,0,1,1};out=evaluate_sampled(in,{7,64,64});
 near(out.flux_watt_per_square_metre,(W(thin[1])-1)/(4*std::numbers::pi_v<W>),"adjacent binary wavelength positivity");
 // Exact binary input witnesses: true positive overlap must never report zero.
 const double s=std::ldexp(1.,-20),e=std::ldexp(1.,-52);
 std::array<double,2> edge_source{s,s*(1+e)},edge_band{s*(1+2*e),2*s};
 in={{edge_source,ones},{edge_band,ones},1,e,1,1};
 auto witness=[&](W width){auto v=evaluate_sampled(in,{1,64,64});
  check(v.admission_status==irred::numerics::Status::ok,"thin witness admission");
  if(v.flux_watt_per_square_metre.availability==Availability::available)
   near(v.flux_watt_per_square_metre,width/(4*std::numbers::pi_v<W>*(1+W(in.redshift))),"exact thin witness positive");
  else check(v.flux_watt_per_square_metre.availability==Availability::failed&&v.flux_watt_per_square_metre.numerical_status==irred::numerics::Status::conditioning_budget_exceeded,"thin witness typed conditioning");};
 witness(std::ldexp(1.L,-124));
 edge_source={s/2,s};edge_band={s,2*s};in.redshift=std::ldexp(1.,-65);witness(std::ldexp(1.L,-85));
 std::cout<<"PASS "<<checks<<" sampled photometry peer controls\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
