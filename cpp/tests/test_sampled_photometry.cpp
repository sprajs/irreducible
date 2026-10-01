// Synthetic exact polynomial controls; rational integration is independent of
// production Bernstein integration. SI/pi ancestry is shared and disclosed.
#include "irred/sampled_photometry.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred;
namespace {
unsigned checks=0;
void check(bool v,const char* why){++checks;if(!v)throw std::runtime_error(why);}
void near(const photometry::Outcome &o,long double ref){
 check(o.availability==photometry::Availability::available && o.value.has_value(),"available");
 check(o.numerical_status==numerics::Status::ok,"numerical ok");
 check(ref<=0 || *o.value>0,"positive control never zero");
 check(std::abs(*o.value-ref)<=2e-12L*std::abs(ref)+1e-300L,"named budget");
}
photometry::SampledResult run(const photometry::SampledInput &p,unsigned mask=7){return photometry::evaluate_sampled(p,{mask,65536,65536});}
}
int main(){try{
 std::array<double,2> wave{1,2},lum{1,1},trans{0.5,0.5};
 photometry::SampledInput p{{wave,lum},{wave,trans},1,0,2,3};
 const long double scale=1/(4*std::numbers::pi_v<long double>),hc=6.62607015e-34L*299792458;
 auto v=run(p);check(v.admission_status==numerics::Status::ok,"admission");
 near(v.flux_watt_per_square_metre,scale);near(v.energy_joule,3*scale);near(v.expected_photons,4.5L*scale/hc);
 lum={1,2};trans={0,1};v=run(p);
 // L=lambda, T=lambda-1 on [1,2]: integral L=3/2,
 // integral L*T=5/6, integral lambda*L*T=17/12.
 near(v.flux_watt_per_square_metre,1.5L*scale);near(v.energy_joule,5*scale);near(v.expected_photons,8.5L*scale/hc);
 // Clip the linear source to a finite band with independent polynomial
 // antiderivatives, keeping transmission explicitly distinct from support.
 std::array<double,2> clipped_wave{1.25,1.75},clipped_trans{0.25,0.75};
 p.passband={clipped_wave,clipped_trans};auto clipped=run(p);
 const long double a=1.25L,b=1.75L;
 const auto energy_antiderivative=[](long double x){return x*x*x/3-x*x/2;};
 const auto photon_antiderivative=[](long double x){return x*x*x*x/4-x*x*x/3;};
 near(clipped.flux_watt_per_square_metre,(b*b-a*a)/2*scale);
 near(clipped.energy_joule,6*(energy_antiderivative(b)-energy_antiderivative(a))*scale);
 near(clipped.expected_photons,6*(photon_antiderivative(b)-photon_antiderivative(a))*scale/hc);
 p.passband={wave,trans};
 std::array<double,3> extra_wave{1,1.5,2},extra_lum{1,1.5,2},extra_trans{0,0.5,1};
 p.spectrum={extra_wave,extra_lum};p.passband={extra_wave,extra_trans};auto refined=run(p);
 near(refined.energy_joule,*v.energy_joule.value);near(refined.expected_photons,*v.expected_photons.value);
 p.spectrum={wave,lum};p.passband={wave,trans};
 for(unsigned field=0;field<4;++field){auto q=p;
  if(field==0)q.luminosity_distance_metre=2;
  if(field==1)q.collecting_area_square_metre=4;
  if(field==2)q.observer_exposure_second=6;
  if(field==3)q.redshift=1;
  auto x=run(q);if(field==3){near(x.energy_joule,0);continue;}
  near(x.flux_watt_per_square_metre,field==0?0.375L*scale:1.5L*scale);
  near(x.energy_joule,(field==0?1.25L:10)*scale);
 }
 // Bolometric support doubles at z=1 while spectral scale halves;
 // observer photon wavelengths double with observer exposure held fixed.
 std::array<double,2> doubled{2,4};p.redshift=1;p.passband={doubled,trans};v=run(p);
 near(v.flux_watt_per_square_metre,1.5L*scale);near(v.energy_joule,5*scale);near(v.expected_photons,17*scale/hc);
 p.redshift=0;p.passband={wave,trans};p.collecting_area_square_metre=0;v=run(p);
 near(v.energy_joule,0);near(v.expected_photons,0);near(v.flux_watt_per_square_metre,1.5L*scale);
 p.collecting_area_square_metre=2;trans={0,0};v=run(p);near(v.energy_joule,0);near(v.flux_watt_per_square_metre,1.5L*scale);
 p.luminosity_distance_metre=std::numeric_limits<double>::quiet_NaN();v=run(p);check(v.admission_status==numerics::Status::nonfinite_input&&!v.energy_joule.value,"validate before zero");p.luminosity_distance_metre=1;
 trans={0.5,0.5};lum={1,1};
 for(unsigned bad=0;bad<6;++bad){auto q=p;auto w=wave,y=lum,t=trans;
  if(bad==0)w[1]=w[0];
  if(bad==1)w[0]=0;
  if(bad==2)y[0]=-1;
  if(bad==3)t[0]=1.01;
  if(bad==4)y[1]=INFINITY;
  q.spectrum={w,y};q.passband={wave,t};if(bad==5)q.spectrum.luminosity_watt_per_metre=std::span(y).first(1);
  auto x=run(q);check(x.admission_status!=numerics::Status::ok&&!x.energy_joule.value,"invalid arrays");
 }
 v=photometry::evaluate_sampled(p,{7,3,2});check(v.admission_status==numerics::Status::work_limit&&!v.energy_joule.value,"knots quota");
 v=photometry::evaluate_sampled(p,{7,4,0});check(v.admission_status==numerics::Status::work_limit&&!v.energy_joule.value,"segments quota");
 v=run(p,2);near(v.energy_joule,3*scale);check(v.expected_photons.availability==photometry::Availability::omitted&&!v.expected_photons.value,"omitted photons");
 // Exact binary redshift endpoint adversaries from rectangular contract.
 const double s=std::ldexp(1.0,-20),e=std::ldexp(1.0,-52);
 wave={s,s*(1+e)};std::array<double,2> band{s*(1+2*e),std::nextafter(s*(1+2*e),INFINITY)};
 p.passband={band,trans};p.redshift=e;v=run(p);
 check(v.energy_joule.availability==photometry::Availability::failed && v.energy_joule.numerical_status==numerics::Status::conditioning_budget_exceeded,"tiny overlap not zero");
 wave={std::nextafter(s,0.0),s};band={s,std::nextafter(s,INFINITY)};p.redshift=std::ldexp(1.0,-65);v=run(p);
 check(v.energy_joule.availability==photometry::Availability::failed,"redshift sub-wide overlap");
 p.redshift=1;band={2*s,3*s};v=run(p);near(v.energy_joule,0);
 // Adjacent binary64 band endpoints inside exact source support remain positive.
 wave={1,2};band={1.5,std::nextafter(1.5,INFINITY)};p.redshift=0;v=run(p);near(v.energy_joule,3*scale*(static_cast<long double>(band[1])-band[0]));
 lum={std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::denorm_min()};p.luminosity_distance_metre=1e100;v=run(p);
 check(v.flux_watt_per_square_metre.availability==photometry::Availability::failed&&!v.flux_watt_per_square_metre.value,"positive underflow");
 wave={1e100,2e100};band={1.2e100,std::nextafter(1.2e100,INFINITY)};lum={1e205,1e205};p.luminosity_distance_metre=1e3;v=run(p);
 check(v.energy_joule.availability==photometry::Availability::available,"energy isolated");
 check(v.expected_photons.availability==photometry::Availability::failed&&v.expected_photons.numerical_status==numerics::Status::overflow,"photon overflow");
 auto only=run(p,2);near(only.energy_joule,*v.energy_joule.value);
 std::cout<<"PASS "<<checks<<" sampled photometry owner controls\n";
}catch(const std::exception&e){std::cerr<<"FAIL check "<<checks<<": "<<e.what()<<'\n';return 1;}}
