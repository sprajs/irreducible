// Original Decimal110/150 direct-charge bisection and joint centered moments.
// Shared physical/SI/atomic facts; no native solver is a numerical oracle.
#include "irred/baryon_abundance_law.hpp"
#include "baryon_abundance_law_peer_facts.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
namespace c=irred::cosmology; namespace f=abundance_law_peer_facts;
using W=long double; using S=irred::numerics::Status;
unsigned checks=0; W maximum_covariance_fraction=0;
void need(bool b,const char* why){++checks;if(!b)throw std::runtime_error(why);}
template<class R,class M,class C>
void compare(const R&r,const M&mu,const C&cov,W mean_budget){
 need(r.status==S::ok&&r.moments,"all joint states admitted");
 const auto&m=*r.moments;const std::size_t n=mu.size();
 need(m.mean.size()==n&&m.covariance.size()==n*n&&r.axes.size()==n,"full ordered moments");
 for(std::size_t i=0;i<n;++i){
  need(std::abs(W(m.mean[i])/mu[i]-1)<=mean_budget,"independent mean comparison");
  need(m.covariance[i*n+i]>0,"resolved positive physical variance");
  for(std::size_t j=0;j<n;++j){
   W scale=std::sqrt(cov[i*n+i]*cov[j*n+j]);
   W fraction=std::abs(W(m.covariance[i*n+j])-cov[i*n+j])/(2e-8L*scale);
   maximum_covariance_fraction=std::max(maximum_covariance_fraction,fraction);
   need(fraction<=1,"independent full covariance comparison");
   need(m.covariance[i*n+j]==m.covariance[j*n+i],"one symmetric triangle");
   need(W(m.covariance_empirical_sensitivity[i*n+j])<=1e-8L*std::sqrt(W(m.covariance[i*n+i])*m.covariance[j*n+j]),"separate numerical covariance allocation");
  }
 }
}
c::BaryonAbundanceLaw owner(){
 const std::array<c::BaryonAbundanceLawRow,2> rows{{{"early-a",f::scale_factors[0]},{"early-b",f::scale_factors[1]}}};
 std::array<c::BaryonAbundanceLawState,3> states{};
 const std::array<const char*,3> ids{"state-a","state-b","state-c"};
 for(std::size_t s=0;s<3;++s){auto&x=f::states[s];states[s]={ids[s],f::weights[s],{x[0],x[1],x[2],x[3],"original joint synthetic physical-density and He4 fraction law","independently supplied neutral effective masses"},f::temperatures[s],"joint supplied synthetic temperatures"};}
 return c::prepare_baryon_abundance_law({rows,states,"original finite masses2/3/5","one joint state controls both epochs and species"});
}
}
int main(){try{
 auto op=owner();need(op.status()==S::ok,"retained joint source");
 auto d=op.evaluate_density();compare(d,f::density_mean,f::density_covariance,2e-15L);
 auto l=op.evaluate_equilibrium();compare(l,f::lte_mean,f::lte_covariance,5e-13L);
 need(d.attempts.size()==3&&l.attempts.size()==3&&l.work.solves==6,"all ordered attempts and cumulative solves");
 need(d.source==l.source&&d.source==op.source(),"one immutable source owner across consumers");
 need(d.moments->covariance[1]>0&&l.moments->covariance[1]<0,"joint density correlation and ionization anticorrelation retained");
 for(unsigned mask=1;mask<64;++mask){auto r=op.evaluate_equilibrium(mask);need(r.status==S::ok&&r.moments,"all independently requested mask moments");
  std::vector<W> mu,cov;std::vector<std::size_t> selected;
  for(std::size_t row=0;row<2;++row)for(unsigned bit=0;bit<6;++bit)if(mask&(1u<<bit)){selected.push_back(row*6+bit);mu.push_back(f::lte_mean[row*6+bit]);}
  for(auto i:selected)for(auto j:selected)cov.push_back(f::lte_covariance[i*12+j]);
  compare(r,mu,cov,5e-13L);
  need(selected.size()==r.axes.size(),"no unrequested axis");}
 auto held=d.source;op={};need(held&&held->states.size()==3&&held->states[2].temperature_kelvin[1]==19500,"returned source outlives preparation");
 std::cout<<"PASS "<<checks<<" independent abundance-law controls; covariance allocation fraction="<<double(maximum_covariance_fraction)<<'\n';
}catch(const std::exception&e){std::cerr<<"FAIL "<<checks<<": "<<e.what()<<'\n';return 1;}}
