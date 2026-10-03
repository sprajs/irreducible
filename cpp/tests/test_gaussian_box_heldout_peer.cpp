// Independent direct 2-D finite-rectangle quadrature of the defining density.
// This synthetic algorithm control is not the released experiment reference.
#include "box_heldout_controls.hpp"
#include <cmath>
#include <iostream>
using namespace irred::statistics;
namespace {
void check(bool v,const char*m){if(!v)throw std::runtime_error(m);}
long double integral(unsigned n,bool joint,long double heldout) {
  constexpr long double lo0=-8,hi0=8,lo1=-9,hi1=9,s=11.L/16;
  const auto h0=(hi0-lo0)/n,h1=(hi1-lo1)/n;long double sum=0;
  // Here rT=(1/4,-1/2), xt=(3/4,-1/8), k=(1/2,1/4),
  // a=(1/4,-3/8), b=(v-3)-k^T*rT=(v-3).
  for(unsigned i=0;i<=n;++i) {
    const auto beta0=lo0+i*h0;const unsigned w0=i==0||i==n?1:i%2?4:2;
    for(unsigned j=0;j<=n;++j) {
      const auto beta1=lo1+j*h1;const unsigned w1=j==0||j==n?1:j%2?4:2;
      auto q=(.25L-beta0)*(.25L-beta0)+(-.5L-beta1)*(-.5L-beta1);
      if(joint){const auto z=(heldout-3)-(.25L*beta0-.375L*beta1);q+=z*z/s;}
      sum+=w0*w1*std::exp(-q/2);
    }
  }
  return sum*h0*h1/9;
}
long double reference(unsigned n,long double value) {
  return std::log(integral(n,true,value)/integral(n,false,value))-
      .5L*std::log(2*std::acos(-1.L)*(11.L/16));
}
}
int main(){try {
  box_heldout_controls::Input in;in.x[4]=.75;in.x[5]=-.125;
  auto g=in.gaussian();auto p=in.policy(g,1,2,2);auto o=in.prepare(std::move(g),p);
  check(o.status()==DensityStatus::finite,"nonzero conditional response prepared");
  const std::vector<double>y={1.25,1.5},v={3,3.375};const std::vector<BoxHeldoutRequest>requests={{0,0},{0,1}};
  const auto out=o.evaluate(y,o.training_row_ids(),v,requests,p);
  for(std::size_t i=0;i<v.size();++i) {
    const auto coarse=reference(128,v[i]),fine=reference(256,v[i]);
    check(std::abs(coarse-fine)<1e-11L,"direct finite-box quadrature refinement");
    check(out.densities[i].density_available,"non-Gaussian finite-box mixture admitted");
    const auto interval=out.densities[i].log_density;
    check(std::max(std::abs(static_cast<long double>(interval.lower)-fine),
        std::abs(static_cast<long double>(interval.upper)-fine))<3e-9L,"both native endpoints agree with independent rectangle algorithm");
  }
  // Change source orientation k while retaining C symmetry; the scalar law
  // changes, so an accidental independence substitution cannot pass this pair.
  in.covariance[2]=in.covariance[6]=-.5;auto other=in.gaussian();auto op=in.policy(other);auto oo=in.prepare(std::move(other),op);
  const std::vector<double>vv={3};const std::vector<BoxHeldoutRequest>rr={{0,0}};
  const auto different=oo.evaluate(y,oo.training_row_ids(),vv,rr,op);
  check(different.densities[0].density_available,"reversed correlation still valid");
  check(std::abs(different.densities[0].log_density.lower-out.densities[0].log_density.lower)>1e-3,"source covariance sign retained");
  std::cout<<"heldout independent synthetic rectangle controls complete\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
