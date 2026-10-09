#include "irred/curved_flrw.hpp"
#include <stdexcept>
#include <cmath>
#include <limits>
#include <utility>
using namespace irred::cosmology;
using S=irred::numerics::Status;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
#define CHECK(expression) require((expression),#expression)
void close(double x,double y,double r=2e-8){CHECK(std::abs(x-y)<=r*std::max(1.,std::abs(y)));}
// Independent composite Simpson in ln(1+z), unlike production adaptive z mesh.
long double radial_reference(CurvedFLRWSpec s,double z,int n=32768){
  long double end=std::log1p(z),step=end/n,sum=0,k=1.L-s.omega_m-s.omega_r-s.omega_lambda;
  for(int j=0;j<=n;++j){long double u=std::exp(j*step),q=s.omega_r*u*u*u*u+s.omega_m*u*u*u+k*u*u+s.omega_lambda;sum+=(j==0||j==n?1:j%2?4:2)*u/std::sqrt(q);}
  return sum*step/3*299792.458L/s.h0_km_s_mpc;
}
int main(){
  const double zs[]={0,.001,.5,2,10};
  CurvedFLRW milne({70,0,0,0});auto m=milne.evaluate(zs);CHECK(m.status==S::ok);
  for(auto r:m.rows){close(r.e,1+r.redshift);close(r.radial_mpc,299792.458/70*std::log1p(r.redshift));close(r.luminosity_mpc,299792.458/70*(r.redshift+r.redshift*r.redshift/2));close(r.luminosity_mpc,r.angular_mpc*std::pow(1+r.redshift,2));}
  CurvedFLRW de_sitter({70,0,0,1});auto ds=de_sitter.evaluate(zs);CHECK(ds.status==S::ok);for(auto r:ds.rows){close(r.e,1);close(r.radial_mpc,299792.458/70*r.redshift);}
  CurvedFLRW dust({70,1,0,0});auto d=dust.evaluate(zs);CHECK(d.status==S::ok);
  for(auto r:d.rows)close(r.radial_mpc,2*299792.458/70*(1-1/std::sqrt(1+r.redshift)));
  CurvedFLRW radiation({70,0,1,0});auto rad=radiation.evaluate(zs);CHECK(rad.status==S::ok);for(auto r:rad.rows)close(r.radial_mpc,299792.458/70*r.redshift/(1+r.redshift));
  for(auto s:{CurvedFLRWSpec{67,.3,.0001,.6},CurvedFLRWSpec{67,.3,.0001,.8},CurvedFLRWSpec{67,.3,.0001,.6999}}){auto b=CurvedFLRW(s).evaluate(zs);CHECK(b.status==S::ok);for(auto r:b.rows){auto refined=radial_reference(s,r.redshift);CHECK(std::abs(refined-radial_reference(s,r.redshift,16384))<1e-10L*std::max(1.L,std::abs(refined)));close(r.radial_mpc,static_cast<double>(refined));CHECK(r.transverse_error_estimate_mpc>=0);}}
  // Flat limit continuity on both signs without catastrophic sqrt(Ok) division.
  double z=3;auto flat=CurvedFLRW({70,.3,0,.7}).evaluate({&z,1});for(double eps:{-1e-12,1e-12}){auto a=CurvedFLRW({70,.3,0,.7+eps}).evaluate({&z,1});CHECK(a.status==S::ok);close(a.rows[0].transverse_mpc,flat.rows[0].transverse_mpc,1e-10);}
  // Positive endpoint does not admit a path with an intervening forbidden band.
  double large=10;auto turn=CurvedFLRW({70,.1,0,1.5}).evaluate({&large,1});CHECK(turn.status==S::outside_domain);
  // Closed Lambda model hits H²=0 before z=1; staticH0=0 is unsupported.
  double beyond=5;CHECK(CurvedFLRW({70,.3,0,1.71345}).evaluate({&beyond,1}).status==S::outside_domain);
  double zz=1;CHECK(CurvedFLRW({70,0,0,2}).evaluate({&zz,1}).status==S::outside_domain);
  CHECK(CurvedFLRW({0,1,0,0}).status()==S::outside_domain);
  CHECK(CurvedFLRW({70,-.1,0,1}).status()==S::outside_domain);CHECK(CurvedFLRW({std::numeric_limits<double>::infinity(),1,0,0}).status()==S::nonfinite_input);
  const double invalid[]={-1,std::numeric_limits<double>::quiet_NaN(),10001};auto bad=milne.evaluate(invalid);CHECK(bad.rows.size()==3);for(auto r:bad.rows)CHECK(r.status!=S::ok);
  auto policy=CurvedFLRWPolicy{};policy.maximum_points=1;CHECK(milne.evaluate(zs,policy).status==S::work_limit);policy={};policy.maximum_native_bytes=1;CHECK(milne.evaluate(zs,policy).status==S::work_limit);
  policy={};policy.maximum_callbacks_per_point=3;auto fail=CurvedFLRW({70,.3,0,.7}).evaluate({&z,1},policy);CHECK(fail.status==S::work_limit);
  policy={};policy.relative_tolerance=1e-18;CHECK(milne.evaluate({&z,1},policy).status!=S::ok);
  CurvedFLRW copy=milne,moved=std::move(copy);CHECK(copy.status()==S::invalid_input);CHECK(moved.status()==S::ok);auto& alias=moved;moved=std::move(alias);CHECK(moved.status()==S::ok);
}
