#include "irred/quantities.hpp"
#include "fixtures/foundations_oracles.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
static unsigned count=0;static void expect(bool b,const char* s){++count;if(!b)throw std::runtime_error(s);}
static Quantity q(double value,Unit unit,Role role=Role::physical_length){return {value,unit,role,Frame::none,role==Role::comoving_distance?LengthConvention::comoving_a0_one:LengthConvention::physical};}
static Target t(Unit unit,Role role=Role::physical_length){return {unit,role,Frame::none,role==Role::comoving_distance?LengthConvention::comoving_a0_one:LengthConvention::physical};}
static double oracle(const char* id){for(auto& o:test_fixtures::foundations_oracles)if(o.id==id)return o.rounded;throw std::runtime_error("missing frozen independent oracle");}
static void near_ulps(double x,double expected,unsigned ulps,const char* label){double lo=expected,hi=expected;for(unsigned i=0;i<ulps;++i){lo=std::nextafter(lo,-INFINITY);hi=std::nextafter(hi,INFINITY);}expect(x>=lo&&x<=hi,label);}
int main(){try{
 // Exact dyadic scaling: multiplication by 1000=125*2^3 preserves these <=10-bit significands.
 for(int e:{-1000,-100,-1,0,500,1000}) {const double x=std::ldexp(1.,e);auto c=convert(q(x,Unit::kilometre),t(Unit::metre));expect(c.status==QuantityStatus::ok,"dyadic scaling finite");expect(c.target.value==std::ldexp(125.,e+3),"independent dyadic result");}
 // At x=1+2^-52 the exact scaled excess is 1000/512=1.953125 ULPs at 1000.
 const double x=std::nextafter(1.,2.);const double rounded=std::nextafter(std::nextafter(1000.,INFINITY),INFINITY);
 auto c=convert(q(x,Unit::kilometre),t(Unit::metre));expect(c.status==QuantityStatus::ok&&c.target.value==rounded,"independent ties-distance rounding");
 near_ulps(parsec_in_metres(),oracle("parsec_in_metres"),1,"Machin parsec one-ULP envelope");
 c=convert(q(1,Unit::megaparsec),t(Unit::metre));near_ulps(c.target.value,oracle("megaparsec_in_metres"),1,"Mpc independent one-ULP envelope");
 c=convert({70,Unit::km_per_s_per_mpc,Role::expansion_rate},{Unit::inverse_second,Role::expansion_rate});near_ulps(c.target.value,oracle("H0_70_inverse_seconds"),2,"inverse-time independent two-ULP envelope");
 const std::array<Role,4> distance_roles{Role::physical_length,Role::comoving_distance,Role::luminosity_distance,Role::angular_diameter_distance};
 for(auto a:distance_roles)for(auto b:distance_roles)if(a!=b)expect(convert(q(1,Unit::metre,a),t(Unit::metre,b)).status==QuantityStatus::role_mismatch,"same dimensions never transform role");
 c=convert(q(std::numeric_limits<double>::max(),Unit::parsec),t(Unit::megaparsec));expect(c.status==QuantityStatus::ok,"ratio-first huge source remains finite");near_ulps(c.target.value,std::numeric_limits<double>::max()/1000000.,2,"huge parsec to Mpc independent scale");
 for(double v:{-0.,0.,std::numeric_limits<double>::min(),std::numeric_limits<double>::max()}) {auto input=q(v,Unit::metre);c=convert(input,t(Unit::metre));expect(c.status==QuantityStatus::ok&&c.target.value==v,"normal/zero identity");expect(std::signbit(c.target.value)==std::signbit(v),"signed zero preserved");}
 auto source=q(17,Unit::kilometre);source.constants="old-unqualified-set";c=convert(source,t(Unit::metre));expect(c.status==QuantityStatus::constant_set_mismatch,"stale identity rejected");expect(c.source.value==17&&c.source.unit==source.unit&&c.source.role==source.role&&c.source.frame==source.frame&&c.source.convention==source.convention&&c.source.constants==source.constants,"failed source metadata untouched");
 for(auto frame:{Frame::heliocentric,Frame::cmb,Frame::model}) {Quantity z{std::nextafter(-1.,0.),Unit::one,Role::redshift,frame};c=convert(z,{Unit::one,Role::redshift,frame});expect(c.status==QuantityStatus::ok&&c.target.value==z.value,"valid redshift immediate boundary");expect(c.source.frame==frame&&c.target.frame==frame,"redshift frame preserved");}
 expect(convert({.1,Unit::one,Role::redshift,static_cast<Frame>(99)},{Unit::one,Role::redshift,static_cast<Frame>(99)}).status==QuantityStatus::missing_convention,"unknown frame rejected even when equal");
 expect(convert(q(1,Unit::metre),{Unit::metre,Role::physical_length,Frame::none,static_cast<LengthConvention>(99)}).status==QuantityStatus::unsupported_transform,"unknown convention cannot transform");
 expect(convert({1,Unit::one,Role::ratio,Frame::none},{Unit::one,Role::ratio,Frame::none,LengthConvention::physical}).status==QuantityStatus::unsupported_transform,"dimensionless cannot inherit length convention");
 std::array<Quantity,3> input{q(-2,Unit::kilometre),q(INFINITY,Unit::kilometre),q(4,Unit::kilometre)};std::array<Conversion,3> output{};expect(convert_batch(input,t(Unit::metre),output)==QuantityStatus::ok,"mixed batch execution completed");expect(output[0].status==QuantityStatus::ok&&output[0].target.value==-2000&&output[1].status==QuantityStatus::nonfinite_input&&output[2].target.value==4000,"mixed per-row failure/order retained");
 std::printf("{\"suite\":\"independent_quantity_hostile\",\"checks\":%u,\"passed\":true,\"scope\":\"bounded F01 unit scaling only\"}\n",count);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
