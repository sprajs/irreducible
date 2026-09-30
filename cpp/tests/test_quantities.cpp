#include "cosmology/quantities.hpp"
#include "fixtures/foundations_oracles.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
using namespace cosmology;
namespace {
int checks=0;double maximum_relative=0;
void check(bool ok,const char* name){++checks;if(!ok)throw std::runtime_error(name);}
double oracle(std::string_view id){for(auto v:test_fixtures::foundations_oracles)if(v.id==id)return v.rounded;throw std::runtime_error("missing oracle");}
void near(double actual,double expected,const char* name){const double r=std::abs((actual-expected)/expected);maximum_relative=std::max(maximum_relative,r);check(r<=1e-13,name);}
Quantity length(double v,Unit u,Role role=Role::physical_length){return {v,u,role,Frame::none,role==Role::comoving_distance?LengthConvention::comoving_a0_one:LengthConvention::physical};}
Target target(Unit u,Role role=Role::physical_length){return {u,role,Frame::none,role==Role::comoving_distance?LengthConvention::comoving_a0_one:LengthConvention::physical};}
}
int main(){try{
 static_assert(!std::is_same_v<Length,Duration>);
 check(speed_of_light_m_per_s==299792458,"exact c");
 near(parsec_in_metres(),oracle("parsec_in_metres"),"independent Machin parsec");
 auto c=convert(length(1,Unit::megaparsec),target(Unit::metre));
 check(c.status==QuantityStatus::ok,"Mpc convert");near(c.target.value,oracle("megaparsec_in_metres"),"Mpc oracle");
 check(c.source.value==1&&c.source.unit==Unit::megaparsec&&c.source.constants==constant_set_id,"source retained");
 c=convert({70,Unit::km_per_s_per_mpc,Role::expansion_rate},{Unit::inverse_second,Role::expansion_rate});
 check(c.status==QuantityStatus::ok,"H0 conversion");near(c.target.value,oracle("H0_70_inverse_seconds"),"inverse time oracle");
 // Independently defined exact decimal scaling, including sign and zero.
 for(double x:{-1e12,-1e-12,0.,1e-12,1.,1e12}){
  c=convert(length(x,Unit::kilometre),target(Unit::metre));check(c.status==QuantityStatus::ok,"finite signed length");
  check(c.target.value==x*1000,"exact km definition");
 }
 for(double x:{1e-12,1.,1e12}){
  c=convert(length(x,Unit::megaparsec),target(Unit::metre));check(c.status==QuantityStatus::ok,"qualified Mpc envelope");
  near(c.target.value,x*oracle("megaparsec_in_metres"),"Mpc envelope independent value");
 }
 c=convert({1,Unit::day,Role::duration},{Unit::second,Role::duration});check(c.status==QuantityStatus::ok&&c.target.value==86400,"exact day scale");
 check(convert(length(0,Unit::metre,Role::luminosity_distance),target(Unit::kilometre,Role::luminosity_distance)).status==QuantityStatus::ok,"zero distance accepted");
 check(convert(length(-1,Unit::metre,Role::luminosity_distance),target(Unit::metre,Role::luminosity_distance)).status==QuantityStatus::invalid_domain,"negative distance rejected");
 check(convert({-1,Unit::one,Role::redshift,Frame::cmb},{Unit::one,Role::redshift,Frame::cmb}).status==QuantityStatus::invalid_domain,"redshift boundary");
 check(convert({0,Unit::inverse_second,Role::expansion_rate},{Unit::inverse_second,Role::expansion_rate}).status==QuantityStatus::invalid_domain,"zero expansion rejected");
 check(convert({.2,Unit::one,Role::redshift},{Unit::one,Role::redshift}).status==QuantityStatus::missing_convention,"missing frame");
 check(convert({.2,Unit::one,Role::redshift,Frame::cmb},{Unit::one,Role::redshift,Frame::heliocentric}).status==QuantityStatus::unsupported_transform,"frame transform rejected");
 check(convert({.2,Unit::one,Role::redshift,Frame::cmb},{Unit::one,Role::ratio,Frame::cmb}).status==QuantityStatus::role_mismatch,"redshift ratio mismatch");
 check(convert(length(1,Unit::metre),target(Unit::metre,Role::luminosity_distance)).status==QuantityStatus::role_mismatch,"distance role mismatch");
 check(convert({1,Unit::metre,Role::physical_length},target(Unit::metre)).status==QuantityStatus::unsupported_transform,"absent length convention");
 check(convert(length(1,Unit::metre),{Unit::second,Role::physical_length}).status==QuantityStatus::dimension_mismatch,"dimension mismatch");
 check(convert({1,static_cast<Unit>(99),Role::ratio},{Unit::one,Role::ratio}).status==QuantityStatus::unknown_unit,"unknown unit");
 auto bad=length(1,Unit::metre);bad.constants="future-constants";
 check(convert(bad,target(Unit::metre)).status==QuantityStatus::constant_set_mismatch,"stale constant identity");
 for(double x:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()})check(convert(length(x,Unit::metre),target(Unit::metre)).status==QuantityStatus::nonfinite_input,"nonfinite rejected");
 check(convert(length(std::numeric_limits<double>::max(),Unit::megaparsec),target(Unit::metre)).status==QuantityStatus::overflow,"overflow");
 check(convert(length(std::numeric_limits<double>::denorm_min(),Unit::metre),target(Unit::kilometre)).status==QuantityStatus::underflow,"subnormal envelope");
 check(convert(length(std::numeric_limits<double>::denorm_min(),Unit::metre),target(Unit::metre)).status==QuantityStatus::underflow,"same-unit subnormal envelope");
 check(convert(length(std::numeric_limits<double>::min(),Unit::metre),target(Unit::metre)).status==QuantityStatus::ok,"smallest normal identity");
 check(convert(length(std::numeric_limits<double>::min(),Unit::kilometre),target(Unit::metre)).status==QuantityStatus::ok,"smallest normal scaled");
 check(convert(length(std::numeric_limits<double>::min(),Unit::metre),target(Unit::kilometre)).status==QuantityStatus::underflow,"normal input subnormal result");
 check(convert(length(std::numeric_limits<double>::max(),Unit::metre),target(Unit::metre)).status==QuantityStatus::ok,"largest finite identity");
 check(convert(length(std::numeric_limits<double>::max(),Unit::metre),target(Unit::kilometre)).status==QuantityStatus::ok,"avoid intermediate overflow");
 std::array<Quantity,3> src{length(1,Unit::kilometre),length(-1,Unit::kilometre),length(std::numeric_limits<double>::infinity(),Unit::kilometre)};
 std::array<Conversion,3> out{};
 check(convert_batch(src,target(Unit::metre),out)==QuantityStatus::ok,"batch completed");
 check(out[0].target.value==1000&&out[1].target.value==-1000&&out[2].status==QuantityStatus::nonfinite_input,"per-row failure retained");
 check(convert_batch({},target(Unit::metre),{})==QuantityStatus::ok,"empty batch");
 check(convert_batch(src,target(Unit::metre),{})==QuantityStatus::invalid_batch,"batch length mismatch");
 std::cout<<"quantity checks="<<checks<<" max_relative_oracle_difference="<<maximum_relative<<'\n';return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
