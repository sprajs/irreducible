#include "cosmology/quantities.hpp"
#include <cmath>
#include <limits>
#include <numbers>
namespace cosmology {
namespace {
struct UnitDefinition { bool known; Dimension dimension; long double scale; };
UnitDefinition definition(Unit u) noexcept {
 const long double pc=648000.L*149597870700.L/std::numbers::pi_v<long double>;
 switch(u) {
 case Unit::one:return {true,Dimension::dimensionless,1};
 case Unit::metre:return {true,Dimension::length,1};
 case Unit::kilometre:return {true,Dimension::length,1000};
 case Unit::second:return {true,Dimension::time,1};
 case Unit::day:return {true,Dimension::time,86400};
 case Unit::inverse_second:return {true,Dimension::inverse_time,1};
 case Unit::parsec:return {true,Dimension::length,pc};
 case Unit::megaparsec:return {true,Dimension::length,pc*1000000};
 case Unit::km_per_s_per_mpc:return {true,Dimension::inverse_time,1000/(pc*1000000)};
 default:return {false,Dimension::dimensionless,0};
 }
}
bool role_dimension(Role r, Dimension& d) noexcept {
 switch(r) {
 case Role::ratio:case Role::redshift:d=Dimension::dimensionless;return true;
 case Role::physical_length:case Role::comoving_distance:case Role::luminosity_distance:case Role::angular_diameter_distance:d=Dimension::length;return true;
 case Role::duration:d=Dimension::time;return true;
 case Role::expansion_rate:d=Dimension::inverse_time;return true;
 default:return false;
 }
}
bool valid_frame(Frame f) noexcept {return f==Frame::heliocentric||f==Frame::cmb||f==Frame::model;}
QuantityStatus validate(Quantity q, UnitDefinition u) noexcept {
 Dimension d{};
 if(!role_dimension(q.role,d)||d!=u.dimension)return QuantityStatus::dimension_mismatch;
 if(q.constants!=constant_set_id)return QuantityStatus::constant_set_mismatch;
 if(q.role==Role::redshift) {
  if(!valid_frame(q.frame))return QuantityStatus::missing_convention;
  if(q.value<=-1)return QuantityStatus::invalid_domain;
 } else if(q.frame!=Frame::none)return QuantityStatus::unsupported_transform;
 const bool length=d==Dimension::length;
 if(length) {
  const auto required=q.role==Role::comoving_distance?LengthConvention::comoving_a0_one:LengthConvention::physical;
  if(q.convention!=required)return QuantityStatus::missing_convention;
 } else if(q.convention!=LengthConvention::none)return QuantityStatus::unsupported_transform;
 if(q.role==Role::comoving_distance||q.role==Role::luminosity_distance||q.role==Role::angular_diameter_distance) {
  if(q.value<0)return QuantityStatus::invalid_domain;
 }
 if(q.role==Role::expansion_rate && q.value<=0)return QuantityStatus::invalid_domain;
 return QuantityStatus::ok;
}
}
double parsec_in_metres() noexcept {return static_cast<double>(definition(Unit::parsec).scale);}
Conversion convert(Quantity q, Target t) noexcept {
 Conversion out;out.source=q;out.target={0,t.unit,t.role,t.frame,t.convention,t.constants};
 auto fail=[&](QuantityStatus s){out.status=s;return out;};
 if(!std::isfinite(q.value))return fail(QuantityStatus::nonfinite_input);
 const auto a=definition(q.unit),b=definition(t.unit);
 if(!a.known||!b.known)return fail(QuantityStatus::unknown_unit);
 if(a.dimension!=b.dimension)return fail(QuantityStatus::dimension_mismatch);
 if(q.role!=t.role)return fail(QuantityStatus::role_mismatch);
 if(q.frame!=t.frame||q.convention!=t.convention)return fail(QuantityStatus::unsupported_transform);
 if(q.constants!=t.constants)return fail(QuantityStatus::constant_set_mismatch);
 if(auto s=validate(q,a);s!=QuantityStatus::ok)return fail(s);
 if(auto s=validate({q.value,t.unit,t.role,t.frame,t.convention,t.constants},b);s!=QuantityStatus::ok)return fail(s);
 // Ratio-first evaluation avoids overflowing an otherwise representable conversion.
 const long double result=static_cast<long double>(q.value)*(a.scale/b.scale);
 if(!std::isfinite(result)||std::abs(result)>std::numeric_limits<double>::max())return fail(QuantityStatus::overflow);
 if(result!=0 && std::abs(result)<std::numeric_limits<double>::min())return fail(QuantityStatus::underflow);
 out.target.value=static_cast<double>(result);out.status=QuantityStatus::ok;return out;
}
QuantityStatus convert_batch(std::span<const Quantity> source,Target target,std::span<Conversion> output) noexcept {
 if(source.size()!=output.size())return QuantityStatus::invalid_batch;
 // Contract: source and output storage are disjoint and remain live for the call.
 for(std::size_t i=0;i<source.size();++i)output[i]=convert(source[i],target);
 return QuantityStatus::ok;
}
} // namespace cosmology
