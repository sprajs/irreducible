// Standalone installed-library linkage contract; no survey/science qualification.
#include <irred/numerics.hpp>
#include <irred/sampled_photometry.hpp>
#include <array>
#include <cmath>
#include <numbers>
int main() {
 const auto result=irred::numerics::log1p_checked(0.5);
 if(result.status!=irred::numerics::Status::ok || std::abs(result.value-0.4054651081081643819780131154643491)>=1e-14) return 1;
 const std::array<double,2> wavelength{1e-6,2e-6},luminosity{1,1},transmission{1,1};
 const irred::photometry::SampledInput input{{wavelength,luminosity},{wavelength,transmission},1,0,1,1};
 const auto sampled=irred::photometry::evaluate_sampled(input,{3,4,1});
 const double expected=1e-6/(4*std::numbers::pi);
 return sampled.admission_status==irred::numerics::Status::ok && sampled.flux_watt_per_square_metre.value && sampled.energy_joule.value && std::abs(*sampled.flux_watt_per_square_metre.value/expected-1)<2e-12 && *sampled.energy_joule.value==*sampled.flux_watt_per_square_metre.value ? 0 : 2;
}
