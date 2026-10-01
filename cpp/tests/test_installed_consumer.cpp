// Standalone installed-library linkage contract; no survey/science qualification.
#include <irred/numerics.hpp>
#include <irred/gaussian_design.hpp>
#include <irred/sampled_photometry.hpp>
#include <irred/early_late.hpp>
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
 if (!(sampled.admission_status==irred::numerics::Status::ok && sampled.flux_watt_per_square_metre.value && sampled.energy_joule.value && std::abs(*sampled.flux_watt_per_square_metre.value/expected-1)<2e-12 && *sampled.energy_joule.value==*sampled.flux_watt_per_square_metre.value)) return 2;
 using namespace irred::statistics;
 const std::array<double,9> covariance{2,1,0,1,2,1,0,1,2};
 const std::array<double,6> design{1,0,1,1,1,2};
 const std::array<double,3> residual{1,0,0};
 Metadata source_metadata;
 source_metadata.ordered_ids={"anchor-a","anchor-b","host"};
 source_metadata.measure="product d(magnitude)";
 source_metadata.ordering_provenance="explicit synthetic order";
 auto gaussian=prepare_gaussian(covariance,MatrixKind::covariance,source_metadata,9,1e-10);
 DesignMetadata metadata;
 metadata.ordered_parameter_ids={"offset","slope"};
 metadata.parameter_units={"magnitude","magnitude"};
 metadata.shared_nuisance_ids={"offset"};
 metadata.residual_unit="magnitude";
 metadata.design_identity="installed synthetic linear control";
 metadata.dependence_identity="declared correlated synthetic covariance";
 auto profile=DesignProfile::prepare(std::move(gaussian),design,source_metadata.ordered_ids,std::move(metadata));
 const auto fitted=profile.evaluate(residual,source_metadata.ordered_ids);
 if (!(fitted.status==DensityStatus::finite && fitted.coefficients.size()==2 && std::abs(fitted.coefficients[0]-1)<2e-12 && std::abs(fitted.coefficients[1]+.5)<2e-12 && std::abs(fitted.quadratic-.25)<1e-10)) return 3;
 using namespace irred::cosmology;
 const SoundHorizonRequest radiation{{70,0,1,0,.1},1059,"installed synthetic supplied drag"};
 const std::array<double,1> redshifts{1};
 const EarlyLatePolicy policy{1e-9,2e-11,1e-11,5e-11,100000,200000,40,4,4*1024*1024,{1e-9,2e-11,100000,40,4,200000,4*1024*1024}};
 const auto geometry=evaluate_early_late(radiation,redshifts,early_late_mask(EarlyLateOutput::dm_mpc)|early_late_mask(EarlyLateOutput::dl_mpc),policy);
 if (geometry.status!=irred::numerics::Status::ok || geometry.rows.size()!=1 || geometry.ruler) return 4;
 const auto &dm=geometry.rows[0].outputs[static_cast<unsigned>(EarlyLateOutput::dm_mpc)], &dl=geometry.rows[0].outputs[static_cast<unsigned>(EarlyLateOutput::dl_mpc)];
 const double exact_dm=299792.458/70/2;
 return dm.status==irred::numerics::Status::ok && dl.status==irred::numerics::Status::ok && dm.value && dl.value && std::abs(*dm.value-exact_dm)<1e-7 && std::abs(*dl.value-2*exact_dm)<1e-7 ? 0 : 5;
}
