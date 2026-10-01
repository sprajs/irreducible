// Standalone installed-library linkage contract; no survey/science qualification.
#include <irred/numerics.hpp>
#include <irred/gaussian_design.hpp>
#include <array>
#include <cmath>
int main() {
 const auto result=irred::numerics::log1p_checked(0.5);
 if(result.status!=irred::numerics::Status::ok || std::abs(result.value-0.4054651081081643819780131154643491)>=1e-14) return 1;
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
 return fitted.status==DensityStatus::finite && fitted.coefficients.size()==2 && std::abs(fitted.coefficients[0]-1)<2e-12 && std::abs(fitted.coefficients[1]+.5)<2e-12 && std::abs(fitted.quadratic-.25)<1e-10 ? 0 : 2;
}
