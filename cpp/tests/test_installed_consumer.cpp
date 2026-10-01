// Standalone installed-library linkage contract; no survey/science qualification.
#include <irred/numerics.hpp>
#include <irred/gaussian_design.hpp>
#include <irred/sampled_photometry.hpp>
#include <irred/calibration_ladder.hpp>
#include <array>
#include <cmath>
#include <numbers>
int installed_ladder() {
 using namespace irred::calibration;
 using namespace irred::statistics;
 Model m;
 m.ordered_host_ids={"A","B"};
 m.magnitude_convention="synthetic common magnitude convention";
 m.metallicity_coordinate_identity="synthetic metallicity dex";
 m.distance_shape_identity="supplied fixed shape at Href70";
 m.calibration_identity="one synthetic shared calibration measurement";
 m.dependence_identity="synthetic independent conditional rows";
 m.conditional_covariance_identity="synthetic .01 diagonal mag squared conditional on delta";
 for(const auto *host:{"A","B"}) m.rows.push_back({RowKind::anchor_modulus,std::string("anchor/")+host,host,"",0,0,0,0});
 const std::array<double,3> periods{.6,1,1.5},metallicities{-.2,.3,-.1};
 for(const auto *host:{"A","B"}) for(unsigned i=0;i<3;++i) {
  const std::string id=std::string("cep/")+host+"/"+std::to_string(i);
  m.rows.push_back({RowKind::cepheid,id,host,id,0,periods[i],metallicities[i],0});
 }
 for(const auto *host:{"A","B"}) m.rows.push_back({RowKind::calibrator_supernova,std::string("cal/")+host,host,std::string("SN/")+host,1,0,0,0});
 m.rows.push_back({RowKind::hubble_supernova,"flow","","SN/flow",1,0,0,36});
 m.rows.push_back({RowKind::calibration_measurement,"delta","","",1,0,0,0});
 Metadata md;for(const auto &row:m.rows) md.ordered_ids.push_back(row.row_id);
 md.measure="product d(mag)";md.source_semantics="synthetic controls";
 md.table_identity="installed analytic ladder control";md.ordering_provenance="explicit synthetic order";
 md.calibration_provenance=m.calibration_identity;md.dependence_provenance=m.dependence_identity;
 md.uncertainty_identity=m.conditional_covariance_identity;
 std::vector<double> c(m.rows.size()*m.rows.size());for(unsigned i=0;i<m.rows.size();++i)c[i*m.rows.size()+i]=.01;
 auto g=prepare_gaussian(c,MatrixKind::covariance,md,c.size(),1e-10);
 auto ladder=Ladder::prepare(std::move(g),std::move(m));
 const std::array<double,12> y{31,32,27.16,26.06,24.48,28.16,27.06,25.48,12.03,13.03,16.53,.03};
 const auto recovered=ladder.fit(y,md.ordered_ids);
 const double h0=70*std::pow(10.,.1);
 return recovered.relative_fit.status==DensityStatus::finite && recovered.relative_fit.coefficients.size()==8 && std::abs(recovered.relative_fit.coefficients[5]+19)<1e-9 && std::abs(recovered.h0_km_s_Mpc-h0)<1e-7 ? 0 : 6;
}
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
 return installed_ladder();
}
