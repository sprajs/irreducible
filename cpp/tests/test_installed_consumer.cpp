// Standalone installed-library linkage contract; no survey/science qualification.
#include <irred/numerics.hpp>
#include <irred/gaussian_design.hpp>
#include <irred/sampled_photometry.hpp>
#include <irred/photometry_calibration.hpp>
#include <irred/calibration_ladder.hpp>
#include <irred/early_late.hpp>
#include <irred/bao_conditional.hpp>
#include <irred/correlated_calibration.hpp>
#include <array>
#include <cmath>
#include <numbers>
int installed_conditional_bao() {
 using namespace irred;
 bao::DensityInput input;
 input.queries={{1,bao::Observable::hubble_over_ruler}};
 input.observed={2.5*std::sqrt(3.)};input.covariance={1};input.ordered_ids={"synthetic DH/rs"};
 input.role=bao::RowRole::synthetic_control;input.covariance_unit=bao::CovarianceUnit::dimensionless_ratio_squared;
 input.table_identity="installed radiation analytic ratio";input.covariance_identity="unit ratio covariance";
 input.ordering_provenance="one ordered row";input.calibration_provenance="synthetic";input.dependence_provenance="one full row";
 const auto source=bao::prepare_density(std::move(input),{1,1,4096,4*1024*1024,1e-10,numerics::Arithmetic::longdouble_cpu_v1});
 const cosmology::SoundHorizonRequest point{{70,0,1,0,1},9,"installed supplied synthetic drag"};
 bao::ConditionalDensityPolicy policy;
 policy.predictions={1e-9,2e-11,1e-11,5e-11,100000,200000,40,4,4*1024*1024,{1e-9,2e-11,100000,40,4,200000,4*1024*1024}};
 policy.maximum_models=1;policy.maximum_queries=1;policy.maximum_string_bytes=4096;
 policy.maximum_native_bytes=8*1024*1024;policy.maximum_total_callbacks=200000;
 policy.maximum_forward_sensitivity=1e-10;policy.requested=static_cast<std::uint32_t>(bao::Output::normalized_density)|static_cast<std::uint32_t>(bao::Output::predictions);
 const auto result=source.evaluate_conditional({&point,1},policy);
 return result.status==statistics::DensityStatus::finite && result.slots.size()==1 && result.slots[0].result &&
  result.slots[0].result->density.status==statistics::DensityStatus::finite &&
  std::abs(result.slots[0].result->density.log_value+.5*std::log(2*std::numbers::pi))<1e-8 ? 0 : 9;
}
int installed_passband_calibration() {
 using namespace irred::photometry;
 const std::array<double,2> wavelength{1e-6,2e-6}, luminosity{1,1}, low{.2,.2}, high{.8,.8};
 const std::array<CalibrationBand,1> bands{{{"synthetic band",wavelength,1,1}}};
 const std::array<std::span<const double>,1> low_arrays{low}, high_arrays{high};
 const std::array<CalibrationState,2> states{{{"low",1,low_arrays},{"high",1,high_arrays}}};
 const CalibrationInput input{{wavelength,luminosity},1,0,bands,states,
  "installed synthetic source","installed synthetic calibration","two equiprobable states","one shared state"};
 CalibrationPolicy policy;policy.requested_outputs=collected_energy;
 const auto result=evaluate_calibration(input,policy);
 const double energy=1e-6/(4*std::numbers::pi);
 return result.status==irred::numerics::Status::ok && result.moments && result.axes.size()==1 &&
  std::abs(result.moments->mean[0]/(.5*energy)-1)<2e-10 &&
  std::abs(result.moments->covariance[0]/(.09*energy*energy)-1)<2e-8 ? 0 : 8;
}
int installed_correlated_calibration() {
 using namespace irred::statistics;
 Metadata md;md.ordered_ids={"a","b"};md.measure="product d(magnitude)";
 md.table_identity="installed synthetic calibration control";md.ordering_provenance="explicit rows";
 CalibrationPrior p;p.ordered_parameter_ids={"zero","colour"};p.parameter_units={"magnitude","magnitude"};
 p.mean={.25,-.5};p.covariance={.25,.125,.125,.5};
 p.prior_identity="installed proper correlated prior";p.response_identity="identity response";
 p.residual_unit="magnitude";p.measure_identity=md.measure;p.calibration_identity="synthetic calibration";
 p.dependence_identity="independent conditional noise";p.noise_independence_declared=true;
 auto g=prepare_gaussian(std::array<double,4>{1,0,0,2},MatrixKind::covariance,md,100,1e-10);
 const auto op=CorrelatedCalibration::prepare(std::move(g),std::array<double,4>{1,0,0,1},md.ordered_ids,std::move(p));
 const auto v=op.evaluate(std::array<double,2>{.25,-.5},md.ordered_ids);
 const double expected=-.5*(std::log(3.109375)+2*std::log(2*std::numbers::pi));
 return v.density.status==DensityStatus::finite && std::abs(v.quadratic)<1e-14 && std::abs(v.density.log_value-expected)<2e-12 ? 0 : 7;
}
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
 using namespace irred::cosmology;
 const SoundHorizonRequest radiation{{70,0,1,0,.1},1059,"installed synthetic supplied drag"};
 const std::array<double,1> redshifts{1};
 const EarlyLatePolicy policy{1e-9,2e-11,1e-11,5e-11,100000,200000,40,4,4*1024*1024,{1e-9,2e-11,100000,40,4,200000,4*1024*1024}};
 const auto geometry=evaluate_early_late(radiation,redshifts,early_late_mask(EarlyLateOutput::dm_mpc)|early_late_mask(EarlyLateOutput::dl_mpc),policy);
 if (geometry.status!=irred::numerics::Status::ok || geometry.rows.size()!=1 || geometry.ruler) return 4;
 const auto &dm=geometry.rows[0].outputs[static_cast<unsigned>(EarlyLateOutput::dm_mpc)], &dl=geometry.rows[0].outputs[static_cast<unsigned>(EarlyLateOutput::dl_mpc)];
 const double exact_dm=299792.458/70/2;
 if (!(dm.status==irred::numerics::Status::ok && dl.status==irred::numerics::Status::ok && dm.value && dl.value && std::abs(*dm.value-exact_dm)<1e-7 && std::abs(*dl.value-2*exact_dm)<1e-7)) return 5;
 if (const auto ladder_status=installed_ladder();ladder_status!=0) return ladder_status;
 if (const auto calibration_status=installed_passband_calibration();calibration_status!=0) return calibration_status;
 if (const auto bao_status=installed_conditional_bao();bao_status!=0) return bao_status;
 return installed_correlated_calibration();
}
