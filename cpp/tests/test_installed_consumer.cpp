#include <irred/gaussian_simulation.hpp>
#include <irred/random.hpp>
// Standalone installed-library linkage contract; no survey/science qualification.
#include <irred/numerics.hpp>
#include <irred/gaussian_design.hpp>
#include <irred/gr_growth.hpp>
#include <irred/hydrogen_equilibrium.hpp>
#include <irred/hydrogen_helium_equilibrium.hpp>
#include <irred/baryon_abundance.hpp>
#include <irred/recombination_drag.hpp>
#include <irred/sis_thin_lens.hpp>
#include <irred/gaussian_posterior.hpp>
#include <irred/gaussian_predictive.hpp>
#include <irred/sampled_photometry.hpp>
#include <irred/temporal_photometry.hpp>
#include <irred/detector_selection.hpp>
#include <irred/photometry_calibration.hpp>
#include <irred/calibration_ladder.hpp>
#include <irred/calibration_predictive.hpp>
#include <irred/early_late.hpp>
#include <irred/bao_conditional.hpp>
#include <irred/correlated_calibration.hpp>
#include <irred/thermal_observables.hpp>
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
int installed_thermal_observables() {
 using namespace irred::cosmology;
 const ThermalObservableRequest request{{70,.0245,.1225,2.7255,1e-5,{}},1059.95,
  "installed synthetic supplied drag","installed explicit physical-density source"};
 const auto prepared=prepare_thermal_observables(request);
 const std::array<double,2> redshifts{0,1};
 const auto result=prepared.evaluate(redshifts,early_late_mask(EarlyLateOutput::e)|
  early_late_mask(EarlyLateOutput::dm_over_rs)|thermal_ruler_mask);
 if (prepared.status()!=irred::numerics::Status::ok || result.status!=irred::numerics::Status::ok ||
  !result.ruler || result.ruler->status!=irred::numerics::Status::ok || !result.ruler->value || result.rows.size()!=2) return 10;
 const auto &present=result.rows[0].outputs[static_cast<unsigned>(EarlyLateOutput::e)];
 const auto &ratio=result.rows[1].outputs[static_cast<unsigned>(EarlyLateOutput::dm_over_rs)];
 return present.value && *present.value==1 && ratio.status==irred::numerics::Status::ok &&
  ratio.value && std::isfinite(*ratio.value) && *ratio.value>0 ? 0 : 10;
}
int installed_gaussian_posterior() {
 using namespace irred::statistics;
 Metadata m; m.ordered_ids={"r0","r1"};m.measure="dr0 dr1";
 m.table_identity="installed synthetic proper posterior";m.ordering_provenance="explicit";
 const std::array<double,4> c{2,.25,.25,1.5},x{1,.5,-.25,1};
 auto g=prepare_gaussian(c,MatrixKind::covariance,m,100,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
 ParameterPrior p{{"zero","colour"},{"mag","mag"},{"zero"},{.5,-.25},{1,.375,.375,.5},
  "synthetic proper prior","explicit X","mag","d(zero) d(colour)","declared independent noise",true};
 auto owner=GaussianPosterior::prepare(std::move(g),x,m.ordered_ids,p);
 const std::array<double,2> r{1.25,-.5},beta{.7,-.2};
 const auto mean=owner.condition(r,m.ordered_ids);
 const auto density=owner.log_density(r,m.ordered_ids,beta,p.ordered_parameter_ids);
 return mean.status==DensityStatus::finite && mean.value.size()==2 &&
  std::abs(mean.value[0]-20604./25511)<2e-12 && std::abs(mean.value[1]+7125./51022)<2e-12 &&
  density.density.status==DensityStatus::finite ? 0 : 11;
}
int installed_gr_growth() {
 using namespace irred::cosmology;
 auto owner=prepare_gr_growth(prepare(LCDM(1),FlatFLRW{}));
 const double a[]{.125,1}; auto result=owner.evaluate(a,growth_d|growth_f);
 return result.rows.size()==2 && result.rows[0].d.value==.125 && result.rows[1].d.value==1 && result.rows[0].f.value==1 && result.rows[1].f.value==1 ? 0 : 12;
}
int installed_sis_thin_lens() {
 using namespace irred::lensing;
 const long double theta=2*std::numbers::pi_v<long double>*std::pow(220/299792.458L,2);
 const SISSource source{{{70,0,1,0,1},1059,"installed synthetic radiation-only geometry"},
  .5,2,220,double(.4L*theta),1,double(.2L*theta),1};
 const auto owner=prepare_sis_thin_lens(source);
 const auto images=owner.predict(lens_positions|lens_magnifications|lens_fluxes|lens_delays);
 const std::array<PixelRectangle,1> rectangles{{{-1e-3,1e-3,-1e-3,1e-3}}};
 const auto pixels=owner.pixels(rectangles);
 if(owner.status()!=irred::numerics::Status::ok || images.status!=irred::numerics::Status::ok ||
  images.image_count!=2 || pixels.status!=irred::numerics::Status::ok || pixels.rows.size()!=1) return 15;
 const auto &positive=images.images[0], &negative=images.images[1];
 return positive.x_radians.value && negative.x_radians.value && positive.flux.value && negative.flux.value &&
  positive.relative_delay_seconds.value==0 && negative.relative_delay_seconds.value && *negative.relative_delay_seconds.value>0 &&
  positive.parity==1 && negative.parity==-1 &&
  std::abs(*positive.x_radians.value/(1.4L*theta)-1)<1e-8 &&
  std::abs(*negative.x_radians.value/(-.6L*theta)-1)<1e-8 &&
  std::abs(*positive.flux.value-3.5)<1e-8 && std::abs(*negative.flux.value-1.5)<1e-8 &&
  pixels.rows[0].flux.value && std::abs(*pixels.rows[0].flux.value-5)<1e-8 ? 0 : 15;
}
int installed_temporal() {
 using namespace irred::photometry;
 const double time[]{0,1},wave[]{1,2},lum[]{1,1,1,1},transmission[]{1,1};
 const TemporalGrid grid[]{ {"constant","synthetic",time,wave,lum} };
 const TemporalBand band[]{ {"optical","fixed",{wave,transmission}} };
 auto owner=prepare_temporal(grid,band);
 const TemporalExposure exposure[]{ {0,0,1,0,1,0,0,1} };
 auto result=evaluate_temporal(owner,exposure);
 return result.rows.size()==1 && result.rows[0].coverage==TemporalCoverage::full && result.rows[0].mean_flux_watt_per_square_metre.value && std::abs(*result.rows[0].mean_flux_watt_per_square_metre.value-1/(4*std::acos(-1.)))<2e-12 && result.rows[0].energy_joule.value && result.rows[0].expected_photons.value ? 0 : 13;
}
int installed_measurement_pipeline() {
 using namespace irred::photometry;
 using namespace irred::detector;
 const double time[]{0,1},wave[]{1,2},lum[]{1e-25,1e-25,1e-25,1e-25},transmission[]{1,1};
 const TemporalGrid grid[]{ {"constant","synthetic detector source",time,wave,lum} };
 const TemporalBand band[]{ {"optical","fixed transmission",{wave,transmission}} };
 const auto owner=prepare_temporal(grid,band);
 const TemporalExposure exposure[]{ {0,0,1,0,1,0,0,1} };
 const auto predicted=evaluate_temporal(owner,exposure);
 if(predicted.rows.size()!=1 || !predicted.rows[0].expected_photons.value) return 14;
 const irred::detector::Input model{PhotonLaw::poisson_arrivals,*predicted.rows[0].expected_photons.value,1,0,0,1,0,1,0};
 const Request request[]{ {model,{1,1}} };
 const auto measured=simulate(request,123,{1,1024*1024});
 if(measured.rows.size()!=1 || measured.rows[0].status!=irred::numerics::Status::ok || !measured.rows[0].poisson_electrons) return 14;
 const auto k=*measured.rows[0].poisson_electrons;
 const Observation observation[]{ {1,k>=1,k>=1 ? std::optional<std::uint32_t>(k) : std::nullopt,{},SelectionMeasure::joint_detection_record} };
 SelectionPolicy policy;policy.maximum_rows=1;policy.maximum_payload_bytes=1024*1024;
 const auto scored=likelihood(model,observation,policy);
 const auto expectation=moments(model);
 if(scored.rows.size()!=1 || scored.rows[0].status!=irred::numerics::Status::ok || !scored.rows[0].log_value || !expectation.poisson_mean_electrons) return 14;
 const double lambda=*expectation.poisson_mean_electrons;
 const double expected=k>=1 ? -lambda+k*std::log(lambda)-std::lgamma(k+1.) : -lambda;
 return std::abs(*scored.rows[0].log_value-expected)<2e-12 ? 0 : 14;
}
int installed_hydrogen_history() {
 namespace a=irred::atomic;
 using namespace irred::cosmology;
 const a::HydrogenState state[]{ {10000,6.769876143152199e20} };
 const auto fractions=a::evaluate_hydrogen_equilibrium(state);
 if(fractions.rows.size()!=1 || !fractions.rows[0].ionized.value || !fractions.rows[0].neutral.value ||
  std::abs(*fractions.rows[0].ionized.value-.5)>1e-12 || std::abs(*fractions.rows[0].neutral.value-.5)>1e-12) return 16;
 const PureHydrogenRequest source{{67.4,.02237,.12,2.7255,1.7e-5,{}},1600,300};
 const auto owner=prepare_pure_hydrogen_history(source);
 const double redshifts[]{300,1000};
 const auto history=owner.evaluate(redshifts,hydrogen_electron_fraction|hydrogen_drag_depth);
 const auto root=owner.conditional_unit_depth_redshift();
 if (!(owner.status()==irred::numerics::Status::ok && history.rows.size()==2 &&
  history.rows[0].drag_depth.value==0 && history.rows[1].electron_fraction.value &&
  *history.rows[1].electron_fraction.value>0 && *history.rows[1].electron_fraction.value<1 &&
  root.value && *root.value>1000 && *root.value<1100)) return 16;
 auto coupled_source=source;
 coupled_source.temperature_model=HydrogenTemperatureModel::evolved_compton_adiabatic;
 const auto coupled=prepare_pure_hydrogen_history(coupled_source);
 if(coupled.status()!=irred::numerics::Status::ok ||
  coupled.model_identity()!=evolved_hydrogen_history_id ||
  coupled.work().temperature_rhs_evaluations==0) return 17;
 const auto before=coupled.work().total();
 const auto optical=coupled.evaluate(redshifts,hydrogen_matter_temperature|
  hydrogen_thomson_depth|hydrogen_thomson_opacity|hydrogen_visibility|hydrogen_survival);
 if(optical.rows.size()!=2 || coupled.work().total()!=before ||
  optical.rows[0].thomson_depth.value!=0 || optical.rows[0].finite_endpoint_survival.value!=1) return 17;
 const auto& middle=optical.rows[1];
 return middle.matter_temperature_kelvin.value && middle.thomson_depth.value &&
  middle.thomson_opacity_per_redshift.value && middle.visibility_per_redshift.value &&
  middle.finite_endpoint_survival.value && *middle.matter_temperature_kelvin.value>0 &&
  *middle.matter_temperature_kelvin.value<source.model.tcmb_kelvin*1001 &&
  *middle.thomson_depth.value>0 && *middle.visibility_per_redshift.value>0 &&
  *middle.finite_endpoint_survival.value>0 && *middle.finite_endpoint_survival.value<1 ? 0 : 17;
}
int installed_relic_hydrogen_history() {
 using namespace irred::cosmology;
 const PureHydrogenRequest source{{67.4,.02237,.12,2.7255,0,
   {{.06,1.95,2},{0,1.95,2},{0,1.95,2}}},1600,300,
   HydrogenTemperatureModel::evolved_compton_adiabatic};
 PureHydrogenPolicy p;p.maximum_total_work=700000000;
 const auto h=prepare_pure_hydrogen_history(source,p);
 if(h.status()!=irred::numerics::Status::ok||h.model_identity()!=evolved_hydrogen_relic_history_id||
   h.work().momentum_callbacks<=4000000) return 20;
 const auto row=h.evaluate(std::array{1000.},127);
 if(row.rows.size()!=1) return 20;
 const auto& r=row.rows[0];
 return r.electron_fraction.value&&r.matter_temperature_kelvin.value&&r.drag_depth.value&&
 r.thomson_depth.value&&r.thomson_opacity_per_redshift.value&&r.visibility_per_redshift.value&&
 r.finite_endpoint_survival.value&&*r.electron_fraction.value>0&&*r.electron_fraction.value<1&&
 *r.matter_temperature_kelvin.value>0&&*r.visibility_per_redshift.value>0&&
 *r.finite_endpoint_survival.value>0&&*r.finite_endpoint_survival.value<1 ? 0 : 20;
}
int installed_gaussian_predictive() {
 using namespace irred::statistics;
 const std::array<double,1> C{1},r{2},R{.5},y{2};
 const std::array<double,2> X{1,0},A{2,0};
 const std::array<std::string,1> train_ids{"training"},future_ids{"future"};
 Metadata source;source.ordered_ids={"training"};source.measure="product d(mag)";
 source.source_semantics="synthetic controls";source.table_identity="installed training";
 source.uncertainty_identity="unit Gaussian noise";source.ordering_provenance="one axis";
 source.calibration_provenance="synthetic calibration";source.dependence_provenance="one full noise law";
 auto training=prepare_gaussian(C,MatrixKind::covariance,source,1024,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
 ParameterPrior prior;prior.ordered_parameter_ids={"beta","unused"};prior.parameter_units={"mag","mag"};
 prior.mean={0,0};prior.covariance={1,0,0,1};prior.prior_identity="installed proper prior";
 prior.design_identity="installed training X";prior.residual_unit="mag";prior.parameter_measure="product d(mag)";
 prior.dependence_identity="prior independent of training noise";prior.noise_independence_declared=true;
 auto posterior=GaussianPosterior::prepare(std::move(training),X,train_ids,prior);
 source.ordered_ids={"future"};source.table_identity="installed future noise";
 auto noise=prepare_gaussian(R,MatrixKind::covariance,source,1024,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
 PredictiveMetadata md;md.ordered_parameter_ids={"beta","unused"};md.parameter_units={"mag","mag"};md.response_units={"mag/mag","mag/mag"};
 md.training_event_ids={"training-event"};md.future_event_ids={"future-event"};md.future_unit="mag";
 md.future_covariance_unit="mag^2";md.future_measure="product d(mag)";md.response_identity="installed future A";
 md.conditioning_identity="installed fixed r";md.conditional_noise_identity="independent future noise";
 md.dependence_identity="future noise independent of prior/training noise";
 md.future_noise_independence_declared=true;md.noise_conditional_on_parameters_declared=true;
 auto prediction=GaussianPredictive::prepare(posterior,r,train_ids,noise,A,md);
 if(prediction.status()!=DensityStatus::finite||prediction.mean().size()!=1||prediction.covariance().size()!=1) return 19;
 auto density=prediction.log_density(y,future_ids);
 return std::abs(prediction.mean()[0]-2)<2e-12&&std::abs(prediction.covariance()[0]-2.5)<2e-12&&
 density.density.status==DensityStatus::finite&&std::abs(density.density.log_value+.5*std::log(5*std::numbers::pi))<2e-11 ? 0 : 19;
}
int installed_hydrogen_helium() {
 const irred::atomic::HydrogenHeliumState source{10000,1e12,1e11};
 const auto result=irred::atomic::evaluate_hydrogen_helium_equilibrium({&source,1});
 if(result.rows.size()!=1 || !result.rows[0].electron_density.value ||
    !result.rows[0].helium_doubly_ionized.value || !result.rows[0].hydrogen_neutral.value ||
    !result.rows[0].helium_singly_ionized.value || !result.rows[0].hydrogen_ionized.value) return 18;
 const auto& row=result.rows[0];
 const double charge=1e12* *row.hydrogen_ionized.value+1e11*( *row.helium_singly_ionized.value+2* *row.helium_doubly_ionized.value);
 return std::abs(*row.electron_density.value/charge-1)<2e-12 ? 0 : 18;
}
int installed_ladder_predictive() {
 using namespace irred::calibration;using namespace irred::statistics;
 Model m;m.ordered_host_ids={"A","B"};m.magnitude_convention="installed synthetic mag";
 m.metallicity_coordinate_identity="installed dex";m.distance_shape_identity="supplied fixed shape";
 m.calibration_identity="shared delta";m.dependence_identity="explicit independent synthetic noise";
 m.conditional_covariance_identity="conditional identity noise";
 m.rows={{RowKind::anchor_modulus,"train-A","A","",0,0,0,0}};
 auto noise=[](const Model&model,std::span<const double> c){Metadata md;
  for(const auto&r:model.rows)md.ordered_ids.push_back(r.row_id);
  md.measure="product d(mag)";md.source_semantics="synthetic controls";
  md.table_identity="installed synthetic covariance";md.ordering_provenance="exact rows";
  md.calibration_provenance=model.calibration_identity;md.dependence_provenance=model.dependence_identity;
  md.uncertainty_identity=model.conditional_covariance_identity;
  return prepare_gaussian(c,MatrixKind::covariance,md,1000,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);};
 auto l=linearize(m);ParameterPrior prior;prior.ordered_parameter_ids=l.metadata.ordered_parameter_ids;
 prior.parameter_units=l.metadata.parameter_units;prior.shared_nuisance_ids=l.metadata.shared_nuisance_ids;
 prior.mean={30,31,0,0,0,0,0,0};prior.covariance.resize(64);
 for(unsigned i=0;i<8;++i)prior.covariance[i*8+i]=1;
 prior.prior_identity="proper original installed Gaussian";prior.design_identity=l.metadata.design_identity;
 prior.residual_unit="mag";prior.parameter_measure="product parameter coordinate measure";
 prior.dependence_identity=l.metadata.dependence_identity;prior.noise_independence_declared=true;
 auto g=noise(m,std::array<double,1>{1});auto p=LadderPosterior::prepare(std::move(g),m,std::move(prior));
 if(p.status()!=DensityStatus::finite)return 23;
 auto future=m;future.rows={{RowKind::anchor_modulus,"future-B","B","",0,0,0,0},
  {RowKind::calibration_measurement,"future-delta","","",1,0,0,0}};
 auto fl=linearize(future);auto n=noise(future,std::array<double,4>{1,0,0,1});PredictiveMetadata md;
 md.ordered_parameter_ids=l.metadata.ordered_parameter_ids;md.parameter_units=l.metadata.parameter_units;
 for(const auto&u:md.parameter_units)md.response_units.push_back("mag/"+u);
 md.training_event_ids=l.predictive_event_ids;md.future_event_ids=fl.predictive_event_ids;
 md.future_unit="mag";md.future_covariance_unit="mag^2";md.future_measure="product d(mag)";
 md.response_identity=fl.metadata.design_identity;md.conditioning_identity="installed training vector";
 md.conditional_noise_identity=future.conditional_covariance_identity;md.dependence_identity=fl.metadata.dependence_identity;
 md.future_noise_independence_declared=true;md.noise_conditional_on_parameters_declared=true;
 auto q=LadderPredictive::prepare(p,std::array<double,1>{32},l.ordered_row_ids,n,future,std::move(md));
 return q.status()==DensityStatus::finite && q.mean().size()==2 && std::abs(q.mean()[0]-31)<2e-12 && q.mean()[1]==0 && std::abs(q.covariance()[0]-2)<2e-12 && std::abs(q.covariance()[3]-2)<2e-12 ? 0:23;
}
int installed_baryon_abundance() {
 using namespace irred::cosmology;
 const auto owner=prepare_baryon_abundance({.02237,.245,1.6735328383153192e-27,6.646479071583153e-27,"installed synthetic density/abundance","explicit synthetic neutral masses"});
 const std::array<BaryonAbundanceLteQuery,1> query{{{1./1601,4500}}};
 const auto result=owner.evaluate_equilibrium(query,"installed supplied matter temperature");
 return owner.status()==irred::numerics::Status::ok && result.status==irred::numerics::Status::ok && result.solves==1 && result.rows.size()==1 && result.rows[0].equilibrium.electron_density.value && *result.rows[0].equilibrium.electron_density.value>0 ? 0:21;
}
int installed_nested_fd() {
 using namespace irred::cosmology;
 ThermalPolicy p;p.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis;
 const auto moments=evaluate_thermal_moments(1,p);
 const auto background=prepare_thermal_background({70,5e-5,2e-5,.05,.25,{{.06,.000168,2}}},p);
 const auto batch=background.evaluate(std::array{.1},thermal_e,p);
 return moments.status==irred::numerics::Status::ok && moments.rho_moment && std::abs(*moments.rho_moment-6.045645184985879)<1e-9 && background.momentum_method()==p.momentum_method && batch.status==irred::numerics::Status::ok && batch.rows.size()==1 && batch.rows[0].e.value && std::abs(*batch.rows[0].e.value-17.398718835753673)<1e-8 ? 0:22;
}
int installed_gaussian_simulation() {
 using namespace irred;
 statistics::Metadata metadata;metadata.ordered_ids={"coordinate-a","coordinate-b"};metadata.measure="product d(mag)";metadata.source_semantics="synthetic controls";metadata.table_identity="installed correlated synthetic generating noise";metadata.ordering_provenance="two declared axes";
 const auto covariance=statistics::prepare_gaussian(std::array<double,4>{1,.5,.5,2},statistics::MatrixKind::covariance,metadata,4,1e-10,numerics::Arithmetic::longdouble_cpu_v1);
 statistics::GeneratingMean mean{{1,2},metadata.ordered_ids,{"mag","mag"},"installed supplied synthetic mean",metadata.measure,"ideal normal/discrete addressed approximation"};
 const std::array<random::Address,2> addresses{{{0x1100,0},{0x1100,1}}};
 statistics::GaussianSimulationPolicy policy;policy.outputs=7;
 const auto r=statistics::GaussianSimulation::simulate(covariance,mean,addresses,17,policy);
 if(r.status!=numerics::Status::ok || r.rows.size()!=2 || r.values.size()!=4 || r.words.size()!=4) return 26;
 for(unsigned i=0;i<2;++i){
  if(r.rows[i].status!=numerics::Status::ok) return 26;
  std::array<long double,2> z{};
  for(unsigned j=0;j<2;++j){const auto w=random::words(17,{0x1100+j,i});if(w!=r.words[2*i+j]) return 26;
   const long double u=(static_cast<long double>(w[1])+.5L)/4294967296.L,v=(static_cast<long double>(w[2])+.5L)/4294967296.L;
   z[j]=std::sqrt(-2*std::log(u))*std::cos(2*std::numbers::pi_v<long double>*v);}
  if(std::abs(r.values[2*i]-(1+z[0]))>2e-12L || std::abs(r.values[2*i+1]-(2+.5L*z[0]+std::sqrt(1.75L)*z[1]))>2e-12L) return 26;
 }
 return 0;
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
 if (const auto thermal_status=installed_thermal_observables();thermal_status!=0) return thermal_status;
 if (const auto posterior_status=installed_gaussian_posterior();posterior_status!=0) return posterior_status;
 if (const auto growth_status=installed_gr_growth();growth_status!=0) return growth_status;
 if (const auto hydrogen_status=installed_hydrogen_history();hydrogen_status!=0) return hydrogen_status;
 if (const auto lens_status=installed_sis_thin_lens();lens_status!=0) return lens_status;
 if (const auto temporal_status=installed_temporal();temporal_status!=0) return temporal_status;
 if (const auto pipeline_status=installed_measurement_pipeline();pipeline_status!=0) return pipeline_status;
 if (const auto relic_status=installed_relic_hydrogen_history();relic_status!=0) return relic_status;

 if (const auto predictive_status=installed_gaussian_predictive();predictive_status!=0) return predictive_status;
 if (const auto mixture_status=installed_hydrogen_helium();mixture_status!=0) return mixture_status;
 if(const auto ladder_future_status=installed_ladder_predictive();ladder_future_status!=0) return ladder_future_status;
 if(const auto abundance_status=installed_baryon_abundance();abundance_status!=0) return abundance_status;
 if(const auto nested_status=installed_nested_fd();nested_status!=0) return nested_status;
 if(const auto simulation_status=installed_gaussian_simulation();simulation_status!=0) return simulation_status;
 return installed_correlated_calibration();
}
