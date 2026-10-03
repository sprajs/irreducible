// Standalone installed-header consumer. Identity arguments are declarations;
// the owning external wrapper must independently attest archive/compiler/binary.
#include "irred/gaussian_box_heldout.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <stdexcept>
namespace box_heldout_controls {
using namespace irred::statistics;
struct Input {
  Metadata rows;
  DesignMetadata design;
  BoxSupport box;
  BoxHeldoutMetadata heldout;
  std::size_t heldout_index=2;
  std::vector<double> covariance={1,0,.5,0,1,.25,.5,.25,1},x={1,0,0,1,.5,.25},offsets={1,2,3};
  Input() {
    rows.ordered_ids={"r0","r1","r2"};rows.measure="product d(u)";
    rows.table_identity="original dyadic correlated control/v1";
    rows.uncertainty_identity="full supplied dyadic3 covariance/v1";
    rows.ordering_provenance="literal r0,r1,r2";rows.calibration_provenance="synthetic fixed offsets";
    rows.dependence_provenance="full C including k=(1/2,1/4)";rows.source_semantics="synthetic controls";
    design={{"original0","original2"},{"u0","u2"},{},"u","dyadic2 active source coordinates/v1",rows.dependence_provenance};
    box={{"original0","original2"},{-8,-9},{8,9},"product d(original0)d(original2)","normalized original finite box/v1",
        "synthetic original1 literal0 point mass outside active2"};
    heldout.source_contract_id="synthetic-original-dyadic/v1";
    heldout.candidate_identity="withhold literal original r2/v1";
    heldout.conditioning_identity="fixed C/X/offsets/box; synthetic conditional law";
    heldout.original_row_lineage="r0,r1,r2 exact original order";
    heldout.event_lineage="synthetic controls; no observational event claim";
    heldout.calibration_dependence_identity=rows.dependence_provenance;
    heldout.heldout_unit="u";heldout.heldout_covariance_unit="u^2";heldout.heldout_measure="d(u)";
    heldout.ordered_original_parameter_ids={"original0","original1","original2"};
    heldout.active_original_parameter_indices={0,2};heldout.fixed_original_parameter_index=1;
    heldout.fixed_original_parameter_value=0;heldout.fixed_coordinate_provenance=box.fixed_coordinate_provenance;
  }
  Gaussian gaussian(MatrixKind kind=MatrixKind::covariance)const {
    return prepare_gaussian(covariance,kind,rows,covariance.size(),1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
  }
  BoxHeldoutPolicy policy(const Gaussian&g,std::size_t trainings=1,std::size_t candidates=1,std::size_t requests=1)const {
    BoxHeldoutPolicy p;p.maximum_training_vectors=trainings;p.maximum_candidate_values=candidates;p.maximum_requests=requests;
    const auto setup=GaussianBoxHeldout::preparation_work_bound(g,design,box,heldout);
    std::size_t identity_bytes=0;for(std::size_t j=0;j<rows.ordered_ids.size();++j)if(j!=heldout_index)identity_bytes+=rows.ordered_ids[j].size();
    const auto evaluation=GaussianBoxHeldout::evaluation_work_bound(rows.ordered_ids.size(),design.ordered_parameter_ids.size(),trainings,candidates,requests,identity_bytes);
    if(!setup||!evaluation)throw std::runtime_error("source work bound unavailable");
    p.maximum_preparation_work_units=*setup;p.maximum_evaluation_work_units=*evaluation;return p;
  }
  GaussianBoxHeldout prepare(Gaussian&&g,BoxHeldoutPolicy p)const {
    return GaussianBoxHeldout::prepare(std::move(g),x,offsets,rows.ordered_ids,design,box,heldout_index,heldout,p);
  }
};
} // namespace box_heldout_controls

#include <cfenv>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>
using namespace irred::statistics;
namespace {
void string(std::string_view s) {
  std::cout<<'"';for(const unsigned char c:s){if(c=='"')std::cout<<"\\\"";
    else if(c=='\\')std::cout<<"\\\\";else if(c<32)std::cout<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(c)<<std::dec<<std::setfill(' ');
    else std::cout<<char(c);}std::cout<<'"';
}
void number(double x,bool available=true) {
  if(!available){std::cout<<"null";return;}if(std::isfinite(x)){std::cout<<x;return;}
  std::cout<<"{\"kind\":\"nonfinite\",\"classification\":";
  string(std::isnan(x)?"nan":std::signbit(x)?"negative-infinity":"positive-infinity");std::cout<<'}';
}
void witness(std::optional<long double>x) {
  if(!x){std::cout<<"null";return;}
  if(!std::isfinite(*x)){number(static_cast<double>(*x));return;}
  // A finite wide witness can exceed binary64. Preserve its decimal as a
  // tagged string rather than asking a JSON binary64 parser to overflow it.
  std::cout<<"{\"kind\":\"finite-long-double\",\"value\":\""<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<*x<<"\"}";
  std::cout<<std::setprecision(std::numeric_limits<double>::max_digits10);
}
void interval(BoxInterval x,bool available=true) {
  if(!available){std::cout<<"null";return;}std::cout<<'[';number(x.lower);std::cout<<',';number(x.upper);std::cout<<']';
}
void count(std::optional<std::size_t>x){if(x)std::cout<<*x;else std::cout<<"null";}
void ids(std::span<const std::string>x){std::cout<<'[';for(std::size_t j=0;j<x.size();++j){if(j)std::cout<<',';string(x[j]);}std::cout<<']';}
void work(const BoxHeldoutWork&w) {
  std::cout<<"{\"factor_attempts\":"<<w.factor_attempts<<",\"factors_completed\":"<<w.factors_completed
    <<",\"whitening_attempts\":"<<w.whitening_attempts<<",\"whitenings_completed\":"<<w.whitenings_completed
    <<",\"training_qr_attempts\":"<<w.training_qr_attempts<<",\"training_qr_completed\":"<<w.training_qr_completed
    <<",\"update_qr_attempts\":"<<w.update_qr_attempts<<",\"update_qr_completed\":"<<w.update_qr_completed
    <<",\"training_attempts\":"<<w.training_attempts<<",\"training_completed\":"<<w.training_completed
    <<",\"joint_attempts\":"<<w.joint_attempts<<",\"joint_completed\":"<<w.joint_completed
    <<",\"variance_screen_attempts\":"<<w.variance_screen_attempts<<",\"variance_screens_completed\":"<<w.variance_screens_completed
    <<",\"tail_bound_attempts\":"<<w.tail_bound_attempts<<",\"tail_bounds_completed\":"<<w.tail_bounds_completed
    <<",\"normalization_attempts\":"<<w.normalization_attempts<<",\"normalizations_completed\":"<<w.normalizations_completed
    <<",\"composition_attempts\":"<<w.composition_attempts<<",\"compositions_completed\":"<<w.compositions_completed
    <<",\"charged_work_units\":"<<w.charged_work_units<<'}';
}
void normalization(const GaussianBoxResult&r) {
  std::cout<<"{\"status\":"<<int(r.status)<<",\"numerical_status\":"<<int(r.numerical_status)
    <<",\"stage\":"<<int(r.stage)<<",\"completion_step\":"<<int(r.completion_step)<<",\"completion_parameter_index\":";
  count(r.completion_parameter_index);std::cout<<",\"availability\":["<<r.gaussian_completion_available<<','<<r.endpoint_margins_available<<','
    <<r.rectangle_enclosure_available<<','<<r.normalization_enclosures_available<<','<<r.quantile_enclosure_available<<','<<r.endpoint_cdf_enclosures_available<<']';
  std::cout<<",\"q_min\":";number(r.minimum_quadratic,r.gaussian_completion_available);
  std::cout<<",\"postcast_profile_q\":";number(r.reported_profile_quadratic,r.gaussian_completion_available);
  std::cout<<",\"logdetH\":";number(r.log_design_precision_determinant,r.gaussian_completion_available);
  std::cout<<",\"logdetC\":";number(r.source_covariance_log_determinant,r.gaussian_completion_available);
  std::cout<<",\"stationarity\":";number(r.profile_stationarity,r.gaussian_completion_available);
  std::cout<<",\"variance_sensitivity\":";number(r.maximum_variance_sensitivity,r.gaussian_completion_available);
  auto scalars=[&](const auto&v){if(!r.gaussian_completion_available){std::cout<<"null";return;}std::cout<<'[';
    for(std::size_t j=0;j<v.size();++j){if(j)std::cout<<',';number(v[j]);}std::cout<<']';};
  std::cout<<",\"mean\":";scalars(r.unboxed_mean);std::cout<<",\"variance\":";scalars(r.unboxed_variance);
  auto margins=[&](const auto&v){if(!r.endpoint_margins_available){std::cout<<"null";return;}std::cout<<'[';
    for(std::size_t j=0;j<v.size();++j){if(j)std::cout<<',';interval(v[j]);}std::cout<<']';};
  std::cout<<",\"lower_margins\":";margins(r.standardized_lower_margins);std::cout<<",\"upper_margins\":";margins(r.standardized_upper_margins);
  std::cout<<",\"excluded_mass_upper\":";number(r.excluded_mass_upper,r.endpoint_margins_available);
  std::cout<<",\"box_probability\":";interval(r.box_probability,r.rectangle_enclosure_available);
  std::cout<<",\"log_box_probability\":";interval(r.log_box_probability,r.rectangle_enclosure_available);
  std::cout<<",\"log_prior_volume\":";interval(r.log_prior_volume,r.endpoint_margins_available);
  std::cout<<",\"log_relative_integral\":";interval(r.log_relative_box_integral,r.normalization_enclosures_available);
  std::cout<<",\"log_prior_normalized_evidence\":";interval(r.log_prior_normalized_relative_evidence,r.normalization_enclosures_available);
  std::cout<<",\"log_observation_normalized_evidence\":";interval(r.log_observation_normalized_evidence,r.normalization_enclosures_available);
  std::cout<<",\"cdf_nodes\":"<<r.cdf_node_evaluations<<",\"cdf_attempts\":"<<r.cdf_evaluations<<",\"bisections\":"<<r.bisections<<'}';
}
bool hex(std::string_view s,std::size_t n){return s.size()==n&&std::all_of(s.begin(),s.end(),[](char c){return(c>='0'&&c<='9')||(c>='a'&&c<='f');});}
void emit(const GaussianBoxHeldout&o,const BoxHeldoutBatch&b,char**argv,std::string_view mode,
    DensityStatus source,irred::numerics::Status source_numerical,std::optional<MatrixKind> source_kind,
    DensityStatus after,const BoxHeldoutPolicy&p) {
  std::cout<<std::boolalpha<<std::setprecision(std::numeric_limits<double>::max_digits10);
  std::cout<<"{\"schema_version\":1,\"declared_wrapper_identity\":{\"engine_head\":";string(argv[1]);
  std::cout<<",\"sdk_manifest_sha256\":";string(argv[2]);std::cout<<",\"compile_manifest_sha256\":";string(argv[3]);
  std::cout<<",\"executable_sha256\":";string(argv[4]);std::cout<<"},\"mode\":";string(mode);
  std::cout<<",\"actual_environment\":{\"rounding_mode\":"<<std::fegetround()<<",\"double_digits\":"<<std::numeric_limits<double>::digits
    <<",\"long_double_digits\":"<<std::numeric_limits<long double>::digits<<",\"long_double_max_exponent\":"<<std::numeric_limits<long double>::max_exponent<<'}';
  std::cout<<",\"method_id\":";string(b.method_id);std::cout<<",\"law_id\":";string(b.law_id);std::cout<<",\"enclosure_scope\":";string(b.enclosure_scope);
  std::cout<<",\"declared_control_inputs\":{\"original_values\":[1.25,1.5,3],\"offsets\":[1,2,3],\"training_pool\":";
  if(mode=="shape-refused")std::cout<<"[1]";else std::cout<<"[1.25,1.5]";
  std::cout<<",\"candidate_pool\":[3],\"requested_pairs\":[[0,0]]}";
  std::cout<<",\"original_input_matrix_kind\":";if(source_kind)std::cout<<int(*source_kind);else std::cout<<"null";
  std::cout<<",\"source_numerical_status\":"<<int(source_numerical);
  std::cout<<",\"source_status\":"<<int(source)<<",\"source_after_prepare_status\":"<<int(after)<<",\"actual_source_order\":";ids(o.original_row_ids());
  std::cout<<",\"actual_training_order\":";ids(o.training_row_ids());const auto&md=o.metadata();
  std::cout<<",\"actual_metadata\":";
  if(!o.preparation().retained_frames_available)std::cout<<"null";else {
  std::cout<<"{\"source_contract_id\":";string(md.source_contract_id);std::cout<<",\"candidate_identity\":";string(md.candidate_identity);
  std::cout<<",\"conditioning_identity\":";string(md.conditioning_identity);
  std::cout<<",\"original_row_lineage\":";string(md.original_row_lineage);std::cout<<",\"event_lineage\":";string(md.event_lineage);
  std::cout<<",\"calibration_dependence_identity\":";string(md.calibration_dependence_identity);
  std::cout<<",\"heldout_unit\":";string(md.heldout_unit);std::cout<<",\"heldout_covariance_unit\":";string(md.heldout_covariance_unit);
  std::cout<<",\"original_parameter_order\":";ids(md.ordered_original_parameter_ids);std::cout<<",\"active_original_axes\":[";
  for(std::size_t j=0;j<md.active_original_parameter_indices.size();++j){if(j)std::cout<<',';std::cout<<md.active_original_parameter_indices[j];}
  std::cout<<"],\"fixed_coordinate_provenance\":";string(md.fixed_coordinate_provenance);
  std::cout<<",\"heldout_measure\":";string(md.heldout_measure);std::cout<<",\"fixed_original_axis\":"<<md.fixed_original_parameter_index<<",\"fixed_value\":";
  number(md.fixed_original_parameter_value);std::cout<<'}';}
  const auto&r=o.preparation();std::cout<<",\"preparation\":{\"status\":"<<int(r.status)<<",\"numerical_status\":"<<int(r.numerical_status)<<",\"stage\":"<<int(r.stage)
    <<",\"availability\":["<<r.cross_whitening_available<<','<<r.schur_available<<','<<r.schur_sqrt_available<<','<<r.appended_response_available<<','
    <<r.training_qr_available<<','<<r.update_qr_available<<','<<r.covariance_determinants_available<<','<<r.retained_frames_available<<']';
  std::cout<<",\"schur_variance\":";number(r.schur_variance,r.schur_available);std::cout<<",\"schur_cancellation_ratio\":";number(r.schur_cancellation_ratio,r.schur_available);
  std::cout<<",\"schur_sensitivity_estimate\":";number(r.schur_sensitivity_estimate,r.schur_available);std::cout<<",\"schur_sqrt_enclosure\":";interval(r.schur_sqrt_enclosure,r.schur_sqrt_available);
  std::cout<<",\"cross_whitening_backward_residual\":";number(r.cross_whitening_backward_residual,r.cross_whitening_available);
  std::cout<<",\"cross_whitening_rounding_estimate\":";number(r.cross_whitening_rounding_estimate,r.cross_whitening_available);
  std::cout<<",\"appended_projection_rounding_estimate\":";number(r.appended_projection_rounding_estimate,r.appended_response_available);
  std::cout<<",\"training_preparation_sensitivity\":";number(r.training_preparation_sensitivity,r.training_qr_available);
  std::cout<<",\"training_triangular_condition\":";number(r.training_triangular_condition,r.training_qr_available);
  std::cout<<",\"training_transpose_condition\":";number(r.training_transpose_condition,r.training_qr_available);
  std::cout<<",\"update_triangular_condition\":";number(r.update_triangular_condition,r.update_qr_available);
  std::cout<<",\"update_transpose_condition\":";number(r.update_transpose_condition,r.update_qr_available);
  std::cout<<",\"full_source_logdetC\":";number(r.full_source_covariance_log_determinant,r.covariance_determinants_available);
  std::cout<<",\"training_logdetC\":";number(r.training_covariance_log_determinant,r.covariance_determinants_available);
  std::cout<<",\"refusal_witness\":";witness(r.refusal_witness);std::cout<<",\"refusal_index\":";count(r.refusal_index);
  std::cout<<",\"peak_payload_bound\":";count(r.peak_payload_bound);std::cout<<",\"retained_payload_bound\":";count(r.retained_payload_bound);
  std::cout<<",\"work\":";work(r.work);std::cout<<"},\"frozen_work_caps\":["<<p.maximum_preparation_work_units<<','<<p.maximum_evaluation_work_units<<']';
  std::cout<<",\"batch_status\":"<<int(b.status)<<",\"batch_numerical_status\":"<<int(b.numerical_status)<<",\"scratch_output_payload_bound\":";count(b.scratch_output_payload_bound);
  std::cout<<",\"training_normalizations\":[";for(std::size_t j=0;j<b.training_normalizations.size();++j){if(j)std::cout<<',';normalization(b.training_normalizations[j]);}
  std::cout<<"],\"training_refusal_witnesses\":[";for(std::size_t j=0;j<b.training_refusal_witnesses.size();++j){if(j)std::cout<<',';witness(b.training_refusal_witnesses[j]);}
  std::cout<<"],\"training_offset_subtraction_rounding_estimates\":[";
  for(std::size_t j=0;j<b.training_offset_subtraction_rounding_estimates.size();++j){if(j)std::cout<<',';
    number(b.training_offset_subtraction_rounding_estimates[j],b.training_offset_diagnostics_available[j]);}
  std::cout<<"],\"training_cross_projection_error_estimates\":[";
  for(std::size_t j=0;j<b.training_cross_projection_error_estimates.size();++j){if(j)std::cout<<',';
    number(b.training_cross_projection_error_estimates[j],b.training_cross_projection_diagnostics_available[j]);}
  std::cout<<"],\"training_offset_diagnostics_available\":[";
  for(std::size_t j=0;j<b.training_offset_diagnostics_available.size();++j){if(j)std::cout<<',';std::cout<<bool(b.training_offset_diagnostics_available[j]);}
  std::cout<<"],\"training_cross_projection_diagnostics_available\":[";
  for(std::size_t j=0;j<b.training_cross_projection_diagnostics_available.size();++j){if(j)std::cout<<',';std::cout<<bool(b.training_cross_projection_diagnostics_available[j]);}
  std::cout<<"],\"actual_request_pairs\":[";
  for(std::size_t j=0;j<b.requests.size();++j){if(j)std::cout<<',';std::cout<<'['<<b.requests[j].training_vector_index<<','<<b.requests[j].candidate_value_index<<']';}
  std::cout<<"],\"densities\":[";for(std::size_t j=0;j<b.densities.size();++j){if(j)std::cout<<',';const auto&d=b.densities[j];
    std::cout<<"{\"status\":"<<int(d.status)<<",\"numerical_status\":"<<int(d.numerical_status)<<",\"stage\":"<<int(d.stage)<<",\"density_available\":"<<d.density_available
      <<",\"log_density\":";interval(d.log_density,d.density_available);std::cout<<",\"constant_available\":"<<d.scalar_conditional_constant_available;
    std::cout<<",\"conditional_constant\":";interval(d.scalar_conditional_log_normalization,d.scalar_conditional_constant_available);
    std::cout<<",\"cross_projection_rounding_estimate\":";number(d.cross_projection_rounding_estimate,d.cross_projection_diagnostics_available);
    std::cout<<",\"offset_subtraction_rounding_estimate\":";number(d.offset_subtraction_rounding_estimate,d.offset_diagnostics_available);
    std::cout<<",\"update_rhs_rounding_estimate\":";number(d.update_rhs_rounding_estimate,d.update_diagnostics_available);
    std::cout<<",\"refusal_witness\":";witness(d.refusal_witness);std::cout<<",\"refusal_index\":";count(d.refusal_index);
    std::cout<<",\"joint_normalization\":";normalization(d.joint_normalization);std::cout<<",\"work\":";work(d.work);std::cout<<'}';}
  std::cout<<"],\"work\":";work(b.work);std::cout<<",\"qualification\":\"synthetic numerical control; physical interpretation unavailable\"}\n";
}
}
int main(int argc,char**argv) {
  if(argc!=6||!hex(argv[1],40)||!hex(argv[2],64)||!hex(argv[3],64)||!hex(argv[4],64))return 1;
  const std::string_view mode=argv[5];if(mode!="accepted"&&mode!="schur-refused"&&mode!="shape-refused"&&mode!="wire-nonfinite-control")return 1;
  box_heldout_controls::Input in;if(mode=="schur-refused")in.covariance[8]=.3125+0x1p-50;
  auto g=in.gaussian();const auto before=g.status();const auto source_numerical=g.numerical_status();const auto source_kind=g.input_matrix_kind();auto p=in.policy(g);auto owner=in.prepare(std::move(g),p);
  BoxHeldoutBatch b;
  if(owner.status()==DensityStatus::finite){const std::vector<double>training=mode=="shape-refused"?std::vector<double>{1}:std::vector<double>{1.25,1.5};
    const std::vector<double>values={3};const std::vector<BoxHeldoutRequest>requests={{0,0}};b=owner.evaluate(training,owner.training_row_ids(),values,requests,p);}
  // Explicit serializer-only witness exercises legal transport of nonfinite
  // diagnostics; it does not claim the engine produced this value.
  if(mode=="wire-nonfinite-control") {std::cout<<"{\"scope\":\"serializer-only-synthetic-witness\",\"witness\":";
    witness(std::numeric_limits<long double>::infinity());std::cout<<",\"absent\":";witness(std::nullopt);std::cout<<"}\n";return 0;}
  emit(owner,b,argv,mode,before,source_numerical,source_kind,g.status(),p);
  if(owner.status()==DensityStatus::finite) {
    const auto same=[](auto a,auto b){return a.size()==b.size()&&std::equal(a.begin(),a.end(),b.begin());};
    const std::vector<std::string> expected_training={"r0","r1"};const auto&md=owner.metadata();
    if(!same(owner.original_row_ids(),in.rows.ordered_ids)||!same(owner.training_row_ids(),expected_training)||
        md.ordered_original_parameter_ids!=in.heldout.ordered_original_parameter_ids||
        md.active_original_parameter_indices!=in.heldout.active_original_parameter_indices||
        md.fixed_original_parameter_index!=in.heldout.fixed_original_parameter_index||
        md.fixed_original_parameter_value!=in.heldout.fixed_original_parameter_value||
        md.fixed_coordinate_provenance!=in.heldout.fixed_coordinate_provenance||
        md.source_contract_id!=in.heldout.source_contract_id||md.candidate_identity!=in.heldout.candidate_identity||
        md.heldout_unit!=in.heldout.heldout_unit||md.heldout_measure!=in.heldout.heldout_measure)return 3;
  }
  if(mode=="accepted")return owner.status()==DensityStatus::finite&&b.status==DensityStatus::finite&&b.densities.size()==1&&b.densities[0].density_available?0:2;
  if(mode=="schur-refused")return owner.status()==DensityStatus::numerical_failure&&g.status()==DensityStatus::finite?0:2;
  return b.status==DensityStatus::invalid_input?0:2;
}
