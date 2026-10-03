// Source controls only until separately admitted compilation/execution.
#include "box_heldout_controls.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
using namespace irred::statistics;
namespace {
void check(bool v,const char*m){if(!v)throw std::runtime_error(m);}
bool contains(BoxInterval v,double x){return v.lower<=x&&x<=v.upper;}
void near(double a,double b){check(std::abs(a-b)<2e-12*(1+std::abs(b)),"analytic completion mismatch");}
BoxHeldoutBatch evaluate(const GaussianBoxHeldout&o,BoxHeldoutPolicy p,
    std::vector<double>training={1.25,1.5},std::vector<double>candidate={3},
    std::vector<BoxHeldoutRequest>requests={{0,0}}) {
  return o.evaluate(training,o.training_row_ids(),candidate,requests,p);
}
void work_bound_ownership_controls() {
  box_heldout_controls::Input base;auto g=base.gaussian();
  const auto expected=GaussianBoxHeldout::preparation_work_bound(g,base.design,base.box,base.heldout);
  check(expected&&*expected>1,"logical preparation work bound available");
  auto reserved=base;auto rows=base.rows;
  auto texts=[](std::vector<std::string>&v){v.reserve(v.size()+32);for(auto&s:v)s.reserve(s.size()+256);};
  auto elements=[](auto&v){v.reserve(v.size()+32);};
  texts(rows.ordered_ids);
  for(auto*s:{&rows.arithmetic_id,&rows.measure,&rows.table_identity,&rows.uncertainty_identity,
      &rows.ordering_provenance,&rows.calibration_provenance,&rows.dependence_provenance,
      &rows.source_semantics,&rows.input_matrix_convention,&rows.treatment})s->reserve(s->size()+256);
  texts(reserved.design.ordered_parameter_ids);texts(reserved.design.parameter_units);texts(reserved.design.shared_nuisance_ids);
  for(auto*s:{&reserved.design.residual_unit,&reserved.design.design_identity,&reserved.design.dependence_identity})s->reserve(s->size()+256);
  texts(reserved.box.ordered_parameter_ids);elements(reserved.box.lower);elements(reserved.box.upper);
  for(auto*s:{&reserved.box.parameter_measure,&reserved.box.prior_identity,&reserved.box.fixed_coordinate_provenance})s->reserve(s->size()+256);
  texts(reserved.heldout.ordered_original_parameter_ids);elements(reserved.heldout.active_original_parameter_indices);
  for(auto*s:{&reserved.heldout.source_contract_id,&reserved.heldout.candidate_identity,&reserved.heldout.conditioning_identity,
      &reserved.heldout.original_row_lineage,&reserved.heldout.event_lineage,&reserved.heldout.calibration_dependence_identity,
      &reserved.heldout.heldout_unit,&reserved.heldout.heldout_covariance_unit,&reserved.heldout.heldout_measure,
      &reserved.heldout.fixed_coordinate_provenance})s->reserve(s->size()+256);
  auto reserved_g=prepare_gaussian(base.covariance,MatrixKind::covariance,std::move(rows),
      base.covariance.size(),1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
  check(reserved_g.status()==DensityStatus::finite,"reserved Gaussian metadata admitted");
  check(reserved_g.metadata().ordered_ids.capacity()>g.metadata().ordered_ids.capacity()&&
      reserved.design.ordered_parameter_ids.capacity()>base.design.ordered_parameter_ids.capacity(),"control owns actual spare metadata capacity");
  const auto payload=GaussianBoxHeldout::preparation_payload_bound(g,base.design,base.box,base.heldout);
  const auto reserved_payload=GaussianBoxHeldout::preparation_payload_bound(reserved_g,reserved.design,reserved.box,reserved.heldout);
  check(payload&&reserved_payload&&*reserved_payload>*payload,"payload still charges actual reserved capacities");
  auto bound_matches=[&](const Gaussian&full,const box_heldout_controls::Input&input){
    check(GaussianBoxHeldout::preparation_work_bound(full,input.design,input.box,input.heldout)==expected,
        "logical work bound invariant under reserve copy and move");
  };
  bound_matches(reserved_g,reserved);
  auto copied=reserved;auto copied_g=reserved_g;bound_matches(copied_g,copied);
  auto moved=std::move(copied);auto moved_g=std::move(copied_g);bound_matches(moved_g,moved);
  auto cap=moved.policy(moved_g);check(cap.maximum_preparation_work_units==*expected,"caller freezes logical bound");
  --cap.maximum_preparation_work_units;
  const auto retained=moved_g.retained_payload_bound();
  check(retained.has_value(),"source retained payload available before work refusal");
  auto refused=moved.prepare(std::move(moved_g),cap);
  std::cout<<"metadata_work_control bound="<<*expected<<" cap="<<cap.maximum_preparation_work_units
      <<" refused_status="<<static_cast<unsigned>(refused.status())
      <<" refused_numerical_status="<<static_cast<unsigned>(refused.preparation().numerical_status)
      <<" refused_charge="<<refused.preparation().work.charged_work_units
      <<" source_status="<<static_cast<unsigned>(moved_g.status())<<'\n';
  check(refused.status()==DensityStatus::numerical_failure&&refused.preparation().numerical_status==irred::numerics::Status::work_limit,
      "reserved copied metadata one-below work bound refuses");
  check(refused.preparation().work.charged_work_units==0&&refused.preparation().work.factor_attempts==0&&
      refused.preparation().work.whitening_attempts==0&&refused.preparation().work.training_qr_attempts==0,
      "one-below work refusal precedes numerical attempts");
  check(moved_g.status()==DensityStatus::finite&&moved_g.input_matrix_kind()==MatrixKind::covariance&&
      moved_g.retained_payload_bound()==retained&&moved_g.metadata().ordered_ids==base.rows.ordered_ids&&
      std::equal(moved_g.covariance().begin(),moved_g.covariance().end(),base.covariance.begin(),base.covariance.end()),
      "one-below work refusal retains complete source identity and covariance");
  auto exact=reserved.policy(reserved_g);auto prepared=reserved.prepare(std::move(reserved_g),exact);
  check(prepared.status()==DensityStatus::finite&&prepared.preparation().work.charged_work_units==*expected,
      "exact logical bound remains sufficient after metadata copies");
}
}
int main(){try {
  work_bound_ownership_controls();
  box_heldout_controls::Input in;auto g=in.gaussian();auto policy=in.policy(g,2,2,3);
  check(g.input_matrix_kind()==MatrixKind::covariance,"authoritative covariance kind");
  auto copy=g;check(copy.input_matrix_kind()==MatrixKind::covariance,"copy preserves authoritative kind");
  auto copy_moved=std::move(copy);check(copy_moved.input_matrix_kind()==MatrixKind::covariance&&!copy.input_matrix_kind(),"kind follows Gaussian move");
  copy_moved=std::move(copy_moved);check(copy_moved.input_matrix_kind()==MatrixKind::covariance,"Gaussian self move retains kind");
  auto owner=in.prepare(std::move(g),policy);
  check(owner.status()==DensityStatus::finite&&g.status()==DensityStatus::invalid_input,"successful preparation consumes source");
  check(!g.input_matrix_kind(),"moved source kind cleared");
  near(owner.preparation().schur_variance,11./16);
  check(owner.preparation().work.factor_attempts==1&&owner.preparation().work.factors_completed==1,"one training factor");
  auto result=evaluate(owner,policy);check(result.status==DensityStatus::finite,"batch shape admitted");
  check(result.output_layout_available,"complete public output layout earned");
  check(result.densities[0].density_available,"analytic correlated heldout admitted");
  // Exact a=0 makes parameter integration cancel for every box and yT. The
  // scalar target is log sqrt(8/(11*pi)); libm here is a convenience check,
  // never the independently pinned original-source reference.
  const auto target=.5*std::log(8/(11*std::acos(-1.)));
  check(contains(result.densities[0].log_density,target),"analytic a=0 conditional law");
  check(result.training_normalizations[0].cdf_node_evaluations==0&&result.densities[0].joint_normalization.cdf_node_evaluations==0,"normalization needs no quantile");
  auto pooled=evaluate(owner,policy,{1.25,1.5,1.5,2},{3,3.25},{{1,0},{0,1},{1,0}});
  check(pooled.training_normalizations.size()==2&&pooled.densities.size()==3,"pooled order preserved");
  check(pooled.work.training_attempts==2&&pooled.work.joint_attempts==3,"training transformed once per pool");
  check(pooled.work.factor_attempts==0,"evaluation never refactors C");
  check(pooled.requests[0].training_vector_index==1&&pooled.requests[1].candidate_value_index==1,"request pairs preserved");
  check(pooled.densities[0].log_density.lower==pooled.densities[2].log_density.lower,"repeat request reproducible");
  const auto retained=owner.retained_payload_bound(),scratch=owner.evaluation_payload_bound(2,3);
  check(retained&&scratch,"actual payload receipts");
  auto elements=policy;elements.box.design.maximum_elements=owner.training_row_ids().size();
  check(evaluate(owner,elements).numerical_status==irred::numerics::Status::work_limit,"tightened cap includes appended full row");
  auto low=policy;low.box.design.maximum_payload_bytes=*retained+*scratch-1;
  check(evaluate(owner,low,{1.25,1.5,1.5,2},{3,3.25},{{1,0},{0,1},{1,0}}).status==DensityStatus::numerical_failure,"simultaneous retained+output cap");
  low=policy;low.maximum_evaluation_work_units=*owner.evaluation_work_bound(1,1,1)-1;
  check(evaluate(owner,low).numerical_status==irred::numerics::Status::work_limit,"one unit below work gate");
  low=policy;low.box.maximum_log_probability_width=1e-12;
  check(evaluate(owner,low).status==DensityStatus::finite,"tightened policy shape remains valid");
  auto invalid=evaluate(owner,policy,{1,2,3},{3});check(invalid.status==DensityStatus::invalid_input,"malformed coarse pool refused");
  invalid=evaluate(owner,policy,{1,2},{3},{{2,0}});check(invalid.status==DensityStatus::invalid_input,"request outside pool refused");
  invalid=evaluate(owner,policy,{1,std::numeric_limits<double>::infinity()},{3});check(invalid.status==DensityStatus::invalid_input,"nonfinite original training refused");
  auto relaxed=policy;relaxed.maximum_conditional_log_density_width*=2;
  check(evaluate(owner,relaxed).status==DensityStatus::invalid_input,"frozen width cannot relax");
  GaussianBoxHeldout moved(std::move(owner));check(moved.status()==DensityStatus::finite&&owner.status()==DensityStatus::invalid_input,"move invalidates source owner");
  check(owner.original_row_ids().empty()&&owner.training_row_ids().empty()&&owner.metadata().source_contract_id.empty(),"moved owner safe accessors");
  moved=std::move(moved);check(evaluate(moved,policy).densities[0].density_available,"self move no-op");
  // Caller string cannot disguise an actual precision preparation.
  in.rows.input_matrix_convention="covariance";
  auto pg=in.gaussian(MatrixKind::precision);check(pg.status()==DensityStatus::finite,"synthetic precision admitted generally");
  auto pp=in.policy(pg);auto refused=in.prepare(std::move(pg),pp);
  check(refused.status()==DensityStatus::incompatible_metadata&&pg.status()==DensityStatus::finite,"precision rejected with rollback");
  auto failed_move=std::move(refused);check(failed_move.status()==DensityStatus::incompatible_metadata&&refused.status()==DensityStatus::invalid_input,"failed receipt moves safely");
  // Training rank can fail despite full joint rank and a proper finite box.
  auto rank=in;rank.x={1,0,0,0,0,1};auto rg=rank.gaussian();auto rp=rank.policy(rg);auto rd=rank.prepare(std::move(rg),rp);
  check(rd.status()==DensityStatus::numerical_failure&&rg.status()==DensityStatus::finite,"training deficient algorithm refuses without dropping axis");
  // Dyadic near-cancellation: positive computed Schur is insufficient.
  auto narrow=in;narrow.covariance[8]=.3125+0x1p-50;auto ng=narrow.gaussian();auto np=narrow.policy(ng);auto nd=narrow.prepare(std::move(ng),np);
  check(nd.status()==DensityStatus::numerical_failure&&ng.status()==DensityStatus::finite,"unresolved Schur refuses with source intact");
  check(nd.preparation().cross_whitening_available&&nd.preparation().schur_available&&!nd.preparation().retained_frames_available,"failure retains earned diagnostic gates");
  auto capg=in.gaussian();auto cap=in.policy(capg);--cap.maximum_preparation_work_units;auto cd=in.prepare(std::move(capg),cap);
  check(cd.preparation().numerical_status==irred::numerics::Status::work_limit&&capg.status()==DensityStatus::finite,"preparation precharge rollback");
  // Independent zero response, tiny variance is still a likelihood row.
  auto zero=in;zero.covariance={1,0,0,0,1,0,0,0,0x1p-40};zero.x={1,0,0,1,0,0};
  auto zg=zero.gaussian();auto zp=zero.policy(zg);auto zo=zero.prepare(std::move(zg),zp);auto zr=evaluate(zo,zp);
  check(zr.densities[0].density_available,"tiny independent zero response retained");
  check(contains(zr.densities[0].log_density,20*std::log(2.)-.5*std::log(2*std::acos(-1.))),"tiny variance determinant enters ratio");
  // Interior row and a nonzero orthogonal training tail. A zero-response row
  // remains in both full and training likelihoods, whose tail must cancel.
  auto interior=in;interior.rows.ordered_ids={"r0","r1","r2","r3"};interior.heldout_index=1;
  interior.heldout.candidate_identity="withhold literal original r1/v1";interior.heldout.original_row_lineage="r0,r1,r2,r3 exact original order";
  interior.rows.ordering_provenance="literal r0,r1,r2,r3; withhold original r1";
  interior.covariance={1,.5,0,0,.5,1,.25,0,0,.25,1,0,0,0,0,1};
  interior.x={1,0,.5,.25,0,1,0,0};interior.offsets={1,3,2,4};
  auto ig=interior.gaussian();auto ip=interior.policy(ig);auto io=interior.prepare(std::move(ig),ip);
  const auto ir=evaluate(io,ip,{1.25,1.5,4.75},{3});
  check(io.training_row_ids().size()==3&&io.training_row_ids()[1]=="r2","interior source order extraction");
  near(ir.training_normalizations[0].minimum_quadratic,.5625);
  near(ir.densities[0].joint_normalization.minimum_quadratic,.5625);
  check(contains(ir.densities[0].log_density,target),"nonzero training tail cancels in predictive ratio");
  auto mapping=in;mapping.heldout.active_original_parameter_indices={2,0};auto mg=mapping.gaussian();auto mp=mapping.policy(mg);
  check(mapping.prepare(std::move(mg),mp).status()==DensityStatus::incompatible_metadata&&mg.status()==DensityStatus::finite,"malformed active mapping refuses with rollback");
  auto nonzero_fixed=in;nonzero_fixed.heldout.fixed_original_parameter_value=1;
  auto nfg=nonzero_fixed.gaussian();auto nfp=nonzero_fixed.policy(nfg);
  check(nonzero_fixed.prepare(std::move(nfg),nfp).status()==DensityStatus::incompatible_metadata&&nfg.status()==DensityStatus::finite,"v1 fixed coordinate is literal zero");
  auto provenance=in;provenance.box.fixed_coordinate_provenance.clear();provenance.heldout.fixed_coordinate_provenance.clear();
  auto fg=provenance.gaussian();auto fp=provenance.policy(fg);
  check(provenance.prepare(std::move(fg),fp).status()==DensityStatus::incompatible_metadata&&fg.status()==DensityStatus::finite,"absent point-mass provenance refuses");
  // Coherent change of all observation coordinates changes the density by its
  // scalar observation Jacobian, while the original parameter box stays fixed.
  auto units=in;for(auto&v:units.covariance)v*=4;for(auto&v:units.x)v*=2;for(auto&v:units.offsets)v*=2;
  units.design.residual_unit=units.heldout.heldout_unit="newu";units.heldout.heldout_covariance_unit="newu^2";units.heldout.heldout_measure="d(newu)";
  auto ug=units.gaussian();auto up=units.policy(ug);auto uo=units.prepare(std::move(ug),up);auto ur=evaluate(uo,up,{2.5,3},{6});
  check(contains(ur.densities[0].log_density,target-std::log(2.)),"native observation-coordinate Jacobian");
  std::cout<<"heldout source controls complete\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
