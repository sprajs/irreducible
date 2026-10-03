#define IRRED_RECOVERY_CAPTURE
#include "gaussian_recovery_controls.hpp"
namespace {
using namespace recovery_controls;namespace f=recovery_peer_facts;
void run(){std::vector<std::string>rows{"training0","training1"},future_rows{"future0","future1"},parameter{"offset","slope"},units{"mag","mag"};
 const std::vector<double>C{1,.25,.25,2},S{2,.5,.5,1},X{1,0,1,1},m{.5,-.25},truth{1.5,-1.25},A{1,2,1,-1},R{.5,.125,.125,.75};
 s::ParameterPrior prior{parameter,units,{"offset"},m,S,"frozen proper prior","fixed offset slope design","mag","product d(mag)","independent original prior and training noise",true};
 auto noise=gaussian(C,metadata(rows));auto original_p=s::GaussianPosterior::prepare(std::move(noise),X,rows,prior);
 need(original_p.status()==s::DensityStatus::finite,"generic posterior preparation");auto ng=gaussian(R,metadata(future_rows));auto pg=gaussian(S,metadata(parameter));
 for(size_t i=0;i<4;++i)near(original_p.covariance()[i],f::generic_V[i]);
 s::PredictiveMetadata md;md.ordered_parameter_ids=parameter;md.parameter_units=units;md.response_units={"mag/mag","mag/mag"};
 md.training_event_ids={"event:train0","event:train1"};md.future_event_ids={"event:future0","event:future1"};md.future_unit="mag";md.future_covariance_unit="mag^2";md.future_measure="product d(mag)";
 md.response_identity="fixed heldout response";md.conditioning_identity="exact generated training vector";md.conditional_noise_identity="explicit independent future noise";md.dependence_identity="explicit independent original prior/training/future noise";md.future_noise_independence_declared=true;md.noise_conditional_on_parameters_declared=true;
 Reference ref{wide(f::generic_K),wide(f::generic_V),wide(f::generic_fixed_sampling),wide(f::generic_fixed_bias),wide(f::generic_W),wide(f::generic_fixed_future_sampling),wide(f::generic_fixed_future_bias),wide(f::generic_posterior_coverage),wide(f::generic_future_coverage)};
 allocation_observation::calls=0;allocation_observation::active=true;
 const auto setup_begin=std::chrono::steady_clock::now();
 auto retained=s::GaussianPredictiveConditioning::prepare(std::move(original_p),ng,A,md);
 const auto setup_seconds=std::chrono::duration<W>(std::chrono::steady_clock::now()-setup_begin).count();
 allocation_observation::active=false;const auto setup_allocations=allocation_observation::calls;
 need(retained.status()==s::DensityStatus::finite,"generic invariant predictive preparation");
 const auto&p=retained.posterior();
 for(size_t i=0;i<4;++i)near(retained.covariance()[i],ref.Wcov[i]);
 std::cout<<std::setprecision(17)<<"{\"kind\":\"retained_predictive_setup\",\"campaign\":\"generic offset slope\",\"seconds\":"<<setup_seconds<<",\"observed_allocations\":"<<setup_allocations<<",\"work_units\":"<<retained.preparation_work_units()<<",\"retained_payload\":"<<*retained.retained_payload_bound()<<",\"maximum_batch_payload\":"<<*retained.batch_payload_bound(chunk)<<",\"same_owner_for_both_ensembles\":true,\"arithmetic\":\"longdouble_cpu_v1\"}\n";
 // Original inputs associated with a generating refusal never enter either arm.
 const std::vector<double>upstream_y{1,2,0,0,3,4},upstream_future{5,6,0,0,7,8};
 const std::array<n::Status,3>upstream_status{n::Status::ok,n::Status::work_limit,n::Status::ok};
 const auto available=available_predictive_inputs(upstream_y,upstream_future,upstream_status,2,2);
 need(available.count==2&&available.position==std::vector<size_t>{0,SIZE_MAX,1}&&available.training==std::vector<double>{1,2,3,4}&&available.future==std::vector<double>{5,6,7,8},"upstream generating refusal retains order and excludes unavailable vectors");
 const auto upstream_batch=retained.evaluate(available.training,rows,available.future,future_rows,available.count);
 need(upstream_batch.rows.size()==2&&upstream_batch.rows[0].status==s::DensityStatus::finite,"no artificial zero vector executed for upstream refusal");
 auto batch_predict=[&](const auto&y,const auto&fy,size_t count){return retained.evaluate(y,rows,fy,future_rows,count);};
 auto condition=[&](const auto&y){return p.condition(y,rows);};auto density=[&](const auto&y,const auto&beta){return p.log_density(y,rows,beta,parameter);};
 auto predict=[&](const auto&y,const auto&fy){auto q=s::GaussianPredictive::prepare(p,y,rows,ng,A,md);FutureEvaluation r;r.status=q.status();r.numerical_status=q.numerical_status();r.work=q.preparation_work_units()+32*future_rows.size()*future_rows.size();if(r.status!=s::DensityStatus::finite)return r;
  r.mean.assign(q.mean().begin(),q.mean().end());r.mean_errors.assign(q.mean_absolute_error_estimates().begin(),q.mean_absolute_error_estimates().end());auto d=q.log_density(fy,future_rows);r.status=d.density.status;r.numerical_status=d.density.numerical_status;r.quadratic=d.quadratic;r.density=d;if(r.status!=s::DensityStatus::finite){r.mean.clear();r.mean_errors.clear();}return r;};
 generation_benchmark(p.source(),generating(p.source(),{0,0},{"mag","mag"},"conditional training noise"));
 campaign("generic offset slope",p.source(),pg,ng,X,A,{0,0},{0,0},m,truth,units,ref,0x210000,condition,density,predict,batch_predict,[](const auto&,const auto&,const auto&){return ProjectionEvaluation{};});
}
}
int main(){try{run();std::cout<<"{\"kind\":\"campaign_completion\",\"status\":\"pass\",\"Philox_independence_proof\":false,\"observational_qualification\":false}\n";}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
