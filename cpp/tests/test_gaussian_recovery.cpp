#include "gaussian_recovery_controls.hpp"
namespace {
using namespace recovery_controls;namespace f=recovery_peer_facts;
void run(){std::vector<std::string>rows{"training0","training1"},future_rows{"future0","future1"},parameter{"offset","slope"},units{"mag","mag"};
 const std::vector<double>C{1,.25,.25,2},S{2,.5,.5,1},X{1,0,1,1},m{.5,-.25},truth{1.5,-1.25},A{1,2,1,-1},R{.5,.125,.125,.75};
 s::ParameterPrior prior{parameter,units,{"offset"},m,S,"frozen proper prior","fixed offset slope design","mag","product d(mag)","independent original prior and training noise",true};
 auto noise=gaussian(C,metadata(rows));auto p=s::GaussianPosterior::prepare(std::move(noise),X,rows,prior);
 need(p.status()==s::DensityStatus::finite,"generic posterior preparation");auto ng=gaussian(R,metadata(future_rows));auto pg=gaussian(S,metadata(parameter));
 for(size_t i=0;i<4;++i)near(p.covariance()[i],f::generic_V[i]);
 s::PredictiveMetadata md;md.ordered_parameter_ids=parameter;md.parameter_units=units;md.response_units={"mag/mag","mag/mag"};
 md.training_event_ids={"event:train0","event:train1"};md.future_event_ids={"event:future0","event:future1"};md.future_unit="mag";md.future_covariance_unit="mag^2";md.future_measure="product d(mag)";
 md.response_identity="fixed heldout response";md.conditioning_identity="exact generated training vector";md.conditional_noise_identity="explicit independent future noise";md.dependence_identity="explicit independent original prior/training/future noise";md.future_noise_independence_declared=true;md.noise_conditional_on_parameters_declared=true;
 Reference ref{wide(f::generic_K),wide(f::generic_V),wide(f::generic_fixed_sampling),wide(f::generic_fixed_bias),wide(f::generic_W),wide(f::generic_fixed_future_sampling),wide(f::generic_fixed_future_bias),wide(f::generic_posterior_coverage),wide(f::generic_future_coverage)};
 auto condition=[&](const auto&y){return p.condition(y,rows);};auto density=[&](const auto&y,const auto&beta){return p.log_density(y,rows,beta,parameter);};
 auto predict=[&](const auto&y,const auto&fy){auto q=s::GaussianPredictive::prepare(p,y,rows,ng,A,md);FutureEvaluation r;r.status=q.status();r.numerical_status=q.numerical_status();r.work=q.preparation_work_units();if(r.status!=s::DensityStatus::finite)return r;
  for(size_t i=0;i<4;++i){near(q.covariance()[i],ref.Wcov[i]);}r.mean.assign(q.mean().begin(),q.mean().end());auto d=q.log_density(fy,future_rows);r.status=d.density.status;r.numerical_status=d.density.numerical_status;r.quadratic=d.quadratic;return r;};
 generation_benchmark(p.source(),generating(p.source(),{0,0},{"mag","mag"},"conditional training noise"));
 campaign("generic offset slope",p.source(),pg,ng,X,A,{0,0},{0,0},m,truth,units,ref,0x210000,condition,density,predict,[](const auto&,const auto&,const auto&){return ProjectionEvaluation{};});
}
}
int main(){try{run();std::cout<<"{\"kind\":\"campaign_completion\",\"status\":\"pass\",\"Philox_independence_proof\":false,\"observational_qualification\":false}\n";}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
