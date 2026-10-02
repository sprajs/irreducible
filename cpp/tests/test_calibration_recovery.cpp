#include "gaussian_recovery_controls.hpp"
#include "irred/calibration_predictive.hpp"
#include "calibration_predictive_peer_facts.hpp"
namespace {
using namespace recovery_controls;using namespace irred::calibration;namespace f=recovery_peer_facts;
Model model(bool future=false){Model m;m.ordered_host_ids={"A","B"};m.magnitude_convention="original common synthetic mag";m.metallicity_coordinate_identity="synthetic dex";m.distance_shape_identity="fixed synthetic shape";m.calibration_identity="one shared delta";m.dependence_identity="explicit independent original prior/training/future noise";m.conditional_covariance_identity="explicit conditional SPD noise";
 if(future)m.rows={{RowKind::cepheid,"future-cep","A","CepF",.5,1.25,.25,0},{RowKind::hubble_supernova,"future-SN","","SNF",.75,0,0,41}};
 else m.rows={{RowKind::anchor_modulus,"anchorA","A","",0,0,0,0},{RowKind::anchor_modulus,"anchorB","B","",0,0,0,0},
 {RowKind::cepheid,"CepA1","A","EA1",-.5,.5,-.25,0},{RowKind::cepheid,"CepA2","A","EA2",1,1.5,.5,0},
 {RowKind::cepheid,"CepB1","B","EB1",.25,.75,.25,0},{RowKind::cepheid,"CepB2","B","EB2",-.75,1.25,-.5,0},
 {RowKind::calibrator_supernova,"calA","A","SNCA",1,0,0,0},{RowKind::calibrator_supernova,"calB","B","SNCB",-.5,0,0,0},
 {RowKind::hubble_supernova,"H1","","SNH1",.25,0,0,40},{RowKind::hubble_supernova,"H2","","SNH2",1.5,0,0,42},{RowKind::calibration_measurement,"delta-measurement","","",1,0,0,0}};
 return m;}
s::Gaussian noise(const Model&m,const std::vector<double>&C){std::vector<std::string>rows;for(auto&r:m.rows)rows.push_back(r.row_id);return gaussian(C,metadata(rows));}
s::ParameterPrior prior(const Linearization&l,std::vector<double>mean,std::vector<double>C){return {l.metadata.ordered_parameter_ids,l.metadata.parameter_units,l.metadata.shared_nuisance_ids,std::move(mean),std::move(C),"frozen narrow fully correlated proper prior",l.metadata.design_identity,"mag","product d(mag) d(mag per log10(day)) d(mag per dex)",l.metadata.dependence_identity,true};}
s::PredictiveMetadata predictive_metadata(const LadderPosterior&p,const Linearization&future){s::PredictiveMetadata md;md.ordered_parameter_ids=p.posterior().prior().ordered_parameter_ids;md.parameter_units=p.posterior().prior().parameter_units;
 for(auto&u:md.parameter_units){md.response_units.push_back("mag/"+u);}md.training_event_ids=p.linearization().predictive_event_ids;md.future_event_ids=future.predictive_event_ids;
 md.future_unit="mag";md.future_covariance_unit="mag^2";md.future_measure="product d(mag)";md.response_identity=future.metadata.design_identity;md.conditioning_identity="exact generated training vector";md.conditional_noise_identity="explicit conditional SPD noise";md.dependence_identity=future.metadata.dependence_identity;md.future_noise_independence_declared=true;md.noise_conditional_on_parameters_declared=true;return md;}
std::vector<double> factor(size_t d,double diagonal){std::vector<double>L(d*d);for(size_t i=0;i<d;++i){L[i*d+i]=diagonal;if(i)L[i*d]=(i%2?-1.:1.)/128;}return L;}
std::vector<double> covariance(const std::vector<double>&L,size_t d){std::vector<double>C(d*d);for(size_t i=0;i<d;++i)for(size_t j=0;j<d;++j)for(size_t k=0;k<d;++k)C[i*d+j]+=L[i*d+k]*L[j*d+k];return C;}
void support(const std::vector<double>&design,const std::vector<double>&offset,const std::vector<double>&m,
 const std::vector<double>&LS,const std::vector<double>&LN){const size_t p=m.size(),rows=offset.size();for(size_t i=0;i<rows;++i)if(offset[i]!=0){
 W mu=offset[i],width=0;for(size_t j=0;j<p;++j)mu+=W(design[i*p+j])*m[j];
 for(size_t k=0;k<p;++k){W col=0;for(size_t j=0;j<p;++j)col+=W(design[i*p+j])*LS[j*p+k];width+=7*std::abs(col);}
 for(size_t k=0;k<rows;++k)width+=7*std::abs(LN[i*rows+k]);
 need(mu-width>16+1e-8L&&mu+width<32-1e-8L,"all correlated finite prior/noise extrema in exact-offset domain");}}
void refusal_probe(){namespace b=ladder_peer_facts;auto m=model(),fm=model(true);auto l=linearize(m),fl=linearize(fm);auto C=doubles(b::C),S=doubles(b::S),R=doubles(b::R),mean=doubles(b::m);
 auto g=noise(m,C);auto pr=prior(l,mean,S);pr.prior_identity="unchanged broad rational prior failure challenge";auto p=LadderPosterior::prepare(std::move(g),m,pr);need(p.status()==s::DensityStatus::finite,"broad immutable fixture preparation");
 auto pg=gaussian(S,metadata(pr.ordered_parameter_ids,pr.parameter_measure)),ng=noise(fm,R);
 auto gm=generating(pg,mean,pr.parameter_units,"unchanged broad prior"),gt=generating(p.posterior().source(),std::vector<double>(11),std::vector<std::string>(11,"mag"),"unchanged broad training noise"),gf=generating(ng,std::vector<double>(2),std::vector<std::string>(2,"mag"),"unchanged broad future noise");
 size_t refused=0,actually_attempted=0;for(auto seed:seeds){reserve_range(seed,0x130000,0,256,8);reserve_range(seed,0x130100,0,256,11);reserve_range(seed,0x130200,0,256,2);
 auto bp=draw(pg,gm,seed,0x130000,0,256),bt=draw(p.posterior().source(),gt,seed,0x130100,0,256),bf=draw(ng,gf,seed,0x130200,0,256);
 for(size_t i=0;i<256;++i){++actually_attempted;std::vector<double>beta,y,fy;
 std::string_view stage="broad_parameter_copy";s::DensityStatus status=s::DensityStatus::numerical_failure;n::Status numerical=n::Status::work_limit;
 try{if(bp.rows[i].status==n::Status::ok)beta.assign(bp.values.begin()+i*8,bp.values.begin()+(i+1)*8);stage="broad_generation";
 numerical=bp.rows[i].status!=n::Status::ok?bp.rows[i].status:bt.rows[i].status!=n::Status::ok?bt.rows[i].status:bf.rows[i].status;
 if(numerical==n::Status::ok){stage="broad_training_observation";y=observation(l.design,l.offsets_mag,beta,std::span(bt.values.data()+i*11,11));
 stage="broad_future_observation";fy=observation(fl.design,fl.offsets_mag,beta,std::span(bf.values.data()+i*2,2));const auto ys=observation_status(y),fs=observation_status(fy);
 if(ys!=n::Status::ok||fs!=n::Status::ok){stage=ys!=n::Status::ok?"broad_training_observation":"broad_future_observation";numerical=ys!=n::Status::ok?ys:fs;}
 else {stage="broad_exact_offset_conditioning";auto c=p.condition(y,l.ordered_row_ids);status=c.status;numerical=c.numerical_status;
 if(c.status==s::DensityStatus::finite){std::cout<<std::setprecision(17)<<"{\"kind\":\"broad_successful_attempt\",\"seed\":\""<<seed<<"\",\"sample\":"<<i<<",\"stream_family\":\"1245184\",\"beta\":";vector_json(beta);
 std::cout<<",\"training\":";vector_json(y);std::cout<<",\"future\":";vector_json(fy);std::cout<<",\"posterior_mean\":";vector_json(c.value);std::cout<<",\"status\":"<<int(c.status)<<",\"numerical_status\":"<<int(c.numerical_status)<<",\"coverage_claim\":false}\n";stage={};}}}
 }catch(const std::bad_alloc&){status=s::DensityStatus::numerical_failure;numerical=n::Status::work_limit;}
 catch(const std::exception&){++refused;failed("unchanged broad rational refusal probe",seed,i,stage,beta,y,fy,bp,bt,bf,i,status,numerical);throw;}
 if(!stage.empty()){++refused;failed("unchanged broad rational refusal probe",seed,i,stage,beta,y,fy,bp,bt,bf,i,status,numerical);}}}
 std::cout<<"{\"kind\":\"broad_admission_inventory\",\"planned_count\":512,\"actually_attempted\":"<<actually_attempted<<",\"attempts\":512,\"refused\":"<<refused<<",\"coverage_claim\":false}\n";
 // Separately declared deterministic cell, never presented as a random draw.
 auto counterexample=doubles(b::y);counterexample[8]=std::nextafter(16.,0.);
 const auto failed_mean=p.condition(counterexample,l.ordered_row_ids);
 const W difference=W(counterexample[8])-l.offsets_mag[8];const double rounded=double(difference);
 std::cout<<std::setprecision(21)<<"{\"kind\":\"deterministic_failed_attempt\",\"operation\":\"Ladder-exact-offset-domain/v1\",\"random_draw\":false,\"attempts\":1,\"refused\":1,\"hits\":0,\"coverage_bounds\":[0,1],\"coverage_claim\":false,\"stage\":\"offset_translation\",\"status\":"<<int(failed_mean.status)<<",\"numerical_status\":"<<int(failed_mean.numerical_status)<<",\"original_y\":";
 vector_json(counterexample);std::cout<<",\"offset\":"<<l.offsets_mag[8]<<",\"wide_difference\":"<<difference<<",\"rounded_difference\":"<<rounded<<",\"remainder\":"<<difference-W(rounded)<<"}\n";
 need(failed_mean.status!=s::DensityStatus::finite&&failed_mean.numerical_status==n::Status::conditioning_budget_exceeded&&difference!=W(rounded),"known rounded offset refusal preserved");
}
void run(){auto m=model(),fm=model(true);auto l=linearize(m),fl=linearize(fm);need(l.status==s::DensityStatus::finite&&fl.status==s::DensityStatus::finite,"compiled row equations");
 for(size_t i=0;i<l.design.size();++i)need(l.design[i]==ladder_peer_facts::X[i],"frozen training equation transcription");
 for(size_t i=0;i<fl.design.size();++i)need(fl.design[i]==ladder_peer_facts::A[i],"frozen future equation transcription");
 std::vector<double>mean{30,31,-4,-3,.25,-19,.25,.125},truth=mean;truth[6]=.5;
 auto LS=factor(8,1./16),LC=factor(11,1./8);std::vector<double>LR{.125,0,.03125,.125};
 auto S=covariance(LS,8),C=covariance(LC,11),R=covariance(LR,2);
 support(l.design,l.offsets_mag,mean,LS,LC);support(fl.design,fl.offsets_mag,mean,LS,LR);
 support(l.design,l.offsets_mag,truth,std::vector<double>(64),LC);support(fl.design,fl.offsets_mag,truth,std::vector<double>(64),LR);
 auto g=noise(m,C);auto pr=prior(l,mean,S);auto p=LadderPosterior::prepare(std::move(g),m,pr);need(p.status()==s::DensityStatus::finite,"narrow proper ladder preparation");
 auto ng=noise(fm,R),pg=gaussian(S,metadata(pr.ordered_parameter_ids,pr.parameter_measure));auto md=predictive_metadata(p,fl);
 for(size_t i=0;i<64;++i)near(p.posterior().covariance()[i],f::ladder_V[i]);
 Reference ref{wide(f::ladder_K),wide(f::ladder_V),wide(f::ladder_fixed_sampling),wide(f::ladder_fixed_bias),wide(f::ladder_W),wide(f::ladder_fixed_future_sampling),wide(f::ladder_fixed_future_bias),wide(f::ladder_posterior_coverage),wide(f::ladder_future_coverage)};
 auto condition=[&](const auto&y){return p.condition(y,l.ordered_row_ids);};auto density=[&](const auto&y,const auto&beta){return p.log_density(y,l.ordered_row_ids,beta,pr.ordered_parameter_ids);};
 auto predict=[&](const auto&y,const auto&fy){auto q=LadderPredictive::prepare(p,y,l.ordered_row_ids,ng,fm,md);FutureEvaluation r;r.status=q.status();r.numerical_status=q.numerical_status();if(r.status!=s::DensityStatus::finite)return r;
 r.work=q.predictive().preparation_work_units();for(size_t i=0;i<4;++i)near(q.covariance()[i],ref.Wcov[i]);r.mean.assign(q.mean().begin(),q.mean().end());auto d=q.log_density(fy,fl.ordered_row_ids);r.status=d.density.status;r.numerical_status=d.density.numerical_status;r.quadratic=d.quadratic;return r;};
 auto projection=[&](const auto&beta,const auto&c,const auto&y){auto h=p.h0_projection(y,l.ordered_row_ids);if(h.status!=s::DensityStatus::finite)return ProjectionEvaluation{h.status,h.numerical_status,false};near(h.eta_mean,c.value[6]);near(h.eta_variance,ref.V[54]);
 const auto width=z95*std::sqrt(ref.V[54]);const auto low=project_h0(m,double(W(c.value[6])-width)),high=project_h0(m,double(W(c.value[6])+width)),actual=project_h0(m,beta[6]);
 if(!std::isfinite(low)||!std::isfinite(high)||!std::isfinite(actual))return ProjectionEvaluation{s::DensityStatus::numerical_failure,n::Status::overflow,false};
 return ProjectionEvaluation{s::DensityStatus::finite,n::Status::ok,(actual>=low&&actual<=high)==(std::abs(W(c.value[6])-beta[6])<=width)};};
 campaign("narrow full correlated proper ladder",p.posterior().source(),pg,ng,l.design,fl.design,l.offsets_mag,fl.offsets_mag,mean,truth,pr.parameter_units,ref,0x110000,condition,density,predict,projection,true);
 refusal_probe();
}
}
int main(){try{run();std::cout<<"{\"kind\":\"campaign_completion\",\"status\":\"pass\",\"observational_qualification\":false}\n";}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
