// Original exact Fraction covariance conditioning and independent prior/noise
// precision/KKT facts, frozen before production. Direct active-block integration
// uses prior-whitened coordinates with no posterior factor oracle.
// moments2e-12*(1+|ref|), density2e-11*(1+|ref|), refinement<=5% allocation.
#include "irred/calibration_predictive.hpp"
#include "calibration_predictive_peer_facts.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numbers>
#include <stdexcept>
namespace {
using namespace irred::calibration;namespace s=irred::statistics;namespace f=ladder_peer_facts;
using W=long double;unsigned checks=0;
void need(bool b,const char*w){++checks;if(!b)throw std::runtime_error(w);}
void near(W a,W b,W tol,const char*w){need(std::abs(a-b)<=tol*(1+std::abs(b)),w);}
template<class A>std::vector<double> doubles(const A&a){return {a.begin(),a.end()};}
Model model(bool future=false){Model m;m.ordered_host_ids={"A","B"};m.magnitude_convention="original common synthetic mag";m.metallicity_coordinate_identity="synthetic dex";m.distance_shape_identity="fixed synthetic shape";m.calibration_identity="one shared delta";m.dependence_identity="conditional synthetic independent prior/noise";m.conditional_covariance_identity="explicit independent conditional SPD noise";
 if(future)m.rows={{RowKind::cepheid,"future-cep","A","CepF",.5,1.25,.25,0},{RowKind::hubble_supernova,"future-SN","","SNF",.75,0,0,41}};
 else m.rows={{RowKind::anchor_modulus,"anchorA","A","",0,0,0,0},{RowKind::anchor_modulus,"anchorB","B","",0,0,0,0},
 {RowKind::cepheid,"CepA1","A","EA1",-.5,.5,-.25,0},{RowKind::cepheid,"CepA2","A","EA2",1,1.5,.5,0},
 {RowKind::cepheid,"CepB1","B","EB1",.25,.75,.25,0},{RowKind::cepheid,"CepB2","B","EB2",-.75,1.25,-.5,0},
 {RowKind::calibrator_supernova,"calA","A","SNCA",1,0,0,0},{RowKind::calibrator_supernova,"calB","B","SNCB",-.5,0,0,0},
 {RowKind::hubble_supernova,"H1","","SNH1",.25,0,0,40},{RowKind::hubble_supernova,"H2","","SNH2",1.5,0,0,42},{RowKind::calibration_measurement,"delta-measurement","","",1,0,0,0}};
 return m;}
s::Gaussian noise(const Model&m,std::span<const double>C){s::Metadata md;for(auto&r:m.rows)md.ordered_ids.push_back(r.row_id);md.measure="product d(mag)";md.source_semantics="synthetic controls";md.table_identity="original Fraction fixture";md.ordering_provenance="explicit exact axes";md.calibration_provenance=m.calibration_identity;md.dependence_provenance=m.dependence_identity;md.uncertainty_identity=m.conditional_covariance_identity;return s::prepare_gaussian(C,s::MatrixKind::covariance,md,1000000,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);}
s::ParameterPrior prior(const Model&m,std::vector<double>mean,std::vector<double>C){auto l=linearize(m);s::ParameterPrior p;p.ordered_parameter_ids=l.metadata.ordered_parameter_ids;p.parameter_units=l.metadata.parameter_units;p.shared_nuisance_ids=l.metadata.shared_nuisance_ids;p.mean=std::move(mean);p.covariance=std::move(C);p.prior_identity="original proper correlated Gaussian";p.design_identity=l.metadata.design_identity;p.residual_unit="mag";p.parameter_measure="product d(mag) d(mag per log10(day)) d(mag per dex)";p.dependence_identity=l.metadata.dependence_identity;p.noise_independence_declared=true;return p;}
s::PredictiveMetadata metadata(const LadderPosterior&p,const Model&m){auto l=linearize(m);s::PredictiveMetadata md;md.ordered_parameter_ids=p.posterior().prior().ordered_parameter_ids;md.parameter_units=p.posterior().prior().parameter_units;for(auto&u:md.parameter_units)md.response_units.push_back("mag/"+u);md.training_event_ids=p.linearization().predictive_event_ids;md.future_event_ids=l.predictive_event_ids;md.future_unit="mag";md.future_covariance_unit="mag^2";md.future_measure="product d(mag)";md.response_identity=l.metadata.design_identity;md.conditioning_identity="exact supplied training vector";md.conditional_noise_identity=m.conditional_covariance_identity;md.dependence_identity=l.metadata.dependence_identity;md.future_noise_independence_declared=true;md.noise_conditional_on_parameters_declared=true;return md;}
W normal2(std::array<W,2>y,std::array<W,2>mu,std::array<W,4>C){W det=C[0]*C[3]-C[1]*C[2],a=y[0]-mu[0],b=y[1]-mu[1];return std::exp(-.5L*((C[3]*a*a-(C[1]+C[2])*a*b+C[0]*b*b)/det))/(2*std::numbers::pi_v<W>*std::sqrt(det));}
void full(){auto m=model(),fm=model(true);auto g=noise(m,doubles(f::C));auto pr=prior(m,doubles(f::m),doubles(f::S));auto p=LadderPosterior::prepare(std::move(g),m,std::move(pr));need(p.status()==s::DensityStatus::finite,"full proper ladder");auto y=doubles(f::y);auto rows=p.linearization().ordered_row_ids;auto mu=p.condition(y,rows);need(mu.status==s::DensityStatus::finite,"full posterior mean");for(size_t i=0;i<8;++i)near(mu.value[i],f::mu[i],2e-12L,"rational covariance-route mean");for(size_t i=0;i<64;++i)near(p.posterior().covariance()[i],f::V[i],2e-12L,"rational independent KKT covariance");for(size_t i=0;i<88;++i)need(p.linearization().design[i]==f::X[i],"single compiled design transcription");
 auto h=p.h0_projection(y,rows);need(h.status==s::DensityStatus::finite,"posterior H0 projection");near(h.median_h0_km_s_Mpc,f::H0_median,2e-12L,"independent lognormal median");near(h.expectation_h0_km_s_Mpc,f::H0_expectation,2e-12L,"independent lognormal expectation");near(h.log_standard_deviation,f::H0_logsd,2e-12L,"independent lognormal spread");
 auto ng=noise(fm,doubles(f::R));auto md=metadata(p,fm);auto q=LadderPredictive::prepare(p,y,rows,ng,fm,md);need(q.status()==s::DensityStatus::finite,"joint ladder future");for(size_t i=0;i<2;++i)near(q.mean()[i],f::future_mean[i],2e-12L,"offset-restored future mean");for(size_t i=0;i<4;++i)near(q.covariance()[i],f::W[i],2e-12L,"shared parameter future covariance");
 std::array<double,2>fy{26.5,21.5};auto den=q.log_density(fy,q.linearization().ordered_row_ids);need(den.density.status==s::DensityStatus::finite,"joint future density");W ref=std::log(normal2({fy[0],fy[1]},{f::future_mean[0],f::future_mean[1]},{f::W[0],f::W[1],f::W[2],f::W[3]}));near(den.density.log_value,ref,2e-11L,"independent cofactor future density");
 auto noiseless=doubles(f::offsets);for(size_t i=0;i<11;++i)for(size_t j=0;j<8;++j)noiseless[i]+=double(f::X[i*8+j]*f::truth[j]);auto fixed=p.condition(noiseless,rows);need(fixed.status==s::DensityStatus::finite,"fixed-truth generating mean");near(fixed.value[6]-f::truth[6],f::fixed_bias[6],2e-12L,"proper posterior fixed-truth mean bias");
 constexpr W z95=1.95996398454005423552L;W width=z95*std::sqrt(f::V[54]),sigma=std::sqrt(f::fixed_sampling[54]),bias=f::fixed_bias[6];auto Phi=[](W x){return std::erfc(-x/std::sqrt(2.L))/2;};W coverage=Phi((width-bias)/sigma)-Phi((-width-bias)/sigma);need(coverage>0&&coverage<1&&std::abs(coverage-.95L)>1e-4L,"fixed-truth eta/H0 coverage differs from nominal posterior interval");
 need(std::abs(f::W[0]-f::fixed_prediction_covariance[0])>1e-4L&&std::abs(f::fixed_bias[6])>1e-4L,"fixed-truth bias and variance differ from posterior law");
 auto overlap=fm;overlap.rows[0].event_id=m.rows[2].event_id;need(LadderPredictive::prepare(p,y,rows,ng,overlap,metadata(p,overlap)).status()!=s::DensityStatus::finite,"same astronomical event cannot be heldout");
 auto untracked=fm;untracked.ordered_host_ids={"B","A"};need(LadderPredictive::prepare(p,y,rows,ng,untracked,metadata(p,untracked)).status()!=s::DensityStatus::finite,"host dictionary identity checked");
}
template<unsigned N>W integrated(){std::array<W,N>x{},w{};for(unsigned i=0;i<N;++i){W z=std::cos(std::numbers::pi_v<W>*(i+.75L)/(N+.5L));for(unsigned k=0;k<80;++k){W a=1,b=z;for(unsigned j=2;j<=N;++j){W n=((2*j-1)*z*b-(j-1)*a)/j;a=b;b=n;}W dp=N*(z*b-a)/(z*z-1),d=b/dp;z-=d;if(std::abs(d)<4*std::numeric_limits<W>::epsilon()){x[i]=z;w[i]=2/((1-z*z)*dp*dp);break;}}}
 W result=0;for(unsigned p=0;p<12;++p)for(unsigned i=0;i<N;++i){W z1=-12.L+2*p+1+x[i];for(unsigned q=0;q<12;++q)for(unsigned j=0;j<N;++j){W z2=-12.L+2*q+1+x[j],b1=30+std::sqrt(2.L)*z1,b2=.25L+.5L/std::sqrt(2.L)*z1+std::sqrt(.875L)*z2;W priorweight=std::exp(-(z1*z1+z2*z2)/2)/(2*std::numbers::pi_v<W>),train=std::exp(-std::pow(32-b1,2)/2)/std::sqrt(2*std::numbers::pi_v<W>);result+=w[i]*w[j]*priorweight*train*normal2({31,.5},{b1,b2},{.5,.125,.125,.75});}}
 // Independent marginal training N(30,3), not posterior normalization.
 return result/(std::exp(-2.L/3)/std::sqrt(6*std::numbers::pi_v<W>));}
template<unsigned N>W ellipse_mass(const LadderPredictive &q){
 std::array<W,N>x{},w{};
 for(unsigned i=0;i<N;++i){W z=std::cos(std::numbers::pi_v<W>*(i+.75L)/(N+.5L));for(unsigned it=0;it<80;++it){W a=1,b=z;for(unsigned j=2;j<=N;++j){W t=((2*j-1)*z*b-(j-1)*a)/j;a=b;b=t;}W dp=N*(z*b-a)/(z*z-1),d=b/dp;z-=d;if(std::abs(d)<4*std::numeric_limits<W>::epsilon()){x[i]=z;w[i]=2/((1-z*z)*dp*dp);break;}}}
 // Exact rational active-block mean/covariance from independent conditioning.
 constexpr W m1=94.L/3,m2=7.L/12,C11=7.L/6,C12=7.L/24,C22=5.L/3;
 W L11=std::sqrt(C11),L21=C12/L11,L22=std::sqrt(C22-L21*L21),jac=L11*L22,mass=0;
 for(unsigned i=0;i<N;++i){W radius=1.5L*(x[i]+1);for(unsigned j=0;j<16;++j){W angle=2*std::numbers::pi_v<W>*j/16,u=radius*std::cos(angle),v=radius*std::sin(angle);std::array<double,2> y{double(m1+L11*u),double(m2+L21*u+L22*v)};auto d=q.log_density(y,q.linearization().ordered_row_ids);need(d.density.status==s::DensityStatus::finite,"joint ellipse density admitted");mass+=w[i]*1.5L*2*std::numbers::pi_v<W>/16*jac*radius*std::exp(W(d.density.log_value));}}
 return mass;
}
void integration(){auto m=model();m.rows={{RowKind::anchor_modulus,"training-only","A","",0,0,0,0}};auto fm=m;fm.rows={{RowKind::anchor_modulus,"future-A","A","",0,0,0,0},{RowKind::calibration_measurement,"future-delta","","",1,0,0,0}};
 std::vector<double>mean{30,31,-4,-3,.25,-19,.5,.25},C(64);for(size_t i=0;i<8;++i)C[i*8+i]=1;C[0]=2;C[7]=C[56]=.5;auto g=noise(m,std::array<double,1>{1});auto p=LadderPosterior::prepare(std::move(g),m,prior(m,mean,C));need(p.status()==s::DensityStatus::finite,"rank deficient proper active block");auto ng=noise(fm,std::array<double,4>{.5,.125,.125,.75});std::array<double,1>y{32};auto q=LadderPredictive::prepare(p,y,p.linearization().ordered_row_ids,ng,fm,metadata(p,fm));need(q.status()==s::DensityStatus::finite,"active block future");std::array<double,2>fy{31,.5};auto r=q.log_density(fy,q.linearization().ordered_row_ids);need(r.density.status==s::DensityStatus::finite,"active block density");W a=integrated<16>(),b=integrated<24>();need(a>0&&b>0&&std::isfinite(a)&&std::isfinite(b),"positive finite independent integration");near(a,b,1e-12L,"direct prior/noise integration refinement <=5%");std::cerr<<std::setprecision(20)<<"active integration coarse="<<a<<" fine="<<b<<" logreference="<<std::log(b)<<" native="<<r.density.log_value<<"\n";near(r.density.log_value,std::log(b),2e-11L,"independent integrated prior*training*future density");W mass16=ellipse_mass<16>(q),mass24=ellipse_mass<24>(q);near(mass16,mass24,1e-12L,"joint ellipse refinement <=5%");near(mass24,1-std::exp(-4.5L),2e-11L,"prior-predictive conditional chi-square2 ellipse mass");}
}
int main(){try{full();integration();std::cout<<"PASS "<<checks<<" independent ladder predictive controls\n";}catch(const std::exception&e){std::cerr<<"FAIL "<<checks<<": "<<e.what()<<'\n';return 1;}}
