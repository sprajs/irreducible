#include "cosmology/numerics.hpp"
#include "fixtures/foundations_oracles.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using namespace cosmology::numerics;
namespace {
int checks=0;double maximum_error=0;
std::array<double,3> integral_errors{};
std::array<std::size_t,3> integral_evaluations{};
void check(bool ok,const char* s){++checks;if(!ok)throw std::runtime_error(s);}
double oracle(std::string_view id){for(auto v:cosmology::test_fixtures::foundations_oracles)if(v.id==id)return v.rounded;throw std::runtime_error("missing oracle");}
void within(double x,double y,double absolute,double relative,const char* s){auto e=std::abs(x-y);maximum_error=std::max(maximum_error,e);check(e<=absolute+relative*std::abs(y),s);}
double exponential(double x,const void*){return std::exp(x);}
double reciprocal(double x,const void*){return 1/(1+x);}
double cubic(double x,const void*){return x*x*x;}
double toy(double x,const void* q){return std::pow(1+x,-1-*static_cast<const double*>(q));}
double peak(double x,const void*){const auto u=(x-.123456)/1e-6;return std::exp(-u*u);}
double nonfinite(double,const void*){return std::numeric_limits<double>::infinity();}
double cancellation_polynomial(double x,const void*){return x*x*x*x*x*x-(1./7.);}
}
int main(){try{
 std::array<double,3> cancellation{1e16,1,-1e16};auto s=compensated_sum(cancellation);check(s.status==Status::ok&&s.value==1,"exact cancellation");
 check(compensated_sum({}).value==0,"empty sum");
 std::array<double,2> overflow{std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};check(compensated_sum(overflow).status==Status::overflow,"sum overflow");
 std::array<double,1> nan{std::numeric_limits<double>::quiet_NaN()};check(compensated_sum(nan).status==Status::nonfinite_input,"sum NaN");
 std::array<double,2> logs{1000,1000};within(log_sum_exp(logs).value,1000+std::log(2.),1e-12,0,"stable logsumexp");check(log_sum_exp({}).status==Status::invalid_input,"empty logsumexp");
 for(auto xs:{"-0.999999999999","-1e-12","1e-12","0.5","10"}){
  // Decimal source representation converted once: oracle is for that decimal input.
  const double x=std::stod(xs);
  // Near -1 reference must account for source binary64 encoding separately.
  if(x>-.9)within(log1p_checked(x).value,oracle(std::string("log1p_")+xs),1e-27,3e-15,"log1p independent Decimal");
  within(expm1_checked(x).value,oracle(std::string("expm1_")+xs),1e-27,3e-15,"expm1 independent Decimal");
 }
 check(log1p_checked(-1).status==Status::outside_domain,"log boundary");check(expm1_checked(1000).status==Status::overflow,"exp overflow");
 within(log_gamma_positive(.5).value,oracle("log_gamma_half"),1e-14,0,"gamma half identity");within(log_gamma_positive(5).value,oracle("log_gamma_5"),1e-14,0,"gamma factorial identity");check(log_gamma_positive(0).status==Status::outside_domain,"gamma domain");
 std::array<double,3> x{0,1,3},y{2,4,8};within(interpolate_linear(x,y,2).value,6,0,0,"analytic linear interpolation");
 check(interpolate_linear(x,y,3).value==8,"last knot");check(interpolate_linear(x,y,-1).status==Status::outside_domain,"no extrapolation");std::array<double,3> duplicate{0,0,3};check(interpolate_linear(duplicate,y,1).status==Status::invalid_input,"duplicate node");
 for(double tol:{1e-8,1e-10,1e-12}){
  IntegrationPolicy p{tol,tol,100000,40};
  const auto er=integrate(exponential,nullptr,0,1,p);check(er.status==Status::ok,"exp integration complete");if(tol==1e-12)within(er.value,oracle("integral_exp_0_1"),1e-12,1e-10,"exp final consumer budget");
  const auto rr=integrate(reciprocal,nullptr,0,1,p);check(rr.status==Status::ok,"rational integration complete");if(tol==1e-12)within(rr.value,oracle("integral_reciprocal_0_1"),1e-12,1e-10,"rational final consumer budget");
  within(integrate(cubic,nullptr,0,1,p).value,.25,1e-12,0,"cubic exact integral");
 }
 for(auto zs:{"0.001","0.5","3"})for(auto qs:{"-1","-0.5","-1e-10","0","1e-10","0.5","1"}){
  const double z=std::stod(zs),q=std::stod(qs);std::size_t level=0;
  for(double tol:{1e-8,1e-10,1e-12}){
   const auto r=integrate(toy,&q,0,z,{tol,tol,100000,40});check(r.status==Status::ok,"toy quadrature completed");
   const double truth=oracle(std::string("toy_integral_z")+zs+"_q"+qs);
   if(tol==1e-12)within(r.value,truth,1e-12,1e-10,"toy final consumer budget");
   integral_errors[level]=std::max(integral_errors[level],std::abs(r.value-truth));
   integral_evaluations[level]+=r.evaluations;
   ++level;
  }
 }
 check(integrate(exponential,nullptr,0,1,{1e-16,1e-16,5,40}).status==Status::work_limit,"quadrature maxwork");check(integrate(nonfinite,nullptr,0,1,{1e-12,1e-12,100,10}).status==Status::nonfinite_input,"nonfinite callback");
 const auto cancellation_integral=integrate(cancellation_polynomial,nullptr,0,1,{1e-12,.1,10000,30});
 check(cancellation_integral.status==Status::ok,"cancelled integral completed");
 const double cancelled_truth=static_cast<double>(1.L/7-static_cast<long double>(1./7.));
 within(cancellation_integral.value,cancelled_truth,1e-12,.1,"final integral relative scale regression");
 check(cancellation_integral.error_estimate<=1e-12+.1*std::abs(cancellation_integral.value),"final estimator budget regression");
 const auto missed=integrate(peak,nullptr,0,1,{1e-12,1e-12,100000,40});
 check(missed.status==Status::ok&&std::abs(missed.value-1.772453850905516e-6)>1e-7,"documented empirical missed peak witness");
 std::array<double,4> matrix{4,1,1,9};const auto f=cholesky(matrix,2,256);check(f.status()==Status::ok,"SPD factor");within(f.log_determinant(),oracle("spd_2x2_logdet"),4e-11,0,"SPD logdet budget");
 std::array<double,2> rhs{2,-3};const auto solution=solve(f,rhs,1e-10);check(solution.status==Status::ok,"SPD solve budget");within(solution.value[0],.6,1e-14,0,"exact solution x");within(solution.value[1],-.4,1e-14,0,"exact solution y");
 const double quadratic=rhs[0]*solution.value[0]+rhs[1]*solution.value[1];within(quadratic,oracle("spd_2x2_quadratic"),4e-11,0,"quadratic budget");
 const double logdensity=-(2*std::log(2*std::acos(-1.))+f.log_determinant()+quadratic)/2;within(logdensity,oracle("spd_2x2_gaussian_logdensity"),1e-10,0,"assembled Gaussian fixture");check(solution.backward_residual<1e-14,"separate residual");
 std::array<double,4> singular{1,1,1,1},nonspd{1,2,2,1},asym{4,1,2,9},near_singular{1,0,0,1e-16};
 check(cholesky(singular,2,256).status()==Status::not_positive_definite,"singular no jitter");check(cholesky(nonspd,2,256).status()==Status::not_positive_definite,"nonSPD no jitter");check(cholesky(asym,2,256).status()==Status::invalid_input,"exact symmetry");check(cholesky(matrix,2,3).status()==Status::work_limit,"matrix resource bound");
 check(solve(cholesky(near_singular,2,256),rhs,1e-10).status==Status::conditioning_budget_exceeded,"forward sensitivity rejects low residual");
 std::cout<<"numerical checks="<<checks<<" maximum_observed_absolute_fixture_error="<<maximum_error<<" missed_peak_limitation=observed\n";
 for(std::size_t i=0;i<3;++i)std::cout<<"toy_refinement_level="<<i<<" maximum_absolute_error="<<integral_errors[i]<<" evaluations="<<integral_evaluations[i]<<'\n';
 std::cout<<"SPD_condition_estimate="<<f.condition_estimate_inf()<<" backward_residual="<<solution.backward_residual<<" estimated_forward_sensitivity="<<solution.estimated_forward_sensitivity<<'\n';return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
