#include "cosmology/numerics.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string_view>
#include <sys/resource.h>
#include <vector>
using namespace cosmology::numerics;
// Synthetic family only: C=I+11^T/64, no released survey covariance or W01 replay.
// Independent rank-one inverse/determinant algebra, exact integer/dyadic sums.
// Native reference logs share standard-library ancestry: no independent libm claim.
// ||C||inf=1+n/64; ||C^-1||inf=1+(n-2)/(64+n), hence
// conditioninf=1+(2*n-2)/64 (n>=2); 50.65625 at1590, 2.96875 at64.
// Budgets below were fixed before runs: synthetic normalized target reporting
// at about seven decimal places; component allocations are conservative policy,
// not survey precision or posterior error. No claimed digits below reference floor.
int main(int argc,char** argv){try{
 if(argc>2||(argc==2&&std::string_view(argv[1])!="--size64"))throw std::runtime_error("expected optional --size64");
 const std::size_t n=argc==2?64:1590;
 std::vector<double> matrix(n*n,1./64),rhs(n);
 std::int64_t total=0,squares=0;
 for(std::size_t i=0;i<n;++i){matrix[i*n+i]+=1;const auto k=static_cast<std::int64_t>(i%11)-5;rhs[i]=k/4.;total+=k;squares+=k*k;}
 // Exact integer reductions and analytic rank-one inverse/determinant, independent of Cholesky.
 const long double sum=static_cast<long double>(total)/4;
 const long double true_quadratic=static_cast<long double>(squares)/16-sum*sum/(64+n);
 const long double true_logdet=std::log(1+static_cast<long double>(n)/64);
 const long double normalization=static_cast<long double>(n)*std::log(2*std::acos(-1.L));
 const long double true_logdensity=-(normalization+true_logdet+true_quadratic)/2;
 const auto start=std::chrono::steady_clock::now();
 const auto factor=cholesky(matrix,n,n*n);
 if(factor.status()!=Status::ok)throw std::runtime_error("factor failure");
 const auto solution=solve(factor,rhs,1e-10);
 if(solution.status!=Status::ok)throw std::runtime_error("solve sensitivity failure");
 const auto finish=std::chrono::steady_clock::now();
 long double quadratic=0,max_component_error=0;
 for(std::size_t i=0;i<n;++i){quadratic+=static_cast<long double>(rhs[i])*solution.value[i];max_component_error=std::max(max_component_error,std::abs(static_cast<long double>(solution.value[i])-(rhs[i]-sum/(64+n))));}
 const long double quadratic_error=std::abs(quadratic-true_quadratic);
 const long double logdet_error=std::abs(static_cast<long double>(factor.log_determinant())-true_logdet);
 const double actual_normalization=n*std::log(2*std::acos(-1.));
 const long double normalization_error=std::abs(actual_normalization-normalization);
 const double logdensity=-(actual_normalization+factor.log_determinant()+static_cast<double>(quadratic))/2;
 const long double logdensity_error=std::abs(logdensity-true_logdensity);
 const bool passed=max_component_error<=1e-10&&quadratic_error<=2e-9&&logdet_error<=2e-9&&normalization_error<=4e-9&&logdensity_error<=1e-8&&solution.backward_residual<=1e-12;
 rusage usage{};getrusage(RUSAGE_SELF,&usage);
 std::printf("{\"n\":%zu,\"passed\":%s,\"maximum_component_error\":%.21Lg,\"quadratic_error\":%.21Lg,\"logdet_error\":%.21Lg,\"normalization_error\":%.21Lg,\"logdensity_error\":%.21Lg,\"backward_residual\":%.17g,\"estimated_forward_sensitivity\":%.17g,\"condition_estimate\":%.17g,\"factor_solve_seconds\":%.9g,\"peak_RSS_KiB\":%ld,\"matrix_bytes\":%zu}\n",n,passed?"true":"false",max_component_error,quadratic_error,logdet_error,normalization_error,logdensity_error,solution.backward_residual,solution.estimated_forward_sensitivity,factor.condition_estimate_inf(),std::chrono::duration<double>(finish-start).count(),usage.ru_maxrss,matrix.size()*sizeof(double));
 return passed?0:1;
 }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
