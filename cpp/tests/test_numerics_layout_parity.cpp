// Compatibility transcript, not an independent scientific oracle. Compile this
// same fixture with baseline and candidate; compare every hex binary64 value,
// status and public diagnostic exactly. Analytic accuracy has separate suites.
#include "irred/numerics.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <limits>
#include <utility>
using namespace irred::numerics;
namespace {
void result(const char*label,const SolveResult&r){
 std::printf("%s status=%u values=%zu",label,(unsigned)r.status,r.value.size());
 for(double x:r.value)std::printf(" %a",x);
 for(double x:{r.backward_residual,r.estimated_forward_sensitivity,
                r.pre_cast_backward_residual,r.pre_cast_sensitivity_estimate,
                r.output_rounding_error_inf,r.output_rounding_error_relative,
                r.post_cast_backward_residual})std::printf(" %a",x);
 std::puts("");
}
void fixture(const char*label,std::span<const double> c,std::size_t n,
             std::span<const double> rhs,Arithmetic arithmetic){
 auto f=cholesky(c,n,n*n,arithmetic);
 std::printf("%s arithmetic=%u factor=%u n=%zu logdet=%a cond=%a\n",label,
 (unsigned)arithmetic,(unsigned)f.status(),f.size(),f.log_determinant(),f.condition_estimate_inf());
 result("ordinary",solve(f,rhs,1e-6));result("tight",solve(f,rhs,1e-30));
 result("zero",solve(f,std::vector<double>(n,0),1e-6));
 result("invalidbudget",solve(f,rhs,0));
 result("wrongsize",solve(f,std::span<const double>{},1e-6));
 auto moved=std::move(f);result("moved",solve(moved,rhs,1e-6));
 result("movedfrom",solve(f,rhs,1e-6));
 auto nonfinite=std::vector<double>(rhs.begin(),rhs.end());
 nonfinite[0]=std::numeric_limits<double>::quiet_NaN();
 result("nonfinite",solve(moved,nonfinite,1e-6));
 if(arithmetic==Arithmetic::longdouble_cpu_v1){
  int before=std::fegetround();std::fesetround(FE_DOWNWARD);
  result("unsupportedround",solve(moved,rhs,1e-6));std::fesetround(before);
 }
}
}
int main(){
 for(auto arithmetic:{Arithmetic::binary64_legacy_v1,Arithmetic::longdouble_cpu_v1}){
  fixture("rational2",std::array<double,4>{4,1,1,9},2,std::array<double,2>{2,-3},arithmetic);
  fixture("rational3",std::array<double,9>{4,1,0,1,3,1,0,1,2},3,std::array<double,3>{1,-2,3},arithmetic);
  std::array<double,256> c{};std::array<double,16> rhs{};
  for(unsigned i=0;i<16;++i){rhs[i]=(int(i%5)-2)/8.;for(unsigned j=0;j<16;++j)c[i*16+j]=(i==j?1.:0.)+((i%3)+1)*((j%3)+1)/64.;}
  fixture("dyadic16",c,16,rhs,arithmetic);
  fixture("near-singular",std::array<double,4>{1,.999999999, .999999999,1},2,std::array<double,2>{1,-1},arithmetic);
  fixture("normal-min",std::array<double,1>{std::numeric_limits<double>::min()},1,std::array<double,1>{std::numeric_limits<double>::min()},arithmetic);
  fixture("normal-max",std::array<double,1>{std::numeric_limits<double>::max()},1,std::array<double,1>{1},arithmetic);
  fixture("output-overflow",std::array<double,1>{std::numeric_limits<double>::min()},1,std::array<double,1>{std::numeric_limits<double>::max()},arithmetic);
  for(auto bad:{std::array<double,4>{1,2,2,1},std::array<double,4>{1,.2,.3,1},std::array<double,4>{1,0,0,std::numeric_limits<double>::infinity()}}){auto f=cholesky(bad,2,4,arithmetic);std::printf("invalid factor=%u\n",(unsigned)f.status());}
  auto capped=cholesky(c,16,255,arithmetic);std::printf("cap factor=%u\n",(unsigned)capped.status());
 }
 result("default",solve(Factorization{},std::array<double,1>{1},1e-6));
 auto unknown=cholesky(std::array<double,1>{1},1,1,static_cast<Arithmetic>(99));
 std::printf("unknown factor=%u id=%s\n",(unsigned)unknown.status(),unknown.arithmetic_id());
}
