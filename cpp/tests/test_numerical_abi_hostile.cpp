// Interface/status parity evidence only: author owns the shared numerical core.
#include "irred/abi.h"
#include "irred/numerics.hpp"
#include "fixtures/foundations_oracles.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>
namespace {
int checks=0;
void check(bool x,const char* name){++checks;if(!x)throw std::runtime_error(name);}
irred_f64_buffer buffer(const double* p,std::size_t n){return {sizeof(irred_f64_buffer),IRRED_ABI_VERSION,2,0,p,n,n*sizeof(double)};}
irred::numerics::ScalarResult direct(std::uint32_t op,std::span<const double> input,std::size_t i){using namespace irred::numerics;switch(op){case 1:return compensated_sum(input);case 2:return log_sum_exp(input);case 3:return log1p_checked(input[i]);case 4:return expm1_checked(input[i]);case 5:return log_gamma_positive(input[i]);default:throw std::runtime_error("invalid direct op");}}
void parity(std::uint32_t op,std::span<const double> input){
 auto b=buffer(input.data(),input.size());irred_result* owner=nullptr;check(irred_numerics_evaluate(op,&b,&owner)==IRRED_OK&&owner,"completed ABI call");
 const double* values=nullptr;const std::uint32_t* statuses=nullptr;const double* errors=nullptr;const std::uint64_t* work=nullptr;std::uint64_t n=0;
 check(irred_result_numerics_view(owner,&values,&statuses,&errors,&work,&n)==IRRED_OK,"numerical view");
 check(n==irred_numerics_output_length(op,input.size()),"result shape");
 for(std::size_t i=0;i<n;++i){const auto d=direct(op,input,i);check(statuses[i]==static_cast<std::uint32_t>(d.status),"direct status parity");if(d.status==irred::numerics::Status::ok){check(values[i]==d.value,"direct finite parity");check(std::isfinite(values[i]),"finite tag finite value");}check(errors[i]==d.error_estimate&&work[i]==d.evaluations,"diagnostic parity");}
 const std::int64_t* wrong=nullptr;std::uint64_t wrongn=123;check(irred_result_view(owner,&wrong,&wrongn)==IRRED_INVALID_INPUT&&!wrong&&wrongn==0,"wrong result view kind");
 check(irred_result_destroy(owner)==IRRED_OK,"same allocator destruction");
}
}
int main(){try{
 static_assert(sizeof(irred_f64_buffer)==40);static_assert(alignof(irred_f64_buffer)==8);
 static_assert(offsetof(irred_f64_buffer,data)==16);static_assert(offsetof(irred_f64_buffer,length)==24);static_assert(offsetof(irred_f64_buffer,byte_length)==32);
 std::array<double,3> cancellation{1e16,1,-1e16};parity(1,cancellation);
 std::array<double,2> shifted{1000,1000};parity(2,shifted);
 std::array<double,5> mixed{0,.5,-1,-2,std::numeric_limits<double>::quiet_NaN()};parity(3,mixed);
 std::array<double,4> exponents{0,1e-12,1000,std::numeric_limits<double>::infinity()};parity(4,exponents);
 std::array<double,4> gamma{.125,.5,100,0};parity(5,gamma);
 std::array<double,2> overflow{std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};parity(1,overflow);
 for(std::uint32_t op=1;op<=5;++op)parity(op,{});
 std::array<double,1> good{.5};auto b=buffer(good.data(),1);
 auto reject=[&](irred_f64_buffer invalid,std::uint32_t expected,std::uint32_t op=3){auto out=reinterpret_cast<irred_result*>(std::uintptr_t(1));check(irred_numerics_evaluate(op,&invalid,&out)==expected&&!out,"invalid descriptor null output");};
 auto bad=b;bad.abi_version++;reject(bad,IRRED_ABI_MISMATCH);
 bad=b;bad.struct_size--;reject(bad,IRRED_INVALID_INPUT);
 bad=b;bad.element_type=1;reject(bad,IRRED_INVALID_INPUT);
 bad=b;bad.reserved=1;reject(bad,IRRED_INVALID_INPUT);
 bad=b;bad.byte_length--;reject(bad,IRRED_INVALID_INPUT);
 bad=b;bad.data=nullptr;reject(bad,IRRED_INVALID_INPUT);
 bad=b;bad.length=IRRED_MAX_BATCH_ELEMENTS+1ULL;bad.byte_length=bad.length*8;reject(bad,IRRED_INVALID_INPUT);
 bad=b;bad.length=std::numeric_limits<std::uint64_t>::max();bad.byte_length=0;reject(bad,IRRED_INVALID_INPUT);
 alignas(double) std::array<std::byte,16> bytes{};bad=b;bad.data=reinterpret_cast<const double*>(bytes.data()+1);reject(bad,IRRED_INVALID_INPUT);
 reject(b,IRRED_INVALID_INPUT,0);reject(b,IRRED_INVALID_INPUT,6);
 irred_result* out=reinterpret_cast<irred_result*>(std::uintptr_t(1));check(irred_numerics_evaluate(3,nullptr,&out)==IRRED_INVALID_INPUT&&!out,"null descriptor");check(irred_numerics_evaluate(3,&b,nullptr)==IRRED_INVALID_INPUT,"null output argument");
 const double* values=reinterpret_cast<const double*>(std::uintptr_t(1));const std::uint32_t* statuses=reinterpret_cast<const std::uint32_t*>(std::uintptr_t(1));const double* errors=values;const std::uint64_t* evaluations=reinterpret_cast<const std::uint64_t*>(std::uintptr_t(1));std::uint64_t n=99;
 check(irred_result_numerics_view(nullptr,&values,&statuses,&errors,&evaluations,&n)==IRRED_INVALID_INPUT&&!values&&!statuses&&!errors&&!evaluations&&n==0,"null view clears outputs");
 check(irred_result_numerics_view(nullptr,nullptr,&statuses,&errors,&evaluations,&n)==IRRED_INVALID_INPUT,"null view argument");
 check(irred_result_destroy(nullptr)==IRRED_OK,"null cleanup");
 // Frozen independent numeric point is convenient bridge evidence, not a new core gate.
 auto x=buffer(good.data(),1);out=nullptr;check(irred_numerics_evaluate(3,&x,&out)==IRRED_OK,"frozen fixture call");check(irred_result_numerics_view(out,&values,&statuses,&errors,&evaluations,&n)==IRRED_OK,"frozen fixture view");
 double expected=0;for(auto f:irred::test_fixtures::foundations_oracles)if(f.id=="log1p_0.5")expected=f.rounded;
 check(statuses[0]==IRRED_NUMERICAL_STATUS_OK&&std::abs(values[0]-expected)<1e-15,"frozen Decimal log1p bridge point");irred_result_destroy(out);
 std::printf("{\"suite\":\"numerical_ABI_hostile\",\"checks\":%d,\"passed\":true,\"ancestry\":\"core author; shared-core direct/ABI parity and contract evidence only\"}\n",checks);return 0;
 }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
