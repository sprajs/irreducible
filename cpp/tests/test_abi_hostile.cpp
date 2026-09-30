#include "irred/abi.h"
#include <array>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <cstdlib>
#include <cstddef>

static unsigned checks=0;
static void require(bool x,const char* name) { ++checks; if(!x) {std::fprintf(stderr,"FAIL %s\n",name); std::exit(1);} }
static cosmo_i64_buffer buffer(const int64_t* p,uint64_t n) {return {sizeof(cosmo_i64_buffer),COSMO_ABI_VERSION,1,0,p,n,n*8};}
int main() {
 static_assert(sizeof(void*)==8,"this evidence case is scoped to the current 64-bit host");
 static_assert(sizeof(cosmo_i64_buffer)==40 && alignof(cosmo_i64_buffer)==8);
 static_assert(offsetof(cosmo_i64_buffer,struct_size)==0 && offsetof(cosmo_i64_buffer,abi_version)==4);
 static_assert(offsetof(cosmo_i64_buffer,element_type)==8 && offsetof(cosmo_i64_buffer,reserved)==12);
 static_assert(offsetof(cosmo_i64_buffer,data)==16 && offsetof(cosmo_i64_buffer,length)==24 && offsetof(cosmo_i64_buffer,byte_length)==32);
 const std::array<int64_t,5> a{0,1,-2,1048576,-1048576},b{0,-1,3,-1048576,1048576}, expected{0,0,1,0,0};
 auto aa=buffer(a.data(),a.size()),bb=buffer(b.data(),b.size());
 cosmo_result* result=nullptr; require(cosmo_add(&aa,&bb,0,&result)==COSMO_OK,"exact arithmetic call");
 const int64_t* values=nullptr;uint64_t n=0; require(cosmo_result_view(result,&values,&n)==COSMO_OK && n==expected.size(),"exact arithmetic shape");
 for(size_t i=0;i<n;i++) require(values[i]==expected[i],"independent integer oracle");
 require(cosmo_result_destroy(result)==COSMO_OK,"normal cleanup");
 const auto reject=[&](cosmo_i64_buffer bad,uint32_t status,const char* name) {
  result=reinterpret_cast<cosmo_result*>(uintptr_t(1));
  require(cosmo_add(&bad,&bb,0,&result)==status && result==nullptr,name);
 };
 auto bad=aa;bad.abi_version++;reject(bad,COSMO_ABI_MISMATCH,"version");
 bad=aa;bad.struct_size--;reject(bad,COSMO_INVALID_INPUT,"short struct");
 bad=aa;bad.reserved=1;reject(bad,COSMO_INVALID_INPUT,"reserved");
 bad=aa;bad.element_type=99;reject(bad,COSMO_INVALID_INPUT,"type");
 bad=aa;bad.byte_length--;reject(bad,COSMO_INVALID_INPUT,"short bytes");
 bad=aa;bad.data=nullptr;reject(bad,COSMO_INVALID_INPUT,"null nonempty");
 alignas(int64_t) std::array<unsigned char,48> raw{};
 bad=aa;bad.data=reinterpret_cast<const int64_t*>(raw.data()+1);reject(bad,COSMO_INVALID_INPUT,"misalignment");
 bad=aa;bad.length=UINT64_MAX;bad.byte_length=UINT64_MAX;reject(bad,COSMO_INVALID_INPUT,"size overflow before access");
 bad=aa;bad.length--;bad.byte_length=bad.length*8;reject(bad,COSMO_INVALID_INPUT,"batch mismatch");
 for(uint32_t fault=1;fault<=2;fault++) {result=reinterpret_cast<cosmo_result*>(uintptr_t(1));require(cosmo_add(&aa,&bb,fault,&result)==(fault==1?COSMO_EXCEPTION:COSMO_ALLOCATION_FAILURE)&&!result,"exception integer fallback");}
 const int64_t max=INT64_MAX,one=1,min=INT64_MIN,negative=-1;
 auto x=buffer(&max,1),y=buffer(&one,1);require(cosmo_add(&x,&y,0,&result)==COSMO_OVERFLOW&&!result,"positive overflow");
 x=buffer(&min,1);y=buffer(&negative,1);require(cosmo_add(&x,&y,0,&result)==COSMO_OVERFLOW&&!result,"negative overflow");
 x=buffer(nullptr,0);y=x;require(cosmo_add(&x,&y,0,&result)==COSMO_OK,"null empty");require(cosmo_result_view(result,&values,&n)==COSMO_OK&&n==0,"empty view");cosmo_result_destroy(result);
 require(cosmo_add(&aa,&bb,0,nullptr)==COSMO_INVALID_INPUT,"null output destination");
 require(cosmo_result_destroy(nullptr)==COSMO_OK,"null destruction");
 std::printf("{\"suite\":\"independent_abi_adversarial\",\"checks\":%u,\"passed\":true,\"scope\":\"synchronous i64 bridge only\"}\n",checks);
}
