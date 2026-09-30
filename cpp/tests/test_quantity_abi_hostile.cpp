#include "irred/abi.h"
#include "irred/quantities.hpp"
#include <array>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <cstdio>
#include <stdexcept>
using namespace irred;
static unsigned checks=0;static void require(bool b,const char* name){++checks;if(!b)throw std::runtime_error(name);}
static irred_quantity_metadata meta(uint32_t u,uint32_t r,uint32_t f=0,uint32_t c=0){return {sizeof(irred_quantity_metadata),IRRED_ABI_VERSION,u,r,f,c,1,0};}
static irred_f64_buffer desc(const double* data,uint64_t n){return {sizeof(irred_f64_buffer),IRRED_ABI_VERSION,2,0,data,n,n*8};}
static Quantity native(double v,const irred_quantity_metadata& m){return {v,static_cast<Unit>(m.unit),static_cast<Role>(m.role),static_cast<Frame>(m.frame),static_cast<LengthConvention>(m.convention)};}
static Target target(const irred_quantity_metadata& m){return {static_cast<Unit>(m.unit),static_cast<Role>(m.role),static_cast<Frame>(m.frame),static_cast<LengthConvention>(m.convention)};}
int main(){try{
 static_assert(sizeof(irred_quantity_metadata)==32&&alignof(irred_quantity_metadata)==4);static_assert(offsetof(irred_quantity_metadata,constant_set)==24&&offsetof(irred_quantity_metadata,reserved)==28);
 static_assert(sizeof(irred_f64_buffer)==40&&alignof(irred_f64_buffer)==8);
 const std::array<double,3> values{2,-3,4};auto d=desc(values.data(),values.size());auto source=meta(IRRED_UNIT_KILOMETRE,IRRED_ROLE_PHYSICAL_LENGTH,0,IRRED_CONVENTION_PHYSICAL);auto dest=meta(IRRED_UNIT_METRE,IRRED_ROLE_PHYSICAL_LENGTH,0,IRRED_CONVENTION_PHYSICAL);
 const auto compare=[&](irred_quantity_metadata a,irred_quantity_metadata b){irred_result* result=nullptr;auto before=a;require(irred_convert_quantities(&d,&a,&b,&result)==IRRED_OK&&result,"conversion call completes");require(std::memcmp(&before,&a,sizeof(a))==0,"borrowed metadata not mutated");const double* out=nullptr;const uint32_t* statuses=nullptr;uint64_t n=0;require(irred_result_f64_view(result,&out,&statuses,&n)==IRRED_OK&&n==values.size()&&out&&statuses,"tagged quantity view");for(std::size_t i=0;i<n;++i){auto direct=convert(native(values[i],a),target(b));require(statuses[i]==static_cast<uint32_t>(direct.status),"direct ABI status parity");if(direct.status==QuantityStatus::ok)require(out[i]==direct.target.value,"direct ABI finite value parity");}const int64_t* wrong=reinterpret_cast<const int64_t*>(uintptr_t(1));n=123;require(irred_result_view(result,&wrong,&n)==IRRED_INVALID_INPUT&&!wrong&&n==0,"wrong result type rejected");irred_result_destroy(result);};
 compare(source,dest);compare(meta(IRRED_UNIT_ONE,IRRED_ROLE_REDSHIFT,IRRED_FRAME_CMB),meta(IRRED_UNIT_ONE,IRRED_ROLE_REDSHIFT,IRRED_FRAME_CMB));compare(source,meta(IRRED_UNIT_METRE,IRRED_ROLE_LUMINOSITY_DISTANCE,0,IRRED_CONVENTION_PHYSICAL));compare(meta(IRRED_UNIT_ONE,IRRED_ROLE_REDSHIFT),meta(IRRED_UNIT_ONE,IRRED_ROLE_REDSHIFT));
 const auto reject=[&](irred_f64_buffer input,irred_quantity_metadata sm,uint32_t expected,const char* name){irred_result* result=reinterpret_cast<irred_result*>(uintptr_t(1));require(irred_convert_quantities(&input,&sm,&dest,&result)==expected&&!result,name);};
 auto bad=d;bad.abi_version++;reject(bad,source,IRRED_ABI_MISMATCH,"f64 stale ABI");bad=d;bad.element_type=1;reject(bad,source,IRRED_INVALID_INPUT,"f64 wrong element type");bad=d;bad.byte_length--;reject(bad,source,IRRED_INVALID_INPUT,"f64 short buffer");bad=d;bad.length=UINT64_MAX;bad.byte_length=UINT64_MAX;reject(bad,source,IRRED_INVALID_INPUT,"f64 count overflow");bad=d;bad.data=nullptr;reject(bad,source,IRRED_INVALID_INPUT,"f64 null nonempty");alignas(double)std::array<unsigned char,32> bytes{};bad=d;bad.data=reinterpret_cast<const double*>(bytes.data()+1);reject(bad,source,IRRED_INVALID_INPUT,"f64 misalignment");auto sm=source;sm.abi_version++;reject(d,sm,IRRED_ABI_MISMATCH,"metadata ABI");sm=source;sm.struct_size--;reject(d,sm,IRRED_INVALID_INPUT,"metadata size");sm=source;sm.reserved=1;reject(d,sm,IRRED_INVALID_INPUT,"metadata reserved");
 irred_result* result=reinterpret_cast<irred_result*>(uintptr_t(1));require(irred_convert_quantities(&d,nullptr,&dest,&result)==IRRED_INVALID_INPUT&&!result,"null source metadata");
 auto empty=desc(nullptr,0);require(irred_convert_quantities(&empty,&source,&dest,&result)==IRRED_OK&&result,"empty conversion");const double* out=nullptr;const uint32_t* statuses=nullptr;uint64_t n=99;require(irred_result_f64_view(result,&out,&statuses,&n)==IRRED_OK&&n==0,"empty tagged view");irred_result_destroy(result);
 std::printf("{\"suite\":\"independent_quantity_ABI_hostile\",\"checks\":%u,\"passed\":true,\"ancestry\":\"shared scientific core; interface parity only\"}\n",checks);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
