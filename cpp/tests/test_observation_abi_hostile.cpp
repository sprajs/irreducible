#include "cosmology/abi.h"
#include "cosmology/observations.hpp"
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <string>
#include <cstdint>
#include <limits>
static std::int64_t fail_after=-1,live=0;
void* operator new(std::size_t n){if(fail_after==0)throw std::bad_alloc();if(fail_after>0)--fail_after;void* p=std::malloc(n?n:1);if(!p)throw std::bad_alloc();++live;return p;}
void operator delete(void* p) noexcept{if(p){--live;std::free(p);}}
void operator delete(void* p,std::size_t) noexcept{::operator delete(p);}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete[](void* p) noexcept{::operator delete(p);}
void operator delete[](void* p,std::size_t) noexcept{::operator delete(p);}
#define CHECK(...) do{if(!(__VA_ARGS__))throw std::runtime_error("observation ABI contract: " #__VA_ARGS__);}while(false)
using namespace cosmology::observations;
static cosmo_bytes bytes(const std::string& s){return {reinterpret_cast<const uint8_t*>(s.data()),s.size()};}
static cosmo_f64_buffer doubles(const double* p,std::size_t n){return {sizeof(cosmo_f64_buffer),COSMO_ABI_VERSION,2,0,p,n,n*sizeof(double)};}
struct Fixture{
 std::string hash=std::string(64,'a'),covhash=std::string(64,'b'),order="supplied source-order declaration",a="measurement-A",b="measurement-B",event="same-event";
 cosmo_bytes ids[2]={bytes(a),bytes(b)},events[2]={bytes(event),bytes(event)};
 double values[2]={2,-3},matrix[4]={4,1.00000003,1,9};uint8_t missing[2]={0,0};uint64_t quality[2]={2,128};
 cosmo_observation_policy policy{2,4,4096};cosmo_observation_descriptor d{};
 Fixture(){d.struct_size=sizeof(d);d.abi_version=COSMO_ABI_VERSION;d.profile=static_cast<uint32_t>(Profile::gaussian_fixture_v1);d.role=static_cast<uint32_t>(Role::synthetic_control);d.unit=static_cast<uint32_t>(Unit::magnitude);d.calibration=static_cast<uint32_t>(Calibration::not_applicable);d.uncertainty=static_cast<uint32_t>(Uncertainty::covariance);d.uncertainty_unit=static_cast<uint32_t>(UncertaintyUnit::magnitude_squared);d.component=static_cast<uint32_t>(Component::total);d.table_sha256=bytes(hash);d.uncertainty_sha256=bytes(covhash);d.ordering_provenance=bytes(order);d.measurement_ids={ids,2,sizeof ids};d.event_ids={events,2,sizeof events};d.uncertainty_axis_ids=d.measurement_ids;d.values=doubles(values,2);d.zhd=d.zcmb=d.zhel=doubles(nullptr,0);d.uncertainty_matrix=doubles(matrix,4);d.missing={missing,2,2};d.quality={quality,2,sizeof quality};}
};
int main(){Fixture f;cosmo_prepared* p=nullptr;uint32_t semantic=99;CHECK(cosmo_prepare_observations(&f.d,&f.policy,&p,&semantic)==COSMO_OK&&semantic==0&&p);cosmo_observation_descriptor v{};f.values[0]=123;CHECK(cosmo_observation_source_view(p,&v)==COSMO_OK);CHECK(v.values.data[0]==2&&v.uncertainty_matrix.data[1]==1.00000003&&v.uncertainty_matrix.data[2]==1);CHECK(v.event_ids.data[0].length==v.event_ids.data[1].length);CHECK(v.quality.data[1]==128);
 cosmo_result* r=nullptr;CHECK(cosmo_observation_select(p,0,&r,&semantic)==COSMO_OK&&semantic==0&&r);const uint8_t* mask=nullptr;const uint64_t* indices=nullptr;uint64_t n=0,m=0;CHECK(cosmo_result_selection_view(r,&mask,&n,&indices,&m)==COSMO_OK&&n==2&&m==2&&mask[0]==1&&mask[1]==1&&indices[0]==0&&indices[1]==1);CHECK(cosmo_result_destroy(r)==COSMO_OK);CHECK(cosmo_observation_select(p,99,&r,&semantic)==COSMO_OK&&semantic==static_cast<uint32_t>(Status::incompatible_semantics)&&!r);CHECK(cosmo_observation_destroy(p)==COSMO_OK);p=nullptr;
 auto bad=f.d;bad.values.abi_version++;CHECK(cosmo_prepare_observations(&bad,&f.policy,&p,&semantic)==COSMO_INVALID_INPUT&&!p);bad=f.d;bad.abi_version++;CHECK(cosmo_prepare_observations(&bad,&f.policy,&p,&semantic)==COSMO_ABI_MISMATCH&&!p);bad=f.d;bad.values.byte_length--;CHECK(cosmo_prepare_observations(&bad,&f.policy,&p,&semantic)==COSMO_INVALID_INPUT&&!p);bad=f.d;bad.values.data=nullptr;CHECK(cosmo_prepare_observations(&bad,&f.policy,&p,&semantic)==COSMO_INVALID_INPUT&&!p);bad=f.d;bad.profile=99;CHECK(cosmo_prepare_observations(&bad,&f.policy,&p,&semantic)==COSMO_OK&&semantic==static_cast<uint32_t>(Status::incompatible_semantics)&&!p);bad=f.d;bad.ordering_provenance={};CHECK(cosmo_prepare_observations(&bad,&f.policy,&p,&semantic)==COSMO_OK&&!p);auto cap=f.policy;cap.maximum_rows=1;CHECK(cosmo_prepare_observations(&f.d,&cap,&p,&semantic)==COSMO_INVALID_INPUT&&!p);cap=f.policy;cap.maximum_string_bytes=1;CHECK(cosmo_prepare_observations(&f.d,&cap,&p,&semantic)==COSMO_INVALID_INPUT&&!p);
 // Fail every allocation in the coarse prepare path; no exception crosses ABI
 // and every temporary/copy/partially built handle must be reclaimed.
 bool succeeded=false;int failures=0;
 for(int point=0;point<200;++point){auto before=live;fail_after=point;auto status=cosmo_prepare_observations(&f.d,&f.policy,&p,&semantic);fail_after=-1;if(status==COSMO_OK){CHECK(p&&semantic==0);CHECK(cosmo_observation_destroy(p)==COSMO_OK);p=nullptr;succeeded=true;CHECK(live==before);break;}CHECK(status==COSMO_ALLOCATION_FAILURE&&!p);CHECK(live==before);++failures;}
 CHECK(succeeded&&failures>10);CHECK(cosmo_prepare_observations(&f.d,&f.policy,&p,&semantic)==COSMO_OK&&p);
 succeeded=false;failures=0;for(int point=0;point<30;++point){auto before=live;fail_after=point;auto status=cosmo_observation_select(p,0,&r,&semantic);fail_after=-1;if(status==COSMO_OK){CHECK(r);CHECK(cosmo_result_destroy(r)==COSMO_OK);r=nullptr;succeeded=true;CHECK(live==before);break;}CHECK(status==COSMO_ALLOCATION_FAILURE&&!r);CHECK(live==before);++failures;}CHECK(succeeded&&failures>=3);CHECK(cosmo_observation_destroy(p)==COSMO_OK);CHECK(cosmo_observation_destroy(nullptr)==COSMO_OK);
}
