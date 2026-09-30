#include "irred/abi.h"
#include "irred/observations.hpp"
#include "irred/observation_enum_checks.inc"
#include "result_internal.hpp"
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <utility>
using namespace irred::observations;
namespace {
template<class T> bool bounded(const T* p,uint64_t n,uint64_t bytes,uint64_t cap) {
 return n<=cap && n<=std::numeric_limits<size_t>::max()/sizeof(T) && bytes==n*sizeof(T) && (!n || (p && reinterpret_cast<uintptr_t>(p)%alignof(T)==0));
}
bool doubles(const irred_f64_buffer& b,uint64_t cap) {
 return b.struct_size==sizeof(b)&&b.abi_version==IRRED_ABI_VERSION&&b.element_type==2&&!b.reserved&&bounded(b.data,b.length,b.byte_length,cap);
}
bool bytes(const irred_bytes& b,uint64_t& left) {
 if(b.length>left || !bounded(b.data,b.length,b.length,left))return false;
 left-=b.length;return true;
}
bool strings(const irred_strings& a,uint64_t cap,uint64_t& left) {
 if(!bounded(a.data,a.length,a.byte_length,cap))return false;
 for(uint64_t i=0;i<a.length;++i)if(!bytes(a.data[i],left))return false;
 return true;
}
std::string copy(const irred_bytes& b) { return b.length ? std::string(reinterpret_cast<const char*>(b.data),static_cast<size_t>(b.length)) : std::string{}; }
std::vector<std::string> copy(const irred_strings& b) {std::vector<std::string> r;r.reserve(b.length);for(uint64_t i=0;i<b.length;++i)r.push_back(copy(b.data[i]));return r;}
template<class T> std::vector<T> copy(const T* p,uint64_t n) {return n?std::vector<T>(p,p+static_cast<size_t>(n)):std::vector<T>{};}
irred_bytes view(const std::string& s) {return {reinterpret_cast<const uint8_t*>(s.data()),s.size()};}
irred_f64_buffer view(const std::vector<double>& v) {return {sizeof(irred_f64_buffer),IRRED_ABI_VERSION,2,0,v.data(),v.size(),v.size()*sizeof(double)};}
irred_u8_buffer view(const std::vector<uint8_t>& v) {return {v.data(),v.size(),v.size()};}
}
struct irred_prepared {
 std::shared_ptr<const Prepared> prepared;
 std::vector<irred_bytes> measurement_ids,event_ids,axis_ids;
 irred_observation_descriptor descriptor{};
 explicit irred_prepared(Prepared p):prepared(std::make_shared<const Prepared>(std::move(p))) {
  const auto& s=prepared->source();
  auto fill=[](const auto& input,auto& output){output.reserve(input.size());for(const auto& x:input)output.push_back(view(x));};
  fill(s.measurement_ids,measurement_ids);fill(s.event_ids,event_ids);fill(s.uncertainty_axis_ids,axis_ids);
  auto strings_view=[](const auto& v){return irred_strings{v.data(),v.size(),v.size()*sizeof(irred_bytes)};};
  auto& d=descriptor;d.struct_size=sizeof(d);d.abi_version=IRRED_ABI_VERSION;
  d.profile=static_cast<uint32_t>(s.profile);d.role=static_cast<uint32_t>(s.role);d.unit=static_cast<uint32_t>(s.unit);d.calibration=static_cast<uint32_t>(s.calibration);d.uncertainty=static_cast<uint32_t>(s.uncertainty);d.uncertainty_unit=static_cast<uint32_t>(s.uncertainty_unit);d.component=static_cast<uint32_t>(s.component);
  d.table_sha256=view(s.table_sha256);d.uncertainty_sha256=view(s.uncertainty_sha256);d.calibration_provenance=view(s.calibration_provenance);d.dependence_provenance=view(s.dependence_provenance);d.quality_dictionary=view(s.quality_dictionary);d.ordering_provenance=view(s.ordering_provenance);
  d.measurement_ids=strings_view(measurement_ids);d.event_ids=strings_view(event_ids);d.uncertainty_axis_ids=strings_view(axis_ids);
  d.values=view(s.values);d.zhd=view(s.zhd);d.zcmb=view(s.zcmb);d.zhel=view(s.zhel);d.uncertainty_matrix=view(s.uncertainty_matrix);
  d.missing=view(s.missing);d.zhd_missing=view(s.zhd_missing);d.zcmb_missing=view(s.zcmb_missing);d.zhel_missing=view(s.zhel_missing);d.source_selection=view(s.source_selection);d.quality={s.quality.data(),s.quality.size(),s.quality.size()*sizeof(uint64_t)};
 }
};
extern "C" uint32_t irred_prepare_observations(const irred_observation_descriptor* d,const irred_observation_policy* p,irred_prepared** out,uint32_t* semantic) {
 if(out)*out=nullptr;
 if(semantic)*semantic=static_cast<uint32_t>(Status::invalid_shape);
 if(!out||!semantic||!d||!p)return IRRED_INVALID_INPUT;
 if(d->abi_version!=IRRED_ABI_VERSION)return IRRED_ABI_MISMATCH;
 if(d->struct_size!=sizeof(*d)||d->reserved||p->maximum_rows>std::numeric_limits<size_t>::max()||p->maximum_matrix_elements>std::numeric_limits<size_t>::max()||p->maximum_string_bytes>std::numeric_limits<size_t>::max())return IRRED_INVALID_INPUT;
 const auto n=d->values.length;
 if(!n||n>p->maximum_rows||n>std::numeric_limits<uint64_t>::max()/n)return IRRED_INVALID_INPUT;
 if(d->uncertainty!=IRRED_OBSERVATION_UNCERTAINTY_NONE && n*n>p->maximum_matrix_elements)return IRRED_INVALID_INPUT;
 if(!doubles(d->values,n)||!doubles(d->zhd,n)||!doubles(d->zcmb,n)||!doubles(d->zhel,n)||!doubles(d->uncertainty_matrix,p->maximum_matrix_elements))return IRRED_INVALID_INPUT;
 for(const auto* b:{&d->missing,&d->zhd_missing,&d->zcmb_missing,&d->zhel_missing,&d->source_selection})if(!bounded(b->data,b->length,b->byte_length,n))return IRRED_INVALID_INPUT;
 if(!bounded(d->quality.data,d->quality.length,d->quality.byte_length,n))return IRRED_INVALID_INPUT;
 uint64_t left=p->maximum_string_bytes;
 for(const auto* b:{&d->table_sha256,&d->uncertainty_sha256,&d->calibration_provenance,&d->dependence_provenance,&d->quality_dictionary,&d->ordering_provenance})if(!bytes(*b,left))return IRRED_INVALID_INPUT;
 for(const auto* a:{&d->measurement_ids,&d->event_ids,&d->uncertainty_axis_ids})if(!strings(*a,n,left))return IRRED_INVALID_INPUT;
 try {
  Input s{};s.profile=static_cast<Profile>(d->profile);s.role=static_cast<Role>(d->role);s.unit=static_cast<Unit>(d->unit);s.calibration=static_cast<Calibration>(d->calibration);s.uncertainty=static_cast<Uncertainty>(d->uncertainty);s.uncertainty_unit=static_cast<UncertaintyUnit>(d->uncertainty_unit);s.component=static_cast<Component>(d->component);
  s.table_sha256=copy(d->table_sha256);s.uncertainty_sha256=copy(d->uncertainty_sha256);s.calibration_provenance=copy(d->calibration_provenance);s.dependence_provenance=copy(d->dependence_provenance);s.quality_dictionary=copy(d->quality_dictionary);s.ordering_provenance=copy(d->ordering_provenance);
  s.measurement_ids=copy(d->measurement_ids);s.event_ids=copy(d->event_ids);s.uncertainty_axis_ids=copy(d->uncertainty_axis_ids);
  s.values=copy(d->values.data,d->values.length);s.zhd=copy(d->zhd.data,d->zhd.length);s.zcmb=copy(d->zcmb.data,d->zcmb.length);s.zhel=copy(d->zhel.data,d->zhel.length);s.uncertainty_matrix=copy(d->uncertainty_matrix.data,d->uncertainty_matrix.length);
  s.missing=copy(d->missing.data,d->missing.length);s.zhd_missing=copy(d->zhd_missing.data,d->zhd_missing.length);s.zcmb_missing=copy(d->zcmb_missing.data,d->zcmb_missing.length);s.zhel_missing=copy(d->zhel_missing.data,d->zhel_missing.length);s.source_selection=copy(d->source_selection.data,d->source_selection.length);s.quality=copy(d->quality.data,d->quality.length);
  auto result=prepare(std::move(s),{static_cast<size_t>(p->maximum_rows),static_cast<size_t>(p->maximum_matrix_elements),static_cast<size_t>(p->maximum_string_bytes)});*semantic=static_cast<uint32_t>(result.status());
  if(result.status()==Status::ok)*out=new irred_prepared(std::move(result));
  return IRRED_OK;
 }catch(const std::bad_alloc&){return IRRED_ALLOCATION_FAILURE;}catch(...){return IRRED_EXCEPTION;}
}
extern "C" uint32_t irred_observation_source_view(const irred_prepared* p,irred_observation_descriptor* d) {if(!d)return IRRED_INVALID_INPUT;*d={};if(!p)return IRRED_INVALID_INPUT;*d=p->descriptor;return IRRED_OK;}
extern "C" uint32_t irred_observation_select(const irred_prepared* p,uint32_t selection,irred_result** out,uint32_t* semantic) {
 if(out)*out=nullptr;
 if(semantic)*semantic=static_cast<uint32_t>(Status::invalid_shape);
 if(!p||!out||!semantic)return IRRED_INVALID_INPUT;
 try {auto selected=p->prepared->select(static_cast<Selection>(selection));*semantic=static_cast<uint32_t>(selected.status);if(selected.status!=Status::ok)return IRRED_OK;
 auto result=new irred_result{};try {result->kind=4;result->selection_mask=std::move(selected.mask);result->selection_indices.assign(selected.source_indices.begin(),selected.source_indices.end());}catch(...){delete result;throw;}*out=result;return IRRED_OK;
 }catch(const std::bad_alloc&){return IRRED_ALLOCATION_FAILURE;}catch(...){return IRRED_EXCEPTION;}
}
extern "C" uint32_t irred_result_selection_view(const irred_result* r,const uint8_t** mask,uint64_t* n,const uint64_t** indices,uint64_t* m) {if(!mask||!n||!indices||!m)return IRRED_INVALID_INPUT;*mask=nullptr;*indices=nullptr;*n=*m=0;if(!r||r->kind!=4)return IRRED_INVALID_INPUT;*mask=r->selection_mask.data();*n=r->selection_mask.size();*indices=r->selection_indices.data();*m=r->selection_indices.size();return IRRED_OK;}
extern "C" uint32_t irred_observation_destroy(irred_prepared* p) {delete p;return IRRED_OK;}

const irred::observations::Prepared *native_observations(const irred_prepared *p) noexcept {
 return p ? p->prepared.get() : nullptr;
}

std::shared_ptr<const irred::observations::Prepared>
shared_native_observations(const irred_prepared *p) noexcept {
 return p ? p->prepared : nullptr;
}
extern "C" uint32_t irred_observation_retained_bytes(const irred_prepared *p,
                                                    uint64_t *bytes) {
 if(bytes)*bytes=0;
 if(!p||!bytes)return IRRED_INVALID_INPUT;
 const auto payload=irred::observations::retained_source_payload_bound(*p->prepared);
 if(!payload)return IRRED_OVERFLOW;
 uint64_t total=*payload;
 auto add=[&](uint64_t n,uint64_t width){if(n>UINT64_MAX/width||n*width>UINT64_MAX-total)return false;total+=n*width;return true;};
 // Conservatively retains the acquisition-view charge until all consumers
 // release the shared source, even if the acquisition handle was released.
 if(!add(1,sizeof(irred::observations::Prepared)+sizeof(irred_prepared))||
    !add(p->measurement_ids.capacity(),sizeof(irred_bytes))||
    !add(p->event_ids.capacity(),sizeof(irred_bytes))||
    !add(p->axis_ids.capacity(),sizeof(irred_bytes)))return IRRED_OVERFLOW;
 *bytes=total;
 return IRRED_OK;
}
extern "C" uint32_t irred_observation_preparation_bytes(const irred_observation_descriptor *d,uint64_t *bytes) {
 if(bytes)*bytes=0;
 if(!d||!bytes)return IRRED_INVALID_INPUT;
 if(d->abi_version!=IRRED_ABI_VERSION)return IRRED_ABI_MISMATCH;
 if(d->struct_size!=sizeof(*d))return IRRED_INVALID_INPUT;
 uint64_t total=sizeof(Prepared)+sizeof(irred_prepared)+4096;
 auto add=[&](uint64_t n,uint64_t width){if(n>UINT64_MAX/width||n*width>UINT64_MAX-total)return false;total+=n*width;return true;};
 uint64_t string_left=16*1024*1024;
 for(const auto* s:{&d->table_sha256,&d->uncertainty_sha256,&d->calibration_provenance,&d->dependence_provenance,&d->quality_dictionary,&d->ordering_provenance}) {
  if(!::bytes(*s,string_left)||!add(s->length+16,2))return IRRED_INVALID_INPUT;
 }
 for(const auto* v:{&d->measurement_ids,&d->event_ids,&d->uncertainty_axis_ids}) {
  if(!strings(*v,10000,string_left)||!add(v->length,2*(sizeof(std::string)+sizeof(irred_bytes))))return IRRED_INVALID_INPUT;
  for(uint64_t i=0;i<v->length;++i)if(!add(v->data[i].length+16,2))return IRRED_OVERFLOW;
 }
 for(const auto* v:{&d->values,&d->zhd,&d->zcmb,&d->zhel,&d->uncertainty_matrix}) {
  if(!doubles(*v,25000000)||!add(v->length,2*sizeof(double)))return IRRED_INVALID_INPUT;
 }
 for(const auto* v:{&d->missing,&d->zhd_missing,&d->zcmb_missing,&d->zhel_missing,&d->source_selection}) {
  if(!bounded(v->data,v->length,v->byte_length,10000)||!add(v->length,2))return IRRED_INVALID_INPUT;
 }
 if(!bounded(d->quality.data,d->quality.length,d->quality.byte_length,10000)||!add(d->quality.length,2*sizeof(uint64_t)))return IRRED_INVALID_INPUT;
 // Conservative supported STL construction capacity envelope. Moves do not
 // duplicate source buffers; caller input, allocator metadata and RSS excluded.
 *bytes=total;return IRRED_OK;
}
