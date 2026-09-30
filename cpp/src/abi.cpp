#include "cosmology/abi.h"
#include "cosmology/arithmetic.hpp"
#include "cosmology/quantities.hpp"
#include "cosmology/quantity_enum_checks.inc"
#include <vector>
#include <limits>
#include <new>
#include <stdexcept>
#include <cstddef>
struct cosmo_result { std::vector<int64_t> values; std::vector<double> f64_values; std::vector<uint32_t> statuses; uint32_t kind=1; };
static uint32_t validate(const cosmo_i64_buffer* b) noexcept {
 if (!b) return COSMO_INVALID_INPUT;
 if (b->abi_version!=COSMO_ABI_VERSION) return COSMO_ABI_MISMATCH;
 if (b->struct_size!=sizeof(*b) || b->reserved || b->element_type!=1) return COSMO_INVALID_INPUT;
 if (b->length>std::numeric_limits<size_t>::max()/sizeof(int64_t) || b->byte_length!=b->length*sizeof(int64_t)) return COSMO_INVALID_INPUT;
 if (b->length && (!b->data || reinterpret_cast<uintptr_t>(b->data)%alignof(int64_t))) return COSMO_INVALID_INPUT;
 return COSMO_OK;
}
extern "C" uint32_t cosmo_add(const cosmo_i64_buffer* a,const cosmo_i64_buffer* b,uint32_t fault,cosmo_result** out) {
 if (!out) return COSMO_INVALID_INPUT;
 *out=nullptr;
 const auto va=validate(a); if(va) return va;
 const auto vb=validate(b); if(vb) return vb;
 if(a->length!=b->length || fault>2) return COSMO_INVALID_INPUT;
 try {
  if(fault==1) throw std::runtime_error("injected");
  if(fault==2) throw std::bad_alloc();
  std::vector<int64_t> values; values.reserve(static_cast<size_t>(a->length));
  for(uint64_t i=0;i<a->length;++i) {
   const auto x=a->data[i],y=b->data[i];
   int64_t sum; if(!cosmology::checked_add(x,y,sum)) return COSMO_OVERFLOW;
   values.push_back(sum);
  }
  *out=new cosmo_result{std::move(values),{},{},1};
  return COSMO_OK;
 } catch(const std::bad_alloc&) {return COSMO_ALLOCATION_FAILURE;} catch(...) {return COSMO_EXCEPTION;}
}
extern "C" uint32_t cosmo_result_view(const cosmo_result* r,const int64_t** data,uint64_t* n) {
 if(!data || !n) return COSMO_INVALID_INPUT;
 *data=nullptr; *n=0;
 if(!r || r->kind!=1) return COSMO_INVALID_INPUT;
 *data=r->values.data(); *n=r->values.size(); return COSMO_OK;
}
extern "C" uint32_t cosmo_result_destroy(cosmo_result* r) {delete r; return COSMO_OK;}

static uint32_t validate_metadata(const cosmo_quantity_metadata* m) noexcept {
 if(!m) return COSMO_INVALID_INPUT;
 if(m->abi_version!=COSMO_ABI_VERSION) return COSMO_ABI_MISMATCH;
 if(m->struct_size!=sizeof(*m) || m->reserved) return COSMO_INVALID_INPUT;
 return COSMO_OK;
}
static std::string_view constants(const cosmo_quantity_metadata& m) noexcept {
 return m.constant_set==COSMO_CONSTANT_SET_SI_IAU_DEFINITIONS_V1 ? cosmology::constant_set_id : std::string_view{};
}
extern "C" uint32_t cosmo_convert_quantities(const cosmo_f64_buffer* values,const cosmo_quantity_metadata* source,const cosmo_quantity_metadata* target,cosmo_result** out) {
 if(!out) return COSMO_INVALID_INPUT;
 *out=nullptr;
 if(!values) return COSMO_INVALID_INPUT;
 if(values->abi_version!=COSMO_ABI_VERSION) return COSMO_ABI_MISMATCH;
 if(values->struct_size!=sizeof(*values) || values->element_type!=2 || values->reserved) return COSMO_INVALID_INPUT;
 if(values->length>std::numeric_limits<size_t>::max()/sizeof(double) || values->byte_length!=values->length*sizeof(double)) return COSMO_INVALID_INPUT;
 if(values->length && (!values->data || reinterpret_cast<uintptr_t>(values->data)%alignof(double))) return COSMO_INVALID_INPUT;
 auto status=validate_metadata(source); if(status) return status;
 status=validate_metadata(target); if(status) return status;
 const cosmology::Target destination{static_cast<cosmology::Unit>(target->unit),static_cast<cosmology::Role>(target->role),static_cast<cosmology::Frame>(target->frame),static_cast<cosmology::LengthConvention>(target->convention),constants(*target)};
 const auto input=[&](double value){return cosmology::Quantity{value,static_cast<cosmology::Unit>(source->unit),static_cast<cosmology::Role>(source->role),static_cast<cosmology::Frame>(source->frame),static_cast<cosmology::LengthConvention>(source->convention),constants(*source)};};
 // Empty arrays still validate scientific metadata, with a normal-domain probe.
 if(values->length==0 && cosmology::convert(input(1.0),destination).status!=cosmology::QuantityStatus::ok) return COSMO_INVALID_INPUT;
 try {
  std::vector<double> converted; std::vector<uint32_t> statuses;
  converted.reserve(static_cast<size_t>(values->length)); statuses.reserve(static_cast<size_t>(values->length));
  for(uint64_t i=0;i<values->length;++i) {
   const auto result=cosmology::convert(input(values->data[i]),destination);
   statuses.push_back(static_cast<uint32_t>(result.status));
   // Failed slots have no public numeric meaning; their status must be consumed.
   converted.push_back(result.status==cosmology::QuantityStatus::ok ? result.target.value : 0.0);
  }
  *out=new cosmo_result{{},std::move(converted),std::move(statuses),2}; return COSMO_OK;
 } catch(const std::bad_alloc&) {return COSMO_ALLOCATION_FAILURE;} catch(...) {return COSMO_EXCEPTION;}
}
extern "C" uint32_t cosmo_result_f64_view(const cosmo_result* r,const double** data,const uint32_t** status,uint64_t* n) {
 if(!data || !status || !n) return COSMO_INVALID_INPUT;
 *data=nullptr; *status=nullptr; *n=0;
 if(!r || r->kind!=2) return COSMO_INVALID_INPUT;
 *data=r->f64_values.data(); *status=r->statuses.data(); *n=r->f64_values.size(); return COSMO_OK;
}
