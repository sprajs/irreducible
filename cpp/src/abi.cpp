#include "cosmology/abi.h"
#include "cosmology/arithmetic.hpp"
#include <vector>
#include <limits>
#include <new>
#include <stdexcept>
#include <cstddef>
struct cosmo_result { std::vector<int64_t> values; };
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
  *out=new cosmo_result{std::move(values)};
  return COSMO_OK;
 } catch(const std::bad_alloc&) {return COSMO_ALLOCATION_FAILURE;} catch(...) {return COSMO_EXCEPTION;}
}
extern "C" uint32_t cosmo_result_view(const cosmo_result* r,const int64_t** data,uint64_t* n) {
 if(!data || !n) return COSMO_INVALID_INPUT;
 *data=nullptr; *n=0;
 if(!r) return COSMO_INVALID_INPUT;
 *data=r->values.data(); *n=r->values.size(); return COSMO_OK;
}
extern "C" uint32_t cosmo_result_destroy(cosmo_result* r) {delete r; return COSMO_OK;}
