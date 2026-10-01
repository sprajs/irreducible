#include "irred/abi.h"
#include "irred/sound_horizon.hpp"
#include "irred/quantities.hpp"
#include "payload_accounting.hpp"
#include <limits>
#include <cmath>
#include <memory>
#include <new>
namespace c = irred::cosmology;
namespace n = irred::numerics;
namespace {
template<class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p)%alignof(T)==0;
}
template<class T> bool header(const T *p) {
  return aligned(p) && p->struct_size==sizeof(T) && p->abi_version==IRRED_ABI_VERSION;
}
irred_bytes bytes(std::string_view s) {
  return {reinterpret_cast<const uint8_t*>(s.data()),s.size()};
}
}
struct irred_sound_horizon_result {
  c::SoundHorizonBatch native;
  std::vector<irred_sound_horizon_row> rows;
};
extern "C" uint32_t irred_sound_horizon_evaluate(
    const irred_sound_horizon_batch *batch,const irred_sound_horizon_policy *policy,
    irred_sound_horizon_result **out) {
  if (!aligned(out)) return IRRED_INVALID_INPUT;
  *out=nullptr;
  if (!aligned(batch)||!aligned(policy)) return IRRED_INVALID_INPUT;
  if (batch->abi_version!=IRRED_ABI_VERSION||policy->abi_version!=IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (!header(batch)||!header(policy)||batch->count>65536 ||
      batch->count>SIZE_MAX/sizeof(irred_sound_horizon_input)||
      batch->byte_length!=batch->count*sizeof(irred_sound_horizon_input)||
      (batch->count&&!aligned(batch->points))||policy->maximum_depth>60 ||
      policy->maximum_callbacks_per_point>SIZE_MAX||policy->maximum_points>SIZE_MAX||
      policy->maximum_total_callbacks>SIZE_MAX||policy->maximum_native_bytes>SIZE_MAX||policy->maximum_native_bytes>1073741824ULL||
      !std::isfinite(policy->absolute_tolerance_mpc)||!std::isfinite(policy->relative_tolerance)||
      policy->absolute_tolerance_mpc<0||policy->relative_tolerance<0||
      (policy->absolute_tolerance_mpc==0&&policy->relative_tolerance==0))
    return IRRED_INVALID_INPUT;
  try {
    size_t origins=0;
    for(size_t i=0;i<batch->count;++i) {
      const auto &p=batch->points[i];
      if(p.abi_version!=IRRED_ABI_VERSION) return IRRED_ABI_MISMATCH;
      if(!header(&p)||p.drag_origin.length>SIZE_MAX ||
         (p.drag_origin.length&&!p.drag_origin.data)||
         !irred::detail::checked_payload_add(origins,p.drag_origin.length,1)||
         !irred::detail::checked_payload_add(origins,1,1)) return IRRED_INVALID_INPUT;
    }
    const auto native_bound=c::sound_horizon_payload_bound(batch->count,origins);
    if(!native_bound) return IRRED_OVERFLOW;
    irred::detail::PayloadAccounting wrapper(sizeof(irred_sound_horizon_result));
    wrapper.add(batch->count,sizeof(c::SoundHorizonRequest));
    wrapper.add(batch->count,sizeof(irred_sound_horizon_row));
    wrapper.add(origins,2); // copied requests + conservative string capacity
    const auto own=wrapper.result();
    size_t combined=own.value_or(0);
    if(!own||!irred::detail::checked_payload_add(combined,1,*native_bound)) return IRRED_OVERFLOW;
    // Minimal diagnostic owner is permitted on quota rejection; no scientific
    // arrays or source snapshot are allocated, and no callback is consumed.
    auto result=std::make_unique<irred_sound_horizon_result>();
    if(batch->count>policy->maximum_points||combined>policy->maximum_native_bytes) {
      result->native.status=n::Status::work_limit;
      *out=result.release();return IRRED_OK;
    }
    std::vector<c::SoundHorizonRequest> requests;
    requests.reserve(batch->count);
    for(size_t i=0;i<batch->count;++i) {
      const auto &p=batch->points[i];
      c::SoundHorizonRequest r{{p.h0_km_s_mpc,p.omega_m,p.omega_r,p.omega_b,p.omega_gamma},p.z_drag,{}};
      if(p.drag_origin.length) r.drag_origin.assign(reinterpret_cast<const char*>(p.drag_origin.data),p.drag_origin.length);
      requests.push_back(std::move(r));
    }
    result->native=c::evaluate_sound_horizon(requests,
        {policy->absolute_tolerance_mpc,policy->relative_tolerance,
         static_cast<size_t>(policy->maximum_callbacks_per_point),
         static_cast<unsigned>(policy->maximum_depth),static_cast<size_t>(policy->maximum_points),
         static_cast<size_t>(policy->maximum_total_callbacks),
         static_cast<size_t>(policy->maximum_native_bytes)-*own});
    result->rows.reserve(result->native.rows.size());
    for(const auto &r:result->native.rows) {
      irred_sound_horizon_row row{};
      row.struct_size=sizeof(row);row.abi_version=IRRED_ABI_VERSION;
      const auto &m=r.source.model;
      row.source={sizeof(irred_sound_horizon_input),IRRED_ABI_VERSION,m.h0_km_s_mpc,
                  m.omega_m,m.omega_r,m.omega_b,m.omega_gamma,r.source.z_drag,bytes(r.source.drag_origin)};
      row.numerical_status=static_cast<uint32_t>(r.status);row.callbacks=r.callbacks;
      if(r.sound_horizon_mpc) {
        row.has_value=1;row.sound_horizon_mpc=*r.sound_horizon_mpc;
        row.error_estimate_mpc=r.error_estimate_mpc;
      }
      result->rows.push_back(row);
    }
    *out=result.release();return IRRED_OK;
  } catch(const std::bad_alloc&) {return IRRED_ALLOCATION_FAILURE;}
    catch(...) {return IRRED_EXCEPTION;}
}
extern "C" uint32_t irred_sound_horizon_result_view(
    const irred_sound_horizon_result *result,irred_sound_horizon_view *out) {
  if(!aligned(out)) return IRRED_INVALID_INPUT;
  *out={};
  if(!aligned(result)) return IRRED_INVALID_INPUT;
  out->struct_size=sizeof(*out);out->abi_version=IRRED_ABI_VERSION;
  out->numerical_status=static_cast<uint32_t>(result->native.status);
  out->rows=result->rows.empty()?nullptr:result->rows.data();out->count=result->rows.size();
  out->callbacks=result->native.callbacks;
  out->model_id=bytes(c::sound_horizon_model_id);out->equation_id=bytes(c::sound_horizon_equation_id);
  out->constant_set_id=bytes(irred::constant_set_id);
  out->coordinate_id=bytes("exact-scale-factor-rescaled-unit-interval");
  out->arithmetic_id=bytes("binary64-callback-wide64-or-more-nearest");return IRRED_OK;
}
extern "C" uint32_t irred_sound_horizon_result_destroy(irred_sound_horizon_result *p) {
  if(p && !aligned(p)) return IRRED_INVALID_INPUT;
  delete p;return IRRED_OK;
}
static_assert(sizeof(irred_sound_horizon_input)==72 && alignof(irred_sound_horizon_input)==8);
static_assert(sizeof(irred_sound_horizon_batch)==32 && alignof(irred_sound_horizon_batch)==8);
static_assert(sizeof(irred_sound_horizon_policy)==64 && alignof(irred_sound_horizon_policy)==8);
static_assert(sizeof(irred_sound_horizon_row)==120 && alignof(irred_sound_horizon_row)==8);
static_assert(sizeof(irred_sound_horizon_view)==120 && alignof(irred_sound_horizon_view)==8);
static_assert(offsetof(irred_sound_horizon_input,drag_origin)==56);
static_assert(offsetof(irred_sound_horizon_policy,maximum_native_bytes)==56);
static_assert(offsetof(irred_sound_horizon_row,callbacks)==112);
static_assert(offsetof(irred_sound_horizon_view,arithmetic_id)==104);
