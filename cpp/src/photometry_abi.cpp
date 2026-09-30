#include "irred/abi.h"
#include "irred/photometry.hpp"
#include "payload_accounting.hpp"
#include <cstddef>
#include <limits>
#include <memory>
#include <new>
#include <vector>
struct irred_photometry_result {
  irred::photometry::Batch native;
  std::vector<irred_photometry_row> rows;
};
namespace {
namespace p = irred::photometry;
namespace n = irred::numerics;
template<class T> bool aligned(const T *value) noexcept {
  return value && reinterpret_cast<std::uintptr_t>(value) % alignof(T) == 0;
}
template<class T> bool header(const T *value) noexcept {
  return aligned(value) && value->struct_size == sizeof(T) && value->abi_version == IRRED_ABI_VERSION;
}
using irred::detail::checked_payload_add;
p::Input input(const irred_photometry_input &x) noexcept {
  return {x.luminosity_watt_per_metre,x.rest_lower_metre,x.rest_upper_metre,
          x.observed_lower_metre,x.observed_upper_metre,x.luminosity_distance_metre,x.redshift,
          x.collecting_area_square_metre,x.optical_transmission,x.observer_exposure_second};
}
irred_photometry_input source(const p::Input &x) noexcept {
  return {sizeof(irred_photometry_input),IRRED_ABI_VERSION,x.luminosity_watt_per_metre,
          x.rest_lower_metre,x.rest_upper_metre,x.observed_lower_metre,x.observed_upper_metre,
          x.luminosity_distance_metre,x.redshift,x.collecting_area_square_metre,
          x.optical_transmission,x.observer_exposure_second};
}
irred_photometry_scalar scalar(const p::Outcome &x) noexcept {
  return {static_cast<std::uint32_t>(x.availability),static_cast<std::uint32_t>(x.numerical_status),x.value.value_or(0)};
}
irred_bytes bytes(std::string_view value) noexcept {
  return {reinterpret_cast<const std::uint8_t *>(value.data()),value.size()};
}
}
extern "C" uint32_t irred_photometry_evaluate(const irred_photometry_batch *batch,
    const irred_photometry_policy *policy, irred_photometry_result **output) {
  if (!aligned(output)) return IRRED_INVALID_INPUT;
  *output=nullptr;
  if (!aligned(batch) || !aligned(policy)) return IRRED_INVALID_INPUT;
  if (batch->abi_version!=IRRED_ABI_VERSION || policy->abi_version!=IRRED_ABI_VERSION) return IRRED_ABI_MISMATCH;
  if (!header(batch) || !header(policy) || policy->reserved || !policy->requested_outputs ||
      (policy->requested_outputs & ~7u) || policy->maximum_rows>65536 ||
      policy->maximum_native_bytes>(1ull<<30) || batch->length>65536 ||
      batch->length>SIZE_MAX/sizeof(irred_photometry_input) ||
      batch->byte_length!=batch->length*sizeof(irred_photometry_input) ||
      (batch->length && !aligned(batch->data))) return IRRED_INVALID_INPUT;
  for (std::uint64_t i=0;i<batch->length;++i) {
    if (batch->data[i].abi_version!=IRRED_ABI_VERSION) return IRRED_ABI_MISMATCH;
    if (!header(batch->data+i)) return IRRED_INVALID_INPUT;
  }
  const auto native_bound=p::output_payload_bound(batch->length);
  std::size_t combined=sizeof(irred_photometry_result)+sizeof(std::vector<p::Input>);
  if (!native_bound || !checked_payload_add(combined,1,*native_bound) ||
      !checked_payload_add(combined,batch->length,sizeof(p::Input)) ||
      !checked_payload_add(combined,batch->length,sizeof(irred_photometry_row))) return IRRED_INVALID_INPUT;
  try {
    auto owned=std::make_unique<irred_photometry_result>();
    if (batch->length>policy->maximum_rows || combined>policy->maximum_native_bytes) {
      owned->native.status=n::Status::work_limit; *output=owned.release();return IRRED_OK;
    }
    std::vector<p::Input> inputs;
    inputs.reserve(batch->length);
    for (std::uint64_t i=0;i<batch->length;++i) inputs.push_back(input(batch->data[i]));
    owned->native=p::evaluate(inputs,{policy->requested_outputs,static_cast<std::size_t>(policy->maximum_rows),*native_bound});
    if (owned->native.status==n::Status::invalid_input) return IRRED_INVALID_INPUT;
    owned->rows.reserve(owned->native.rows.size());
    for (const auto &row:owned->native.rows)
      owned->rows.push_back({sizeof(irred_photometry_row),IRRED_ABI_VERSION,source(row.source),
          static_cast<std::uint32_t>(row.admission_status),0,scalar(row.flux_watt_per_square_metre),
          scalar(row.energy_joule),scalar(row.expected_photons)});
    *output=owned.release();return IRRED_OK;
  } catch(const std::bad_alloc &) {return IRRED_ALLOCATION_FAILURE;}
    catch(...) {return IRRED_EXCEPTION;}
}
extern "C" uint32_t irred_photometry_result_view(const irred_photometry_result *result,
    irred_photometry_view *output) {
  if (!aligned(output)) return IRRED_INVALID_INPUT;
  *output={};
  if (!aligned(result)) return IRRED_INVALID_INPUT;
  *output={sizeof(*output),IRRED_ABI_VERSION,static_cast<std::uint32_t>(result->native.status),0,
      result->rows.data(),result->rows.size(),bytes(p::model_id),bytes(p::constants_id),
      bytes("binary64_storage_longdouble_intermediate"),
      bytes("isotropic_luminosity_distance_standard_redshift")};
  return IRRED_OK;
}
extern "C" uint32_t irred_photometry_result_destroy(irred_photometry_result *result) {
  if (result && !aligned(result)) return IRRED_INVALID_INPUT;
  delete result; return IRRED_OK;
}
static_assert(sizeof(irred_photometry_input)==88);
static_assert(sizeof(irred_photometry_batch)==32);
static_assert(sizeof(irred_photometry_policy)==32);
static_assert(sizeof(irred_photometry_scalar)==16);
static_assert(sizeof(irred_photometry_row)==152);
static_assert(sizeof(irred_photometry_view)==96);
