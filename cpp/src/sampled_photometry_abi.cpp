#include "irred/abi.h"
#include "irred/sampled_photometry.hpp"
#include "payload_accounting.hpp"
#include <cstring>
#include <memory>
#include <new>
#include <vector>
namespace {
namespace p = irred::photometry;
namespace n = irred::numerics;
struct Curve {
  std::vector<std::uint8_t> id;
  std::vector<double> wavelength, values;
};
template <class T> bool aligned(const T *x) noexcept {
  return x && reinterpret_cast<std::uintptr_t>(x) % alignof(T) == 0;
}
template <class T> bool header(const T *x) noexcept {
  return aligned(x) && x->struct_size == sizeof(T) &&
         x->abi_version == IRRED_ABI_VERSION;
}
template <class T>
bool array(const T *x, std::uint64_t count, std::uint64_t bytes,
           std::uint64_t cap) noexcept {
  return count <= cap && count <= SIZE_MAX / sizeof(T) &&
         bytes == count * sizeof(T) && (!count || aligned(x));
}
irred_bytes bytes(std::string_view s) noexcept {
  return {reinterpret_cast<const std::uint8_t *>(s.data()), s.size()};
}
irred_photometry_scalar scalar(const p::Outcome &x, n::Status admission,
                               bool requested) noexcept {
  if (requested && x.availability == p::Availability::omitted &&
      admission != n::Status::ok)
    return {2, static_cast<unsigned>(admission), 0};
  return {static_cast<unsigned>(x.availability),
          static_cast<unsigned>(x.numerical_status), x.value.value_or(0)};
}
} // namespace
struct irred_sampled_photometry_result {
  n::Status status = n::Status::ok;
  std::vector<Curve> spectra, passbands;
  std::vector<irred_sampled_photometry_curve> spectrum_views, passband_views;
  std::vector<irred_sampled_photometry_row> rows;
};
namespace {
void copy_pool(std::vector<Curve> &owners,
               std::vector<irred_sampled_photometry_curve> &views,
               const irred_sampled_photometry_curve *source,
               std::size_t count) {
  owners.resize(count);
  views.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    auto &o = owners[i];
    const auto &s = source[i];
    o.id.assign(s.id.data, s.id.data + s.id.length);
    if (s.length) {
      o.wavelength.assign(s.wavelength_metre, s.wavelength_metre + s.length);
      o.values.assign(s.values, s.values + s.length);
    }
    views[i] = {sizeof(irred_sampled_photometry_curve),
                IRRED_ABI_VERSION,
                {o.id.data(), o.id.size()},
                o.wavelength.data(),
                o.values.data(),
                s.length,
                s.byte_length};
  }
}
} // namespace
extern "C" uint32_t
irred_sampled_photometry_evaluate(const irred_sampled_photometry_batch *b,
                                  const irred_sampled_photometry_policy *q,
                                  irred_sampled_photometry_result **output) {
  if (!aligned(output))
    return IRRED_INVALID_INPUT;
  *output = nullptr;
  if (!aligned(b) || !aligned(q))
    return IRRED_INVALID_INPUT;
  if (b->abi_version != IRRED_ABI_VERSION ||
      q->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (!header(b) || !header(q) || q->reserved || !q->requested_outputs ||
      (q->requested_outputs & ~7u) || q->maximum_rows > 65536 ||
      q->maximum_native_bytes > (1ull << 30) ||
      q->maximum_total_samples > 65536 || q->maximum_samples > 65536 ||
      q->maximum_segments > 65536 ||
      !array(b->spectra, b->spectrum_count, b->spectrum_byte_length, 4096) ||
      !array(b->passbands, b->passband_count, b->passband_byte_length, 4096) ||
      !array(b->exposures, b->exposure_count, b->exposure_byte_length, 65536))
    return IRRED_INVALID_INPUT;
  std::size_t total = 0, payload = sizeof(irred_sampled_photometry_result);
  using irred::detail::checked_payload_add;
  if (!checked_payload_add(payload, b->exposure_count,
                           sizeof(irred_sampled_photometry_row)))
    return IRRED_INVALID_INPUT;
  // Validate every descriptor/count before scanning any sample values or
  // allocating.
  for (unsigned pool = 0; pool < 2; ++pool) {
    const auto *c = pool ? b->passbands : b->spectra;
    const auto count = pool ? b->passband_count : b->spectrum_count;
    if (!checked_payload_add(payload, count,
                             sizeof(Curve) +
                                 sizeof(irred_sampled_photometry_curve)))
      return IRRED_INVALID_INPUT;
    for (std::uint64_t i = 0; i < count; ++i) {
      if (c[i].abi_version != IRRED_ABI_VERSION)
        return IRRED_ABI_MISMATCH;
      if (!header(c + i) || !c[i].id.length || c[i].id.length > 256 ||
          !c[i].id.data ||
          !array(c[i].wavelength_metre, c[i].length, c[i].byte_length, 65536) ||
          (c[i].length && !aligned(c[i].values)) || c[i].length > 65536 - total)
        return IRRED_INVALID_INPUT;
      total += c[i].length;
      if (!checked_payload_add(payload, c[i].id.length, 1) ||
          !checked_payload_add(payload, c[i].length, 2 * sizeof(double)))
        return IRRED_INVALID_INPUT;
    }
  }
  for (std::uint64_t i = 0; i < b->exposure_count; ++i) {
    if (b->exposures[i].abi_version != IRRED_ABI_VERSION)
      return IRRED_ABI_MISMATCH;
    if (!header(b->exposures + i))
      return IRRED_INVALID_INPUT;
  }
  // Bounded identifier scans; no temporary hash/matrix allocation.
  for (unsigned pool = 0; pool < 2; ++pool) {
    const auto *c = pool ? b->passbands : b->spectra;
    const auto count = pool ? b->passband_count : b->spectrum_count;
    for (std::uint64_t i = 0; i < count; ++i)
      for (std::uint64_t j = 0; j < i; ++j)
        if (c[i].id.length == c[j].id.length &&
            !std::memcmp(c[i].id.data, c[j].id.data, c[i].id.length))
          return IRRED_INVALID_INPUT;
  }
  try {
    auto owner = std::make_unique<irred_sampled_photometry_result>();
    if (b->exposure_count > q->maximum_rows ||
        total > q->maximum_total_samples || payload > q->maximum_native_bytes) {
      owner->status = n::Status::work_limit;
      *output = owner.release();
      return IRRED_OK;
    }
    copy_pool(owner->spectra, owner->spectrum_views, b->spectra,
              b->spectrum_count);
    copy_pool(owner->passbands, owner->passband_views, b->passbands,
              b->passband_count);
    owner->rows.resize(b->exposure_count);
    for (std::uint64_t i = 0; i < b->exposure_count; ++i) {
      const auto &e = b->exposures[i];
      p::SampledResult r;
      if (e.spectrum_index < b->spectrum_count &&
          e.passband_index < b->passband_count) {
        const auto &s = owner->spectra[e.spectrum_index],
                   &t = owner->passbands[e.passband_index];
        r = p::evaluate_sampled(
            {{s.wavelength, s.values},
             {t.wavelength, t.values},
             e.luminosity_distance_metre,
             e.redshift,
             e.collecting_area_square_metre,
             e.observer_exposure_second},
            {q->requested_outputs, static_cast<std::size_t>(q->maximum_samples),
             static_cast<std::size_t>(q->maximum_segments)});
      }
      owner->rows[i] = {
          sizeof(irred_sampled_photometry_row),
          IRRED_ABI_VERSION,
          e,
          static_cast<unsigned>(r.admission_status),
          0,
          scalar(r.flux_watt_per_square_metre, r.admission_status,
                 q->requested_outputs & 1),
          scalar(r.energy_joule, r.admission_status, q->requested_outputs & 2),
          scalar(r.expected_photons, r.admission_status,
                 q->requested_outputs & 4)};
    }
    *output = owner.release();
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t
irred_sampled_photometry_result_view(const irred_sampled_photometry_result *r,
                                     irred_sampled_photometry_view *output) {
  if (!aligned(output))
    return IRRED_INVALID_INPUT;
  *output = {};
  if (!aligned(r))
    return IRRED_INVALID_INPUT;
  *output = {sizeof(*output),
             IRRED_ABI_VERSION,
             static_cast<unsigned>(r->status),
             0,
             r->spectrum_views.data(),
             r->spectrum_views.size(),
             r->passband_views.data(),
             r->passband_views.size(),
             r->rows.data(),
             r->rows.size(),
             bytes(p::sampled_model_id),
             bytes(p::constants_id),
             bytes("analytic_piecewise_linear_merged_bernstein"),
             bytes("binary64_storage_longdouble_intermediate"),
             bytes("isotropic_luminosity_distance_standard_redshift")};
  return IRRED_OK;
}
extern "C" uint32_t
irred_sampled_photometry_result_destroy(irred_sampled_photometry_result *r) {
  if (r && !aligned(r))
    return IRRED_INVALID_INPUT;
  delete r;
  return IRRED_OK;
}
