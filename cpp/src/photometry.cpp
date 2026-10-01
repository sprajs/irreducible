#include "irred/photometry.hpp"
#include "photometry_detail.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
namespace irred::photometry {
namespace {
using namespace detail;
bool requested(std::uint32_t mask, std::uint32_t bit) noexcept { return (mask & bit) != 0; }
numerics::Status admission(const Input &p) noexcept {
  for (double x : {p.luminosity_watt_per_metre, p.rest_lower_metre, p.rest_upper_metre,
                   p.observed_lower_metre, p.observed_upper_metre, p.luminosity_distance_metre,
                   p.redshift, p.collecting_area_square_metre, p.optical_transmission,
                   p.observer_exposure_second})
    if (!std::isfinite(x)) return numerics::Status::nonfinite_input;
  if (p.luminosity_watt_per_metre < 0 || p.rest_lower_metre <= 0 ||
      p.rest_upper_metre <= p.rest_lower_metre || p.observed_lower_metre <= 0 ||
      p.observed_upper_metre <= p.observed_lower_metre || p.luminosity_distance_metre <= 0 ||
      p.redshift < 0 || p.collecting_area_square_metre < 0 ||
      p.optical_transmission < 0 || p.optical_transmission > 1 || p.observer_exposure_second < 0)
    return numerics::Status::outside_domain;
  return numerics::Status::ok;
}
void evaluate_row(Row &row, std::uint32_t mask) noexcept {
  const Input &p = row.source;
  Outcome *groups[] = {&row.flux_watt_per_square_metre, &row.energy_joule, &row.expected_photons};
  const std::uint32_t bits[] = {incident_flux, collected_energy, transmitted_photons};
  row.admission_status = admission(p);
  if (row.admission_status != numerics::Status::ok) {
    for (unsigned i = 0; i < 3; ++i) if (requested(mask, bits[i])) fail(*groups[i], row.admission_status);
    return;
  }
  const bool zero_source = p.luminosity_watt_per_metre == 0;
  const bool zero_collection = zero_source || p.collecting_area_square_metre == 0 ||
      p.optical_transmission == 0 || p.observer_exposure_second == 0;
  if (requested(mask, incident_flux) && zero_source) store(*groups[0], 0);
  if (requested(mask, collected_energy) && zero_collection) store(*groups[1], 0);
  if (requested(mask, transmitted_photons) && zero_collection) store(*groups[2], 0);
  bool needs_geometry = false;
  for (unsigned i = 0; i < 3; ++i)
    needs_geometry |= requested(mask, bits[i]) && groups[i]->availability == Availability::omitted;
  if (!needs_geometry) return;
  const Interval redshift_factor = add(1, static_cast<Wide>(p.redshift));
  const Interval left = scale(redshift_factor, p.rest_lower_metre);
  const Interval right = scale(redshift_factor, p.rest_upper_metre);
  const Interval a{std::max(static_cast<Wide>(p.observed_lower_metre), left.lower),
                   std::max(static_cast<Wide>(p.observed_lower_metre), left.upper)};
  const Interval b{std::min(static_cast<Wide>(p.observed_upper_metre), right.lower),
                   std::min(static_cast<Wide>(p.observed_upper_metre), right.upper)};
  if (b.upper <= a.lower) {
    for (unsigned i = 0; i < 3; ++i)
      if (requested(mask, bits[i]) && groups[i]->availability == Availability::omitted) store(*groups[i], 0);
    return;
  }
  const Wide width_lower = add(b.lower, -a.upper).lower;
  const Wide width_upper = add(b.upper, -a.lower).upper;
  const Wide r = 1 + static_cast<Wide>(p.redshift);
  const Wide central_a = std::max(static_cast<Wide>(p.observed_lower_metre), r * p.rest_lower_metre);
  const Wide central_b = std::min(static_cast<Wide>(p.observed_upper_metre), r * p.rest_upper_metre);
  const Wide width = central_b - central_a;
  const Wide conservative_rounding = 128 * std::numeric_limits<Wide>::epsilon() +
                                    4 * std::numeric_limits<double>::epsilon();
  const Wide edge_uncertainty = width_lower > 0
      ? (width_upper - width_lower) / width_lower +
        (redshift_factor.upper - redshift_factor.lower) / redshift_factor.lower
      : std::numeric_limits<Wide>::infinity();
  if (!(width_lower > 0) || !(width > 0) ||
      edge_uncertainty + conservative_rounding > 2e-12L) {
    for (unsigned i = 0; i < 3; ++i)
      if (requested(mask, bits[i]) && groups[i]->availability == Availability::omitted)
        fail(*groups[i], numerics::Status::conditioning_budget_exceeded);
    return;
  }
  const Wide distance = p.luminosity_distance_metre;
  const Wide spectral_flux = detail::spectral_flux(p.luminosity_watt_per_metre, distance, r);
  const Wide band_flux = spectral_flux * width;
  if (requested(mask, incident_flux)) store(*groups[0], band_flux);
  const Wide collection = static_cast<Wide>(p.collecting_area_square_metre) *
                         p.optical_transmission * p.observer_exposure_second;
  if (requested(mask, collected_energy) && !zero_collection) store(*groups[1], collection * band_flux);
  if (requested(mask, transmitted_photons) && !zero_collection)
    store(*groups[2], collection * band_flux * (central_a + central_b) /
          (2 * planck_constant_joule_second * static_cast<Wide>(speed_of_light_m_per_s)));
}
} // namespace
std::optional<std::size_t> output_payload_bound(std::size_t rows) noexcept {
  if (rows > (std::numeric_limits<std::size_t>::max() - sizeof(Batch)) / sizeof(Row)) return {};
  return sizeof(Batch) + rows * sizeof(Row);
}
Batch evaluate(std::span<const Input> inputs, Policy policy) {
  Batch batch;
  if (policy.requested_outputs == 0 || (policy.requested_outputs & ~7u) ||
      policy.maximum_rows > 65536 || inputs.size() > 65536 ||
      std::numeric_limits<Wide>::digits < 64 || std::numeric_limits<Wide>::max_exponent < 16384 ||
      std::fegetround() != FE_TONEAREST) return batch;
  const auto bytes = output_payload_bound(inputs.size());
  if (!bytes) return batch;
  if (inputs.size() > policy.maximum_rows || *bytes > policy.maximum_output_bytes) {
    batch.status = numerics::Status::work_limit; return batch;
  }
  batch.rows.resize(inputs.size());
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    batch.rows[i].source = inputs[i];
    evaluate_row(batch.rows[i], policy.requested_outputs);
  }
  batch.status = numerics::Status::ok;
  return batch;
}
} // namespace irred::photometry
