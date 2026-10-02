#include "irred/temporal_photometry.hpp"
#include "payload_accounting.hpp"
#include "photometry_detail.hpp"
#include <algorithm>
#include <cfenv>
#include <utility>
namespace irred::photometry {
namespace {
using S = numerics::Status;
using detail::Interval;
using detail::Wide;
constexpr std::size_t hard_knots = 65536, hard_values = 4000000,
                      hard_rows = 65536, hard_pools = 4096,
                      hard_bytes = 1ull << 30;
struct OwnedGrid {
  std::vector<char> id, origin;
  std::vector<double> time, wavelength, luminosity;
};
struct OwnedBand {
  std::vector<char> id, origin;
  std::vector<double> wavelength, transmission;
};
std::string_view text(const std::vector<char> &s) noexcept {
  return {s.data(), s.size()};
}
bool text_valid(std::string_view s, std::size_t max) noexcept {
  return !s.empty() && s.size() <= max;
}
bool increase(std::span<const double> a, bool positive) noexcept {
  if (a.size() < 2)
    return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (!std::isfinite(a[i]) || (positive && !(a[i] > 0)) ||
        (i && !(a[i] > a[i - 1])))
      return false;
  return true;
}
Interval sum(Interval a, Interval b) noexcept {
  return {detail::add(a.lower, b.lower).lower,
          detail::add(a.upper, b.upper).upper};
}
Interval subtract(Interval a, Interval b) noexcept {
  return {detail::add(a.lower, -b.upper).lower,
          detail::add(a.upper, -b.lower).upper};
}
Interval product(Interval a, Interval b) noexcept {
  return {detail::multiply(a.lower, b.lower).lower,
          detail::multiply(a.upper, b.upper).upper};
}
Interval divide(Interval a, Interval b) noexcept {
  return {a.lower == 0 ? 0 : std::nextafter(a.lower / b.upper, -INFINITY),
          a.upper == 0 ? 0 : std::nextafter(a.upper / b.lower, INFINITY)};
}
bool tight(Interval a) noexcept {
  return a.lower >= 0 && std::isfinite(a.upper) &&
         (a.upper == 0 ||
          (a.lower > 0 && (a.upper - a.lower) / a.lower <= 1e-13L));
}
Interval mapped(double epoch, double time, Interval r) noexcept {
  const auto a = detail::multiply(r.lower, time),
             b = detail::multiply(r.upper, time);
  return {detail::add(epoch, std::min(a.lower, b.lower)).lower,
          detail::add(epoch, std::max(a.upper, b.upper)).upper};
}
Interval interpolate(Interval at, Interval lo, Interval hi, double a,
                     double b) noexcept {
  const auto width = subtract(hi, lo), offset = subtract(at, lo);
  auto fraction = divide(offset, width);
  fraction.lower = std::max(Wide{0}, fraction.lower);
  fraction.upper = std::min(Wide{1}, fraction.upper);
  const Interval reverse{
      std::max(Wide{0}, detail::add(1, -fraction.upper).lower),
      std::min(Wide{1}, detail::add(1, -fraction.lower).upper)};
  return sum(detail::scale(reverse, a), detail::scale(fraction, b));
}
void fail_row(TemporalRow &row, S status, std::uint32_t mask) noexcept {
  row.admission_status = status;
  Outcome *groups[] = {&row.mean_flux_watt_per_square_metre, &row.energy_joule,
                       &row.expected_photons};
  for (unsigned i = 0; i < 3; ++i)
    if (mask & (1u << i))
      detail::fail(*groups[i], status);
}
S scalar_admission(const TemporalExposure &e) noexcept {
  for (double x :
       {e.luminosity_distance_metre, e.redshift, e.collecting_area_square_metre,
        e.source_epoch_observer_second, e.observer_lower_second,
        e.observer_upper_second})
    if (!std::isfinite(x))
      return S::nonfinite_input;
  if (!(e.luminosity_distance_metre > 0) || e.redshift < 0 ||
      e.collecting_area_square_metre < 0 ||
      !(e.observer_upper_second > e.observer_lower_second))
    return S::outside_domain;
  return S::ok;
}
S duration(Interval x, double &value) noexcept {
  if (!tight(x) || !(x.lower > 0))
    return S::conditioning_budget_exceeded;
  Outcome out;
  detail::store(out, x.lower + (x.upper - x.lower) / 2);
  if (!out.value)
    return out.numerical_status;
  const Wide error = std::abs(Wide(*out.value) - x.lower) / x.lower;
  if (error + (x.upper - x.lower) / x.lower > 1e-13L)
    return S::conditioning_budget_exceeded;
  value = *out.value;
  return S::ok;
}
} // namespace
struct TemporalStorage {
  std::vector<OwnedGrid> grids;
  std::vector<OwnedBand> bands;
  std::size_t payload = 0;
};
PreparedTemporal::PreparedTemporal(PreparedTemporal &&other) noexcept
    : status_(other.status_), storage_(std::move(other.storage_)) {
  other.status_ = S::invalid_input;
}
PreparedTemporal &
PreparedTemporal::operator=(PreparedTemporal &&other) noexcept {
  if (this != &other) {
    status_ = other.status_;
    storage_ = std::move(other.storage_);
    other.status_ = S::invalid_input;
  }
  return *this;
}
std::size_t PreparedTemporal::grid_count() const noexcept {
  return storage_ ? storage_->grids.size() : 0;
}
std::size_t PreparedTemporal::band_count() const noexcept {
  return storage_ ? storage_->bands.size() : 0;
}
std::optional<TemporalGrid>
PreparedTemporal::grid(std::size_t i) const noexcept {
  if (!storage_ || i >= storage_->grids.size())
    return {};
  const auto &g = storage_->grids[i];
  return TemporalGrid{text(g.id), text(g.origin), g.time, g.wavelength,
                      g.luminosity};
}
std::optional<TemporalBand>
PreparedTemporal::band(std::size_t i) const noexcept {
  if (!storage_ || i >= storage_->bands.size())
    return {};
  const auto &b = storage_->bands[i];
  return TemporalBand{
      text(b.id), text(b.origin), {b.wavelength, b.transmission}};
}
std::optional<std::size_t>
PreparedTemporal::retained_payload_bytes() const noexcept {
  return storage_ ? std::optional(storage_->payload) : std::nullopt;
}
PreparedTemporal prepare_temporal(std::span<const TemporalGrid> grids,
                                  std::span<const TemporalBand> bands,
                                  TemporalPreparationPolicy q) {
  PreparedTemporal result;
  auto fail = [&](S s) {
    result.status_ = s;
    return std::move(result);
  };
  if (q.maximum_grids > hard_pools || q.maximum_bands > hard_pools ||
      q.maximum_total_knots > hard_knots ||
      q.maximum_grid_values > hard_values ||
      q.maximum_native_bytes > hard_bytes)
    return fail(S::invalid_input);
  if (grids.size() > hard_pools || bands.size() > hard_pools ||
      grids.size() > q.maximum_grids || bands.size() > q.maximum_bands)
    return fail(S::work_limit);
  if (grids.empty() || bands.empty())
    return fail(S::invalid_input);
  std::size_t knots = 0, values = 0,
              bytes = sizeof(PreparedTemporal) + sizeof(TemporalStorage);
  using irred::detail::checked_payload_add;
  if (!checked_payload_add(bytes, grids.size(), sizeof(OwnedGrid)) ||
      !checked_payload_add(bytes, bands.size(), sizeof(OwnedBand)))
    return fail(S::work_limit);
  // Descriptor/count admission precedes all scientific array scans and copies.
  for (const auto &g : grids) {
    if (!text_valid(g.id, 256) || !text_valid(g.source_origin, 1024))
      return fail(S::invalid_input);
    const auto nt = g.rest_time_second.size(),
               nw = g.rest_wavelength_metre.size();
    if (nt > hard_knots || nw > hard_knots || nt > hard_knots - knots)
      return fail(S::work_limit);
    knots += nt;
    if (nw > hard_knots - knots)
      return fail(S::work_limit);
    knots += nw;
    if (nw && nt > hard_values / nw)
      return fail(S::work_limit);
    const auto count = nt * nw;
    if (count != g.luminosity_watt_per_metre.size())
      return fail(S::invalid_input);
    if (count > hard_values - values)
      return fail(S::work_limit);
    values += count;
    if (!checked_payload_add(bytes, g.id.size() + g.source_origin.size(), 1) ||
        !checked_payload_add(bytes, nt + nw + count, sizeof(double)))
      return fail(S::work_limit);
  }
  for (const auto &b : bands) {
    if (!text_valid(b.id, 256) || !text_valid(b.calibration_origin, 1024))
      return fail(S::invalid_input);
    const auto n = b.passband.wavelength_metre.size();
    if (n > hard_knots - knots)
      return fail(S::work_limit);
    knots += n;
    if (n != b.passband.optical_transmission.size())
      return fail(S::invalid_input);
    if (!checked_payload_add(bytes, b.id.size() + b.calibration_origin.size(),
                             1) ||
        !checked_payload_add(bytes, n, 2 * sizeof(double)))
      return fail(S::work_limit);
  }
  if (knots > q.maximum_total_knots || values > q.maximum_grid_values ||
      bytes > q.maximum_native_bytes)
    return fail(S::work_limit);
  for (std::size_t i = 0; i < grids.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j)
      if (grids[i].id == grids[j].id)
        return fail(S::invalid_input);
    const auto &g = grids[i];
    for (double v : g.rest_time_second)
      if (!std::isfinite(v))
        return fail(S::nonfinite_input);
    for (double v : g.rest_wavelength_metre)
      if (!std::isfinite(v))
        return fail(S::nonfinite_input);
    if (!increase(g.rest_time_second, false) ||
        !increase(g.rest_wavelength_metre, true))
      return fail(S::outside_domain);
    for (double v : g.luminosity_watt_per_metre) {
      if (!std::isfinite(v))
        return fail(S::nonfinite_input);
      if (v < 0)
        return fail(S::outside_domain);
    }
  }
  for (std::size_t i = 0; i < bands.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j)
      if (bands[i].id == bands[j].id)
        return fail(S::invalid_input);
    const auto &b = bands[i].passband;
    for (double v : b.wavelength_metre)
      if (!std::isfinite(v))
        return fail(S::nonfinite_input);
    if (!increase(b.wavelength_metre, true))
      return fail(S::outside_domain);
    for (double v : b.optical_transmission) {
      if (!std::isfinite(v))
        return fail(S::nonfinite_input);
      if (v < 0 || v > 1)
        return fail(S::outside_domain);
    }
  }
  auto storage = std::make_shared<TemporalStorage>();
  storage->payload = bytes;
  storage->grids.resize(grids.size());
  storage->bands.resize(bands.size());
  for (std::size_t i = 0; i < grids.size(); ++i) {
    const auto &g = grids[i];
    auto &o = storage->grids[i];
    o.id.assign(g.id.begin(), g.id.end());
    o.origin.assign(g.source_origin.begin(), g.source_origin.end());
    o.time.assign(g.rest_time_second.begin(), g.rest_time_second.end());
    o.wavelength.assign(g.rest_wavelength_metre.begin(),
                        g.rest_wavelength_metre.end());
    o.luminosity.assign(g.luminosity_watt_per_metre.begin(),
                        g.luminosity_watt_per_metre.end());
  }
  for (std::size_t i = 0; i < bands.size(); ++i) {
    const auto &b = bands[i];
    auto &o = storage->bands[i];
    o.id.assign(b.id.begin(), b.id.end());
    o.origin.assign(b.calibration_origin.begin(), b.calibration_origin.end());
    o.wavelength.assign(b.passband.wavelength_metre.begin(),
                        b.passband.wavelength_metre.end());
    o.transmission.assign(b.passband.optical_transmission.begin(),
                          b.passband.optical_transmission.end());
  }
  result.storage_ = std::move(storage);
  result.status_ = S::ok;
  return result;
}
TemporalBatch evaluate_temporal(const PreparedTemporal &prepared,
                                std::span<const TemporalExposure> exposures,
                                TemporalPolicy q) {
  TemporalBatch out;
  if (!q.requested_outputs || (q.requested_outputs & ~7u) ||
      q.maximum_rows > hard_rows || q.maximum_native_bytes > hard_bytes ||
      q.maximum_segment_work > hard_values)
    return out;
  if (prepared.status_ != S::ok || !prepared.storage_) {
    out.status = prepared.status_;
    return out;
  }
  if (std::numeric_limits<Wide>::digits < 64 ||
      std::numeric_limits<Wide>::max_exponent < 16384 ||
      std::fegetround() != FE_TONEAREST)
    return out;
  if (exposures.size() > hard_rows || exposures.size() > q.maximum_rows) {
    out.status = S::work_limit;
    return out;
  }
  const auto &storage = *prepared.storage_;
  std::size_t max_time = 0, max_wave = 0;
  // Prepared pools may serve many consumers. Unreferenced grids do not impose
  // scratch on this batch; invalid references need only their failed row.
  for (const auto &e : exposures) {
    if (e.grid_index >= storage.grids.size() ||
        e.band_index >= storage.bands.size())
      continue;
    const auto &g = storage.grids[e.grid_index];
    max_time = std::max(max_time, g.time.size());
    max_wave = std::max(max_wave, g.wavelength.size());
  }
  std::size_t payload = sizeof(TemporalBatch) +
                        2 * sizeof(std::vector<Interval>) +
                        sizeof(std::vector<double>);
  using irred::detail::checked_payload_add;
  if (!checked_payload_add(payload, exposures.size(), sizeof(TemporalRow)) ||
      !checked_payload_add(payload, max_time, sizeof(Interval)) ||
      !checked_payload_add(payload, max_wave,
                           sizeof(Interval) + sizeof(double)) ||
      payload > q.maximum_native_bytes) {
    out.status = S::work_limit;
    return out;
  }
  std::vector<Interval> knots(max_time), integrated(max_wave);
  std::vector<double> mean(max_wave);
  out.rows.resize(exposures.size());
  out.status = S::ok;
  for (std::size_t index = 0; index < exposures.size(); ++index) {
    auto &row = out.rows[index];
    const auto &e = exposures[index];
    row.source = e;
    if (e.grid_index >= storage.grids.size() ||
        e.band_index >= storage.bands.size()) {
      fail_row(row, S::invalid_input, q.requested_outputs);
      continue;
    }
    const auto &g = storage.grids[e.grid_index];
    const auto &band = storage.bands[e.band_index];
    auto admission = scalar_admission(e);
    if (admission != S::ok) {
      fail_row(row, admission, q.requested_outputs);
      continue;
    }
    const auto nt = g.time.size(), nw = g.wavelength.size(),
               nb = band.wavelength.size();
    if (nw + nb > hard_knots) {
      fail_row(row, S::work_limit, q.requested_outputs);
      continue;
    }
    const std::size_t work = (nt - 1) * nw + nw + nb - 3;
    if (work > q.maximum_segment_work - out.segment_work) {
      fail_row(row, S::work_limit, q.requested_outputs);
      continue;
    }
    row.segment_work = work;
    out.segment_work += work;
    const Interval full = detail::add(Wide(e.observer_upper_second),
                                      -Wide(e.observer_lower_second));
    double observer_duration = 0;
    admission = duration(full, observer_duration);
    if (admission != S::ok) {
      fail_row(row, admission, q.requested_outputs);
      continue;
    }
    row.observer_duration_second = full.lower + (full.upper - full.lower) / 2;
    const auto r = detail::add(1, e.redshift);
    for (std::size_t i = 0; i < nt; ++i)
      knots[i] = mapped(e.source_epoch_observer_second, g.time[i], r);
    const Interval exposure_lo{e.observer_lower_second,
                               e.observer_lower_second},
        exposure_hi{e.observer_upper_second, e.observer_upper_second};
    const Interval covered_lo{std::max(exposure_lo.lower, knots[0].lower),
                              std::max(exposure_lo.upper, knots[0].upper)},
        covered_hi{std::min(exposure_hi.lower, knots[nt - 1].lower),
                   std::min(exposure_hi.upper, knots[nt - 1].upper)};
    if (covered_hi.upper <= covered_lo.lower) {
      row.coverage = TemporalCoverage::no_overlap;
      row.admission_status = S::ok;
      Outcome *groups[] = {&row.mean_flux_watt_per_square_metre,
                           &row.energy_joule, &row.expected_photons};
      for (unsigned i = 0; i < 3; ++i)
        if (q.requested_outputs & (1u << i))
          detail::store(*groups[i], 0);
      continue;
    }
    const auto covered = subtract(covered_hi, covered_lo);
    if (!tight(covered) || !(covered.lower > 0)) {
      fail_row(row, S::conditioning_budget_exceeded, q.requested_outputs);
      continue;
    }
    row.covered_observer_second =
        covered.lower + (covered.upper - covered.lower) / 2;
    row.covered_fraction =
        row.covered_observer_second / row.observer_duration_second;
    if (exposure_lo.lower >= knots[0].upper &&
        exposure_hi.upper <= knots[nt - 1].lower) {
      row.coverage = TemporalCoverage::full;
    } else if (exposure_lo.lower < knots[0].lower ||
               exposure_hi.upper > knots[nt - 1].upper) {
      row.coverage = TemporalCoverage::partial;
    } else {
      // Unresolved support-edge inclusion cannot be labelled full or partial.
      fail_row(row, S::conditioning_budget_exceeded, q.requested_outputs);
      continue;
    }
    std::fill(integrated.begin(), integrated.begin() + nw, Interval{});
    bool temporal_ok = true;
    for (std::size_t i = 0; i + 1 < nt; ++i) {
      const Interval a{std::max(exposure_lo.lower, knots[i].lower),
                       std::max(exposure_lo.upper, knots[i].upper)},
          b{std::min(exposure_hi.lower, knots[i + 1].lower),
            std::min(exposure_hi.upper, knots[i + 1].upper)};
      if (b.upper <= a.lower)
        continue;
      const auto width = subtract(b, a),
                 time_width = subtract(knots[i + 1], knots[i]);
      if (!tight(width) || !(width.lower > 0) || !tight(time_width) ||
          !(time_width.lower > 0)) {
        temporal_ok = false;
        break;
      }
      for (std::size_t j = 0; j < nw; ++j) {
        const auto l0 = interpolate(a, knots[i], knots[i + 1],
                                    g.luminosity[i * nw + j],
                                    g.luminosity[(i + 1) * nw + j]),
                   l1 = interpolate(b, knots[i], knots[i + 1],
                                    g.luminosity[i * nw + j],
                                    g.luminosity[(i + 1) * nw + j]);
        integrated[j] =
            sum(integrated[j], product(width, divide(sum(l0, l1), {2, 2})));
      }
    }
    if (!temporal_ok) {
      fail_row(row, S::conditioning_budget_exceeded, q.requested_outputs);
      continue;
    }
    for (std::size_t j = 0; j < nw; ++j) {
      const auto average = divide(integrated[j], full);
      if (!tight(average)) {
        admission = S::conditioning_budget_exceeded;
        break;
      }
      Outcome v;
      detail::store(v, average.lower + (average.upper - average.lower) / 2);
      if (!v.value) {
        admission = v.numerical_status;
        break;
      }
      if (average.lower > 0 &&
          std::abs(Wide(*v.value) - average.lower) / average.lower +
                  (average.upper - average.lower) / average.lower >
              1e-13L) {
        admission = S::conditioning_budget_exceeded;
        break;
      }
      mean[j] = *v.value;
    }
    if (admission != S::ok) {
      fail_row(row, admission, q.requested_outputs);
      continue;
    }
    const auto result =
        evaluate_sampled({{g.wavelength, {mean.data(), nw}},
                          {band.wavelength, band.transmission},
                          e.luminosity_distance_metre,
                          e.redshift,
                          e.collecting_area_square_metre,
                          observer_duration},
                         {q.requested_outputs, nw + nb, nw + nb - 3});
    row.admission_status = result.admission_status;
    row.mean_flux_watt_per_square_metre = result.flux_watt_per_square_metre;
    row.energy_joule = result.energy_joule;
    row.expected_photons = result.expected_photons;
    if (row.admission_status != S::ok)
      fail_row(row, row.admission_status, q.requested_outputs);
  }
  return out;
}
} // namespace irred::photometry
