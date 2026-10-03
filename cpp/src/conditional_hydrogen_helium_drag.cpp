#include "irred/conditional_hydrogen_helium_drag.hpp"
#include "hydrogen_helium_cell.hpp"
#include "hydrogen_helium_rates.hpp"
#include "payload_accounting.hpp"
#include "thermal_ruler.hpp"
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <utility>

namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
constexpr std::size_t hard_work = 4000000, hard_bytes = 32 * 1024 * 1024;
constexpr W arithmetic = 4096 * std::numeric_limits<W>::epsilon();
constexpr W projection_floor = 64 * std::numeric_limits<double>::epsilon();
bool valid_policy(const ConditionalHydrogenHeliumDragPolicy &p) {
  return std::fegetround() == FE_TONEAREST && std::numeric_limits<W>::digits >= 64 &&
      std::numeric_limits<W>::max_exponent >= 16384 && p.maximum_total_work <= hard_work &&
      p.maximum_native_bytes <= hard_bytes && p.maximum_intervals <= 4 &&
      p.maximum_root_endpoints <= 8 && p.maximum_endpoint_trials <= 128 &&
      p.maximum_bracket_proposals <= 8 && p.maximum_outer_callbacks_per_integral <= 200000 &&
      p.maximum_depth <= 30 && std::isfinite(p.total_ruler_absolute_tolerance_mpc) &&
      p.total_ruler_absolute_tolerance_mpc > 0 && p.total_ruler_absolute_tolerance_mpc <= 1e-3 &&
      p.total_ruler_relative_tolerance == 0 && std::isfinite(p.capacity_absolute_tolerance) &&
      p.capacity_absolute_tolerance > 0 && p.capacity_absolute_tolerance <= 1e-6 &&
      std::isfinite(p.capacity_relative_tolerance) && p.capacity_relative_tolerance >= 0 &&
      p.capacity_relative_tolerance <= 1e-8;
}
bool add(std::size_t &sum, std::size_t n) noexcept {
  return irred::detail::checked_payload_add(sum, n, 1);
}
W outward(W x) noexcept { return std::nextafter(x, std::numeric_limits<W>::infinity()); }
double diagnostic(W x) noexcept {
  const auto d = static_cast<double>(x);
  return W(d) < x ? std::nextafter(d, INFINITY) : d;
}
// Explicit bounded nonrecursive scratch and simultaneous return headers. Heap
// payload is counted separately exactly once after moves; no NRVO is assumed.
std::size_t scratch_envelope() noexcept {
  return sizeof(ConditionalHydrogenHeliumDrag) + sizeof(ConditionalDragBatch) +
      2 * sizeof(HydrogenHeliumHistory) + sizeof(HydrogenHeliumHistoryRequest) +
      sizeof(ThermalPhysicalMapping) + 4 * sizeof(detail::HHeDragCell) +
      16 * sizeof(detail::HHePositivePolynomial) + 2 * sizeof(detail::HydrogenHeliumCell) +
      4 * sizeof(detail::ThermalRulerIntegral) + sizeof(detail::ThermalRulerWorkBudget);
}
}

std::optional<std::size_t> ConditionalDragWork::checked_total() const noexcept {
  std::size_t sum = 0;
  for (auto n : {source_map, background, momentum, charge, rhs, cell, prefix,
                primitive, root, bound, capacity_outer, ruler_outer})
    if (!add(sum, n)) return {};
  return sum;
}
namespace detail {
struct ConditionalDragAbundanceAccess {
  static W rho(const BaryonAbundance &a) noexcept { return a.rho0_; }
  static W rho_error(const BaryonAbundance &a) noexcept { return a.rho_error_; }
  static std::optional<std::size_t> payload(const BaryonAbundance &a) noexcept { return a.payload(); }
};
struct ConditionalDragAccess {
  struct Evaluation {
    const ConditionalHydrogenHeliumDrag &owner;
    ConditionalDragWork work;
    std::size_t spent;
    S status = S::ok;
    bool charge(std::size_t &category, std::size_t n = 1) noexcept {
      if (spent > owner.policy_.maximum_total_work || n > owner.policy_.maximum_total_work - spent ||
          category > SIZE_MAX - n) { status = S::work_limit; return false; }
      spent += n; category += n; return true;
    }
  };
  static bool charge(ConditionalHydrogenHeliumDrag &o, std::size_t &field, std::size_t n = 1) noexcept {
    if (o.spent_ > o.policy_.maximum_total_work || n > o.policy_.maximum_total_work - o.spent_ || field > SIZE_MAX - n) {
      o.status_ = S::work_limit; return false;
    }
    o.spent_ += n; field += n; return true;
  }
  static ThermalBaryonLoading loading(const ConditionalHydrogenHeliumDrag &o) noexcept {
    return {S::ok, o.loading_, o.loading_error_, o.loading_numerator_, o.loading_denominator_};
  }
  static HHeDragCell cell(const ConditionalHydrogenHeliumDrag &o, std::size_t ascending) noexcept {
    const auto n = HydrogenHeliumCellAccess::count(o.history_);
    const auto c = HydrogenHeliumCellAccess::cell(o.history_, n - 1 - ascending);
    return hydrogen_helium_drag_cell(c, *o.snapshot_->nuclei_today.hydrogen_nuclei.value,
        *o.snapshot_->nuclei_today.helium_nuclei.value, o.loading_, o.loading_error_);
  }
  static std::optional<std::size_t> payload(const ConditionalHydrogenHeliumDrag &o) noexcept {
    irred::detail::PayloadAccounting p(sizeof(o));
    p.embedded(ConditionalDragAbundanceAccess::payload(o.abundance_), sizeof(BaryonAbundance));
    p.embedded(HydrogenHeliumCellAccess::retained_payload(o.history_), sizeof(HydrogenHeliumHistory));
    p.vector(o.prefix_);
    if (o.source_) {
      p.string(o.source_->source_origin); p.string(o.source_->mass_origin);
      p.vector(o.source_->model.species);
    }
    if (o.snapshot_) p.vector(o.snapshot_->fixed_mapped_thermal_source.species);
    return p.result();
  }
  static bool record_peak(ConditionalHydrogenHeliumDrag &o, std::size_t extra) noexcept {
    const auto retained = payload(o);
    if (!retained) { o.status_ = S::work_limit; return false; }
    std::size_t peak = *retained;
    if (!add(peak, extra) || peak > o.policy_.maximum_native_bytes) { o.status_ = S::work_limit; return false; }
    o.peak_payload_ = std::max(o.peak_payload_, peak); return true;
  }
  struct Depth { S status = S::invalid_input; W value = 0, error = 0; };
  static Depth depth(Evaluation &e, W z) {
    const auto &o = e.owner;
    if (e.status != S::ok) return {e.status, 0, 0};
    if (z < o.source_->late_redshift || z > o.source_->initial_redshift) return {S::outside_domain, 0, 0};
    if (z == o.source_->late_redshift) return {S::ok, 0, 0};
    if (z == o.source_->initial_redshift) return {S::ok, o.prefix_.back().value, o.prefix_.back().error};
    const auto n = o.prefix_.size() - 1;
    std::size_t low = 0, high = n;
    while (high - low > 1) {
      const auto mid = low + (high - low) / 2;
      const auto low_redshift = HydrogenHeliumCellAccess::redshift(o.history_, n - mid);
      if (z < low_redshift) high = mid; else low = mid;
    }
    if (!e.charge(e.work.cell) || !e.charge(e.work.primitive)) return {e.status, 0, 0};
    const auto c = cell(o, low);
    const auto q = c.integrate(c.low, z);
    if (q.status != S::ok) return {q.status, 0, 0};
    const W value = o.prefix_[low].value + q.value,
        error = o.prefix_[low].error + q.error + arithmetic * value;
    if (!std::isnormal(value) || !std::isnormal(error)) return {S::conditioning_budget_exceeded, 0, 0};
    return {S::ok, value, error};
  }
  struct BracketBounds { S status = S::invalid_input; W EK = 0, m = 0, L = 0; };
  static BracketBounds bounds(Evaluation &e, W a, W b) {
    BracketBounds out;
    if (!e.charge(e.work.bound)) { out.status = e.status; return out; }
    const auto upper = depth(e, b);
    if (upper.status != S::ok) { out.status = upper.status; return out; }
    out.EK = outward(upper.error + arithmetic * upper.value);
    const auto n = e.owner.prefix_.size() - 1;
    out.m = std::numeric_limits<W>::infinity();
    // Only touched cells are constructed. Index bounds are binary-searched
    // below rather than storing a duplicate cell vector.
    std::size_t lo = 0, hi = n;
    while (lo < hi) {
      const auto mid = lo + (hi - lo) / 2;
      const auto high_redshift = HydrogenHeliumCellAccess::redshift(e.owner.history_, n - 1 - mid);
      if (high_redshift < a) lo = mid + 1; else hi = mid;
    }
    for (std::size_t i = lo; i < n; ++i) {
      const auto low_redshift = HydrogenHeliumCellAccess::redshift(e.owner.history_, n - i);
      if (low_redshift > b) break;
      if (!e.charge(e.work.cell)) { out.status = e.status; return out; }
      const auto c = cell(e.owner, i);
      if (c.status != S::ok) { out.status = c.status; return out; }
      out.m = std::min(out.m, c.lower(std::max(a, c.low), std::min(b, c.high)));
    }
    if (!(out.m > 0) || !std::isnormal(out.m) || !std::isnormal(out.EK)) {
      out.status = S::conditioning_budget_exceeded; return out;
    }
    ThermalRulerWorkBudget budget(e.spent, e.owner.policy_.maximum_total_work, e.work.background,
        e.work.momentum, e.owner.policy_.maximum_outer_callbacks_per_integral,
        e.owner.policy_.history.thermal.maximum_total_callbacks);
    const auto L = thermal_sound_distance_bound(*e.owner.history_.background(),
        loading(e.owner), 1 / (1 + b), 1 / (1 + a),
        e.owner.policy_.history.thermal, budget);
    out.status = L.status; out.L = L.upper_mpc_per_redshift; return out;
  }

  struct CapacityContext {
    const ThermalBackground &background;
    ThermalPolicy thermal;
    ThermalRulerWorkBudget &work;
    W low, bg_relative = 0, cast_relative = 0;
    S status = S::ok;
  };
  static double capacity_integrand(double t, const void *ptr) {
    auto &c = *const_cast<CapacityContext *>(static_cast<const CapacityContext *>(ptr));
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    if (c.status != S::ok) return nan;
    if (!c.work.charge_outer()) { c.status = S::work_limit; return nan; }
    auto p = c.thermal; p.maximum_total_callbacks = std::min(p.maximum_total_callbacks, c.work.momentum_remaining());
    const W a = c.low + (1 - c.low) * t;
    const auto q = c.background.scaled_expansion(a, p);
    if (!c.work.charge_momentum(q.callbacks)) { c.status = S::work_limit; return nan; }
    c.status = q.status;
    if (c.status != S::ok) return nan;
    const W lower = q.a4_e2 - q.error_estimate;
    if (!(lower > 0)) { c.status = S::conditioning_budget_exceeded; return nan; }
    c.bg_relative = std::max(c.bg_relative, q.error_estimate /
        (std::sqrt(lower) * (std::sqrt(q.a4_e2) + std::sqrt(lower))) + arithmetic);
    const W value = 1 / (a * a * a * std::sqrt(q.a4_e2));
    const double rounded = static_cast<double>(value);
    if (!(rounded > 0) || !physical_representable(value)) { c.status = S::outside_domain; return nan; }
    c.cast_relative = std::max(c.cast_relative, std::abs(value - rounded) / value);
    return rounded;
  }
  static ConditionalDragCapacity capacity(Evaluation &e) {
    ConditionalDragCapacity out;
    const auto &o = e.owner;
    const W aL = 1 / (1 + W(o.source_->late_redshift)),
        H0 = W(o.source_->model.h0_km_s_mpc) * 1000 / megaparsec_in_metres_wide(),
        n = *o.snapshot_->nuclei_today.hydrogen_nuclei.value + 2 * W(*o.snapshot_->nuclei_today.helium_nuclei.value),
        scale = W(speed_of_light_m_per_s) * hhe_thomson * n * (1 - aL) / (o.loading_ * H0);
    ThermalRulerWorkBudget budget(e.spent, o.policy_.maximum_total_work, e.work.capacity_outer,
        e.work.momentum, o.policy_.maximum_outer_callbacks_per_integral,
        o.policy_.history.thermal.maximum_total_callbacks);
    if (budget.remaining() < 3) { out.status = S::work_limit; return out; }
    CapacityContext c{*o.history_.background(), o.policy_.history.thermal, budget, aL};
    const double absolute = static_cast<double>(W(o.policy_.capacity_absolute_tolerance) / (8 * scale));
    const auto q = numerics::integrate(capacity_integrand, &c, 0, 1,
        {absolute, o.policy_.capacity_relative_tolerance / 8, budget.remaining(), o.policy_.maximum_depth});
    out.status = c.status == S::ok ? q.status : c.status;
    if (out.status != S::ok) return out;
    if (!(c.cast_relative < 1)) { out.status = S::conditioning_budget_exceeded; return out; }
    const W value = scale * q.value, quad = scale * q.error_estimate,
        dependency = (std::abs(value) + quad) * (c.bg_relative + c.cast_relative) / (1 - c.cast_relative),
        E0 = quad + dependency + projection_floor * value + std::abs(value - W(static_cast<double>(value))),
        loading_relative = o.loading_error_ / (o.loading_ - o.loading_error_),
        loading = (std::abs(value) + E0) * loading_relative,
        error = E0 + loading;
    if (!physical_representable(value) || !std::isnormal(error) ||
        (q.error_estimate > 0 && !std::isnormal(quad)) ||
        ((c.bg_relative > 0 || c.cast_relative > 0) && !std::isnormal(dependency)) ||
        (loading_relative > 0 && !std::isnormal(loading)) ||
        error > o.policy_.capacity_absolute_tolerance + o.policy_.capacity_relative_tolerance * value) {
      out.status = S::conditioning_budget_exceeded; return out;
    }
    out.upper_depth_capacity = static_cast<double>(value); out.numerical_estimate = diagnostic(error); return out;
  }

  static ConditionalDragEndpoint endpoint(Evaluation &e, double D);
};
} // namespace detail

std::optional<std::size_t> ConditionalHydrogenHeliumDrag::payload() const noexcept { return detail::ConditionalDragAccess::payload(*this); }
S ConditionalHydrogenHeliumDrag::status() const noexcept { return status_; }
const ConditionalHydrogenHeliumDragRequest *ConditionalHydrogenHeliumDrag::source() const noexcept { return source_ ? &*source_ : nullptr; }
const ConditionalDragSourceSnapshot *ConditionalHydrogenHeliumDrag::source_snapshot() const noexcept {
  return snapshot_ && snapshot_->fixed_mapped_thermal_source.h0_km_s_mpc > 0 ? &*snapshot_ : nullptr;
}
const HydrogenHeliumHistory *ConditionalHydrogenHeliumDrag::history() const noexcept { return source_ ? &history_ : nullptr; }
ConditionalDragWork ConditionalHydrogenHeliumDrag::preparation_work() const noexcept { return work_; }
const ConditionalDragBudgetDiagnostics *ConditionalHydrogenHeliumDrag::budget_diagnostics() const noexcept { return budget_ ? &*budget_ : nullptr; }

ConditionalHydrogenHeliumDrag::ConditionalHydrogenHeliumDrag(ConditionalHydrogenHeliumDrag &&o) noexcept { *this = std::move(o); }
ConditionalHydrogenHeliumDrag &ConditionalHydrogenHeliumDrag::operator=(ConditionalHydrogenHeliumDrag &&o) noexcept {
  if (this != &o) {
    status_ = o.status_; source_ = std::move(o.source_); abundance_ = std::move(o.abundance_);
    history_ = std::move(o.history_); snapshot_ = std::move(o.snapshot_); budget_ = std::move(o.budget_);
    policy_ = o.policy_; work_ = o.work_; capacity_ = o.capacity_;
    loading_ = o.loading_; loading_error_ = o.loading_error_;
    loading_numerator_ = o.loading_numerator_; loading_denominator_ = o.loading_denominator_; prefix_ = std::move(o.prefix_);
    spent_ = o.spent_; peak_payload_ = o.peak_payload_;
    o.status_ = S::invalid_input; o.source_.reset(); o.snapshot_.reset(); o.budget_.reset();
    o.work_ = {}; o.capacity_ = {}; o.spent_ = o.peak_payload_ = 0;
    o.loading_ = o.loading_error_ = o.loading_numerator_ = o.loading_denominator_ = 0; o.prefix_.clear();
  }
  return *this;
}

ConditionalHydrogenHeliumDrag prepare_conditional_hydrogen_helium_drag(
    const ConditionalHydrogenHeliumDragRequest &r, ConditionalHydrogenHeliumDragPolicy p) {
  using A = detail::ConditionalDragAccess;
  ConditionalHydrogenHeliumDrag out;
  if (!valid_policy(p) || r.source_origin.empty() || r.mass_origin.empty()) return out;
  if (r.source_origin.size() > 65536 || r.mass_origin.size() > 65536) { out.status_ = S::work_limit; return out; }
  // This separately admitted consumer owns the empty explicit-species profile.
  // Existing standalone H/He and thermal domains are unchanged.
  if (!r.model.species.empty()) { out.status_ = S::outside_domain; return out; }
  irred::detail::PayloadAccounting initial(sizeof(out) + scratch_envelope());
  initial.add(r.source_origin.capacity() + 1, 4); initial.add(r.mass_origin.capacity() + 1, 3);
  const auto initial_bound = initial.result();
  if (!initial_bound || *initial_bound > p.maximum_native_bytes) { out.status_ = S::work_limit; return out; }
  out.policy_ = p; out.peak_payload_ = *initial_bound;
  out.budget_ = ConditionalDragBudgetDiagnostics{p.maximum_total_work, p.maximum_native_bytes,
      p.history.maximum_total_work, p.history.maximum_native_bytes, {}, {}, {}, {}, false};
  try {
    out.source_ = r;
    if (!A::charge(out, out.work_.source_map)) return out;
    {
      BaryonAbundanceSource abundance{r.model.physical_baryon_density, r.helium4_neutral_mass_fraction,
          r.hydrogen1_neutral_effective_mass_kg, r.helium4_neutral_effective_mass_kg, r.source_origin, r.mass_origin};
      if (!A::record_peak(out, scratch_envelope() + abundance.source_origin.capacity() + abundance.mass_origin.capacity() + 2)) return out;
      out.abundance_ = prepare_baryon_abundance(abundance, {1, p.maximum_native_bytes - *out.payload()});
      if (out.abundance_.status() != S::ok) { out.status_ = out.abundance_.status(); return out; }
    }
    if (!A::charge(out, out.work_.source_map)) return out;
    {
      const double a = 1;
      if (!A::record_peak(out, scratch_envelope() + sizeof(BaryonAbundanceBatch) + sizeof(BaryonAbundanceRow))) return out;
      const auto row = out.abundance_.evaluate({&a, 1}, {1, p.maximum_native_bytes - *out.payload()});
      if (row.status != S::ok || row.rows.size() != 1) { out.status_ = row.status == S::ok ? S::invalid_input : row.status; return out; }
      out.snapshot_.emplace(ConditionalDragSourceSnapshot{});
      out.snapshot_->nuclei_today = row.rows[0];
    }
    const auto &n = out.snapshot_->nuclei_today;
    if (n.admission_status != S::ok || n.hydrogen_nuclei.status != S::ok || n.helium_nuclei.status != S::ok ||
        !n.hydrogen_nuclei.value || !n.helium_nuclei.value || !(*n.hydrogen_nuclei.value > 0) || !(*n.helium_nuclei.value > 0)) {
      out.status_ = S::outside_domain; return out;
    }
    out.snapshot_->emitted_nuclei_binary64_bits = {std::bit_cast<std::uint64_t>(*n.hydrogen_nuclei.value),
                                                 std::bit_cast<std::uint64_t>(*n.helium_nuclei.value)};
    const W rho = detail::ConditionalDragAbundanceAccess::rho(out.abundance_),
        reconstructed = W(r.hydrogen1_neutral_effective_mass_kg) * *n.hydrogen_nuclei.value +
                        W(r.helium4_neutral_effective_mass_kg) * *n.helium_nuclei.value;
    out.snapshot_->relative_mass_reconstruction_residual = static_cast<double>((reconstructed - rho) / rho);
    out.snapshot_->reconstruction_arithmetic_estimate = diagnostic(arithmetic +
        detail::ConditionalDragAbundanceAccess::rho_error(out.abundance_) +
        n.hydrogen_nuclei.relative_arithmetic_estimate + n.helium_nuclei.relative_arithmetic_estimate);
    HydrogenHeliumHistoryRequest hs{r.model, *n.hydrogen_nuclei.value, *n.helium_nuclei.value,
        r.source_origin, r.initial_redshift, r.late_redshift};
    auto child = p.history;
    // Requested values remain immutable; tighten by actual parent live/work.
    std::size_t parent_live = *out.payload();
    if (!add(parent_live, scratch_envelope()) || !add(parent_live, hs.nuclei_origin.capacity() + 1) ||
        parent_live > p.maximum_native_bytes || out.spent_ > p.maximum_total_work) {
      out.status_ = S::work_limit; return out;
    }
    child.maximum_total_work = std::min(child.maximum_total_work, p.maximum_total_work - out.spent_);
    child.maximum_native_bytes = std::min(child.maximum_native_bytes, p.maximum_native_bytes - parent_live);
    out.budget_->earned_work_before_history = out.spent_;
    out.budget_->parent_live_bytes_before_history = parent_live;
    out.budget_->served_history_work = child.maximum_total_work;
    out.budget_->served_history_native_bytes = child.maximum_native_bytes;
    const auto history_bound = child.base_intervals <= 16384 ? hydrogen_helium_history_payload_bound(
        4 * child.base_intervals, 0, 0, hs.nuclei_origin.capacity()) : std::nullopt;
    if (!history_bound || *history_bound > child.maximum_native_bytes) { out.status_ = S::work_limit; return out; }
    out.peak_payload_ = std::max(out.peak_payload_, parent_live + *history_bound);
    out.budget_->history_preparation_attempted = true;
    out.history_ = prepare_hydrogen_helium_history(hs, child);
    const auto hwork = out.history_.work();
    // Import each bounded field once, even when the child refused. The public
    // HHe output total() is not trusted as a freshness or transfer predicate.
    if (!A::charge(out, out.work_.background, hwork.background_evaluations) ||
        !A::charge(out, out.work_.momentum, hwork.momentum_callbacks) ||
        !A::charge(out, out.work_.charge, hwork.initial_charge_evaluations) ||
        !A::charge(out, out.work_.rhs, hwork.rhs_evaluations)) return out;
    if (const auto *bg = out.history_.background(); bg && out.history_.thermal_mapping_witnesses()) {
      out.snapshot_->fixed_mapped_thermal_source = bg->source();
      out.snapshot_->thermal_source_map = *out.history_.thermal_mapping_witnesses();
      const auto &photon = out.snapshot_->thermal_source_map[0];
      out.snapshot_->relative_photon_reconstruction_residual = static_cast<double>(
          (W(photon.emitted_value) - photon.wide_value) / photon.wide_value);
    }
    if (out.history_.status() != S::ok) { out.status_ = out.history_.status(); return out; }
    const auto loading = detail::thermal_baryon_loading(*out.history_.background());
    if (loading.status != S::ok || !(loading.ratio_today > loading.arithmetic_estimate)) {
      out.status_ = loading.status == S::ok ? S::conditioning_budget_exceeded : loading.status; return out;
    }
    out.loading_ = loading.ratio_today; out.loading_error_ = loading.arithmetic_estimate;
    out.loading_numerator_ = loading.numerator; out.loading_denominator_ = loading.denominator;
    out.snapshot_->baryon_loading = static_cast<double>(out.loading_);
    out.snapshot_->loading_arithmetic_estimate = diagnostic(out.loading_error_);
    const auto count = detail::HydrogenHeliumCellAccess::count(out.history_);
    std::size_t prefix_bytes = 0;
    if (count == SIZE_MAX || !irred::detail::checked_payload_add(prefix_bytes, count + 1, sizeof(ConditionalHydrogenHeliumDrag::Prefix)) ||
        !A::record_peak(out, scratch_envelope() + prefix_bytes)) { out.status_ = S::work_limit; return out; }
    out.prefix_.resize(count + 1);
    W sum = 0, compensation = 0, error = 0, error_compensation = 0;
    for (std::size_t i = 0; i < count; ++i) {
      if (!A::charge(out, out.work_.cell) || !A::charge(out, out.work_.primitive) || !A::charge(out, out.work_.prefix)) return out;
      const auto c = A::cell(out, i); const auto q = c.integrate(c.low, c.high);
      if (q.status != S::ok) { out.status_ = q.status; return out; }
      const W y = q.value - compensation, next = sum + y;
      compensation = (next - sum) - y; sum = next;
      const W ey = q.error - error_compensation, enext = error + ey;
      error_compensation = (enext - error) - ey; error = enext;
      out.prefix_[i + 1] = {sum, outward(error + arithmetic * sum)};
      if (!std::isnormal(sum) || !std::isnormal(out.prefix_[i + 1].error) || !(sum > out.prefix_[i].value)) {
        out.status_ = S::conditioning_budget_exceeded; return out;
      }
    }
    A::Evaluation capacity_work{out, out.work_, out.spent_};
    out.capacity_ = A::capacity(capacity_work);
    out.work_ = capacity_work.work; out.spent_ = capacity_work.spent;
    out.status_ = out.capacity_.status;
    if (!A::record_peak(out, scratch_envelope())) return out;
  } catch (const std::bad_alloc &) { out.status_ = S::work_limit; }
  return out;
}

ConditionalDragEndpoint detail::ConditionalDragAccess::endpoint(Evaluation &e, double D) {
  ConditionalDragEndpoint out;
  const auto &o = e.owner;
  const W target = 1 - W(D), ET = (D == 1 ? 0 : arithmetic * (1 + std::abs(target)));
  const W J = o.prefix_.back().value, EJ = o.prefix_.back().error;
  auto fail = [&](S s) { out.root_status = out.ruler_status = s; return out; };
  if (target > J + EJ) return fail(S::outside_domain);
  if (target != 0 && target + ET >= J - EJ) return fail(S::conditioning_budget_exceeded);
  W z = o.source_->late_redshift, inherited = 0, locator = 0, target_ruler = 0, Ez = 0;
  unsigned trials = 0;
  auto trial = [&]() {
    if (trials >= o.policy_.maximum_endpoint_trials) { e.status = S::work_limit; return false; }
    ++trials; return e.charge(e.work.root);
  };
  if (target != 0) {
    W low = o.source_->late_redshift, high = o.source_->initial_redshift;
    for (unsigned i = 0; i < 64; ++i) {
      if (!trial()) return fail(e.status);
      const W mid = (low + high) / 2;
      const auto q = depth(e, mid);
      if (q.status != S::ok) return fail(q.status);
      if (q.value < target) low = mid; else high = mid;
      if (high - low <= 64 * std::numeric_limits<double>::epsilon() * std::max(1.L, std::abs(mid))) break;
    }
    z = (low + high) / 2;
    const W rounded_z = static_cast<double>(z);
    const W radius = std::max(std::abs(rounded_z - low), std::abs(high - rounded_z));
    // Initial local bounds propose padding only. Every actual proposal below
    // recomputes EK,m,L on the entire expanded bracket, including all shifts.
    const auto initial = bounds(e, low, high);
    if (initial.status != S::ok) return fail(initial.status);
    W padding = outward(2 * (initial.EK + ET) / initial.m);
    bool closed = false;
    for (unsigned proposal = 0; proposal < o.policy_.maximum_bracket_proposals; ++proposal) {
      if (!trial()) return fail(e.status);
      const W a = std::nextafter(low - padding, -std::numeric_limits<W>::infinity()),
          b = outward(high + padding);
      if (!std::isfinite(padding) || a < o.source_->late_redshift || b > o.source_->initial_redshift)
        return fail(S::conditioning_budget_exceeded);
      const auto bound = bounds(e, a, b);
      if (bound.status != S::ok) return fail(bound.status);
      const W rho = outward((bound.EK + ET) / bound.m);
      if (!trial()) return fail(e.status);
      const auto ka = depth(e, a);
      if (!trial()) return fail(e.status);
      const auto kb = depth(e, b);
      if (ka.status != S::ok || kb.status != S::ok) return fail(ka.status == S::ok ? kb.status : ka.status);
      if (low - rho >= a && high + rho <= b && rounded_z >= a && rounded_z <= b &&
          ka.value + bound.EK <= target - ET && kb.value - bound.EK >= target + ET) {
        inherited = bound.L * bound.EK / bound.m;
        locator = bound.L * radius;
        target_ruler = bound.L * ET / bound.m;
        Ez = radius + rho; z = rounded_z; closed = true; break;
      }
      padding = outward(2 * std::max(padding, rho));
    }
    if (!closed) return fail(S::conditioning_budget_exceeded);
  }
  out.root_status = S::ok; out.redshift = static_cast<double>(z);
  out.redshift_numerical_estimate = diagnostic(Ez);
  ThermalRulerWorkBudget budget(e.spent, o.policy_.maximum_total_work, e.work.ruler_outer,
      e.work.momentum, o.policy_.maximum_outer_callbacks_per_integral,
      o.policy_.history.thermal.maximum_total_callbacks);
  const W total = o.policy_.total_ruler_absolute_tolerance_mpc;
  const ThermalRulerAllowance allowance{o.policy_.history.thermal, static_cast<double>(.3L * total), 0,
      o.policy_.maximum_outer_callbacks_per_integral, o.policy_.maximum_depth};
  const auto ruler = integrate_thermal_ruler(*o.history_.background(), loading(o),
      1 / (1 + z), allowance, budget);
  out.ruler_status = ruler.status;
  if (ruler.status != S::ok) return out;
  const double value = static_cast<double>(ruler.value_mpc);
  const W projection = std::abs(ruler.value_mpc - W(value)) + projection_floor * ruler.value_mpc,
      extra = target_ruler + projection,
      combined = inherited + locator + ruler.error_estimate_mpc + extra;
  out.total_ruler_numerical_estimate_mpc = diagnostic(combined);
  if (!(value > 0) || !std::isfinite(value) || !std::isnormal(combined) || inherited > .5L * total ||
      locator > .1L * total || ruler.error_estimate_mpc > .3L * total || extra > .1L * total || combined > total) {
    out.ruler_status = S::conditioning_budget_exceeded; return out;
  }
  out.comoving_ruler_mpc = value; return out;
}

ConditionalDragBatch ConditionalHydrogenHeliumDrag::evaluate(std::span<const ConditionalDragInterval> intervals) const {
  using A = detail::ConditionalDragAccess;
  ConditionalDragBatch out;
  out.work = work_; out.late_charge_capacity = capacity_;
  out.retained_payload_bytes = payload().value_or(SIZE_MAX); out.peak_payload_bytes = peak_payload_;
  if (!valid_policy(policy_) || intervals.size() > policy_.maximum_intervals ||
      intervals.size() > policy_.maximum_root_endpoints / 2) { out.status = S::work_limit; return out; }
  if (status_ != S::ok) { out.status = status_; return out; }
  irred::detail::PayloadAccounting bytes(out.retained_payload_bytes + sizeof(out) + scratch_envelope());
  // Two return/copy headers, one heap row vector. All original interval strings
  // are copied once into their ordered rows, including refused rows.
  bytes.add(intervals.size(), sizeof(ConditionalDragIntervalRow));
  for (const auto &i : intervals) { bytes.string(i.id); bytes.string(i.origin); }
  const auto peak = bytes.result();
  if (!peak || *peak > policy_.maximum_native_bytes) { out.status = S::work_limit; return out; }
  out.peak_payload_bytes = std::max(out.peak_payload_bytes, *peak);
  A::Evaluation e{*this, work_, spent_};
  try {
    out.rows.reserve(intervals.size());
    for (const auto &i : intervals) {
      out.rows.push_back({i, {}}); auto &row = out.rows.back();
      S input = S::ok;
      if (!std::isfinite(i.lower_depth) || !std::isfinite(i.upper_depth)) input = S::nonfinite_input;
      else if (i.id.empty() || i.origin.empty() || i.lower_depth < 0 || i.upper_depth > 1 || i.upper_depth < i.lower_depth)
        input = S::outside_domain;
      else if (!capacity_.upper_depth_capacity || i.upper_depth > W(*capacity_.upper_depth_capacity) + capacity_.numerical_estimate)
        input = S::outside_domain;
      else if (i.upper_depth > W(*capacity_.upper_depth_capacity) - capacity_.numerical_estimate)
        input = S::conditioning_budget_exceeded;
      if (input != S::ok) {
        for (auto &endpoint : row.endpoints) endpoint.root_status = endpoint.ruler_status = input;
        continue;
      }
      row.endpoints[0] = A::endpoint(e, i.upper_depth);
      row.endpoints[1] = A::endpoint(e, i.lower_depth);
      if (row.endpoints[0].root_status != S::ok || row.endpoints[1].root_status != S::ok) {
        const auto cause = row.endpoints[0].root_status == S::ok ? row.endpoints[1].root_status : row.endpoints[0].root_status;
        for (auto &endpoint : row.endpoints) {
          endpoint.root_status = endpoint.ruler_status = cause;
          endpoint.redshift.reset(); endpoint.comoving_ruler_mpc.reset();
        }
      }
    }
    out.status = S::ok;
  } catch (const std::bad_alloc &) { out.status = S::work_limit; }
  out.work = e.work; return out;
}
} // namespace irred::cosmology
