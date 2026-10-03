#pragma once
// Local proof caller. Only installed public irred headers are used; this is
// neither a new SDK owner nor a reusable arbitrary-weight law framework.
#include "abundance_history_cohort_source.hpp"
#include <irred/baryon_abundance.hpp>
#include <irred/hydrogen_helium_history.hpp>
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
namespace abundance_history_cohort {
namespace c = irred::cosmology;
using S = irred::numerics::Status;
using W = long double;
inline bool profile() noexcept {
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0) || !defined(__GLIBCXX__)
  return false;
#else
  return std::fegetround() == FE_TONEAREST &&
      std::numeric_limits<W>::digits >= 64 &&
      std::numeric_limits<W>::max_exponent >= 16384;
#endif
}
inline bool add(std::size_t &total, std::size_t n) noexcept {
  if (n > std::numeric_limits<std::size_t>::max() - total) return false;
  total += n; return true;
}
struct Bytes {
  std::size_t value = 0;
  bool valid = true;
  void count(std::size_t n, std::size_t width = 1) noexcept {
    if (!valid || (n && width > (SIZE_MAX - value) / n)) { valid = false; return; }
    value += n * width;
  }
  void copy_string(const std::string &s, std::size_t copies = 1) noexcept {
    if (s.size() == SIZE_MAX) { valid = false; return; }
    count(copies, std::max(std::size_t(32), s.size() + 1));
  }
  void actual_string(const std::string &s) noexcept {
    if (s.capacity() == SIZE_MAX) { valid = false; return; }
    count(s.capacity() + 1);
  }
  std::optional<std::size_t> result() const noexcept {
    return valid ? std::optional(value) : std::nullopt;
  }
};
struct State {
  std::string id;
  double mass = relative_mass;
  c::ThermalPhysicalModel model{};
  double helium_fraction = 0, hydrogen_mass_kg = 0, helium_mass_kg = 0;
  std::string source_origin, mass_origin, photon_temperature_origin;
};
struct Query { std::string id; double redshift = 0; };
struct Input {
  std::array<State, states> support;
  std::array<Query, rows> queries;
  std::string distribution_origin, dependence_origin;
  std::string role = std::string(abundance_history_cohort::role);
  std::string temperature_role = std::string(abundance_history_cohort::temperature_role);
  double initial_redshift = abundance_history_cohort::initial_redshift;
  double late_redshift = abundance_history_cohort::late_redshift;
  unsigned requested_outputs = outputs;
};
inline Input original_input() {
  Input in;
  in.distribution_origin = distribution_origin; in.dependence_origin = dependence_origin;
  const std::array<double, states> h0{h0_a, h0_b}, T{tcmb_a, tcmb_b}, Y{y_a, y_b};
  const std::array<std::string_view, states> origin{source_origin_a, source_origin_b},
      photon{photon_origin_a, photon_origin_b};
  for (std::size_t s = 0; s < states; ++s) {
    auto &x = in.support[s]; x.id = state_ids[s];
    x.model = {h0[s], omega_b, omega_cdm, T[s], omega_other, {}};
    x.helium_fraction = Y[s]; x.hydrogen_mass_kg = neutral_hydrogen_mass_kg;
    x.helium_mass_kg = neutral_helium_mass_kg;
    x.source_origin = origin[s]; x.mass_origin = mass_origin;
    x.photon_temperature_origin = photon[s];
  }
  const std::array<double, rows> z{initial_redshift, interior_redshift, late_redshift};
  for (std::size_t r = 0; r < rows; ++r) in.queries[r] = {std::string(row_ids[r]), z[r]};
  return in;
}
struct Policy {
  std::size_t maximum_preparation_work = preparation_work_cap,
      maximum_evaluation_work = evaluation_work_cap,
      maximum_two_evaluation_work = two_evaluation_work_cap,
      maximum_whole_live_bytes = whole_live_payload_cap;
  double covariance_error_fraction = covariance_resolution;
  c::HydrogenHeliumHistoryPolicy history;
};
inline bool valid_policy(const Policy &p) noexcept {
  const auto &h = p.history;
  return p.maximum_preparation_work <= preparation_work_cap &&
      p.maximum_evaluation_work <= evaluation_work_cap &&
      p.maximum_two_evaluation_work <= two_evaluation_work_cap &&
      p.maximum_whole_live_bytes <= whole_live_payload_cap &&
      std::isfinite(p.covariance_error_fraction) && p.covariance_error_fraction > 0 &&
      p.covariance_error_fraction <= covariance_resolution &&
      (h.base_intervals == default_base_intervals || h.base_intervals == 16384) &&
      h.maximum_fine_intervals >= 4 * h.base_intervals &&
      h.maximum_fine_intervals <= maximum_fine_intervals &&
      h.maximum_total_work <= child_work_cap && h.maximum_native_bytes <= child_payload_cap &&
      h.absolute_fraction_tolerance >= 0 && h.absolute_fraction_tolerance <= 1e-8 &&
      h.relative_fraction_tolerance >= 0 && h.relative_fraction_tolerance <= 2e-6 &&
      h.absolute_temperature_tolerance_kelvin >= 0 && h.absolute_temperature_tolerance_kelvin <= 2e-4 &&
      h.relative_temperature_tolerance >= 0 && h.relative_temperature_tolerance <= 2e-6 &&
      h.absolute_opacity_tolerance >= 0 && h.absolute_opacity_tolerance <= 2e-9 &&
      h.relative_opacity_tolerance >= 0 && h.relative_opacity_tolerance <= 3e-6 &&
      (h.absolute_fraction_tolerance > 0 || h.relative_fraction_tolerance > 0) &&
      (h.absolute_temperature_tolerance_kelvin > 0 || h.relative_temperature_tolerance > 0) &&
      (h.absolute_opacity_tolerance > 0 || h.relative_opacity_tolerance > 0) &&
      h.thermal.momentum_method == c::ThermalMomentumMethod::direct_adaptive &&
      h.thermal.absolute_tolerance >= 0 && h.thermal.absolute_tolerance <= 1e-12 &&
      h.thermal.relative_tolerance >= 0 && h.thermal.relative_tolerance <= 2e-12 &&
      (h.thermal.absolute_tolerance > 0 || h.thermal.relative_tolerance > 0) &&
      h.thermal.maximum_callbacks_per_evaluation <= 200000 &&
      h.thermal.maximum_total_callbacks <= child_work_cap && h.thermal.maximum_depth <= 30 &&
      h.thermal.maximum_points <= 4096 && h.thermal.maximum_species <= 16 &&
      h.thermal.maximum_native_bytes <= 16 * 1024 * 1024;
}
struct Source {
  Input input;
  Policy requested_policy;
  std::array<double, states> probabilities{probability, probability};
  Source(const Input &in, const Policy &p) : input(in), requested_policy(p) {}
};
struct NativeSnapshot {
  std::optional<c::BaryonAbundanceRow> nuclei;
  std::array<std::uint64_t, 2> emitted_nuclei_bits{};
  std::optional<c::ThermalFlatModel> emitted_thermal;
  std::optional<std::array<c::ThermalScalarMapWitness, 4>> thermal_witnesses;
  std::optional<double> maximum_log_excluded_heiii_activity;
  bool native_history_source_present = false;
};
struct PreparationWork {
  std::size_t complete_intakes = 0, abundance_prepares = 0, abundance_maps = 0,
      background = 0, momentum = 0, charge = 0, rhs = 0;
  bool overflowed = false;
  std::optional<std::size_t> checked_total() const noexcept {
    if (overflowed) return {};
    std::size_t n = 0;
    for (auto v : {complete_intakes, abundance_prepares, abundance_maps,
                  background, momentum, charge, rhs}) if (!add(n, v)) return {};
    return n;
  }
};
struct EvaluationWork {
  std::size_t dispatches = 0, queried_rows = 0, collections = 0, witnesses = 0,
      means = 0, differences = 0, products = 0, projections = 0,
      serialization_fields = 0, serialization_wide_nibbles = 0,
      serialization_wide_decompositions = 0;
  std::optional<std::size_t> checked_total() const noexcept {
    std::size_t n = 0;
    for (auto v : {dispatches, queried_rows, collections, witnesses, means,
                  differences, products, projections, serialization_fields,
                  serialization_wide_nibbles, serialization_wide_decompositions}) if (!add(n, v)) return {};
    return n;
  }
};
struct Budget {
  std::size_t requested_work = 0, requested_bytes = 0;
  std::optional<std::size_t> served_work, served_bytes;
  bool history_attempted = false;
};
struct Slot {
  S status = S::invalid_input;
  c::BaryonAbundance abundance;
  c::BaryonAbundanceBatch abundance_attempt;
  c::HydrogenHeliumHistory history;
  NativeSnapshot snapshot;
  Budget budget;
};
struct Moments {
  std::array<double, axes> mean{}, mean_error{}, mean_cast_and_round{};
  std::array<double, covariance_cells> covariance{}, covariance_error{}, covariance_cast_and_round{};
  std::array<double, axes> difference_round{};
};
struct ReductionScratch {
  std::array<W, axes> mean{}, difference{}, mean_error{}, difference_error{},
      mean_round{}, difference_round{};
};
struct Evaluation {
  S status = S::invalid_input;
  std::shared_ptr<const Source> source;
  std::array<NativeSnapshot, states> snapshots;
  std::array<S, states> state_status{S::invalid_input, S::invalid_input};
  std::array<c::HydrogenHeliumHistoryBatch, states> attempts;
  EvaluationWork work;
  std::optional<Moments> moments;
  std::size_t conservative_whole_live_bytes = 0, known_result_payload_bytes = 0;
  Evaluation() = default;
  Evaluation(const Evaluation &) = delete;
  Evaluation &operator=(const Evaluation &) = delete;
  Evaluation(Evaluation &&) noexcept = default;
  Evaluation &operator=(Evaluation &&) noexcept = default;
};
inline W up(W x) noexcept {
  return x == 0 ? 0 : std::nextafter(x, std::numeric_limits<W>::infinity());
}
inline W positive_add(W x, W y) noexcept {
  if (x == 0 && y == 0) return 0;
  const W value = x + y;
  return value > 0 && std::isnormal(value) ? up(value) : std::numeric_limits<W>::infinity();
}
inline W positive_product(W x, W y) noexcept {
  if (x == 0 || y == 0) return 0;
  const W value = x * y;
  return value > 0 && std::isnormal(value) ? up(value) : std::numeric_limits<W>::infinity();
}
inline W round_allowance(W result) noexcept {
  if (result == 0) return 0;
  const W unit = std::numeric_limits<W>::epsilon() / 2;
  const W denominator = std::nextafter(1 - unit, W(0));
  return up(positive_product(std::abs(result), unit) / denominator);
}
inline bool wide_coordinate(W x, bool required_nonzero = false) noexcept {
  return std::isfinite(x) && (x == 0 ? !required_nonzero : std::isnormal(x));
}
inline std::optional<double> outward(W x) noexcept {
  if (!std::isfinite(x) || x < 0) return {};
  if (x == 0) return 0;
  double y = static_cast<double>(x);
  if (!std::isfinite(y)) return {};
  if (y == 0) y = std::nextafter(0.0, INFINITY);
  if (!(y > 0)) return {};
  if (W(y) < x) y = std::nextafter(y, INFINITY);
  if (!std::isfinite(y)) return {};
  return y;
}
inline S project(W value, W error, W budget, bool positive,
                 double &out, double &diagnostic, double &cast_and_round,
                 W operation_round) noexcept {
  if (!wide_coordinate(value, positive) || !std::isfinite(error) || error < 0 ||
      !std::isfinite(budget) || budget < 0) return S::conditioning_budget_exceeded;
  const double y = static_cast<double>(value);
  if (!std::isfinite(y) || (value != 0 && y == 0) || (positive && !(y > 0)))
    return S::conditioning_budget_exceeded;
  const W cast_loss = std::abs(value - W(y));
  const auto e = outward(positive_add(error, cast_loss));
  const auto a = outward(positive_add(operation_round, cast_loss));
  if (!e || !a || W(*e) > budget) return S::conditioning_budget_exceeded;
  out = y; diagnostic = *e; cast_and_round = *a; return S::ok;
}
// Independent toy controls feed this fixed two-half reduction; no physics is
// implemented here. Errors are absolute deterministic parent diagnostics.
inline S reduce(const std::array<std::array<W, axes>, states> &f,
                const std::array<std::array<W, axes>, states> &e,
                const std::array<bool, axes> &constant,
                double covariance_fraction, EvaluationWork &work,
                std::size_t work_cap, std::optional<Moments> &output) noexcept {
  output.reset();
  if (!profile() || !std::isfinite(covariance_fraction) || covariance_fraction <= 0 ||
      covariance_fraction > covariance_resolution) return S::invalid_input;
  auto charge = [&](std::size_t &field) {
    const auto used = work.checked_total();
    if (!used || *used >= work_cap || field == SIZE_MAX) return false;
    ++field; return true;
  };
  ReductionScratch w;
  Moments m;
  for (std::size_t i = 0; i < axes; ++i) {
    if (!charge(work.means)) return S::work_limit;
    if (!wide_coordinate(f[0][i], true) || !wide_coordinate(f[1][i], true) ||
        !std::isfinite(e[0][i]) || !std::isfinite(e[1][i]) || e[0][i] < 0 || e[1][i] < 0)
      return S::conditioning_budget_exceeded;
    const W sum = f[0][i] + f[1][i];
    if (!wide_coordinate(sum, true)) return S::conditioning_budget_exceeded;
    w.mean[i] = sum / 2; w.mean_round[i] = round_allowance(sum) / 2;
    w.mean_error[i] = positive_add(positive_add(e[0][i], e[1][i]) / 2, w.mean_round[i]);
    const W absolute = i % 2 == 0 ? W(2e-4) : W(2e-9),
        relative = i % 2 == 0 ? W(2e-6) : W(3e-6),
        relative_part = std::nextafter(relative * std::abs(w.mean[i]), W(0)),
        budget = std::nextafter(absolute + relative_part, W(0));
    if (!charge(work.projections)) return S::work_limit;
    if (auto s = project(w.mean[i], w.mean_error[i], budget, true,
                         m.mean[i], m.mean_error[i], m.mean_cast_and_round[i],
                         w.mean_round[i]); s != S::ok) return s;
    if (!charge(work.differences)) return S::work_limit;
    w.difference[i] = constant[i] ? 0 : f[0][i] - f[1][i];
    w.difference_round[i] = constant[i] ? 0 : round_allowance(w.difference[i]);
    w.difference_error[i] = constant[i] ? 0 : positive_add(
        positive_add(e[0][i], e[1][i]), w.difference_round[i]);
    if (!constant[i] && !wide_coordinate(w.difference[i], true))
      return S::conditioning_budget_exceeded;
    const auto ar = outward(w.difference_round[i]);
    if (!ar) return S::conditioning_budget_exceeded;
    m.difference_round[i] = *ar;
  }
  for (std::size_t i = 0; i < axes; ++i) for (std::size_t j = 0; j <= i; ++j) {
    const std::size_t ij = i * axes + j, ji = j * axes + i;
    if (!charge(work.products) || !charge(work.projections)) return S::work_limit;
    if (constant[i] || constant[j]) continue;
    const W product = w.difference[i] * w.difference[j];
    if (!wide_coordinate(product, true)) return S::conditioning_budget_exceeded;
    const W value = product / 4, arithmetic = round_allowance(product) / 4;
    // Exact two-state geometric variance identity; inward allowance prevents
    // product rounding from manufacturing allocation slack. No sqrt/libm.
    const W scale = std::nextafter(std::abs(value) - arithmetic, W(0));
    if (!(scale > 0) || !wide_coordinate(scale)) return S::conditioning_budget_exceeded;
    const W inherited = positive_add(positive_add(
        positive_product(std::abs(w.difference[i]), w.difference_error[j]),
        positive_product(std::abs(w.difference[j]), w.difference_error[i])),
        positive_product(w.difference_error[i], w.difference_error[j])) / 4;
    const W error = positive_add(inherited, arithmetic);
    const W budget = std::nextafter(W(covariance_fraction) * scale, W(0));
    if (auto s = project(value, error, budget, i == j,
                         m.covariance[ij], m.covariance_error[ij],
                         m.covariance_cast_and_round[ij], arithmetic); s != S::ok) return s;
    m.covariance[ji] = m.covariance[ij]; m.covariance_error[ji] = m.covariance_error[ij];
    m.covariance_cast_and_round[ji] = m.covariance_cast_and_round[ij];
  }
  output = m; return S::ok;
}
inline bool same_state(const State &a, const State &b) noexcept {
  const auto &x = a.model, &y = b.model;
  if (!x.species.empty() || !y.species.empty()) return false;
  return bits(x.h0_km_s_mpc) == bits(y.h0_km_s_mpc) &&
      bits(x.physical_baryon_density) == bits(y.physical_baryon_density) &&
      bits(x.physical_cdm_density) == bits(y.physical_cdm_density) &&
      bits(x.tcmb_kelvin) == bits(y.tcmb_kelvin) &&
      bits(x.physical_massless_nonphoton_density) == bits(y.physical_massless_nonphoton_density) &&
      bits(a.helium_fraction) == bits(b.helium_fraction) &&
      bits(a.hydrogen_mass_kg) == bits(b.hydrogen_mass_kg) && bits(a.helium_mass_kg) == bits(b.helium_mass_kg);
}
inline bool constant_axis(const Input &in, std::size_t axis) noexcept {
  if (same_state(in.support[0], in.support[1])) return true;
  return axis % 2 == 0 && bits(in.queries[axis / 2].redshift) == bits(in.initial_redshift) &&
      bits(in.support[0].model.tcmb_kelvin) == bits(in.support[1].model.tcmb_kelvin);
}
class Cohort {
public:
  Cohort() = default;
  Cohort(const Cohort &) = delete;
  Cohort &operator=(const Cohort &) = delete;
  Cohort(Cohort &&other) noexcept { *this = std::move(other); }
  Cohort &operator=(Cohort &&other) noexcept {
    if (this != &other) {
      status_ = other.status_; source_ = std::move(other.source_);
      slots_ = std::move(other.slots_); work_ = other.work_; policy_ = other.policy_;
      preparation_bound_ = other.preparation_bound_; known_payload_ = other.known_payload_;
      spent_ = other.spent_; evaluation_count_ = other.evaluation_count_;
      serialization_recorded_ = other.serialization_recorded_;
      other.status_ = S::invalid_input; other.work_ = {}; other.preparation_bound_ = 0;
      other.known_payload_ = 0; other.spent_ = 0; other.evaluation_count_ = 0;
      other.serialization_recorded_ = false;
      for (auto &slot : other.slots_) { slot.status = S::invalid_input; slot.snapshot = {}; slot.budget = {}; }
    }
    return *this;
  }
  S status() const noexcept { return status_; }
  const std::shared_ptr<const Source> &source() const noexcept { return source_; }
  const Slot &slot(std::size_t s) const { return slots_.at(s); }
  PreparationWork preparation_work() const noexcept { return work_; }
  std::size_t conservative_whole_live_bytes() const noexcept { return preparation_bound_; }
  std::size_t known_retained_payload_bytes() const noexcept { return known_payload_; }
  std::size_t evaluation_count() const noexcept { return evaluation_count_; }
  Evaluation evaluate();
  std::size_t serialization_budget(const Evaluation &) const noexcept;
  S record_serialization(Evaluation &, std::size_t fields, std::size_t nibbles,
                         std::size_t decompositions) noexcept;
  static std::optional<std::size_t> whole_live_bound(const Input &, const Policy &) noexcept;
private:
  S status_ = S::invalid_input;
  std::shared_ptr<const Source> source_;
  std::array<Slot, states> slots_;
  PreparationWork work_;
  Policy policy_;
  std::size_t preparation_bound_ = 0, known_payload_ = 0, spent_ = 0, evaluation_count_ = 0;
  bool serialization_recorded_ = false;
  bool charge(std::size_t &category, std::size_t n = 1) noexcept {
    const auto total = work_.checked_total();
    if (!total || *total > policy_.maximum_preparation_work ||
        n > policy_.maximum_preparation_work - *total || !add(category, n)) return false;
    return true;
  }
  std::optional<std::size_t> known_payload() const noexcept;
  friend Cohort prepare(const Input &, Policy);
};
static_assert(!std::is_copy_constructible_v<Cohort> && !std::is_copy_assignable_v<Cohort>);
static_assert(std::is_nothrow_move_constructible_v<Cohort> && std::is_nothrow_move_assignable_v<Cohort>);
inline std::optional<std::size_t> Cohort::whole_live_bound(const Input &in, const Policy &p) noexcept {
  Bytes b;
  b.count(2, sizeof(Cohort)); b.count(1, sizeof(Source));
  b.count(3, sizeof(Evaluation)); // prior receipt, local return, destination
  b.count(1, sizeof(ReductionScratch));
  b.count(2, sizeof(std::array<std::array<W, axes>, states>)); // collected values/errors
  b.count(1, sizeof(std::array<double, rows>));
  b.count(1, sizeof(std::array<bool, axes>));
  b.count(2, sizeof(Moments)); // unelided reducer result/optional transfer
  b.count(maximum_wire_bytes + 128); // fixed output buffer and writer counters
  b.count(maximum_wire_scratch_bytes); // fixed prefix/trailer/index arrays
  b.count(1, sizeof(c::BaryonAbundanceSource));
  b.count(2, sizeof(c::BaryonAbundance)); // unelided native prepare/return headers
  b.count(2, sizeof(c::BaryonAbundanceBatch));
  b.count(2, sizeof(c::HydrogenHeliumHistoryRequest));
  b.count(2, sizeof(c::HydrogenHeliumHistory)); // child return/assignment headers
  for (const auto &s : {&in.distribution_origin, &in.dependence_origin, &in.role, &in.temperature_role}) b.copy_string(*s);
  for (const auto &r : in.queries) b.copy_string(r.id);
  for (const auto &s : in.support) {
    b.copy_string(s.id); b.copy_string(s.source_origin, 5);
    b.copy_string(s.mass_origin, 3); b.copy_string(s.photon_temperature_origin);
    b.count(s.model.species.size(), sizeof(c::ThermalPhysicalSpecies));
    const auto origin = std::max(std::size_t(32), s.source_origin.size() + 1);
    const auto h = c::hydrogen_helium_history_payload_bound(
        p.history.maximum_fine_intervals, rows, s.model.species.size(), origin);
    if (!h) return {};
    b.count(*h);
    b.count(rows * 2, sizeof(c::HydrogenHeliumHistoryRow)); // both query receipts
    b.count(1, sizeof(c::BaryonAbundanceRow));
  }
  return b.result();
}
inline std::optional<std::size_t> Cohort::known_payload() const noexcept {
  Bytes b; b.count(1, sizeof(Cohort));
  if (source_) {
    b.count(1, sizeof(Source));
    const auto &in = source_->input;
    for (const auto &s : {&in.distribution_origin, &in.dependence_origin, &in.role, &in.temperature_role}) b.actual_string(*s);
    for (const auto &r : in.queries) b.actual_string(r.id);
    for (const auto &s : in.support) {
      b.actual_string(s.id); b.actual_string(s.source_origin); b.actual_string(s.mass_origin);
      b.actual_string(s.photon_temperature_origin);
      b.count(s.model.species.capacity(), sizeof(c::ThermalPhysicalSpecies));
    }
  }
  for (const auto &s : slots_) {
    b.count(s.abundance_attempt.rows.capacity(), sizeof(c::BaryonAbundanceRow));
    if (s.abundance.source()) {
      b.actual_string(s.abundance.source()->source_origin); b.actual_string(s.abundance.source()->mass_origin);
    }
    if (const auto *h = s.history.source()) {
      b.actual_string(h->nuclei_origin); b.count(h->model.species.capacity(), sizeof(c::ThermalPhysicalSpecies));
    }
    if (s.snapshot.emitted_thermal) b.count(s.snapshot.emitted_thermal->species.capacity(), sizeof(c::ThermalSpecies));
  }
  return b.result(); // Opaque nodes excluded here; full public envelopes remain reserved.
}
inline bool structural(const Input &in) noexcept {
  auto text = [](const std::string &s) { return !s.empty() && s.size() <= maximum_origin_bytes; };
  if (!text(in.distribution_origin) || !text(in.dependence_origin) || !text(in.role) ||
      !text(in.temperature_role) || in.requested_outputs != outputs ||
      bits(in.initial_redshift) != bits(initial_redshift) || bits(in.late_redshift) != bits(late_redshift)) return false;
  for (std::size_t s = 0; s < states; ++s) {
    const auto &x = in.support[s];
    if (!text(x.id) || !text(x.source_origin) || !text(x.mass_origin) ||
        !text(x.photon_temperature_origin) || bits(x.mass) != bits(relative_mass) ||
        x.model.species.size() > 16) return false;
    if (s && x.id == in.support[0].id) return false;
  }
  for (std::size_t r = 0; r < rows; ++r) {
    if (!text(in.queries[r].id)) return false;
    for (std::size_t j = 0; j < r; ++j) if (in.queries[r].id == in.queries[j].id) return false;
  }
  return true;
}
inline bool copied_string_fits(const std::string &s) noexcept {
  return s.capacity() != SIZE_MAX &&
      s.capacity() + 1 <= std::max(std::size_t(32), s.size() + 1);
}
inline bool source_copy_fits(const Input &in) noexcept {
  for (const auto &s : {&in.distribution_origin, &in.dependence_origin, &in.role, &in.temperature_role})
    if (!copied_string_fits(*s)) return false;
  for (const auto &r : in.queries) if (!copied_string_fits(r.id)) return false;
  for (const auto &s : in.support)
    for (const auto &text : {&s.id, &s.source_origin, &s.mass_origin, &s.photon_temperature_origin})
      if (!copied_string_fits(*text)) return false;
  for (const auto &s : in.support) if (s.model.species.capacity() != s.model.species.size()) return false;
  return true;
}
inline Cohort prepare(const Input &in, Policy p = {}) {
  Cohort out; out.policy_ = p;
  if (!profile() || !valid_policy(p) || !structural(in)) return out;
  const auto bound = Cohort::whole_live_bound(in, p);
  if (!bound || *bound > p.maximum_whole_live_bytes) { out.status_ = S::work_limit; return out; }
  out.preparation_bound_ = *bound;
  // Failure to acquire the source throws bad_alloc. Once acquired, all later
  // refusals retain the independent caller source even if native source absent.
  out.source_ = std::make_shared<Source>(in, p);
  if (!source_copy_fits(out.source_->input)) {
    out.status_ = S::work_limit;
    for (auto &slot : out.slots_) slot.status = S::work_limit;
    return out;
  }
  out.status_ = S::ok;
  auto preserve = [&](S s) { if (out.status_ == S::ok && s != S::ok) out.status_ = s; };
  for (std::size_t s = 0; s < states; ++s) {
    auto &slot = out.slots_[s]; const auto &x = out.source_->input.support[s];
    slot.budget.requested_work = p.history.maximum_total_work;
    slot.budget.requested_bytes = p.history.maximum_native_bytes;
    if (!out.charge(out.work_.complete_intakes)) { slot.status = S::work_limit; preserve(slot.status); continue; }
    if (in.role != abundance_history_cohort::role || in.temperature_role != abundance_history_cohort::temperature_role ||
        !x.model.species.empty()) { slot.status = S::outside_domain; preserve(slot.status); continue; }
    try {
      {
        c::BaryonAbundanceSource abundance{x.model.physical_baryon_density, x.helium_fraction,
            x.hydrogen_mass_kg, x.helium_mass_kg, x.source_origin, x.mass_origin};
        if (!out.charge(out.work_.abundance_prepares)) { slot.status = S::work_limit; preserve(slot.status); continue; }
        slot.abundance = c::prepare_baryon_abundance(abundance, {1, p.history.maximum_native_bytes});
      }
      if (slot.abundance.source() &&
          (!copied_string_fits(slot.abundance.source()->source_origin) ||
           !copied_string_fits(slot.abundance.source()->mass_origin))) {
        slot.status = S::work_limit; preserve(slot.status); continue;
      }
      if (slot.abundance.status() != S::ok) { slot.status = slot.abundance.status(); preserve(slot.status); continue; }
      if (!out.charge(out.work_.abundance_maps)) { slot.status = S::work_limit; preserve(slot.status); continue; }
      const double a = 1;
      slot.abundance_attempt = slot.abundance.evaluate({&a, 1}, {1, p.history.maximum_native_bytes});
      if (slot.abundance_attempt.status != S::ok || slot.abundance_attempt.rows.size() != 1) {
        slot.status = slot.abundance_attempt.status == S::ok ? S::invalid_input : slot.abundance_attempt.status;
        preserve(slot.status); continue;
      }
      slot.snapshot.nuclei = slot.abundance_attempt.rows.front();
      const auto &n = *slot.snapshot.nuclei;
      S mapping = n.admission_status;
      if (mapping == S::ok && (n.hydrogen_nuclei.status != S::ok || !n.hydrogen_nuclei.value))
        mapping = n.hydrogen_nuclei.status == S::ok ? S::invalid_input : n.hydrogen_nuclei.status;
      if (mapping == S::ok && (n.helium_nuclei.status != S::ok || !n.helium_nuclei.value))
        mapping = n.helium_nuclei.status == S::ok ? S::invalid_input : n.helium_nuclei.status;
      if (mapping == S::ok && (!(*n.hydrogen_nuclei.value > 0) || !(*n.helium_nuclei.value > 0))) mapping = S::outside_domain;
      if (mapping != S::ok) { slot.status = mapping; preserve(slot.status); continue; }
      slot.snapshot.emitted_nuclei_bits = {bits(*n.hydrogen_nuclei.value), bits(*n.helium_nuclei.value)};
      c::HydrogenHeliumHistoryRequest request{x.model, *n.hydrogen_nuclei.value, *n.helium_nuclei.value,
          x.source_origin, in.initial_redshift, in.late_redshift};
      const auto used = out.work_.checked_total(), known = out.known_payload();
      if (!used || !known || *used > p.maximum_preparation_work || *known > p.maximum_whole_live_bytes) {
        slot.status = S::work_limit; preserve(slot.status); continue;
      }
      auto child = p.history;
      child.maximum_total_work = std::min(child.maximum_total_work, p.maximum_preparation_work - *used);
      // Parent reservation includes the other complete opaque history envelope
      // and every caller/request/result lifetime phase. This child's entire
      // envelope is the only released reservation at its dispatch.
      const auto child_envelope = c::hydrogen_helium_history_payload_bound(
          p.history.maximum_fine_intervals, rows, 0,
          std::max(std::size_t(32), x.source_origin.size() + 1));
      if (!child_envelope || *child_envelope > *bound) {
        slot.status = S::work_limit; preserve(slot.status); continue;
      }
      const std::size_t parent_reservation = *bound - *child_envelope;
      if (parent_reservation > p.maximum_whole_live_bytes) {
        slot.status = S::work_limit; preserve(slot.status); continue;
      }
      child.maximum_native_bytes = std::min(child.maximum_native_bytes,
          p.maximum_whole_live_bytes - parent_reservation);
      slot.budget.served_work = child.maximum_total_work; slot.budget.served_bytes = child.maximum_native_bytes;
      slot.budget.history_attempted = true;
      slot.history = c::prepare_hydrogen_helium_history(request, child);
      const auto hw = slot.history.work();
      // Import every actual category independently, even if a parent contract
      // violation exceeds its cap. Never truncate a failed child's work ledger.
      const bool b = add(out.work_.background, hw.background_evaluations),
          m = add(out.work_.momentum, hw.momentum_callbacks),
          q = add(out.work_.charge, hw.initial_charge_evaluations),
          r = add(out.work_.rhs, hw.rhs_evaluations);
      out.work_.overflowed = out.work_.overflowed || !(b && m && q && r);
      if (const auto *w = slot.history.thermal_mapping_witnesses()) slot.snapshot.thermal_witnesses = *w;
      slot.snapshot.native_history_source_present = slot.history.source() != nullptr;
      if (const auto *b = slot.history.background()) slot.snapshot.emitted_thermal = b->source();
      slot.snapshot.maximum_log_excluded_heiii_activity = slot.history.maximum_log_excluded_heiii_activity();
      slot.status = slot.history.status(); preserve(slot.status);
      const auto actual_work = out.work_.checked_total();
      if (!actual_work || *actual_work > p.maximum_preparation_work) {
        slot.status = S::work_limit; preserve(slot.status);
      }
      if (const auto *native = slot.history.source(); native && !copied_string_fits(native->nuclei_origin)) {
        slot.status = S::work_limit; preserve(slot.status);
      }
    } catch (const std::bad_alloc &) { slot.status = S::work_limit; preserve(slot.status); }
  }
  const auto known = out.known_payload(), work = out.work_.checked_total();
  if (!known || *known > *bound || !work) preserve(S::work_limit);
  else { out.known_payload_ = *known; out.spent_ = *work; }
  return out;
}
inline Evaluation Cohort::evaluate() {
  Evaluation out; out.source = source_; out.conservative_whole_live_bytes = preparation_bound_;
  for (std::size_t s = 0; s < states; ++s) {
    out.snapshots[s] = slots_[s].snapshot;
    out.state_status[s] = slots_[s].status;
  }
  if (!source_) { out.status = status_; return out; }
  if (!profile()) return out;
  if (evaluation_count_ >= 2) { out.status = S::work_limit; return out; }
  ++evaluation_count_;
  auto charge = [&](std::size_t &field, std::size_t n = 1) {
    const auto used = out.work.checked_total();
    if (!used || *used > policy_.maximum_evaluation_work || n > policy_.maximum_evaluation_work - *used ||
        spent_ > policy_.maximum_two_evaluation_work || n > policy_.maximum_two_evaluation_work - spent_ ||
        !add(field, n) || !add(spent_, n)) return false;
    return true;
  };
  S failure = S::ok;
  auto preserve = [&](S s) { if (failure == S::ok && s != S::ok) failure = s; };
  std::array<double, rows> z{};
  std::array<std::array<W, axes>, states> values{}, errors{};
  std::array<bool, axes> constant{};
  for (std::size_t r = 0; r < rows; ++r) z[r] = source_->input.queries[r].redshift;
  for (std::size_t s = 0; s < states; ++s) {
    if (!charge(out.work.dispatches)) { preserve(S::work_limit); continue; }
    if (slots_[s].status != S::ok) { preserve(slots_[s].status); continue; }
    // Charge the requested row attempts before the pooled native dispatch.
    // Allocation refusal may leave no native rows, but still made this request.
    if (!charge(out.work.queried_rows, rows)) { preserve(S::work_limit); continue; }
    out.attempts[s] = slots_[s].history.evaluate(z, outputs, rows, policy_.history.maximum_native_bytes);
    const auto &attempt = out.attempts[s];
    preserve(attempt.status);
    if (attempt.status != S::ok || attempt.rows.size() != rows) {
      if (attempt.status == S::ok) preserve(S::invalid_input);
      continue;
    }
    for (std::size_t r = 0; r < rows; ++r) {
      const auto &row = attempt.rows[r];
      const std::array<const c::HydrogenHeliumHistoryValue *, 2> groups{
          &row.matter_temperature_kelvin, &row.thomson_opacity_per_redshift};
      for (std::size_t k = 0; k < 2; ++k) {
        if (!charge(out.work.collections)) { preserve(S::work_limit); continue; }
        const auto &v = *groups[k];
        if (v.status != S::ok || !v.value) {
          preserve(v.status == S::ok ? S::invalid_input : v.status); continue;
        }
        values[s][2 * r + k] = *v.value; errors[s][2 * r + k] = v.absolute_error_estimate;
      }
    }
  }
  Bytes result; result.count(1, sizeof(Evaluation));
  for (const auto &a : out.attempts) result.count(a.rows.capacity(), sizeof(c::HydrogenHeliumHistoryRow));
  for (const auto &s : out.snapshots) if (s.emitted_thermal) result.count(s.emitted_thermal->species.capacity(), sizeof(c::ThermalSpecies));
  const auto actual = result.result();
  if (!actual || *actual > preparation_bound_) preserve(S::work_limit);
  else out.known_result_payload_bytes = *actual;
  if (failure == S::ok) {
    for (std::size_t i = 0; i < axes; ++i) {
      if (!charge(out.work.witnesses)) { preserve(S::work_limit); break; }
      constant[i] = constant_axis(source_->input, i);
    }
    if (failure == S::ok) {
      const auto before = out.work.checked_total();
      const auto remaining_scope = spent_ <= policy_.maximum_two_evaluation_work
          ? policy_.maximum_two_evaluation_work - spent_ : 0;
      std::size_t served = before.value_or(SIZE_MAX);
      if (!before || !add(served, remaining_scope)) preserve(S::work_limit);
      else {
        served = std::min(served, policy_.maximum_evaluation_work);
        const auto reduction_status = reduce(values, errors, constant, policy_.covariance_error_fraction,
                                             out.work, served, out.moments);
        const auto after = out.work.checked_total();
        if (!after || *after < *before || !add(spent_, *after - *before)) {
          out.moments.reset(); preserve(S::work_limit);
        } else preserve(reduction_status);
      }
    }
  }
  if (failure == S::ok && status_ != S::ok) { out.moments.reset(); failure = status_; }
  out.status = failure;
  return out;
}
inline std::size_t Cohort::serialization_budget(const Evaluation &receipt) const noexcept {
  const auto work = receipt.work.checked_total();
  if (serialization_recorded_ || !work || *work > policy_.maximum_evaluation_work ||
      spent_ > policy_.maximum_two_evaluation_work) return 0;
  return std::min(policy_.maximum_evaluation_work - *work,
                  policy_.maximum_two_evaluation_work - spent_);
}
inline S Cohort::record_serialization(Evaluation &receipt, std::size_t fields,
                                      std::size_t nibbles, std::size_t decompositions) noexcept {
  std::size_t n = fields;
  const auto budget = serialization_budget(receipt);
  if (receipt.source != source_ || serialization_recorded_ ||
      !add(n, nibbles) || !add(n, decompositions) || n > budget)
    return S::work_limit;
  serialization_recorded_ = true;
  receipt.work.serialization_fields = fields;
  receipt.work.serialization_wide_nibbles = nibbles;
  receipt.work.serialization_wide_decompositions = decompositions;
  if (!add(spent_, n)) return S::work_limit;
  return S::ok;
}
} // namespace abundance_history_cohort
