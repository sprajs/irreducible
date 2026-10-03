#pragma once
// Fixed proof-output bridge only. Integer/double bit fields and exact binary
// wide significands preserve inputs and diagnostics without decimal roundtrip.
#include "abundance_history_cohort.hpp"
#include <cstdio>
namespace abundance_history_cohort {
struct Wire {
  std::array<char, maximum_wire_bytes> bytes{};
  std::size_t used = 0, fields = 0, wide_nibbles = 0, wide_decompositions = 0,
      work_cap = evaluation_work_cap;
  bool failed = false, finished = false;
  static constexpr std::string_view hex = "0123456789abcdef";
  static constexpr std::size_t trailer_reservation = 256;
  std::size_t work() const noexcept { return fields + wide_nibbles + wide_decompositions; }
  void limit_work(std::size_t cap) noexcept {
    work_cap = std::min(cap, evaluation_work_cap);
    if (work() + 5 > work_cap) failed = true;
  }
  bool begin(std::string_view prefix, std::string_view key, std::size_t payload,
             std::size_t extra_work = 0) noexcept {
    if (failed || finished) return false;
    std::size_t n = prefix.size();
    if (!add(n, key.size()) || !add(n, payload) || !add(n, 2) ||
        fields >= maximum_wire_fields - 5 || work() > work_cap ||
        extra_work + 6 > work_cap - work() || n > maximum_wire_bytes - trailer_reservation - used) {
      failed = true; return false;
    }
    ++fields; append(prefix); append(key); append('='); return true;
  }
  void append(char c) noexcept { bytes[used++] = c; }
  void append(std::string_view value) noexcept { for (char c : value) append(c); }
  void hex64(std::uint64_t value) noexcept {
    for (unsigned n = 16; n-- > 0;) append(hex[(value >> (4 * n)) & 15]);
  }
  void integer(std::string_view prefix, std::string_view key, std::uint64_t value) noexcept {
    if (!begin(prefix, key, 16)) return;
    hex64(value); append('\n');
  }
  void scalar(std::string_view prefix, std::string_view key, double value) noexcept {
    integer(prefix, key, bits(value));
  }
  void text(std::string_view prefix, std::string_view key, std::string_view value) noexcept {
    if (value.size() > (SIZE_MAX - 17) / 2) { failed = true; return; }
    if (!begin(prefix, key, 17 + 2 * value.size())) return;
    hex64(value.size()); append(':');
    for (unsigned char c : value) { append(hex[c >> 4]); append(hex[c & 15]); }
    append('\n');
  }
  void wide(std::string_view prefix, std::string_view key, W value) noexcept {
    // Binary64-style raw storage would serialize x87 padding. Instead emit
    // exact sign/class/exponent/64 significand bits under this wire profile.
    if (std::numeric_limits<W>::radix != 2 || std::numeric_limits<W>::digits != 64) {
      failed = true; return;
    }
    if (!begin(prefix, key, 1 + 1 + 1 + 1 + 1 + 16 + 1 + 16, 17)) return;
    append(std::signbit(value) ? '1' : '0'); append(':');
    if (!std::isfinite(value)) {
      append(std::isnan(value) ? 'n' : 'i'); append(':');
      hex64(0); append(':'); hex64(0); append('\n'); return;
    }
    append('f'); append(':');
    int exponent = 0;
    W fraction = std::frexp(std::abs(value), &exponent); ++wide_decompositions;
    append(exponent < 0 ? '1' : '0');
    hex64(exponent < 0 ? std::uint64_t(-std::int64_t(exponent)) : std::uint64_t(exponent));
    append(':');
    for (unsigned n = 0; n < 16; ++n) {
      fraction *= 16;
      const unsigned digit = static_cast<unsigned>(fraction);
      fraction -= digit; ++wide_nibbles; append(hex[digit]);
    }
    append('\n');
    if (fraction != 0) failed = true;
  }
  void finish() noexcept {
    if (finished) return;
    // Five fixed trailer fields are reserved before every attempted body field.
    const auto status = failed ? S::work_limit : S::ok;
    constexpr std::array<std::string_view, 5> keys{
        "wire.fields", "wire.wide_nibbles", "wire.wide_decompositions", "wire.bytes", "wire.status"};
    std::size_t trailer = 0;
    for (auto key : keys) trailer += key.size() + 18;
    const std::array<std::uint64_t, 5> values{fields + 5, wide_nibbles, wide_decompositions,
                                            used + trailer, std::uint64_t(status)};
    for (std::size_t i = 0; i < keys.size(); ++i) {
      append(keys[i]); append('='); hex64(values[i]); append('\n'); ++fields;
    }
    finished = true;
  }
  bool publish() noexcept {
    finish();
    return std::fwrite(bytes.data(), 1, used, stdout) == used && std::fflush(stdout) == 0;
  }
};
static_assert(sizeof(Wire) <= maximum_wire_bytes + 128);
inline std::uint64_t wire_status(S s) noexcept { return static_cast<std::uint64_t>(s); }
inline void wire_input(Wire &w, const Input &in, const Policy &p) noexcept {
  w.text("", "wire.identity", "fixed-original-two-half-history-output/v1");
  w.text("", "model", model_id); w.text("", "method", method_id);
  w.text("", "abundance_assets", c::baryon_abundance_assets_id);
  w.text("", "history_model", c::hydrogen_helium_history_model_id);
  w.text("", "history_method", c::hydrogen_helium_history_method_id);
  w.text("", "role", in.role); w.text("", "temperature_role", in.temperature_role);
  w.text("", "distribution_origin", in.distribution_origin);
  w.text("", "dependence_origin", in.dependence_origin);
  w.scalar("", "initial_z_bits", in.initial_redshift); w.scalar("", "late_z_bits", in.late_redshift);
  w.integer("", "output_mask", in.requested_outputs);
  for (std::size_t s = 0; s < states; ++s) {
    const std::string_view prefix = s ? "B." : "A.";
    const auto &x = in.support[s]; const auto &m = x.model;
    w.text(prefix, "id", x.id); w.text(prefix, "source_origin", x.source_origin);
    w.text(prefix, "mass_origin", x.mass_origin); w.text(prefix, "photon_origin", x.photon_temperature_origin);
    w.scalar(prefix, "mass_bits", x.mass); w.scalar(prefix, "probability_bits", probability);
    w.scalar(prefix, "H0_bits", m.h0_km_s_mpc); w.scalar(prefix, "omega_b_bits", m.physical_baryon_density);
    w.scalar(prefix, "omega_cdm_bits", m.physical_cdm_density); w.scalar(prefix, "Tcmb_bits", m.tcmb_kelvin);
    w.scalar(prefix, "omega_other_bits", m.physical_massless_nonphoton_density);
    w.scalar(prefix, "Y_bits", x.helium_fraction); w.scalar(prefix, "mH_bits", x.hydrogen_mass_kg);
    w.scalar(prefix, "mHe_bits", x.helium_mass_kg); w.integer(prefix, "species_count", m.species.size());
  }
  for (std::size_t r = 0; r < rows; ++r) {
    const std::array<std::string_view, rows> prefix{"row0.", "row1.", "row2."};
    w.text(prefix[r], "id", in.queries[r].id); w.scalar(prefix[r], "z_bits", in.queries[r].redshift);
  }
  w.integer("policy.", "preparation_work", p.maximum_preparation_work);
  w.integer("policy.", "evaluation_work", p.maximum_evaluation_work);
  w.integer("policy.", "two_evaluation_work", p.maximum_two_evaluation_work);
  w.integer("policy.", "whole_live_bytes", p.maximum_whole_live_bytes);
  w.scalar("policy.", "covariance_fraction_bits", p.covariance_error_fraction);
  const auto &h = p.history;
  w.integer("history.", "N", h.base_intervals); w.integer("history.", "maximum_fine", h.maximum_fine_intervals);
  w.integer("history.", "requested_work", h.maximum_total_work); w.integer("history.", "requested_bytes", h.maximum_native_bytes);
  w.scalar("history.", "fraction_atol_bits", h.absolute_fraction_tolerance);
  w.scalar("history.", "fraction_rtol_bits", h.relative_fraction_tolerance);
  w.scalar("history.", "temperature_atol_bits", h.absolute_temperature_tolerance_kelvin);
  w.scalar("history.", "temperature_rtol_bits", h.relative_temperature_tolerance);
  w.scalar("history.", "opacity_atol_bits", h.absolute_opacity_tolerance);
  w.scalar("history.", "opacity_rtol_bits", h.relative_opacity_tolerance);
  w.scalar("thermal.", "atol_bits", h.thermal.absolute_tolerance);
  w.scalar("thermal.", "rtol_bits", h.thermal.relative_tolerance);
  w.integer("thermal.", "per_evaluation_callbacks", h.thermal.maximum_callbacks_per_evaluation);
  w.integer("thermal.", "total_callbacks", h.thermal.maximum_total_callbacks);
  w.integer("thermal.", "depth", h.thermal.maximum_depth); w.integer("thermal.", "points", h.thermal.maximum_points);
  w.integer("thermal.", "species", h.thermal.maximum_species); w.integer("thermal.", "bytes", h.thermal.maximum_native_bytes);
  w.integer("thermal.", "method", static_cast<std::uint64_t>(h.thermal.momentum_method));
}
inline void wire_work(Wire &w, std::string_view prefix, const EvaluationWork &x) noexcept {
  w.integer(prefix, "dispatches", x.dispatches); w.integer(prefix, "queried_rows", x.queried_rows);
  w.integer(prefix, "collections", x.collections); w.integer(prefix, "witnesses", x.witnesses);
  w.integer(prefix, "means", x.means); w.integer(prefix, "differences", x.differences);
  w.integer(prefix, "products", x.products); w.integer(prefix, "projections", x.projections);
}
inline void wire_owner(Wire &w, const Cohort &owner) noexcept {
  w.integer("owner.", "status", wire_status(owner.status()));
  w.integer("owner.", "source_present", bool(owner.source()));
  w.integer("owner.", "whole_live_bound", owner.conservative_whole_live_bytes());
  w.integer("owner.", "known_retained_bytes", owner.known_retained_payload_bytes());
  const auto x = owner.preparation_work();
  w.integer("prep.", "intakes", x.complete_intakes); w.integer("prep.", "abundance_prepares", x.abundance_prepares);
  w.integer("prep.", "abundance_maps", x.abundance_maps); w.integer("prep.", "background", x.background);
  w.integer("prep.", "momentum", x.momentum); w.integer("prep.", "charge", x.charge);
  w.integer("prep.", "rhs", x.rhs); w.integer("prep.", "overflowed", x.overflowed);
  for (std::size_t s = 0; s < states; ++s) {
    const std::string_view prefix = s ? "B.native." : "A.native.";
    const auto &slot = owner.slot(s); const auto &n = slot.snapshot;
    w.integer(prefix, "status", wire_status(slot.status));
    w.integer(prefix, "abundance_owner_status", wire_status(slot.abundance.status()));
    w.integer(prefix, "abundance_batch_status", wire_status(slot.abundance_attempt.status));
    w.integer(prefix, "history_owner_status", wire_status(slot.history.status()));
    const auto hw = slot.history.work();
    w.integer(prefix, "background_work", hw.background_evaluations);
    w.integer(prefix, "momentum_work", hw.momentum_callbacks);
    w.integer(prefix, "charge_work", hw.initial_charge_evaluations);
    w.integer(prefix, "rhs_work", hw.rhs_evaluations);
    w.integer(prefix, "history_attempted", slot.budget.history_attempted);
    w.integer(prefix, "requested_work", slot.budget.requested_work);
    w.integer(prefix, "requested_bytes", slot.budget.requested_bytes);
    w.integer(prefix, "served_work_present", slot.budget.served_work.has_value());
    w.integer(prefix, "served_work", slot.budget.served_work.value_or(0));
    w.integer(prefix, "served_bytes_present", slot.budget.served_bytes.has_value());
    w.integer(prefix, "served_bytes", slot.budget.served_bytes.value_or(0));
    w.integer(prefix, "source_present", n.native_history_source_present);
    w.integer(prefix, "nuclei_present", n.nuclei.has_value());
    if (n.nuclei) {
      w.integer(prefix, "map_admission", wire_status(n.nuclei->admission_status));
      w.integer(prefix, "nH_status", wire_status(n.nuclei->hydrogen_nuclei.status));
      w.integer(prefix, "nHe_status", wire_status(n.nuclei->helium_nuclei.status));
      w.integer(prefix, "nH_present", n.nuclei->hydrogen_nuclei.value.has_value());
      w.integer(prefix, "nHe_present", n.nuclei->helium_nuclei.value.has_value());
      if (n.nuclei->hydrogen_nuclei.value) w.scalar(prefix, "nH_bits", *n.nuclei->hydrogen_nuclei.value);
      if (n.nuclei->helium_nuclei.value) w.scalar(prefix, "nHe_bits", *n.nuclei->helium_nuclei.value);
      w.scalar(prefix, "nH_map_error_bits", n.nuclei->hydrogen_nuclei.relative_arithmetic_estimate);
      w.scalar(prefix, "nHe_map_error_bits", n.nuclei->helium_nuclei.relative_arithmetic_estimate);
    }
    w.integer(prefix, "emitted_thermal_present", n.emitted_thermal.has_value());
    if (n.emitted_thermal) {
      const auto &t = *n.emitted_thermal;
      w.scalar(prefix, "emitted_H0_bits", t.h0_km_s_mpc); w.scalar(prefix, "Omega_gamma_bits", t.omega_gamma);
      w.scalar(prefix, "Omega_b_bits", t.omega_b); w.scalar(prefix, "Omega_cdm_bits", t.omega_cdm);
      w.scalar(prefix, "Omega_other_bits", t.omega_massless_nonphoton); w.integer(prefix, "emitted_species_count", t.species.size());
    }
    w.integer(prefix, "thermal_witnesses_present", n.thermal_witnesses.has_value());
    if (n.thermal_witnesses) {
      const std::array<std::string_view, 4> ap{s ? "B.map.gamma." : "A.map.gamma.", s ? "B.map.b." : "A.map.b.",
          s ? "B.map.cdm." : "A.map.cdm.", s ? "B.map.other." : "A.map.other."};
      for (std::size_t k = 0; k < 4; ++k) {
        const auto &q = (*n.thermal_witnesses)[k];
        w.wide(ap[k], "wide", q.wide_value); w.wide(ap[k], "operation", q.wide_operation_estimate);
        w.wide(ap[k], "cast_loss", q.measured_absolute_cast_loss); w.scalar(ap[k], "emitted_bits", q.emitted_value);
      }
    }
    w.integer(prefix, "maximum_Heiii_log_present", n.maximum_log_excluded_heiii_activity.has_value());
    if (n.maximum_log_excluded_heiii_activity) w.scalar(prefix, "maximum_Heiii_log_bits", *n.maximum_log_excluded_heiii_activity);
  }
}
inline void wire_evaluations(Wire &w, const Evaluation &first, const Evaluation &last) noexcept {
  w.integer("first.", "status", wire_status(first.status)); wire_work(w, "first.work.", first.work);
  w.integer("last.", "status", wire_status(last.status)); wire_work(w, "last.work.", last.work);
  w.integer("last.", "known_result_bytes", last.known_result_payload_bytes);
  for (std::size_t s = 0; s < states; ++s) {
    const std::string_view prefix = s ? "B.query." : "A.query.";
    w.integer(prefix, "prepared_status", wire_status(last.state_status[s]));
    w.integer(prefix, "batch_status", wire_status(last.attempts[s].status));
    w.integer(prefix, "mask", last.attempts[s].requested_outputs);
    w.integer(prefix, "actual_rows", last.attempts[s].rows.size());
    const std::array<std::string_view, rows> rp{s ? "B.row0." : "A.row0.", s ? "B.row1." : "A.row1.", s ? "B.row2." : "A.row2."};
    for (std::size_t r = 0; r < last.attempts[s].rows.size() && r < rows; ++r) {
      const auto &row = last.attempts[s].rows[r]; w.scalar(rp[r], "z_bits", row.redshift);
      const std::array<const c::HydrogenHeliumHistoryValue *, 2> values{&row.matter_temperature_kelvin, &row.thomson_opacity_per_redshift};
      const std::array<std::string_view, 2> names{"Tm.", "opacity."};
      // Prefix components use a fixed local array, not allocation or formatting.
      for (std::size_t k = 0; k < 2; ++k) {
        std::array<char, 32> key{}; std::size_t len = 0;
        for (char ch : rp[r]) key[len++] = ch;
        for (char ch : names[k]) key[len++] = ch;
        const std::string_view p{key.data(), len}; const auto &v = *values[k];
        w.integer(p, "status", wire_status(v.status)); w.integer(p, "present", v.value.has_value());
        if (v.value) w.scalar(p, "value_bits", *v.value);
        w.scalar(p, "error_bits", v.absolute_error_estimate);
      }
      w.integer(rp[r], "omitted_H_present", row.hydrogen_ionized_fraction.value.has_value());
      w.integer(rp[r], "omitted_He_present", row.helium_singly_ionized_fraction.value.has_value());
      w.integer(rp[r], "omitted_ne_present", row.electron_number_density_per_cubic_metre.value.has_value());
    }
  }
  w.integer("last.", "moments_present", last.moments.has_value());
  w.integer("", "physical_joint_law_qualified", 0); w.integer("", "complete_reference_qualified", 0);
  if (!last.moments) return;
  const auto &m = *last.moments;
  for (std::size_t i = 0; i < axes; ++i) {
    // Repeated keys form ordered arrays; an independent parser preserves every
    // occurrence rather than using a map that silently overwrites entries.
    w.integer("mean.", "axis", i); w.scalar("mean.", "value_bits", m.mean[i]);
    w.scalar("mean.", "error_bits", m.mean_error[i]); w.scalar("mean.", "arithmetic_bits", m.mean_cast_and_round[i]);
    w.scalar("difference.", "arithmetic_bits", m.difference_round[i]);
  }
  for (std::size_t ij = 0; ij < covariance_cells; ++ij) {
    w.integer("covariance.", "cell", ij); w.scalar("covariance.", "value_bits", m.covariance[ij]);
    w.scalar("covariance.", "error_bits", m.covariance_error[ij]);
    w.scalar("covariance.", "arithmetic_bits", m.covariance_cast_and_round[ij]);
  }
}
} // namespace abundance_history_cohort
