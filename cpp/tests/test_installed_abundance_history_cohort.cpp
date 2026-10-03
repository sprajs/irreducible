// Independent installed-header/static-library consumer; local proof helper
// includes public irred headers only. No private node, loading or cell access.
#include "abundance_history_cohort_wire.hpp"
#include <iostream>
#include <stdexcept>
namespace h = abundance_history_cohort;
int main(int argc, char **) {
  h::Wire wire;
  bool publish_attempted = false;
  try {
    if (argc != 1) throw std::runtime_error("proof caller accepts no input overrides");
    auto input = h::original_input();
    const h::Policy policy;
    h::wire_input(wire, input, policy); // originally offered bits precede acquisition
    auto owner = h::prepare(input, policy);
    h::wire_owner(wire, owner); // preserve every earned/refused preparation prefix
    if (owner.status() != h::S::ok || !owner.source()) throw std::runtime_error("fixed synthetic histories refused");
    const auto work = owner.preparation_work().checked_total();
    if (!work || *work > h::preparation_work_cap ||
        owner.conservative_whole_live_bytes() > h::whole_live_payload_cap)
      throw std::runtime_error("bounded preparation receipt");
    input.support[0].model.tcmb_kelvin = 0;
    input.support[0].source_origin = "changed borrowed caller input";
    input.queries[0].redshift = 0;
    auto retained = std::move(owner);
    if (owner.source() || owner.status() != h::S::invalid_input)
      throw std::runtime_error("move invalidates source owner");
    retained = std::move(retained);
    auto first = retained.evaluate();
    auto final = retained.evaluate();
    wire.limit_work(retained.serialization_budget(final));
    h::wire_evaluations(wire, first, final);
    wire.finish();
    if (wire.failed || retained.record_serialization(final, wire.fields, wire.wide_nibbles,
                                                    wire.wide_decompositions) != h::S::ok)
      throw std::runtime_error("fixed output serialization work/byte cap");
    if (first.status != h::S::ok || !first.moments || first.source != retained.source())
      throw std::runtime_error("complete original two-half moments");
    if (final.status != h::S::ok || !final.moments ||
        retained.preparation_work().checked_total() != work)
      throw std::runtime_error("repeated pooled evaluation adds no trajectory work");
    if (wire.fields != h::successful_wire_fields || wire.wide_nibbles != h::successful_wire_nibbles ||
        wire.wide_decompositions != h::successful_wire_decompositions || wire.work() != h::successful_wire_work)
      throw std::runtime_error("frozen complete successful output graph");
    for (std::size_t s = 0; s < h::states; ++s) {
      const auto &slot = retained.slot(s);
      if (!slot.budget.history_attempted || !slot.budget.served_work || !slot.budget.served_bytes ||
          slot.budget.requested_work != h::child_work_cap || slot.budget.requested_bytes != h::child_payload_cap ||
          !final.snapshots[s].nuclei || !final.snapshots[s].emitted_thermal ||
          !final.snapshots[s].thermal_witnesses || !final.snapshots[s].native_history_source_present)
        throw std::runtime_error("original and emitted source/budget roles retained");
      for (std::size_t r = 0; r < h::rows; ++r) {
        const auto &row = final.attempts[s].rows.at(r);
        if (h::bits(row.redshift) != h::bits(final.source->input.queries[r].redshift) ||
            !row.matter_temperature_kelvin.value || !row.thomson_opacity_per_redshift.value ||
            row.hydrogen_ionized_fraction.value || row.helium_singly_ionized_fraction.value ||
            row.electron_number_density_per_cubic_metre.value)
          throw std::runtime_error("literal ordered mask24 outputs");
      }
    }
    const auto first_work = first.work.checked_total(), final_work = final.work.checked_total();
    std::size_t scope_work = *work;
    if (!first_work || !final_work || *first_work > h::evaluation_work_cap || *final_work > h::evaluation_work_cap ||
        !h::add(scope_work, *first_work) || !h::add(scope_work, *final_work) || scope_work > h::two_evaluation_work_cap)
      throw std::runtime_error("caller logical evaluation and serialization cap");
    const auto source = final.source;
    retained = h::Cohort{};
    if (h::bits(source->input.support[0].model.tcmb_kelvin) != h::bits(h::tcmb_a) ||
        source->input.support[0].source_origin != h::source_origin_a ||
        source->probabilities[0] != .5 || source->probabilities[1] != .5 || !final.moments)
      throw std::runtime_error("receipt outlives inputs and history owners");
    publish_attempted = true;
    if (!wire.publish()) throw std::runtime_error("fixed output publication failure");
  } catch (const std::exception &e) {
    wire.failed = true;
    if (!publish_attempted) wire.publish();
    std::cerr << "FAIL synthetic_history_cohort: " << e.what() << '\n'; return 1;
  }
}
