// Independent installed-header/static-library consumer; local proof helper
// includes public irred headers only. No private node, loading or cell access.
#include "abundance_history_cohort.hpp"
#include <iostream>
#include <stdexcept>
namespace h = abundance_history_cohort;
int main(int argc, char **) {
  try {
    if (argc != 1) throw std::runtime_error("proof caller accepts no input overrides");
    auto input = h::original_input();
    auto owner = h::prepare(input);
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
    if (first.status != h::S::ok || !first.moments || first.source != retained.source())
      throw std::runtime_error("complete original two-half moments");
    auto final = retained.evaluate();
    if (final.status != h::S::ok || !final.moments ||
        retained.preparation_work().checked_total() != work)
      throw std::runtime_error("repeated pooled evaluation adds no trajectory work");
    const auto first_work = first.work.checked_total(), final_work = final.work.checked_total();
    std::size_t scope_work = *work;
    if (!first_work || !final_work || *first_work > h::evaluation_work_cap || *final_work > h::evaluation_work_cap ||
        !h::add(scope_work, *first_work) || !h::add(scope_work, *final_work) || scope_work > h::two_evaluation_work_cap)
      throw std::runtime_error("caller logical evaluation cap");
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
    const auto source = final.source;
    retained = h::Cohort{};
    if (h::bits(source->input.support[0].model.tcmb_kelvin) != h::bits(h::tcmb_a) ||
        source->input.support[0].source_origin != h::source_origin_a ||
        source->probabilities[0] != .5 || source->probabilities[1] != .5 || !final.moments)
      throw std::runtime_error("receipt outlives inputs and history owners");
    std::cout << "synthetic_history_cohort=accepted source_role=exact_emitted_working_input"
              << " physical_law=unqualified complete_reference=withheld preparation_work=" << *work
              << " whole_live_bound=" << final.conservative_whole_live_bytes << '\n';
  } catch (const std::exception &e) {
    std::cerr << "FAIL synthetic_history_cohort: " << e.what() << '\n'; return 1;
  }
}
