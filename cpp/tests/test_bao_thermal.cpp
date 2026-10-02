#include "test_bao_thermal_fixture.hpp"
#include <cfenv>
#include <iostream>
#include <limits>
using namespace thermal_bao_test;
int main() {
  try {
    const W nominal[]{20, 20, 20};
    auto d = prepared(nominal);
    auto req = massive();
    auto p = policy();
    // Strong synthetic residuals may exceed strict density projection; request
    // ratios alone for resource/ownership controls, then center observations.
    p.requested = 2;
    auto ratio = evaluate(d, req, p);
    need(ratio.slots.size() == 1 && ratio.slots[0].predictions.size() == 3,
         "thermal predictions admitted");
    std::vector<W> center(ratio.slots[0].predictions.begin(),
                          ratio.slots[0].predictions.end());
    d = prepared(center);
    const auto *matrix = d.source().covariance.data();
    const auto *queries = d.source().queries.data();
    auto full = evaluate(d, req);
    auto &s = full.slots[0];
    need(s.result &&
             s.density_state.availability == cosmology::Availability::available,
         "strict thermal density accepted");
    need(d.source().covariance.data() == matrix &&
             d.source().queries.data() == queries,
         "retained immutable observation unchanged");
    need(s.preparation_status == S::ok && s.preparation_callbacks > 0 &&
             s.outer_callbacks > 0 &&
             s.momentum_callbacks > s.preparation_callbacks,
         "preparation outer nested work visible");
    need(s.callbacks == s.outer_callbacks + s.momentum_callbacks &&
             full.callbacks == s.callbacks &&
             full.preparation_callbacks == s.preparation_callbacks,
         "combined callbacks charged once");
    auto provider = cosmology::prepare_thermal_observables(req);
    const double zs[]{1, 1, 1};
    auto predicted = provider.evaluate(
        zs,
        cosmology::early_late_mask(cosmology::EarlyLateOutput::dm_over_rs) |
            cosmology::early_late_mask(cosmology::EarlyLateOutput::dh_over_rs) |
            cosmology::early_late_mask(cosmology::EarlyLateOutput::dv_over_rs));
    need(s.preparation_callbacks ==
                 provider.background().preparation_callbacks() &&
             s.outer_callbacks == predicted.outer_callbacks &&
             s.momentum_callbacks ==
                 predicted.momentum_callbacks + s.preparation_callbacks,
         "one mapping/preparation state per model");
    p = policy();
    p.maximum_projection_log_density_error = 1e-30;
    auto refused = evaluate(d, req, p);
    need(!refused.slots[0].result && refused.slots[0].predictions.size() == 3 &&
             refused.slots[0].residuals.size() == 3,
         "projection refusal preserves outputs");
    need(refused.slots[0].predictions_state.availability ==
                 cosmology::Availability::available &&
             refused.slots[0].residuals_state.availability ==
                 cosmology::Availability::available &&
             refused.slots[0].density_state.numerical_status ==
                 S::conditioning_budget_exceeded,
         "independent requested states");
    for (unsigned mask : {1u, 2u, 4u}) {
      p = policy();
      p.requested = mask;
      auto only = evaluate(d, req, p);
      auto &o = only.slots[0];
      need(bool(o.result) == bool(mask & 1) &&
               o.predictions.empty() == !(mask & 2) &&
               o.residuals.empty() == !(mask & 4),
           "output omission");
      need((mask & 1) || o.density_state.availability ==
                             cosmology::Availability::not_requested,
           "density state omission");
    }
    p = policy();
    auto bad = req;
    bad.model.physical_cdm_density = -1;
    cosmology::ThermalObservableRequest mixed[]{req, bad, req};
    auto mix = d.evaluate_thermal(mixed, p);
    need(mix.status == statistics::DensityStatus::finite &&
             mix.slots.size() == 3 && mix.slots[0].result &&
             !mix.slots[1].result && mix.slots[2].result,
         "mixed independently failed mapping");
    need(mix.slots[1].source.model.physical_cdm_density == -1 &&
             mix.slots[1].preparation_callbacks == 0,
         "failed source remains inspectable");
    p.maximum_total_callbacks = s.callbacks;
    auto limited = d.evaluate_thermal(mixed, p);
    need(limited.callbacks <= p.maximum_total_callbacks &&
             limited.slots[0].result && !limited.slots[2].result,
         "whole batch includes normalization and failed work");
    p = policy();
    p.predictions.thermal.maximum_callbacks_per_evaluation = 32;
    auto prepfailed = evaluate(d, req, p);
    need(prepfailed.slots[0].preparation_callbacks > 0 &&
             prepfailed.slots[0].outer_callbacks == 0 &&
             !prepfailed.slots[0].result &&
             prepfailed.callbacks <= p.maximum_total_callbacks,
         "failed preparation charged");
    p = policy();
    p.maximum_total_callbacks = 0;
    auto zerowork = evaluate(d, req, p);
    need(zerowork.callbacks == 0 && !zerowork.slots[0].result,
         "zero callbacks forbid normalization");
    for (unsigned quota = 0; quota < 4; ++quota) {
      p = policy();
      if (quota == 0)
        p.maximum_models = 0;
      if (quota == 1)
        p.maximum_queries = 2;
      if (quota == 2)
        p.maximum_string_bytes = 1;
      if (quota == 3)
        p.maximum_native_bytes = 1;
      auto r = evaluate(d, req, p);
      need(r.slots.empty() && r.callbacks == 0 &&
               r.numerical_status == S::work_limit,
           "quota admission before copy/work");
    }
    p = policy();
    p.predictions.thermal.maximum_species = 0;
    auto species_limit = evaluate(d, req, p);
    need(species_limit.callbacks == 0 &&
             species_limit.slots[0].preparation_status == S::work_limit &&
             !species_limit.slots[0].result,
         "species quota forbids preparation work");
    auto too_many = req;
    too_many.model.species.resize(17, req.model.species[0]);
    auto structural = evaluate(d, too_many);
    need(structural.slots.empty() && structural.callbacks == 0 &&
             structural.numerical_status == S::outside_domain,
         "unsupported species count rejected before source copies");
    for (unsigned mask : {0u, 8u}) {
      p = policy();
      p.requested = mask;
      need(evaluate(d, req, p).slots.empty(), "invalid output mask");
    }
    size_t origins = 0;
    for (auto *o : {&req.drag_origin, &req.source_origin})
      origins += std::max(size_t(32), o->size() + 1);
    auto bytes = bao::thermal_density_payload_bound(1, 3, origins, 1, 1, 7);
    need(bytes.has_value(), "payload bound available");
    p = policy();
    p.maximum_native_bytes = *bytes - 1;
    need(evaluate(d, req, p).slots.empty(), "payload threshold refused");
    p.maximum_native_bytes = *bytes;
    need(evaluate(d, req, p).slots[0].result.has_value(),
         "payload threshold accepted");
    need(!bao::thermal_density_payload_bound(SIZE_MAX, 3, 64, 1, 1, 7),
         "payload overflow");
    req.model.species.clear();
    req.drag_origin.clear();
    need(s.source.model.species.size() == 1 && !s.source.drag_origin.empty(),
         "returned source owns independent payload");
    auto copy = d;
    auto moved = std::move(d);
    auto *self = &moved;
    moved = std::move(*self);
    req = massive();
    need(evaluate(d, req).slots.empty() &&
             evaluate(moved, req).slots[0].result &&
             evaluate(copy, req).slots[0].result,
         "retained source copy/move lifetime");
    p = policy();
    p.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    need(evaluate(moved, req, p).status ==
             statistics::DensityStatus::incompatible_metadata,
         "retained Gaussian arithmetic required");
    const auto rounding = std::fegetround();
    need(std::fesetround(FE_DOWNWARD) == 0, "set rounding");
    auto rounded = evaluate(moved, req);
    need(std::fesetround(rounding) == 0, "restore rounding");
    need(!rounded.slots[0].result && rounded.callbacks == 0,
         "unsupported arithmetic cannot expose ratios");
    std::cout << "PASS " << checks << " thermal BAO owner controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
