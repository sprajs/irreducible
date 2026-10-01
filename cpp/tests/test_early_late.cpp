#include "irred/early_late.hpp"
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
int main() {
  unsigned checks = 0;
  auto require = [&](bool b, const char *s) {
    ++checks;
    if (!b)
      throw std::runtime_error(s);
  };
  EarlyLatePolicy p{
      1e-9,        2e-11,
      1e-11,       5e-11,
      1000000,     4000000,
      40,          64,
      1024 * 1024, {1e-9, 2e-11, 1000000, 40, 1, 1000000, 1024 * 1024}};
  const unsigned all = (1u << early_late_output_count) - 1;
  require(early_late_mask(EarlyLateOutput::count) == 0,
          "sentinel output has no mask");
  require(early_late_mask(static_cast<EarlyLateOutput>(32u)) == 0,
          "invalid output tag avoids undefined shift");
  const double zs[]{0, 1e-12, .1, 1, 10, 1000};
  for (double matter : {0., .75}) {
    SoundHorizonRequest source{
        {70, matter, 1 - matter, 0, .1}, 1000, "analytic control"};
    auto batch = evaluate_early_late(source, zs, all, p);
    require(batch.status == S::ok && batch.rows.size() == 6, "bounded rows");
    for (const auto &row : batch.rows) {
      for (const auto &v : row.outputs)
        require(v.status == S::ok && v.value.has_value(),
                "requested output accepted");
      const long double a = 1.L / (1.L + row.redshift), r = 1.L - matter;
      // Direct elementary antiderivative of da/sqrt(r+m*a), rationalized.
      const long double expected =
          299792.458L / 70 * 2 * (1 - a) /
          (std::sqrt(r + matter) + std::sqrt(r + matter * a));
      const auto dm = *row.outputs[3].value;
      require(std::abs(dm - expected) <= 1e-9L + 2e-11L * std::abs(expected),
              "analytic distance budget");
      require(*row.outputs[4].value == 0 ||
                  std::abs(*row.outputs[4].value / (1 + row.redshift) - dm) <=
                      1e-9,
              "duality");
    }
    const auto old = batch;
    source.model.h0_km_s_mpc = 140;
    batch = evaluate_early_late(source, zs, all, p);
    for (unsigned j = 0; j < 6; ++j)
      for (unsigned i = 6; i < 9; ++i)
        require(std::abs(*batch.rows[j].outputs[i].value -
                         *old.rows[j].outputs[i].value) <=
                    1e-11 + 5e-11 * std::abs(*old.rows[j].outputs[i].value),
                "H0 ratio cancellation");
    require(batch.ruler && batch.ruler->callbacks > 0, "one ruler");
    std::size_t callbacks = batch.ruler->callbacks;
    for (auto &row : batch.rows)
      callbacks += row.callbacks;
    require(callbacks == batch.callbacks, "single ruler callback accounting");
  }
  SoundHorizonRequest source{{70, .3, 1e-4, .05, 5e-5}, 1059, "supplied epoch"};
  const double invalid[]{-1, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()};
  auto batch = evaluate_early_late(source, invalid, all, p);
  require(batch.rows[0].outputs[0].status == S::outside_domain &&
              batch.rows[1].outputs[0].status == S::nonfinite_input,
          "invalid redshift");
  source.model.omega_r = 0;
  batch = evaluate_early_late(source, zs, all, p);
  require(batch.callbacks == 0 &&
              batch.rows[0].outputs[0].status == S::outside_domain,
          "zero radiation rejected");
  source.model = {70, .3, 1e-4, .05, 5e-5};
  auto q = p;
  q.maximum_total_callbacks = 0;
  const unsigned eonly = early_late_mask(EarlyLateOutput::e);
  batch = evaluate_early_late(source, zs, eonly, q);
  require(batch.callbacks == 0 && !batch.ruler &&
              batch.rows[1].outputs[0].status == S::ok,
          "state requests no quadrature");
  q = p;
  q.maximum_native_bytes = 0;
  require(evaluate_early_late(source, zs, all, q).status == S::work_limit,
          "preallocation quota");
  q = p;
  q.sound.maximum_total_callbacks = 0;
  batch = evaluate_early_late(source, zs, all, q);
  require(batch.rows[2].outputs[3].status == S::ok &&
              batch.rows[2].outputs[6].status == S::work_limit,
          "distance survives ruler failure");
  source.model.h0_km_s_mpc = std::numeric_limits<double>::min();
  batch = evaluate_early_late(source, zs, all, p);
  require(batch.rows[2].outputs[3].status == S::outside_domain,
          "overflow rejected");
  source.model = {70, 0, 1, 0, 1};
  const double hostile[]{std::numeric_limits<double>::denorm_min(),
                         std::numeric_limits<double>::max()};
  batch = evaluate_early_late(source, hostile, all, p);
  require(batch.rows[0].outputs[3].status == S::outside_domain,
          "positive distance underflow rejected");
  require(batch.rows[1].outputs[0].status == S::outside_domain,
          "highz E overflow rejected");
  for (bool reverse : {false, true}) {
    const double tiny = std::numeric_limits<double>::denorm_min();
    source.model = {70, reverse ? tiny : 1., reverse ? 1. : tiny, 0,
                    reverse ? 1. : tiny};
    batch = evaluate_early_late(source, zs, all, p);
    require(batch.callbacks == 0 &&
                batch.rows[0].outputs[0].status == S::outside_domain,
            "exact closure excess rejected");
  }
  source.model = {70, .3, 1e-4, .31, 5e-5};
  batch = evaluate_early_late(source, zs, all, p);
  require(batch.callbacks == 0 &&
              batch.rows[0].outputs[0].status == S::outside_domain,
          "baryon subset rejected");
  source.model = {70, .3, 1e-4, .05, 2e-4};
  batch = evaluate_early_late(source, zs, all, p);
  require(batch.callbacks == 0 &&
              batch.rows[0].outputs[0].status == S::outside_domain,
          "photon subset rejected");
  source.model = {70, .3, 1e-4, .05, 5e-5};
  q = p;
  q.absolute_tolerance_ratio = 1e-30;
  q.relative_tolerance_ratio = 1e-30;
  batch = evaluate_early_late(source, zs, all, q);
  require(batch.rows[2].outputs[3].status == S::ok &&
              batch.rows[2].outputs[6].status ==
                  S::conditioning_budget_exceeded,
          "unattainable ratio budget separate");
  q = p;
  q.maximum_total_callbacks = 3;
  batch = evaluate_early_late(source, zs, all, q);
  require(batch.callbacks <= 3, "combined callback quota never exceeded");
  q = p;
  q.maximum_points = 0;
  require(evaluate_early_late(source, zs, all, q).rows.empty(),
          "point quota before row allocation");
  require(evaluate_early_late(source, zs, 1u << early_late_output_count, p)
                  .status == S::invalid_input,
          "unknown mask rejected");
  const int mode = std::fegetround();
  std::fesetround(FE_UPWARD);
  batch = evaluate_early_late(source, zs, all, p);
  std::fesetround(mode);
  require(batch.status == S::invalid_input && batch.rows.empty(),
          "unsupported arithmetic mode rejected");
  std::cout << "early_late owner checks=" << checks << '\n';
}
