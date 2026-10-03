// Synthetic supplied spectrum/window mean; no observational likelihood.
#include <irred/windowed_linear_power.hpp>
#include <iostream>
#include <utility>

int main() {
  namespace p = irred::windowed_linear_power;
  p::Source source;
  source.identity = {"synthetic/affine-P-at-z0.32/v1",
                     "synthetic/signed-cross-mixing-window/v1",
                     "supplied-matter-contrast",
                     "supplied-scale-independent-linear-f",
                     "three-columns/two-P0-then-two-P2-rows",
                     "fixed-synthetic-operator/not-release-calibration"};
  source.spectrum = {.32, {.01, .03, .06, .10, .20},
                         {1050, 1150, 1300, 1500, 2000}};
  source.window.coordinates = p::Coordinates::fixed_h_reference;
  source.window.h_reference = .6711;
  source.window.theory_k = {.04, .08, .12};
  source.window.output0_k = {.06, .10};
  source.window.output2_k = {.06, .10};
  source.window.output0_ids = {10, 11};
  source.window.output2_ids = {20, 21};
  source.window.w00 = {.5, .5, 0, 0, .25, .75};
  source.window.w02 = {.02, -.01, 0, 0, .01, -.02};
  source.window.w20 = {.03, -.02, 0, 0, -.01, .02};
  source.window.w22 = {.5, .5, 0, 0, .25, .75};
  constexpr std::size_t payload = 16*1024*1024;
  auto owner = p::prepare(std::move(source), {2000000, payload});
  if (owner.status() != p::Status::ok) return 1;
  const p::ModelPoint model{.32, 2, .8, 1.03, .97};
  const p::EvaluationPolicy policy{1e-7, 1e-6, 8000000, 2000000,
                                    4096, 20, payload, false};
  const auto result = owner.evaluate(model, policy);
  if (result.status != p::Status::ok || result.mean.size() != 4) return 2;
  for (const auto& row : result.mean)
    std::cout << row.multipole << ' ' << row.source_row_id << ' '
              << row.source_k << ' ' << row.power.value << ' '
              << row.power.combined_estimate << '\n';
  // Means and diagnostics are in (Mpc/h_ref)^3; this is no observed score.
  return 0;
}
