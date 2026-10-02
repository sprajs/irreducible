// Frozen original direct-SI RK4 reference precedes owner implementation.
// Constants ancestry is explicitly shared CODATA2022/ASD5.12; equations
// describe the same conditional model, algorithms/implementation are
// independent.
#include "irred/recombination_drag.hpp"
#include "recombination_drag_reference.hpp"
#include <algorithm>
#include <cfenv>
#include <cstdint>
#include <iostream>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace rr = recombination_reference;
namespace {
unsigned checks = 0;
long double maximum_x_fraction = 0, maximum_tau_fraction = 0,
            maximum_refinement_fraction = 0;
void need(bool b, const char *label) {
  ++checks;
  if (!b)
    throw std::runtime_error(label);
}
void near(long double got, long double ref, long double allowance,
          const char *label) {
  if (std::abs(got - ref) > allowance)
    std::cerr << label << " got=" << (double)got << " ref=" << (double)ref
              << " difference=" << (double)std::abs(got - ref)
              << " allowance=" << (double)allowance << '\n';
  need(std::abs(got - ref) <= allowance, label);
}
long double scalar(const HydrogenHistoryValue &v) {
  need(v.status == S::ok && v.value, "accepted history coordinate");
  return *v.value;
}
PureHydrogenRequest request(const rr::Model &m = {}) {
  return {{(double)m.H0,
           (double)m.b,
           (double)m.cdm,
           (double)m.T0,
           (double)m.other,
           {}},
          1600,
          300};
}
} // namespace
int main() {
  try {
    PureHydrogenHistory retained;
    for (rr::Model model : {rr::Model{}, rr::Model{60, .018, .08, 2.7, 0},
                            rr::Model{80, .029, .149, 2.75, 3e-5}}) {
      auto coarse = rr::integrate(model, 1600, 300, .01),
           fine = rr::integrate(model, 1600, 300, .005);
      need(fine.stability < .2L,
           "resolved explicit reference stiffness margin");
      auto owner = prepare_pure_hydrogen_history(request(model));
      if (owner.status() != S::ok)
        std::cerr << "preparation status=" << (int)owner.status()
                  << " work=" << owner.work().total()
                  << " rhs=" << owner.work().rhs_evaluations << "\n";
      need(owner.status() == S::ok && owner.source() && owner.background() &&
               owner.background()->status() == S::ok,
           "owned physical thermal history");
      double zs[]{1500, 1400, 1300, 1200, 1100, 1000, 900, 800, 600, 300};
      auto predicted = owner.evaluate(zs, 3);
      need(predicted.status == S::ok && predicted.rows.size() == 10,
           "ordered requested history");
      for (unsigned i = 0; i < 10; ++i) {
        auto expected = rr::sample(fine, zs[i]),
             old = rr::sample(coarse, zs[i]);
        long double xa = 1e-8L + 1e-6L * expected.x,
                    ta = 2e-7L + 1e-6L * expected.tau;
        maximum_refinement_fraction = std::max(
            {maximum_refinement_fraction, std::abs(expected.x - old.x) / xa,
             std::abs(expected.tau - old.tau) / ta});
        near(expected.x, old.x, .05L * xa,
             "independent xe reference refinement");
        near(expected.tau, old.tau, .05L * ta,
             "independent drag reference refinement");
        auto &row = predicted.rows[i];
        long double x = scalar(row.electron_fraction),
                    tau = scalar(row.drag_depth);
        maximum_x_fraction =
            std::max(maximum_x_fraction, std::abs(x - expected.x) / xa);
        maximum_tau_fraction =
            std::max(maximum_tau_fraction, std::abs(tau - expected.tau) / ta);
        near(x, expected.x, xa, "independent non-equilibrium xe");
        near(tau, expected.tau, ta, "independent conditional optical depth");
        need(x > 0 && x < 1 && tau >= 0, "history positivity");
        near(x, expected.x,
             row.electron_fraction.absolute_error_estimate + .05L * xa,
             "xe empirical diagnostic covers named difference");
        near(tau, expected.tau,
             row.drag_depth.absolute_error_estimate + .05L * ta,
             "tau empirical diagnostic covers named difference");
        if (i > 0)
          need(x < scalar(predicted.rows[i - 1].electron_fraction) &&
                   tau < scalar(predicted.rows[i - 1].drag_depth),
               "ordered physical monotonic control");
      }
      near(scalar(owner.conditional_unit_depth_redshift()), fine.root, .002,
           "conditional positive root");
      near(coarse.root, fine.root, .05L * .002, "root reference refinement");
      need(owner.work().background_evaluations == 32769 &&
               owner.work().equilibrium_solves == 1 &&
               owner.work().rhs_evaluations > 0 &&
               owner.work().momentum_callbacks == 0 &&
               owner.work().total() <= 4000000,
           "complete preparation work counts");
      auto points = owner.evaluate(std::array<double, 1>{1000}, 1);
      need(scalar(points.rows[0].electron_fraction) >
               rr::saha(model, 1000) + 1e-8L +
                   1e-6L * scalar(points.rows[0].electron_fraction),
           "nonequilibrium differs from Saha after boundary");
      // Original ODE finite-difference residual with a separate directSI rate.
      for (double z : {1400., 1200., 1000., 800.}) {
        double neighbors[]{z - 1, z + 1, z - .5, z + .5, z};
        auto values = owner.evaluate(neighbors, 1);
        long double d1 = (scalar(values.rows[1].electron_fraction) -
                          scalar(values.rows[0].electron_fraction)) /
                         2,
                    d2 = scalar(values.rows[3].electron_fraction) -
                         scalar(values.rows[2].electron_fraction);
        long double expected =
            rr::rhs(rr::coefficients(model, z),
                    {scalar(values.rows[4].electron_fraction), 0})[0];
        near(d2, expected, 2e-6, "independent ODE residual");
        near(d1, d2, 2e-6, "ODE residual interval refinement");
      }
      retained = std::move(owner);
    }
    // Ground-state detailed balance, from direct SI equations; no solver
    // readback.
    for (long double z : {1500.L, 1600.L}) {
      rr::Model m;
      auto q = rr::coefficients(m, z);
      auto x = rr::saha(m, z);
      near(q.D * x * x, q.E * (1 - x), 2e-12L * std::abs(q.D * x * x),
           "atomic detailed balance limit");
    }
    // An explicit massless FD population has an independent radiation density.
    auto massless_source = request();
    massless_source.model.species.push_back({0, 1.95, 2});
    rr::Model radiation;
    radiation.other += rr::photon_physical(radiation) * 7.L / 8 *
                       std::pow(1.95L / radiation.T0, 4);
    auto radiation_old = rr::integrate(radiation, 1600, 300, .01),
         radiation_fine = rr::integrate(radiation, 1600, 300, .005);
    auto massless = prepare_pure_hydrogen_history(massless_source);
    need(massless.status() == S::ok && massless.work().momentum_callbacks == 0,
         "massless FD retained analytic background");
    double radiation_z[]{1200, 1000};
    auto radiation_rows = massless.evaluate(radiation_z, 3);
    for (unsigned i = 0; i < 2; ++i) {
      auto r = rr::sample(radiation_fine, radiation_z[i]),
           older = rr::sample(radiation_old, radiation_z[i]);
      auto xa = 1e-8L + 1e-6L * r.x, ta = 2e-7L + 1e-6L * r.tau;
      near(r.x, older.x, .05L * xa, "massless reference xe refinement");
      near(r.tau, older.tau, .05L * ta, "massless reference tau refinement");
      near(scalar(radiation_rows.rows[i].electron_fraction), r.x, xa,
           "independent massless FD evolution");
      near(scalar(radiation_rows.rows[i].drag_depth), r.tau, ta,
           "independent massless FD drag");
    }
    // Declared endpoints remain exact for nonbinary mesh spacing.
    for (std::size_t base : {8191u, 8003u}) {
      auto endpoint_source = request();
      endpoint_source.initial_redshift = 1600.1;
      endpoint_source.late_redshift = 300.3;
      PureHydrogenPolicy endpoint_policy;
      endpoint_policy.base_intervals = base;
      auto endpoint =
          prepare_pure_hydrogen_history(endpoint_source, endpoint_policy);
      double endpoint_z[]{300.3, 1600.1, 300.3};
      auto e = endpoint.evaluate(endpoint_z, 3);
      need(endpoint.status() == S::ok && e.rows.size() == 3,
           "nonbinary mesh domain");
      need(scalar(e.rows[0].drag_depth) == 0 &&
               scalar(e.rows[2].drag_depth) == 0,
           "exact late endpoint optical depth");
      near(scalar(e.rows[1].electron_fraction), rr::saha({}, 1600.1L),
           2e-12L * rr::saha({}, 1600.1L), "independent initial Saha boundary");
    }
    // Frozen baseline reference facts, independently obtained before
    // production.
    rr::Model baseline;
    auto baseline_owner = prepare_pure_hydrogen_history(request());
    double fixed_z[]{1200, 1000, 300};
    auto fixed = baseline_owner.evaluate(fixed_z, 3);
    near(scalar(fixed.rows[0].electron_fraction), .290840819123972106879L,
         1e-8L + 1e-6L * .290840819123972106879L, "frozen original xe1200");
    near(scalar(fixed.rows[1].drag_depth), .472937715689268879660L,
         2e-7L + 1e-6L * .472937715689268879660L, "frozen original tau1000");
    near(scalar(baseline_owner.conditional_unit_depth_redshift()),
         1053.40958016048116241L, .002, "frozen original conditional root");
    for (double initial : {1550., 1650.}) {
      auto source = request();
      source.initial_redshift = initial;
      auto changed = prepare_pure_hydrogen_history(source);
      double zs[]{1400, 1200, 1000, 300};
      auto values = changed.evaluate(zs, 1);
      auto independent = rr::integrate(baseline, initial, 300, .005);
      for (unsigned i = 0; i < 4; ++i) {
        auto reference = rr::sample(independent, zs[i]);
        near(scalar(values.rows[i].electron_fraction), reference.x,
             1e-8L + 1e-6L * reference.x,
             "initial boundary conditional evolution");
        auto standard =
            rr::sample(rr::integrate(baseline, 1600, 300, .01), zs[i]);
        near(reference.x, standard.x, .05L * (1e-8L + 1e-6L * reference.x),
             "initial boundary reference sensitivity gate");
      }
    }
    auto late = request();
    late.late_redshift = 600;
    auto truncated = prepare_pure_hydrogen_history(late);
    auto late_ref = rr::integrate(baseline, 1600, 600, .005);
    near(scalar(truncated.conditional_unit_depth_redshift()), late_ref.root,
         .002, "endpoint defined conditional root");
    need(std::abs(scalar(truncated.conditional_unit_depth_redshift()) -
                  scalar(baseline_owner.conditional_unit_depth_redshift())) >
             .1,
         "late endpoint dependence retained");
    PureHydrogenPolicy p;
    p.base_intervals = 16384;
    auto refined =
        prepare_pure_hydrogen_history(request(), p).evaluate(fixed_z, 3);
    for (unsigned i = 0; i < 3; ++i) {
      near(scalar(refined.rows[i].electron_fraction),
           scalar(fixed.rows[i].electron_fraction),
           1e-8 + 1e-6 * scalar(fixed.rows[i].electron_fraction),
           "owner mesh refinement xe");
      near(scalar(refined.rows[i].drag_depth), scalar(fixed.rows[i].drag_depth),
           2e-7 + 1e-6 * scalar(fixed.rows[i].drag_depth),
           "owner mesh refinement tau");
    }
    double ordered_z[]{1000, 1200, 1000};
    auto ordered = baseline_owner.evaluate(ordered_z, 3);
    need(scalar(ordered.rows[0].electron_fraction) ==
                 scalar(ordered.rows[2].electron_fraction) &&
             scalar(ordered.rows[1].electron_fraction) ==
                 scalar(fixed.rows[0].electron_fraction),
         "query axis order preserved");
    double mixed_z[]{1000, -1, std::numeric_limits<double>::quiet_NaN(), 1700,
                     1000};
    auto mixed = baseline_owner.evaluate(mixed_z, 3);
    need(mixed.rows[0].electron_fraction.value &&
             !mixed.rows[1].electron_fraction.value &&
             mixed.rows[2].electron_fraction.status == S::nonfinite_input &&
             !mixed.rows[3].electron_fraction.value &&
             mixed.rows[4].electron_fraction.value,
         "mixed history failures retain valid rows");
    for (unsigned mask : {1u, 2u}) {
      auto subset = baseline_owner.evaluate(std::span(fixed_z, 1), mask);
      need(bool(subset.rows[0].electron_fraction.value) == bool(mask & 1) &&
               bool(subset.rows[0].drag_depth.value) == bool(mask & 2),
           "requested history groups only");
    }
    need(baseline_owner.evaluate(fixed_z, 0).rows.empty() &&
             baseline_owner.evaluate(fixed_z, 128).rows.empty(),
         "mask validation");
    need(baseline_owner.evaluate(fixed_z, 3, 0).rows.empty() &&
             baseline_owner.evaluate(fixed_z, 3, 4096, 1).rows.empty(),
         "query resource admission");
    need(!pure_hydrogen_payload_bound(SIZE_MAX, 0, 0) &&
             !pure_hydrogen_payload_bound(1, SIZE_MAX, 0),
         "payload overflow");
    p = {};
    p.maximum_total_work = 0;
    need(prepare_pure_hydrogen_history(request(), p).status() == S::work_limit,
         "zero combined work refusal");
    p = {};
    p.maximum_total_work = 40000;
    auto exhausted = prepare_pure_hydrogen_history(request(), p);
    need(exhausted.status() == S::work_limit &&
             exhausted.work().total() <= 40000,
         "nonlinear work ceiling exact");
    p = {};
    p.maximum_native_bytes = 1;
    need(prepare_pure_hydrogen_history(request(), p).status() == S::work_limit,
         "preparation payload preflight");
    p = {};
    p.base_intervals = SIZE_MAX;
    need(prepare_pure_hydrogen_history(request(), p).status() == S::work_limit,
         "grid multiplication guard");
    p = {};
    p.maximum_output_points = 1;
    auto bounded = prepare_pure_hydrogen_history(request(), p);
    need(bounded.evaluate(fixed_z, 3).rows.empty(), "retained query ceiling");
    p = {};
    p.maximum_native_bytes = *pure_hydrogen_payload_bound(32768, 0, 0);
    auto retained_cap = prepare_pure_hydrogen_history(request(), p);
    need(retained_cap.status() == S::ok &&
             retained_cap.evaluate(fixed_z, 3).rows.empty(),
         "retained payload ceiling applies to output");
    p = {};
    p.absolute_root_tolerance = 1e-30;
    auto root_refusal = prepare_pure_hydrogen_history(request(), p);
    need(root_refusal.status() == S::ok &&
             !root_refusal.conditional_unit_depth_redshift().value &&
             root_refusal.evaluate(fixed_z, 3).rows[0].electron_fraction.value,
         "conditional root refusal preserves history");
    p = {};
    p.absolute_x_tolerance = 1e-30;
    p.relative_x_tolerance = 1e-30;
    auto tight = prepare_pure_hydrogen_history(request(), p)
                     .evaluate(std::span(fixed_z, 1), 3);
    need(!tight.rows[0].electron_fraction.value &&
             tight.rows[0].drag_depth.value,
         "group quality refusal preserves tau");
    for (unsigned choice = 0; choice < 6; ++choice) {
      auto invalid = request();
      if (choice == 0)
        invalid.model.species.push_back({.31, 1.9, 2});
      if (choice == 1)
        invalid.late_redshift = 0;
      if (choice == 2)
        invalid.initial_redshift = 2000;
      if (choice == 3)
        invalid.model.physical_baryon_density = 0;
      if (choice == 4)
        invalid.model.tcmb_kelvin = 4;
      if (choice == 5)
        invalid.model.h0_km_s_mpc = std::numeric_limits<double>::infinity();
      need(prepare_pure_hydrogen_history(invalid).status() != S::ok,
           "unsupported source closure or domain");
    }
    auto copy = baseline_owner;
    auto moved = std::move(baseline_owner);
    auto *self = &moved;
    moved = std::move(*self);
    need(copy.evaluate(fixed_z, 1).rows[0].electron_fraction.value &&
             moved.evaluate(fixed_z, 1).rows[0].electron_fraction.value &&
             baseline_owner.evaluate(fixed_z, 1).rows.empty(),
         "owned source copy/move lifetime");
    need(prepare_pure_hydrogen_history(request()).source()->model.h0_km_s_mpc ==
             67.4,
         "temporary source retained");
    int rounding = std::fegetround();
    need(std::fesetround(FE_DOWNWARD) == 0, "set hostile rounding");
    auto rejected = moved.evaluate(fixed_z, 1);
    need(std::fesetround(rounding) == 0, "restore rounding");
    need(rejected.rows.empty(), "arithmetic refusal");
    std::cout << "PASS " << checks
              << " pureH conditional history controls; maximum xe fraction="
              << (double)maximum_x_fraction
              << " tau fraction=" << (double)maximum_tau_fraction
              << " reference fraction=" << (double)maximum_refinement_fraction
              << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
