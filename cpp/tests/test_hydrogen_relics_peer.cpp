#include "hydrogen_relic_reference.hpp"
#ifndef IRRED_REFERENCE_ONLY
#include "irred/recombination_drag.hpp"
#endif
#include <iomanip>
#include <iostream>
#include <string>
namespace {
namespace ref = hydrogen_relic_reference;
namespace old = hydrogen_thermal_reference;
using ref::W;
unsigned checks = 0;
void require(bool ok, const char *what) {
  ++checks;
  if (!ok)
    throw std::runtime_error(what);
}
W budget(unsigned group, W value) {
  constexpr W absolute[]{1e-8L, 2e-7L, 2e-5L, 2e-7L, 2e-9L, 2e-9L, 0};
  constexpr W relative[]{1e-6L, 1e-6L, 2e-7L, 1e-6L, 2e-6L, 3e-6L, 1e-5L};
  return absolute[group] + relative[group] * std::abs(value);
}
std::array<W, 7> coordinates(old::Value a) {
  return {a.x,       a.drag,       a.temperature, a.thomson,
          a.opacity, a.visibility, a.survival};
}
void momentum() {
  ref::MomentumRule<16> coarse;
  ref::MomentumRule<32> fine;
  for (W y : {0.L, 1e-4L, .001L, .003L, .01L, .1L, 1.L, 10.L, 100.L, 2000.L}) {
    W a = coarse.moment(y), b = fine.moment(y), delta = std::abs(a / b - 1),
      tail = ref::tail(y) / b;
    std::cout << "FD y=" << y << " GL16=" << a << " GL32=" << b
              << " refinement_relative=" << delta << " tail_relative=" << tail
              << '\n';
    require(b > 0 && a > 0 && delta <= 1e-11L, "FD16/32 reference refinement");
    require(tail <= 1e-11L, "FD relative omitted-tail envelope");
    if (y == 0)
      require(std::abs(b / (7 * std::pow(std::numbers::pi_v<W>, 4) / 120) -
                       1) <= 1e-14L,
              "independent massless FD moment");
  }
}
struct Case {
  const char *name;
  ref::Model model;
  double initial, late;
  bool evolved = true;
  bool fixed_density_control = false;
};
#ifndef IRRED_REFERENCE_ONLY
void fixed_density_h0(const irred::cosmology::ThermalBackground &background,
                      irred::cosmology::ThermalPhysicalModel model,
                      ref::Model source) {
  using namespace irred::cosmology;
  // Change only H0: all physical densities, temperatures and species stay fixed.
  const W h1 = source.H0 / 100;
  ref::Background reference1(source, 32);
  model.h0_km_s_mpc = 80;
  source.H0 = 80;
  const W h2 = source.H0 / 100, expected = (h2 - h1) * (h2 + h1);
  auto mapping = map_thermal_physical_model(model);
  require(mapping.status == irred::numerics::Status::ok && mapping.model,
          "fixed-density alternate physical mapping");
  auto alternate = prepare_thermal_background(*mapping.model);
  require(alternate.status() == irred::numerics::Status::ok,
          "fixed-density alternate background only");
  std::vector<double> scale;
  for (double z : {0., 1., 100., 300., 600., 1000., 1600.})
    scale.push_back(1 / (1 + z));
  auto first = background.evaluate(scale, thermal_h),
       second = alternate.evaluate(scale, thermal_h);
  require(first.status == irred::numerics::Status::ok &&
              second.status == irred::numerics::Status::ok &&
              first.rows.size() == scale.size() &&
              second.rows.size() == scale.size(),
          "fixed-density ordered H batches");
  ref::Background reference2(source, 32);
  for (std::size_t i = 0; i < scale.size(); ++i) {
    const auto &a = first.rows[i].h_km_s_mpc,
               &b = second.rows[i].h_km_s_mpc;
    require(a.status == irred::numerics::Status::ok && a.value &&
                b.status == irred::numerics::Status::ok && b.value,
            "fixed-density independently available H");
    const W z = 1 / W(scale[i]) - 1,
            ar = reference1.hubble(z) * old::mpc() / 100000,
            br = reference2.hubble(z) * old::mpc() / 100000,
            av = W(*a.value) / 100, bv = W(*b.value) / 100;
    require(std::abs(av / ar - 1) <= 2e-10L &&
                std::abs(bv / br - 1) <= 2e-10L,
            "both fixed-density native H versus original direct SI FD");
    const W actual = (bv - av) * (bv + av),
            reference = (br - ar) * (br + ar),
            allowance = 32 * std::numeric_limits<double>::epsilon() *
                        (av * av + bv * bv),
            reference_allowance = 32 * std::numeric_limits<W>::epsilon() *
                                  (ar * ar + br * br);
    std::cout << "FIXED_DENSITY_H0 z=" << z << " expected=" << expected
              << " native=" << actual << " reference=" << reference
              << " native_error=" << std::abs(actual - expected)
              << " native_allowance=" << allowance
              << " reference_error=" << std::abs(reference - expected)
              << " reference_allowance=" << reference_allowance
              << " condition=" << (av * av + bv * bv) / expected << '\n';
    // Stored binary64 H is ill-conditioned for the tiny Lambda difference at
    // high z. This is an explicit rounding control, not a relative error floor.
    require(actual > 0 && reference > 0 &&
                std::abs(actual - expected) <= allowance &&
                std::abs(reference - expected) <= reference_allowance,
            "positive fixed-density squared H law with declared rounding");
  }
}
#endif
void physical(const Case &test) {
  std::vector<W> z{test.initial, test.initial - .001, test.initial - .01,
                   test.initial - .1};
  for (double value : {1500., 1400., 1300., 1200., 1100., 1000., 900., 800.,
                       700., 600., 500., 400.})
    if (value < test.initial - .1 && value > test.late)
      z.push_back(value);
  z.push_back(test.late);
  ref::Background h16(test.model, 16), h32(test.model, 32);
  W max_h_refinement = 0;
  for (W at : {0.L, 1.L, 100.L, 300.L, 600.L, 1000.L, W(test.initial)}) {
    W delta = std::abs(h16.hubble(at) / h32.hubble(at) - 1);
    max_h_refinement = std::max(max_h_refinement, delta);
    require(delta <= 1e-11L, "direct SI H16/32 refinement");
  }
  require(std::abs(h32.hubble(0) / (test.model.H0 * 1000 / old::mpc()) - 1) <=
              1e-16L,
          "source-defined present closure");
  auto coarse = ref::integrate(test.model, test.initial, test.late, z, 1, 16,
                               test.evolved);
  auto medium = ref::integrate(test.model, test.initial, test.late, z, 2, 16,
                               test.evolved);
  auto fine = ref::integrate(test.model, test.initial, test.late, z, 4, 32,
                             test.evolved);
  const unsigned count = test.evolved ? 7 : 3;
  W maximum_reference = 0, boundary_reference = 0;
  for (std::size_t i = 0; i < z.size(); ++i) {
    auto a = coordinates(medium.query(i)), b = coordinates(fine.query(i));
    for (unsigned group = 0; group < count; ++group) {
      W fraction = std::abs(a[group] - b[group]) / budget(group, b[group]);
      maximum_reference = std::max(maximum_reference, fraction);
      if (i < 4)
        boundary_reference = std::max(boundary_reference, fraction);
      if (!(fraction <= .05L))
        std::cerr << "REFINEMENT case=" << test.name << " z=" << z[i]
                  << " group=" << group << " fraction=" << fraction << '\n';
      require(fraction <= .05L,
              "combined FD/ODE refinement within5% allocation");
    }
    require(b[0] > 0 && b[0] < 1 && b[2] > 0,
            "reference history query positivity");
    if (!test.evolved)
      require(b[2] == test.model.T0 * (1 + z[i]),
              "prescribed temperature is a distinct fixed closure");
  }
  for (const auto *r : {&coarse, &medium, &fine}) {
    require(r->stats.residual <= 8e-17L, "actual converged stage residual");
    require(r->stats.min_x > 0 && r->stats.min_neutral > 0 &&
                r->stats.min_theta > 0,
            "reference stage positivity");
    require(r->stats.maximum_initial_step_stiffness <= .125L * (1 + 1e-15L),
            "initial boundary-layer stiffness resolved");
  }
  require(std::abs(medium.root - fine.root) <= .0001L,
          "reference conditional root refinement");
  W gl_mass_error = 0;
  if (test.evolved) {
    auto quad = ref::opacity_gl(fine);
    require(std::abs(quad.thomson - fine.nodes.back().thomson) <=
                .05L * budget(3, quad.thomson),
            "reference independent GL Thomson depth");
    require(std::abs(quad.drag - fine.nodes.back().drag) <=
                .05L * budget(1, quad.drag),
            "reference independent GL drag depth");
    require(std::abs(quad.survival - fine.query(0).survival) <=
                .05L * budget(6, quad.survival),
            "reference independent GL survivor");
    gl_mass_error = std::abs(quad.mass + quad.survival - 1);
    require(gl_mass_error <= 1e-8L,
            "reference GL finite visible mass plus survivor");
  }
  std::cout << "REFERENCE " << test.name << " steps=" << coarse.stats.steps
            << ',' << medium.stats.steps << ',' << fine.stats.steps
            << " rhs=" << fine.stats.rhs << " Newton=" << fine.stats.newton
            << " residual=" << fine.stats.residual
            << " minimum_stage_x=" << fine.stats.min_x
            << " minimum_stage_neutral=" << fine.stats.min_neutral
            << " minimum_stage_theta=" << fine.stats.min_theta
            << " initial_hA=" << fine.stats.maximum_initial_step_stiffness
            << " H_refinement=" << max_h_refinement
            << " max_refinement_fraction=" << maximum_reference
            << " boundary_fraction=" << boundary_reference
            << " root_delta=" << std::abs(medium.root - fine.root)
            << " GL_mass_error=" << gl_mass_error
            << " H_queries=" << fine.physics.background.h_queries
            << " momentum_terms=" << fine.physics.background.momentum_terms
            << '\n';
#ifndef IRRED_REFERENCE_ONLY
  using namespace irred::cosmology;
  ThermalPhysicalModel model{
      double(test.model.H0), double(test.model.baryon), double(test.model.cdm),
      double(test.model.T0), double(test.model.other),  {}};
  for (auto a : test.model.species)
    model.species.push_back(
        {double(a.mass), double(a.temperature), double(a.weight)});
  PureHydrogenPolicy policy;
  policy.maximum_total_work = 700000000;
  auto owner = prepare_pure_hydrogen_history(
      {model, test.initial, test.late,
       test.evolved ? HydrogenTemperatureModel::evolved_compton_adiabatic
                    : HydrogenTemperatureModel::prescribed_radiation},
      policy);
  require(owner.status() == irred::numerics::Status::ok,
          "native explicitly budgeted positive-mass history");
  require(owner.work().total() <= policy.maximum_total_work &&
              owner.work().momentum_callbacks > 0,
          "actual native massive momentum work");
  const auto preparation_work = owner.work().total();
  std::vector<double> redshifts;
  for (W value : z)
    redshifts.push_back(double(value));
  auto native = owner.evaluate(redshifts, test.evolved ? 127 : 7);
  require(native.status == irred::numerics::Status::ok &&
              native.rows.size() == z.size(),
          "ordered native history query");
  W maximum_native = 0;
  for (std::size_t i = 0; i < z.size(); ++i) {
    auto reference = coordinates(fine.query(i));
    const auto &row = native.rows[i];
    const HydrogenHistoryValue *groups[]{&row.electron_fraction,
                                         &row.drag_depth,
                                         &row.matter_temperature_kelvin,
                                         &row.thomson_depth,
                                         &row.thomson_opacity_per_redshift,
                                         &row.visibility_per_redshift,
                                         &row.finite_endpoint_survival};
    for (unsigned group = 0; group < count; ++group) {
      require(groups[group]->status == irred::numerics::Status::ok &&
                  groups[group]->value.has_value(),
              "native requested group independently available");
      W fraction = std::abs(*groups[group]->value - reference[group]) /
                   budget(group, reference[group]);
      maximum_native = std::max(maximum_native, fraction);
      if (!(fraction <= 1))
        std::cerr << "NATIVE case=" << test.name << " z=" << z[i]
                  << " group=" << group << " fraction=" << fraction << '\n';
      require(fraction <= 1, "native independent massive FD Radau comparison");
    }
  }
  auto root = owner.conditional_unit_depth_redshift();
  require(root.status == irred::numerics::Status::ok && root.value &&
              std::abs(*root.value - fine.root) <= .002L,
          "native conditional root independent comparison");
  std::vector<double> scale;
  for (double value : {0., 1., 100., 300., 600., 1000., test.initial})
    scale.push_back(1 / (1 + value));
  auto H = owner.background()->evaluate(scale, thermal_h);
  W maximum_H = 0;
  require(H.status == irred::numerics::Status::ok &&
              H.rows.size() == scale.size(),
          "native retained-background H batch");
  for (std::size_t i = 0; i < scale.size(); ++i) {
    require(H.rows[i].h_km_s_mpc.status == irred::numerics::Status::ok &&
                H.rows[i].h_km_s_mpc.value.has_value(),
            "native retained-background H available");
    W expected = fine.physics.hubble(1 / W(scale[i]) - 1),
      actual = *H.rows[i].h_km_s_mpc.value * 1000 / old::mpc();
    W error = std::abs(actual / expected - 1);
    maximum_H = std::max(maximum_H, error);
    require(error <= 2e-10L, "independent direct-SI massive H");
  }
  std::cout << "NATIVE " << test.name
            << " maximum_allocation_fraction=" << maximum_native
            << " maximum_H_relative=" << maximum_H
            << " work=" << owner.work().total()
            << " momentum=" << owner.work().momentum_callbacks << '\n';
  if (test.fixed_density_control)
    fixed_density_h0(*owner.background(), model, test.model);
  if (test.evolved) {
    old::Gauss8 gauss;
    W survivor = *native.rows.front().finite_endpoint_survival.value,
      previous = 0;
    for (unsigned panels : {128u, 256u, 512u}) {
      W mass = 0, width = (test.initial - test.late) / W(panels);
      std::vector<double> points;
      std::vector<W> weights;
      points.reserve(1024);
      weights.reserve(1024);
      auto consume = [&] {
        auto values = owner.evaluate(points, hydrogen_visibility);
        require(values.status == irred::numerics::Status::ok &&
                    values.rows.size() == points.size(),
                "bounded native visibility query");
        for (std::size_t i = 0; i < points.size(); ++i) {
          const auto &g = values.rows[i].visibility_per_redshift;
          require(g.status == irred::numerics::Status::ok && g.value &&
                      *g.value > 0,
                  "native quadrature positive visibility");
          mass += weights[i] * *g.value;
        }
        points.clear();
        weights.clear();
      };
      for (unsigned i = 0; i < panels; ++i)
        for (unsigned j = 0; j < 8; ++j) {
          points.push_back(
              double(test.late + (i + .5L) * width + width / 2 * gauss.x[j]));
          weights.push_back(width / 2 * gauss.w[j]);
          if (points.size() == 1024)
            consume();
        }
      if (!points.empty())
        consume();
      W error = std::abs(mass + survivor - 1),
        refinement = panels == 128 ? 0 : std::abs(mass - previous);
      std::cout << "NATIVE_MASS " << test.name << " panels=" << panels
                << " mass=" << mass << " survivor=" << survivor
                << " conservation_error=" << error
                << " refinement=" << refinement << '\n';
      require(error <= 2e-7L,
              "actual native finite visible mass plus surviving mass");
      if (panels == 512)
        require(refinement <= 1e-8L,
                "actual native mass quadrature refinement");
      previous = mass;
    }
  }
  require(owner.work().total() == preparation_work,
          "ordered queries retain prepared background/evolution work");
#endif
}
} // namespace
int main() {
  try {
    std::cout << std::setprecision(21) << std::unitbuf;
    std::cerr << std::setprecision(21);
    momentum();
    ref::Model fid{67.4, .02237,
                   .12,  2.7255,
                   0,    {{.06, 1.95, 2}, {0, 1.95, 2}, {0, 1.95, 2}}};
    ref::Model high{80, .03, .15, 2.75, 0, {{.3, 2, 2}, {0, 2, 2}, {0, 2, 2}}};
    ref::Model low{60, .015, .08, 2.70, 0, {{.001, 1.8, 1}, {0, 1.8, 1}}};
    physical({"fiducial-.06-two-massless", fid, 1600, 300, true, true});
    physical({"heavy-.3-high-physical", high, 1600, 300});
    physical({"light-.001-low-physical", low, 1550, 600});
    physical({"light-.001-prescribed1650", low, 1650, 300, false});
    std::cout << "PASS " << checks
              << " massive-relic hydrogen peer controls; original "
                 "SI/GL16-32/Radau, common constants and long-double/libm\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL after " << checks << " controls: " << e.what() << '\n';
    return 1;
  }
}
