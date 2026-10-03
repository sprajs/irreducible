// Original synthetic controls for the separately identified massless-FD
// transfer. MB47/48/23 are shared source ancestry, not independent cosmology.
// Closed characteristic/Bessel and constraint solutions below are mathematical
// witnesses. Numerical admission, installed-consumer acceptance and observation
// are separate gates. Source-only creation did not execute these controls.
#include "massless_fd_peer.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using W = long double;
using C = std::complex<W>;
using S = irred::numerics::Status;
namespace peer = massless_fd_peer;
unsigned checks = 0;
constexpr C imaginary{0, 1};
constexpr W arithmetic = 4096 * std::numeric_limits<W>::epsilon();

void need(bool condition, const char *name) {
  ++checks;
  if (!condition)
    throw std::runtime_error(name);
}

void near(C actual, C expected, W absolute, W relative, const char *name) {
  if (!(std::abs(actual - expected) <=
        absolute + relative * std::abs(expected))) {
    std::cerr << name << ": " << actual << " vs " << expected << '\n';
    need(false, name);
  }
  ++checks;
}

void near(W actual, W expected, W absolute, W relative, const char *name) {
  near(C(actual), C(expected), absolute, relative, name);
}

W legendre(unsigned ell, W mu) {
  switch (ell) {
  case 0:
    return 1;
  case 1:
    return mu;
  case 2:
    return (3 * mu * mu - 1) / 2;
  default:
    throw std::runtime_error(
        "control only defines original low Legendre modes");
  }
}

template <class F> C mean(const peer::Rule<W> &rule, F function) {
  C sum{};
  for (std::size_t j = 0; j < rule.mu.size(); ++j)
    sum += rule.weight[j] * function(j, rule.mu[j]);
  return sum / rule.normalization;
}

// Independent elementary low-order spherical Bessel formulae. Small-phase
// expansions avoid a cancelled oracle; their first omitted terms are bounded
// below the arithmetic allocation for |phase|<.01.
W bessel(unsigned ell, W phase) {
  const W q = phase * phase;
  if (std::abs(phase) < .01L) {
    if (ell == 0)
      return 1 - q / 6 + q * q / 120 - q * q * q / 5040;
    if (ell == 1)
      return phase * (1.L / 3 - q / 30 + q * q / 840 - q * q * q / 45360);
    return q * (1.L / 15 - q / 210 + q * q / 7560 - q * q * q / 498960);
  }
  if (ell == 0)
    return std::sin(phase) / phase;
  if (ell == 1)
    return (std::sin(phase) / phase - std::cos(phase)) / phase;
  return ((3 / q - 1) * std::sin(phase) - 3 * std::cos(phase) / phase) / phase;
}

C characteristic_remainder(W t) {
  // exp(-it)=1-it+remainder. The even part is never obtained as cos(t)-1.
  const W half_sine = std::sin(t / 2), q = t * t;
  const W odd =
      std::abs(t) < .01L
          ? t * q * (1.L / 6 - q / 120 + q * q / 5040 - q * q * q / 362880)
          : t - std::sin(t);
  return {-2 * half_sine * half_sine, odd};
}

peer::Rule<W> discrete_witness_rule() {
  // A positive paired synthetic quadrature deliberately has m2 != 1/3.
  // It exercises declared discretization forcing, not solve_case admission.
  peer::Rule<W> rule;
  rule.mu = {-.8L, -.3L, .3L, .8L};
  rule.weight = {.25L, .25L, .25L, .25L};
  peer::inspect_rule(rule);
  need(rule.status == S::ok,
       "positive synthetic paired rule admitted as algebra control");
  near(rule.m2, (.8L * .8L + .3L * .3L) / 2, arithmetic, arithmetic,
       "independent stored second moment");
  near(rule.m4, (std::pow(.8L, 4) + std::pow(.3L, 4)) / 2, arithmetic,
       arithmetic, "independent stored fourth moment");
  return rule;
}

peer::Angular<W> centered_polynomial(const peer::Rule<W> &rule) {
  peer::Angular<W> angular;
  angular.A = 1.3L;
  angular.B = -.4L;
  for (W mu : rule.mu)
    angular.Q.emplace_back(.7L * (mu * mu - rule.m2),
                           -.2L * (mu * mu * mu - rule.m4 / rule.m2 * mu));
  return angular;
}

C brightness(const peer::Angular<W> &angular, std::size_t j, W mu, W k) {
  return angular.A + imaginary * k * mu * angular.B + k * k * angular.Q[j];
}

void characteristic_controls() {
  for (unsigned nodes : {64U, 128U, 256U}) {
    peer::Budget budget;
    const auto rule = peer::make_rule<W>(nodes, budget);
    need(rule.status == S::ok, "original angular rule construction");
    near(mean(rule, [](std::size_t, W) { return C(1); }), C(1), arithmetic,
         arithmetic, "normalized half mean");
    near(rule.m2, 1.L / 3, arithmetic, arithmetic,
         "actual GL low polynomial quadrature allocation");
    near(rule.m4, 1.L / 5, arithmetic, arithmetic,
         "actual GL fourth polynomial quadrature allocation");
    for (W phase : {0.L, 1e-7L, .3L, 2.L, 20.L}) {
      constexpr W amplitude = -2;
      peer::Angular<W> angular;
      angular.A = amplitude;
      angular.B = -amplitude * phase; // k=1: -ik mu eta initial transport
      for (W mu : rule.mu)
        angular.Q.push_back(amplitude * characteristic_remainder(mu * phase));
      std::vector<C> original;
      for (std::size_t j = 0; j < rule.mu.size(); ++j) {
        original.push_back(brightness(angular, j, rule.mu[j], 1));
        near(original.back(),
             amplitude * std::exp(-imaginary * rule.mu[j] * phase),
             arithmetic * std::max(1.L, phase), arithmetic,
             "characteristic phase and low-basis reconstruction");
      }
      const auto translation = peer::recenter(rule, angular, 1.L, &budget);
      near(translation.maximum_node_change, 0,
           arithmetic * std::max(1.L, phase), 0,
           "characteristic centering preserves nodes");
      const auto readout = peer::hybrid(rule, angular, 1.L, 1.L, 0.L, 0.L);
      near(readout.delta, amplitude * bessel(0, phase), 2e-14L, 2e-14L,
           "characteristic monopole spherical Bessel");
      near(4 * readout.theta / 3, amplitude * bessel(1, phase), 2e-14L, 2e-14L,
           "characteristic dipole sign and normalization");
      if (phase == 0)
        near(2 * readout.sigma, 0, 0, 0, "zero-phase quadrupole exact zero");
      else
        near(2 * readout.sigma / (phase * phase),
             amplitude * bessel(2, phase) / (phase * phase), 2e-14L, 2e-14L,
             "characteristic quadrupole without raw O(1) cancellation");
    }
  }
}

void nodal_transport_controls() {
  const auto rule = discrete_witness_rule();
  auto angular = centered_polynomial(rule);
  constexpr W k = .6L, metric_source = .25L;
  peer::Budget budget;
  const auto derivative =
      peer::angular_rhs(rule, angular, k, metric_source, &budget);
  near(mean(rule, [&](std::size_t j, W) { return angular.Q[j]; }), C(0),
       arithmetic, 0, "synthetic discrete centered mean");
  near(mean(rule, [&](std::size_t j, W mu) { return mu * angular.Q[j]; }), C(0),
       arithmetic, 0, "synthetic discrete centered dipole");
  near(mean(rule, [&](std::size_t j, W) { return derivative.Q[j]; }), C(0),
       arithmetic, 0, "actual-m2 RHS preserves discrete mean");
  near(mean(rule, [&](std::size_t j, W mu) { return mu * derivative.Q[j]; }),
       C(0), arithmetic, 0, "actual-m2 RHS preserves discrete dipole");
  for (std::size_t j = 0; j < rule.mu.size(); ++j) {
    const W mu = rule.mu[j];
    near(brightness(derivative, j, mu, k),
         metric_source - imaginary * k * mu * brightness(angular, j, mu, k),
         arithmetic, arithmetic, "original MB48 at each stored node");
  }
  const auto hybrid = peer::hybrid(rule, angular, k, 1.L, -.3L, .8L);
  const auto mQ =
      mean(rule, [&](std::size_t j, W mu) { return mu * mu * angular.Q[j]; });
  const W phi_prime = .02L;
  const W psi_prime = metric_source / 4 - phi_prime;
  const C density_prime = derivative.A - 4 * psi_prime;
  near(density_prime - (-4 * hybrid.theta / 3 + 4 * phi_prime),
       k * k * (rule.m2 - 1.L / 3) * angular.B, arithmetic, arithmetic,
       "declared actual-m2 continuum continuity defect");
  const C theta_prime = -k * k * derivative.B / 4.L;
  near(theta_prime - k * k * (hybrid.delta / 4 - hybrid.sigma + hybrid.psi),
       k * k * k * k * (1 / rule.m2 - 3) * mQ / 4.L, arithmetic, arithmetic,
       "declared actual-m2 Euler shear defect");
  const C actual_shear_moment_prime = mean(rule, [&](std::size_t j, W mu) {
    return legendre(2, mu) * derivative.Q[j];
  });
  const C expected_shear_moment_prime =
      1.5L * angular.B * (rule.m4 - rule.m2 * rule.m2) -
      imaginary * k * mean(rule, [&](std::size_t j, W mu) {
        return mu * legendre(2, mu) * angular.Q[j];
      });
  near(actual_shear_moment_prime, expected_shear_moment_prime, arithmetic,
       arithmetic, "actual differentiated shear uses stored fourth moment");
}

void translation_controls() {
  const auto rule = discrete_witness_rule();
  auto angular = centered_polynomial(rule);
  constexpr W k = 1e-8L;
  const C alpha{.13L, 0}, beta{0, -.07L};
  for (std::size_t j = 0; j < rule.mu.size(); ++j)
    angular.Q[j] += alpha + beta * rule.mu[j];
  std::vector<C> original;
  for (std::size_t j = 0; j < rule.mu.size(); ++j)
    original.push_back(brightness(angular, j, rule.mu[j], k));
  const auto before = peer::hybrid(rule, angular, k, 1.L, -.3L, .8L);
  peer::Budget budget;
  const auto translation = peer::recenter(rule, angular, k, &budget);
  near(translation.alpha, alpha, arithmetic, arithmetic,
       "actual Gram recovers constant translation");
  near(translation.beta, beta, arithmetic, arithmetic,
       "actual Gram recovers dipole translation");
  for (std::size_t j = 0; j < rule.mu.size(); ++j)
    near(brightness(angular, j, rule.mu[j], k), original[j], arithmetic,
         arithmetic, "Gram translation preserves each original D node");
  near(mean(rule, [&](std::size_t j, W) { return angular.Q[j]; }), C(0),
       arithmetic, 0, "Gram-centered residual mean");
  near(mean(rule, [&](std::size_t j, W mu) { return mu * angular.Q[j]; }), C(0),
       arithmetic, 0, "Gram-centered residual dipole");
  const auto after = peer::hybrid(rule, angular, k, 1.L, -.3L, .8L);
  const C p2_mean =
      mean(rule, [](std::size_t, W mu) { return C(legendre(2, mu)); });
  const C mu_p2_mean =
      mean(rule, [](std::size_t, W mu) { return C(mu * legendre(2, mu)); });
  near((after.sigma - before.sigma) / (k * k),
       (alpha * p2_mean + beta * mu_p2_mean) / 2.L, arithmetic, arithmetic,
       "node preservation does not silently erase hybrid shear translation");
  near((after.theta - before.theta) / (k * k * k), imaginary * beta / 4.L,
       arithmetic / k, arithmetic, "hybrid velocity translation retained");

  // A deliberately redundant represented D=0 must be translated, not merely
  // have its Q mean/dipole overwritten. Its canonical low moments vanish.
  peer::Angular<W> redundant;
  redundant.A = -k * k * alpha.real();
  redundant.B = imaginary * k * beta;
  for (W mu : rule.mu)
    redundant.Q.push_back(alpha + beta * mu);
  for (std::size_t j = 0; j < rule.mu.size(); ++j)
    near(brightness(redundant, j, rule.mu[j], k), C(0), arithmetic * k * k, 0,
         "redundant constant/dipole represents zero D");
  peer::recenter(rule, redundant, k, &budget);
  const auto zero = peer::hybrid(rule, redundant, k, 1.L, 0.L, .8L);
  near(zero.delta, 0, arithmetic, 0,
       "redundant D canonical density within absolute arithmetic floor");
  near(zero.theta / (k * k * k), 0, arithmetic, 0,
       "redundant D canonical velocity remains zero");
  near(zero.sigma / (k * k), 0, arithmetic, 0,
       "redundant D canonical shear remains zero");
}

void initializer_controls() {
  const auto rule = discrete_witness_rule();
  constexpr W k = .002L, H = .8L, eta = 1.25L, p = -10.L / 19;
  peer::Epoch epoch;
  epoch.status = S::ok;
  epoch.hcal = H;
  epoch.fr = .9999L - 1e-12L;
  epoch.fc = .0001L;
  epoch.fl = 1e-12L;
  epoch.x2 = (k / H) * (k / H);
  epoch.g = -1 + epoch.fc / 2 + 2 * epoch.fl;
  epoch.closure_defect = epoch.fc + epoch.fr + epoch.fl - 1;
  epoch.background_identity_defect =
      epoch.g - 1 + 1.5L * epoch.fc + 2 * epoch.fr;
  peer::Budget budget;
  const auto initial = peer::initialize(rule, epoch, k, eta, 0.L, &budget);
  need(budget.status == S::ok,
       "synthetic projected initializer control bounded");
  const auto &ideal = initial.ideal;
  const W F = epoch.fc + 4 * epoch.fr / 3;
  near(ideal.phi, 7 * p / 5, arithmetic, arithmetic,
       "ideal source retains original spatial potential");
  near(ideal.delta_c, -3 * p / 2, arithmetic, arithmetic,
       "ideal source retains signed Newtonian CDM density");
  near(ideal.delta_r, -2 * p, arithmetic, arithmetic,
       "ideal source retains radiation leading density");
  near(ideal.Delta, -epoch.x2 * ideal.phi / (1.5L * F), arithmetic * epoch.x2,
       arithmetic, "ideal projected Hamiltonian correction");
  near(ideal.Vc, p / 2 + ideal.Delta / 3, arithmetic, arithmetic,
       "ideal CDM velocity uses actual projected correction");
  near(ideal.Vr, ideal.Vc, 0, 0,
       "ideal projected relative velocity remains zero");
  near(ideal.Z, -ideal.psi + 1.5L * F * ideal.Vc, arithmetic, arithmetic,
       "ideal trace derivative ancestry is original momentum initializer");
  need(initial.Z == ideal.Z,
       "finite Gram/hybrid readout does not reset original projected Z");

  const W q_amplitude = -2 * p * eta * eta / 3;
  const W p2_mean = ((3 * rule.m2) - 1) / 2;
  const W p2_variance = rule.squared_p2 - p2_mean * p2_mean;
  near(initial.translation.alpha, C(q_amplitude * p2_mean), arithmetic,
       arithmetic, "sampled source quadrupole actual mean is translated");
  near(initial.actual.sigma / (k * k), -q_amplitude * p2_variance / 2,
       arithmetic, arithmetic,
       "initialized hybrid quadrupole actual normalization");
  near(initial.actual.psi,
       ideal.phi + 3 * epoch.fr * H * H * q_amplitude * p2_variance, arithmetic,
       arithmetic, "initialized hybrid lapse follows actual shear");
  near(initial.actual.delta, initial.angular.A - 4 * initial.actual.psi,
       arithmetic, arithmetic,
       "initialized hybrid density follows actual lapse");
  near(initial.E,
       ideal.Z + initial.actual.psi + epoch.x2 * ideal.phi / 3 +
           (epoch.fc * ideal.delta_c + epoch.fr * initial.actual.delta) / 2,
       arithmetic, arithmetic, "initial original energy defect is retained");
  near(initial.W,
       ideal.Z + initial.actual.psi -
           1.5L * (epoch.fc * ideal.Vc + 4 * epoch.fr * initial.actual.Vr / 3),
       arithmetic, arithmetic, "initial original momentum defect is retained");
  need(std::abs(initial.E) > 1e-6L && std::abs(initial.W) > 1e-6L,
       "synthetic finite-quadrature initial defects are not projected away");
  for (std::size_t j = 0; j < rule.mu.size(); ++j) {
    const W mu = rule.mu[j];
    const C ideal_D = ideal.delta_r + 4 * ideal.psi +
                      imaginary * k * mu * (-4 * ideal.Vr / H) +
                      k * k * q_amplitude * legendre(2, mu);
    near(brightness(initial.angular, j, mu, k), ideal_D, arithmetic, arithmetic,
         "initial centering preserves original projected D nodes");
  }

  // This pure-radiation equation limit is a mathematical fixture, outside the
  // positive-CDM physical admission. It tests the selected 7/19 identity.
  epoch.fc = epoch.fl = epoch.closure_defect = 0;
  epoch.fr = 1;
  epoch.g = -1;
  epoch.background_identity_defect = 0;
  const auto radiation = peer::initialize(rule, epoch, k, eta, 0.L, &budget);
  near(radiation.ideal.Delta / epoch.x2, 7.L / 19, arithmetic, arithmetic,
       "projected source coefficient 7/19 retained");
  need(std::abs(radiation.ideal.Delta / epoch.x2 - .25L) > .1L,
       "regular growing 1/4 control is not silently substituted");
}

std::filesystem::path evidence_directory() {
  const char *supplied = std::getenv("IRRED_MASSLESS_FD_PEER_EVIDENCE_DIR");
  const auto root = supplied ? std::filesystem::path(supplied)
                             : std::filesystem::temp_directory_path();
  std::filesystem::create_directories(root);
  const auto stamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
  for (unsigned attempt = 0; attempt < 32; ++attempt) {
    const auto path =
        root / ("irred-massless-fd-peer-" + std::to_string(stamp) + "-" +
                std::to_string(attempt));
    if (std::filesystem::create_directory(path))
      return path;
  }
  throw std::runtime_error(
      "cannot create fresh retained peer evidence directory");
}

void resource_controls(const std::filesystem::path &directory) {
  {
    peer::State<W> source, destination;
    source.v = {1, 2, 3, 4, 5, 6, 7};
    source.R = {{.1L, .2L}, {.3L, .4L}, {.5L, .6L}};
    destination.v = {-1, -2, -3, -4, -5, -6, -7};
    destination.R = {{-.1L, -.2L}};
    const auto original = destination;
    peer::Budget refused;
    refused.maximum_scalar_updates = peer::active_scalars(source) - 1;
    need(!peer::copy_state(destination, source, refused) &&
             refused.status == S::work_limit && refused.scalar_updates == 0 &&
             destination.v == original.v && destination.R == original.R,
         "real active-state baseline copy refuses before any destination "
         "mutation");
    peer::Budget exact;
    exact.maximum_scalar_updates = peer::active_scalars(source);
    need(peer::copy_state(destination, source, exact) &&
             exact.scalar_updates == source.v.size() + 2 * source.R.size() &&
             destination.v == source.v && destination.R == source.R,
         "real active-state copy counts complex components at exact cap");
  }
  {
    peer::State<W> state;
    state.v = {-1, -2, -3, -4, -5, -6, -7};
    const auto original = state.v;
    bool invoked = false;
    peer::Budget budget;
    budget.maximum_scalar_updates = 6;
    need(
        !peer::assign_initial_scalars(state, budget,
                                      [&](auto &v) {
                                        invoked = true;
                                        v.fill(0);
                                      }) &&
            !invoked && state.v == original && budget.scalar_updates == 0 &&
            budget.status == S::work_limit,
        "real initializer assignment cap refuses before its writer is invoked");
    peer::Budget exact;
    exact.maximum_scalar_updates = 7;
    need(peer::assign_initial_scalars(
             state, exact, [&](auto &v) { v = {1, 2, 3, 4, 5, 6, 7}; }) &&
             exact.scalar_updates == 7 && state.v[6] == 7,
         "real initializer charges its seven active scalar destinations");
  }
  {
    peer::State<W> state;
    state.R = {{.1L, .2L}};
    const auto original = state.R;
    peer::Budget budget;
    budget.maximum_scalar_updates = 1;
    need(!peer::resize_residual(state, 2, budget) && state.R == original &&
             budget.status == S::work_limit && budget.scalar_updates == 0,
         "real residual resize refuses before initializing an extra complex "
         "node");
    peer::Budget exact;
    exact.maximum_scalar_updates = 2;
    need(peer::resize_residual(state, 2, exact) && exact.scalar_updates == 2 &&
             state.R[0] == original[0] && state.R[1] == C(0),
         "residual resize charges newly initialized real and imaginary "
         "components");
  }
  {
    peer::Budget budget;
    budget.maximum_rhs = 0;
    need(!budget.stage() && budget.status == S::work_limit && budget.rhs == 0 &&
             budget.attempted_stages == 1,
         "failed RHS attempt retains count and refuses without cap overshoot");
    need(
        !budget.stage() && budget.attempted_stages == 2,
        "already-refused budget does not erase later attempted stage identity");
  }
  {
    peer::Budget budget;
    budget.maximum_scalar_updates = 0;
    const auto rule = peer::make_rule<W>(64, budget);
    need(budget.status == S::work_limit && rule.status != S::ok &&
             budget.scalar_updates == 0,
         "quadrature scalar cap refuses without a partial accepted rule");
  }
  {
    peer::Budget budget;
    budget.maximum_payload_bytes = 0;
    const auto rule = peer::make_rule<W>(64, budget);
    need(rule.status != S::ok, "angular payload zero withholds rule");
  }
  {
    peer::Budget budget;
    budget.maximum_ledger_bytes = 3;
    const auto path = directory / "ledger-byte-refusal.txt";
    peer::Ledger ledger(path.string());
    need(!ledger.write("stage\n", budget) && budget.status == S::work_limit &&
             budget.logging_refusals == 1 && budget.ledger_bytes == 0,
         "encoded byte cap refuses before append and retains refusal counter");
    need(
        std::filesystem::file_size(path) == 0,
        "ledger byte refusal preserves existing file without truncated record");
  }
  {
    std::size_t campaign_bytes = 0;
    peer::Budget budget;
    budget.campaign_bytes = &campaign_bytes;
    budget.maximum_campaign_bytes = 3;
    peer::Ledger ledger((directory / "campaign-byte-refusal.txt").string());
    need(!ledger.write("stage\n", budget) && budget.status == S::work_limit &&
             campaign_bytes == 0 && budget.ledger_bytes == 0,
         "campaign evidence byte cap is separate from per-case cap");
  }
  {
    peer::Budget budget;
    peer::Ledger ledger((directory / "record-size-refusal.txt").string());
    need(!ledger.write(std::string(8193, 'x'), budget) &&
             budget.status == S::work_limit && budget.ledger_bytes == 0,
         "encoded record scratch cap refuses before append");
  }
}

// The callback binds the independently assembled peer trace/constraint RHS.
// Exact curves exercise its signs and normalization rather than integrating
// a second copy of the same residual ODE to make an expected fixture.
template <class ResidualRHS>
void constraint_controls(ResidualRHS residual_rhs) {
  constexpr W initial_energy = .03L, initial_momentum = -.07L;
  for (W n : {0.L, .25L, 2.L, 9.L}) {
    const W decaying = std::exp(-n), growing = std::exp(n);
    {
      const W energy = initial_energy * growing;
      const W momentum = initial_momentum * decaying;
      const auto derivative =
          residual_rhs(-1.L, 0.L, 0.L, energy, momentum, 0.L, 0.L);
      near(derivative.first, energy, arithmetic, arithmetic,
           "undamped radiation Hamiltonian growing mode");
      near(derivative.second, -momentum, arithmetic, arithmetic,
           "undamped zero-k momentum decaying mode");
    }
    {
      const W energy = initial_energy * decaying;
      const W momentum = (initial_momentum - 2 * initial_energy * n) * decaying;
      const auto derivative =
          residual_rhs(-1.L, 0.L, 2.L, energy, momentum, 0.L, 0.L);
      near(derivative.first, -initial_energy * decaying, arithmetic, arithmetic,
           "kappa2 zero-k triangular analytic energy");
      near(derivative.second,
           (-initial_momentum - 2 * initial_energy + 2 * initial_energy * n) *
               decaying,
           arithmetic, arithmetic,
           "kappa2 zero-k repeated-rate analytic momentum");
    }
    {
      constexpr W energy_force = .004L, momentum_force = -.006L;
      const W energy =
          energy_force + (initial_energy - energy_force) * decaying;
      const W momentum = initial_momentum * decaying +
                         (momentum_force - 2 * energy_force) * (1 - decaying) -
                         2 * (initial_energy - energy_force) * n * decaying;
      const auto derivative = residual_rhs(-1.L, 0.L, 2.L, energy, momentum,
                                           energy_force, momentum_force);
      near(derivative.first, -(initial_energy - energy_force) * decaying,
           arithmetic, arithmetic, "injected constant energy forcing retained");
      near(derivative.second,
           (-initial_momentum + momentum_force - 2 * energy_force -
            2 * (initial_energy - energy_force) * (1 - n)) *
               decaying,
           arithmetic, arithmetic,
           "injected constant momentum forcing retained");
    }
    {
      // g=0, kappa=2 and x^2=3/8 give a repeated eigenvalue -5/2.
      // exp(A n)=exp(-5n/2)[I+n(A+5I/2)] supplies a finite-k exact control.
      constexpr W x2 = 3.L / 8;
      const W energy_slope = -initial_energy / 2 + initial_momentum / 8;
      const W momentum_slope = -2 * initial_energy + initial_momentum / 2;
      const W common = std::exp(-2.5L * n);
      const W energy = (initial_energy + n * energy_slope) * common;
      const W momentum = (initial_momentum + n * momentum_slope) * common;
      const auto derivative =
          residual_rhs(0.L, std::sqrt(x2), 2.L, energy, momentum, 0.L, 0.L);
      near(derivative.first,
           (energy_slope - 2.5L * (initial_energy + n * energy_slope)) * common,
           arithmetic, arithmetic, "finite-k analytic Hamiltonian coupling");
      near(derivative.second,
           (momentum_slope - 2.5L * (initial_momentum + n * momentum_slope)) *
               common,
           arithmetic, arithmetic, "finite-k analytic momentum coupling");
    }
  }

  // g and x are time-varying coefficient witnesses, not frozen-eigenvalue
  // assertions. Their transformations explicitly include x_N=-g*x.
  for (W g : {-1.L, -.2L, .8L, 1.L}) {
    for (W x : {1e-12L, .2L, 4.L}) {
      for (W kappa : {1.5L, 2.L, 3.L}) {
        constexpr W energy = .012L, momentum = -.023L;
        const auto derivative =
            residual_rhs(g, x, kappa, energy, momentum, 0.L, 0.L);
        const W scale = std::sqrt(3 * kappa);
        const W y = scale * energy / x;
        const W y_prime = scale * (derivative.first + g * energy) / x;
        const W first_energy_prime =
            2 * y * y_prime + 2 * momentum * derivative.second;
        const W first_expected =
            -2 * (g + 1 + kappa) * y * y - 2 * (g + 2) * momentum * momentum;
        near(first_energy_prime, first_expected, arithmetic, arithmetic,
             "nonautonomous y norm includes derivative of x");
        need(first_energy_prime < 0,
             "kappa positive finite-k homogeneous contraction");
        const W t = x * momentum / scale;
        const W t_prime = x * (derivative.second - g * momentum) / scale;
        const W second_energy_prime =
            2 * energy * derivative.first + 2 * t * t_prime;
        const W second_expected =
            -2 * (2 * g + 1 + kappa) * energy * energy - 4 * (g + 1) * t * t;
        near(second_energy_prime, second_expected, arithmetic, arithmetic,
             "regular t norm includes derivative of x");
        need(second_energy_prime <= 0,
             "kappa at least one regular homogeneous norm");
      }
    }
  }
  const auto marginal = residual_rhs(-1.L, 0.L, 1.L, .03L, -.03L, 0.L, 0.L);
  near(marginal.first, 0, arithmetic, 0,
       "kappa1 zero-k energy can remain constant");
  const auto insufficient = residual_rhs(-1.L, 0.L, .5L, .03L, 0.L, 0.L, 0.L);
  need(insufficient.first > 0,
       "positive kappa below one does not cure zero-k energy growth");
}

template <class CDMReadout, class CDMRHS>
void centered_cdm_controls(CDMReadout readout, CDMRHS rhs, W g, W x) {
  constexpr W p = -10.L / 19, phi_star = 7 * p / 5;
  for (const auto state : std::array<std::array<W, 3>, 3>{
           {{{.02L, -.03L, .04L}}, {{-.1L, .2L, -.3L}}, {{0, 1e-30L, 0}}}}) {
    const W h = state[0], w = state[1], e = state[2];
    near(readout(h, w, e), 3 * w + e + 3 * h,
         arithmetic * std::max({std::abs(h), std::abs(w), std::abs(e)}),
         arithmetic, "centered CDM readout retains tiny comoving response");
    if (std::abs(w) > 1e-20L) {
      const W delta_c = 3 + e + 3 * (phi_star + h);
      const W velocity = p / 2 + w;
      near(readout(h, w, e), delta_c + 3 * velocity, arithmetic, arithmetic,
           "centered CDM readout agrees with primitive identity above floor");
    }
    constexpr W v = .01L, z = -.004L;
    const auto derivative = rhs(g, h, v, w, x, z);
    const W velocity = p / 2 + w, psi = p + h - v;
    near(derivative[0], (g - 1) * velocity + psi, arithmetic, arithmetic,
         "CDM centered velocity is primitive conformal Euler transform");
    near(derivative[1], -x * x * velocity, arithmetic, arithmetic,
         "CDM centered density is primitive continuity transform");
    near(derivative[2], z, arithmetic, arithmetic,
         "CDM centered spatial potential retains trace derivative");
  }
}

void eds_equation_controls() {
  // Formal radiation-free EdS equations are a mathematical limit, outside the
  // admitted positive explicit-FD owner. No zero-species solve is substituted.
  const auto rule = discrete_witness_rule();
  constexpr W phi = -3.L / 5, velocity = 2 * phi / 3;
  constexpr W p = -10.L / 19, phi_star = 7 * p / 5;
  peer::Epoch epoch;
  epoch.status = S::ok;
  epoch.hcal = 1;
  epoch.fc = 1;
  epoch.fr = epoch.fl = epoch.closure_defect = 0;
  epoch.g = -.5L;
  epoch.background_identity_defect = 0;
  for (W x : {0.L, .2L, 4.L}) {
    epoch.x2 = x * x;
    peer::State<W> state;
    state.v[1] = phi - phi_star;
    state.v[3] = velocity - p / 2;
    state.v[4] = -2 * epoch.x2 * phi / 3;
    state.R.resize(rule.mu.size());
    const auto readout = peer::readout(state, rule, epoch);
    near(readout.phi, phi, arithmetic, arithmetic,
         "formal EdS constant spatial potential for unit asymptotic zeta");
    near(readout.psi, phi, arithmetic, arithmetic,
         "formal EdS no anisotropic stress lapse");
    near(readout.s, 0, arithmetic, 0,
         "formal EdS radiation-free slip vanishes");
    near(readout.E, 0, arithmetic, 0,
         "formal EdS original Hamiltonian constraint");
    near(readout.W, 0, arithmetic, 0,
         "formal EdS original momentum constraint");
    const W delta_c = -2 * phi - 2 * epoch.x2 * phi / 3;
    near(peer::centered_cdm(state), -2 * epoch.x2 * phi / 3, arithmetic,
         arithmetic, "formal EdS growing comoving CDM identity");
    near(delta_c + 3 * velocity, peer::centered_cdm(state), arithmetic,
         arithmetic, "formal EdS Newtonian and comoving CDM distinction");
    const W h_N = 0, Z_N = 0, psi_N = 0;
    const W w_N = (epoch.g - 1) * state.v[3] + state.v[1] - readout.v +
                  p * (epoch.g + 1) / 2;
    const W e_N = -epoch.x2 * velocity;
    near(w_N, 0, arithmetic, 0,
         "formal EdS centered Euler retains constant growing velocity");
    near(e_N + 3 * h_N, -2 * epoch.x2 * phi / 3, arithmetic, arithmetic,
         "formal EdS continuity includes x2_N=x2");
    near(Z_N + (epoch.g + 2) * state.v[2] + psi_N +
             (2 * epoch.g + 1) * readout.psi + epoch.x2 * readout.s / 3,
         0, arithmetic, 0,
         "formal EdS original Einstein trace has constant potential");
    const auto constraint =
        peer::constraint_rhs(epoch.g, x, 2.L, readout.E, readout.W);
    near(constraint[0], 0, arithmetic, 0,
         "formal EdS damping preserves zero energy constraint");
    near(constraint[1], 0, arithmetic, 0,
         "formal EdS damping preserves zero momentum constraint");
  }
}

void full_trace_controls(const irred::cosmology::ThermalBackground &background,
                         const std::filesystem::path &directory,
                         std::size_t &campaign_bytes) {
  peer::Budget budget;
  budget.campaign_bytes = &campaign_bytes;
  const auto rule = peer::make_rule<W>(64, budget);
  need(rule.status == S::ok, "full trace algebra control angular rule");
  irred::cosmology::ThermalPolicy thermal;
  thermal.momentum_method = background.momentum_method();
  need(budget.charge(budget.background_queries, 1,
                     budget.maximum_background_queries),
       "full trace control radiation query charged");
  const auto radiation = background.scaled_expansion(0, thermal);
  budget.momentum_callbacks += radiation.callbacks;
  need(radiation.status == S::ok && radiation.a4_e2 > 0,
       "full trace control shares retained thermal massless normalization");
  constexpr W scale_factor = 1e-4L;
  auto epoch = peer::shared_epoch(background, scale_factor, 0, radiation.a4_e2,
                                  radiation.error_estimate, 0, budget);
  need(epoch.status == S::ok, "full trace algebra control retained epoch");
  const W wave = .3L * epoch.hcal;
  epoch = peer::shared_epoch(background, scale_factor, wave, radiation.a4_e2,
                             radiation.error_estimate, 0, budget);
  need(epoch.status == S::ok,
       "full trace algebra control retained nonzero phase");
  peer::Ledger ledger(
      (directory / "injected-full-trace-controls.txt").string());
  peer::Controls control;
  control.nodes = 64;
  control.run_id = "synthetic-full-trace-injection";
  peer::Diagnostics diagnostics;

  const W p = -10.L / 19, p2_mean = (3 * rule.m2 - 1) / 2;
  const W p2_variance = rule.squared_p2 - p2_mean * p2_mean;
  auto state_with_slip = [&](W h, W w, W v, W z) {
    peer::State<W> state;
    state.v[0] = std::log(scale_factor);
    state.v[1] = h;
    state.v[2] = z;
    state.v[3] = w;
    const W q2 = -(2 * p / 5 + v) / (3 * epoch.fr);
    const W coefficient = q2 / p2_variance;
    need(budget.scalars(7 + 2 * rule.mu.size()),
         "injected full-state assembly charged");
    for (W mu : rule.mu)
      state.R.emplace_back(coefficient * (legendre(2, mu) - p2_mean) +
                               2 * p * legendre(2, mu) / 3,
                           0);
    std::ostringstream input;
    input << std::setprecision(std::numeric_limits<W>::max_digits10)
          << "control-input a=" << scale_factor << " k=" << wave << " h=" << h
          << " w=" << w << " v=" << v << " Z=" << z
          << " e=0 c=0 beta=0 U_basis=mean-centered-P2 coefficient="
          << coefficient << " m2=" << rule.m2 << " m4=" << rule.m4
          << " H0=" << background.source().h0_km_s_mpc
          << " cdm=" << background.source().omega_cdm
          << " radiation=" << radiation.a4_e2
          << " radiation_error=" << radiation.error_estimate << '\n';
    need(ledger.write(input.str(), budget),
         "injected original input record retained");
    return state;
  };

  centered_cdm_controls(
      [](W h, W w, W e) {
        peer::State<W> state;
        state.v[1] = h;
        state.v[3] = w;
        state.v[4] = e;
        return peer::centered_cdm(state);
      },
      [&](W, W h, W v, W w, W, W z) {
        const auto state = state_with_slip(h, w, v, z);
        peer::State<W> derivative;
        need(peer::rhs(background, radiation.a4_e2, radiation.error_estimate,
                       wave, 1.L, state, derivative, rule, control, budget,
                       ledger, diagnostics) == S::ok,
             "independent full RHS available for CDM coordinate identity");
        return std::array<W, 3>{derivative.v[3] / epoch.hcal,
                                derivative.v[4] / epoch.hcal,
                                derivative.v[1] / epoch.hcal};
      },
      epoch.g, std::sqrt(epoch.x2));

  // Reconstruct ORIGINAL Einstein residual derivatives from the primitive
  // species laws and the actual angular derivative. This route does not call
  // constraint_rhs to manufacture an expected propagation result.
  auto state = state_with_slip(.02L, -.03L, .01L, -.004L);
  state.v[4] = .04L;
  state.v[5] = -.02L;
  state.v[6] = .08L;
  need(ledger.write("control-update e=.04 c=-.02 beta=.08\n", budget),
       "additional injected original coordinates retained");
  const auto readout = peer::readout(state, rule, epoch, &budget);
  peer::State<W> derivative;
  need(peer::rhs(background, radiation.a4_e2, radiation.error_estimate, wave,
                 1.L, state, derivative, rule, control, budget, ledger,
                 diagnostics) == S::ok,
       "full angular/trace derivatives computed for injected control");
  const W H = epoch.hcal, g = epoch.g, fr = epoch.fr, fc = epoch.fc;
  const W Z = state.v[2], phi = readout.phi, psi = readout.psi;
  const W delta_c = 3 + state.v[4] + 3 * phi;
  const W delta_r = readout.delta_r;
  const W Vc = p / 2 + state.v[3], Vr = -readout.b / 4;
  const W fc_N = -(2 * g + 1) * fc, fr_N = -2 * (g + 1) * fr;
  const W delta_c_N = 3 * Z + derivative.v[4] / H;
  const W delta_r_N = 4 * Z + derivative.v[5] / H;
  const W Vc_N = derivative.v[3] / H, Vr_N = -derivative.v[6] / (4 * H);
  const W Z_N = derivative.v[2] / H;
  const W q2_N = mean(rule, [&](std::size_t j, W mu) {
                   return legendre(2, mu) * derivative.R[j] / H;
                 }).real();
  const W s_N = -3 * (fr_N * readout.q2 + fr * q2_N);
  const W psi_N = Z - s_N;
  const W actual_E_N =
      Z_N + psi_N + epoch.x2 * (Z - 2 * g * phi) / 3 +
      (fc_N * delta_c + fc * delta_c_N + fr_N * delta_r + fr * delta_r_N) / 2;
  const W actual_W_N =
      Z_N + psi_N -
      1.5L * (fc_N * Vc + fc * Vc_N + 4 * (fr_N * Vr + fr * Vr_N) / 3);
  const W sigma = -epoch.x2 * readout.q2 / 2;
  const W fE = fr * epoch.x2 * readout.b * (rule.m2 - 1.L / 3) / 2;
  const W fW = -2 * fr * sigma * (1 - 1 / (3 * rule.m2));
  const W background_defect = g - 1 + 1.5L * fc + 2 * fr;
  near(actual_E_N,
       -(2 * g + 3) * readout.E + epoch.x2 * readout.W / 3 + fE +
           background_defect * Z,
       arithmetic, arithmetic,
       "full trace energy propagation retains all forcing");
  near(actual_W_N,
       -2 * readout.E - (g + 2) * readout.W + fW - background_defect * psi,
       arithmetic, arithmetic,
       "full trace momentum propagation retains all forcing");
  near(Z_N + (g + 2) * Z + psi_N + (2 * g + 1) * psi +
           epoch.x2 * readout.s / 3 - fr * delta_r / 2,
       -2 * readout.E, arithmetic, arithmetic,
       "original Einstein trace records actual damping correction");
  need(budget.rhs == 4 && budget.attempted_stages == 4 &&
           budget.ledger_bytes > 0,
       "injected controls retain original stage and resource lineage");
}

struct ReferencePoint {
  std::array<W, 3> value{}, error{};
};

ReferencePoint
native_comparison(const irred::cosmology::ThermalBackground &background,
                  const std::filesystem::path &directory,
                  std::size_t &campaign_bytes) {
  using namespace irred::cosmology;
  ReferencePoint shared_reference;
  // Exact declared request tuple. This is a synthetic source, not photons,
  // standard Planck species, ordinary sigma8, or a primary CMB calculation.
  need(background.status() == S::ok,
       "exact retained massless FD background admitted");
  constexpr double initial = 1e-14, target = 1e-4;
  const auto native = prepare_massless_fd_transfer(background, initial);
  need(native.status() == S::ok && native.background() != nullptr,
       "new native collisionless owner preparation");
  need(native.background()->source().h0_km_s_mpc == 70 &&
           native.background()->source().omega_gamma == 0 &&
           native.background()->source().omega_massless_nonphoton == 0 &&
           native.background()->source().omega_b == 0 &&
           native.background()->source().omega_cdm == .3 &&
           native.background()->source().species.size() == 1 &&
           native.background()->source().species[0].mass_ev == 0 &&
           native.background()->source().species[0].temperature_today_ev ==
               .0002 &&
           native.background()->source().species[0].statistical_weight == 2,
       "new native owner retains exact explicit FD source");
  const std::array<double, 2> waves{1e-7, .01};
  constexpr unsigned outputs = massless_fd_comoving_cdm |
                               massless_fd_spatial_potential |
                               massless_fd_lapse_potential;
  const MasslessFDTransferPolicy policy;
  need(policy.absolute_tolerance == 1e-7 && policy.relative_tolerance == 1e-4 &&
           policy.maximum_log_step == .04 && policy.maximum_phase_step == .04 &&
           policy.maximum_constraint_residual == 1e-7 &&
           policy.maximum_points == 16 &&
           policy.maximum_rhs_per_point == 1000000 &&
           policy.maximum_total_rhs == 4000000 &&
           policy.maximum_scalar_updates_per_point == 128000000 &&
           policy.maximum_total_scalar_updates == 512000000 &&
           policy.maximum_background_queries == 4000000 &&
           policy.maximum_age_quadrature_evaluations == 200000 &&
           policy.maximum_native_bytes == 16 * 1024 * 1024,
       "native defaults retain exact accepted request and resource contract");
  peer::Budget native_evidence_budget;
  native_evidence_budget.campaign_bytes = &campaign_bytes;
  peer::Ledger native_evidence(
      (directory / "native-request-and-return.txt").string());
  auto retain = [&](const std::ostringstream &record) {
    need(native_evidence.write(record.str(), native_evidence_budget),
         "native request or returned witness retained and flushed");
  };
  auto work_record = [](std::ostream &record,
                        const MasslessFDTransferWork &work) {
    record << " rhs=" << work.rhs << " scalars=" << work.scalar_updates
           << " background=" << work.background_queries
           << " age_evaluations=" << work.age_quadrature_evaluations
           << " momentum_callbacks=" << work.momentum_callbacks;
  };
  const auto &source = background.source();
  std::ostringstream request;
  request << std::setprecision(std::numeric_limits<W>::max_digits10)
          << "native-request status=started completion=unobserved"
          << " incomplete_disposition=interrupted actual_counters=unavailable"
          << " final_counts=unavailable model=" << massless_fd_transfer_model_id
          << " method=" << massless_fd_transfer_method_id
          << " arithmetic=" << massless_fd_transfer_arithmetic_id
          << " H0_km_s_Mpc=" << source.h0_km_s_mpc
          << " omega_gamma=" << source.omega_gamma
          << " omega_massless_nonphoton=" << source.omega_massless_nonphoton
          << " omega_b=" << source.omega_b << " omega_cdm=" << source.omega_cdm
          << " species_count=" << source.species.size()
          << " species0_mass_eV=" << source.species[0].mass_ev
          << " species0_T0_eV=" << source.species[0].temperature_today_ev
          << " species0_weight=" << source.species[0].statistical_weight
          << " initial_a=" << initial << " target_a=" << target
          << " ordered_k_Mpc_inverse=" << waves[0] << ',' << waves[1]
          << " output_mask=" << outputs
          << " absolute=" << policy.absolute_tolerance
          << " relative=" << policy.relative_tolerance
          << " maximum_log_step=" << policy.maximum_log_step
          << " maximum_phase_step=" << policy.maximum_phase_step
          << " maximum_constraint_residual="
          << policy.maximum_constraint_residual
          << " maximum_points=" << policy.maximum_points
          << " maximum_rhs_per_point=" << policy.maximum_rhs_per_point
          << " maximum_total_rhs=" << policy.maximum_total_rhs
          << " maximum_scalars_per_point="
          << policy.maximum_scalar_updates_per_point
          << " maximum_total_scalars=" << policy.maximum_total_scalar_updates
          << " maximum_background_queries=" << policy.maximum_background_queries
          << " maximum_age_evaluations="
          << policy.maximum_age_quadrature_evaluations
          << " maximum_native_bytes=" << policy.maximum_native_bytes
          << " thermal_method="
          << thermal_momentum_method_id(policy.thermal.momentum_method)
          << " thermal_absolute=" << policy.thermal.absolute_tolerance
          << " thermal_relative=" << policy.thermal.relative_tolerance
          << " thermal_maximum_callbacks="
          << policy.thermal.maximum_callbacks_per_evaluation
          << " thermal_maximum_total_callbacks="
          << policy.thermal.maximum_total_callbacks
          << " thermal_maximum_depth=" << policy.thermal.maximum_depth
          << " thermal_maximum_points=" << policy.thermal.maximum_points
          << " thermal_maximum_species=" << policy.thermal.maximum_species
          << " thermal_maximum_native_bytes="
          << policy.thermal.maximum_native_bytes << " evidence_case_cap="
          << native_evidence_budget.maximum_ledger_bytes
          << " evidence_campaign_cap="
          << native_evidence_budget.maximum_campaign_bytes
          << " preparation_callbacks=" << background.preparation_callbacks()
          << '\n';
  retain(request);
  const auto got = native.evaluate(waves, target, outputs, policy);
  std::ostringstream returned;
  returned
      << "native-return completion=returned actual_counters=available status="
      << static_cast<int>(got.status) << " rows=" << got.rows.size()
      << " target_a=" << std::setprecision(std::numeric_limits<W>::max_digits10)
      << got.scale_factor << " output_mask=" << got.requested_outputs;
  work_record(returned, got.work);
  returned << '\n';
  retain(returned);
  // Preserve failed values and all attempted trial states before any admission
  // assertion. A watchdog interruption leaves only the flushed started record;
  // it does not fabricate work counters for an evaluate call that never
  // returned.
  for (std::size_t row_index = 0; row_index < got.rows.size(); ++row_index) {
    const auto &row = got.rows[row_index];
    std::ostringstream row_record;
    row_record << std::setprecision(std::numeric_limits<W>::max_digits10)
               << "native-row row=" << row_index
               << " k=" << row.wavenumber_mpc_inverse
               << " constraint=" << row.maximum_constraint_residual
               << " closure=" << row.maximum_closure_defect
               << " background_identity="
               << row.maximum_background_identity_defect
               << " eta=" << row.conformal_age_mpc
               << " eta_error=" << row.conformal_age_error_mpc
               << " attempts_recorded=" << row.attempts_recorded;
    work_record(row_record, row.work);
    row_record << '\n';
    retain(row_record);
    const std::array<const MasslessFDTransferValue *, 3> values{
        &row.comoving_cdm, &row.spatial_potential, &row.lapse_potential};
    for (std::size_t field = 0; field < values.size(); ++field) {
      const auto &value = *values[field];
      std::ostringstream record;
      record << std::setprecision(std::numeric_limits<W>::max_digits10)
             << "native-field row=" << row_index << " field=" << field
             << " status=" << static_cast<int>(value.status)
             << " available=" << bool(value.value);
      if (value.value)
        record << " value=" << *value.value;
      record << " error=" << value.absolute_error_estimate
             << " time=" << value.time_refinement
             << " initial=" << value.initial_refinement
             << " hierarchy=" << value.hierarchy_refinement
             << " background_age=" << value.background_age_sensitivity
             << " arithmetic_cast=" << value.arithmetic_cast_sensitivity
             << '\n';
      retain(record);
    }
    for (unsigned trial = 0;
         trial <
         std::min<std::size_t>(row.attempts_recorded, row.attempts.size());
         ++trial) {
      const auto &attempt = row.attempts[trial];
      const auto &initial_state = attempt.initial;
      std::ostringstream record;
      record << std::setprecision(std::numeric_limits<W>::max_digits10)
             << "native-attempt row=" << row_index << " trial=" << trial
             << " status=" << static_cast<int>(attempt.status)
             << " hierarchy_l=" << attempt.hierarchy_l
             << " time_divisor=" << attempt.time_divisor
             << " control_kind=" << attempt.control_kind
             << " control_sign=" << attempt.control_sign
             << " endpoint_available=" << attempt.endpoint_available
             << " lapse_available=" << attempt.lapse_available
             << " endpoint_scale_factor=" << attempt.endpoint_scale_factor
             << " metric_epoch_scale_factor="
             << attempt.metric_epoch_scale_factor
             << " maximum_stage_phase_bound="
             << attempt.maximum_stage_phase_bound
             << " maximum_phase_increment_upper="
             << attempt.maximum_phase_increment_upper;
      if (attempt.endpoint_available)
        record << " Delta=" << attempt.endpoint[0]
               << " phi=" << attempt.endpoint[1]
               << " endpoint_scaled_shear=" << attempt.endpoint_scaled_shear;
      else
        record << " Delta=unavailable phi=unavailable "
                  "endpoint_scaled_shear=unavailable";
      if (attempt.endpoint_available && attempt.lapse_available)
        record << " psi=" << attempt.endpoint[2];
      else
        record << " psi=unavailable";
      record << " eta=" << attempt.endpoint_eta_mpc
             << " constraint=" << attempt.maximum_constraint_residual
             << " initial_status=" << static_cast<int>(initial_state.status);
      // An early refusal leaves the default invalid initializer unpopulated.
      // A later refusal may retain a populated but unadmitted source witness.
      if (initial_state.status == S::invalid_input)
        record << " initial_fields=unavailable";
      else
        record << " initial_fields=available"
               << " initial_a=" << initial_state.scale_factor
               << " initial_eta=" << initial_state.eta_mpc
               << " initial_eta_error=" << initial_state.eta_error_mpc
               << " initial_H=" << initial_state.hcal_mpc_inverse
               << " leading_phi=" << initial_state.leading_phi
               << " leading_psi=" << initial_state.leading_psi
               << " leading_delta_c=" << initial_state.leading_delta_c
               << " leading_delta_r=" << initial_state.leading_delta_r
               << " leading_theta=" << initial_state.leading_theta
               << " leading_sigma=" << initial_state.leading_sigma
               << " projected_phi=" << initial_state.projected_phi
               << " projected_psi=" << initial_state.projected_psi
               << " projected_Vc=" << initial_state.projected_vc
               << " projected_Delta=" << initial_state.projected_delta
               << " projected_Z=" << initial_state.projected_phi_n
               << " scaled_shear=" << initial_state.scaled_shear
               << " lapse_correction=" << initial_state.lapse_correction
               << " velocity_correction=" << initial_state.velocity_correction
               << " hamiltonian_before=" << initial_state.hamiltonian_before
               << " hamiltonian_after=" << initial_state.hamiltonian_after
               << " radiation_fraction=" << initial_state.radiation_fraction
               << " initial_closure=" << initial_state.closure_defect;
      work_record(record, attempt.work);
      record << '\n';
      retain(record);
    }
  }
  need(got.status == S::ok && got.rows.size() == waves.size() &&
           got.requested_outputs == outputs && got.scale_factor == target,
       "native exact ordered SDK request admitted");
  need(got.work.rhs > 0 && got.work.scalar_updates > 0 &&
           got.work.background_queries > 0 &&
           got.work.age_quadrature_evaluations > 0,
       "native all-attempt setup and calculation counters retained");
  need(got.work.rhs <= 4000000 && got.work.scalar_updates <= 512000000 &&
           got.work.background_queries <= 4000000 &&
           got.work.age_quadrature_evaluations <= 200000,
       "native exact batch caps remain bounded");

  for (std::size_t row_index = 0; row_index < waves.size(); ++row_index) {
    const auto &row = got.rows[row_index];
    need(row.wavenumber_mpc_inverse == waves[row_index],
         "native preserves original source axis order");
    need(row.work.rhs <= policy.maximum_rhs_per_point &&
             row.work.scalar_updates <=
                 policy.maximum_scalar_updates_per_point &&
             row.attempts_recorded == row.attempts.size(),
         "native original per-point caps and all trial identities retained");
    for (const auto &attempt : row.attempts) {
      need(attempt.status == S::ok && attempt.endpoint_available &&
               attempt.lapse_available &&
               attempt.endpoint_scale_factor ==
                   attempt.metric_epoch_scale_factor &&
               attempt.initial.status == S::ok,
           "native accepted row retains every actual trial and source "
           "initializer");
      need(
          std::isfinite(attempt.maximum_stage_phase_bound) &&
              attempt.maximum_stage_phase_bound >= 0 &&
              attempt.maximum_stage_phase_bound <=
                  W(policy.maximum_phase_step) &&
              std::isfinite(attempt.maximum_phase_increment_upper) &&
              attempt.maximum_phase_increment_upper >= 0 &&
              attempt.maximum_phase_increment_upper <=
                  W(policy.maximum_phase_step),
          "native accepted trial obeys original stage and attained phase caps");
    }
    peer::Budget budget;
    budget.campaign_bytes = &campaign_bytes;
    const auto reference = peer::campaign(
        background, waves[row_index], initial, target,
        (directory / ("comparison-row-" + std::to_string(row_index) + ".txt"))
            .string(),
        budget);
    std::ostringstream peer_return;
    peer_return
        << std::setprecision(std::numeric_limits<W>::max_digits10)
        << "peer-return row=" << row_index
        << " status=" << static_cast<int>(reference.status)
        << " admitted_fields="
        << (reference.status == S::ok ? "available" : "withheld")
        << " diagnostics_partial=" << (reference.status == S::ok ? 0 : 1)
        << " budget_status=" << static_cast<int>(budget.status)
        << " rhs=" << budget.rhs << " scalars=" << budget.scalar_updates
        << " background=" << budget.background_queries
        << " age_evaluations=" << budget.age_quadrature_evaluations
        << " momentum_callbacks=" << budget.momentum_callbacks
        << " attempts=" << budget.attempted_stages
        << " ledger_bytes=" << budget.ledger_bytes
        << " logging_refusals=" << budget.logging_refusals
        << " maximum_E=" << reference.diagnostics.maximum_E
        << " maximum_W=" << reference.diagnostics.maximum_W
        << " maximum_trace=" << reference.diagnostics.maximum_trace
        << " maximum_correction=" << reference.diagnostics.maximum_correction
        << " maximum_slip=" << reference.diagnostics.maximum_slip
        << " maximum_conservation="
        << reference.diagnostics.maximum_conservation
        << " maximum_recenter=" << reference.diagnostics.maximum_recenter
        << " maximum_subspace=" << reference.diagnostics.maximum_subspace_defect
        << " eta=" << reference.eta << " eta_error=" << reference.eta_error;
    for (std::size_t field = 0; field < reference.value.size(); ++field) {
      if (reference.status == S::ok)
        peer_return << " value" << field << '=' << reference.value[field];
      peer_return << " diagnostic_error" << field << '='
                  << reference.error[field];
      for (std::size_t term = 0; term < reference.contributions[field].size();
           ++term)
        peer_return << " diagnostic_contribution" << field << '_' << term << '='
                    << reference.contributions[field][term];
    }
    peer_return << '\n';
    retain(peer_return);
    need(reference.status == S::ok && budget.status == S::ok,
         "independent angular trace campaign meets original allocation");
    if (row_index == 1) {
      shared_reference.value = reference.value;
      shared_reference.error = reference.error;
    }
    need(budget.rhs > 0 && budget.scalar_updates > 0 &&
             budget.background_queries > 0 && budget.ledger_bytes > 0 &&
             budget.attempted_stages > 0,
         "independent peer attempts and original trace evidence retained");
    need(budget.rhs <= budget.maximum_rhs &&
             budget.scalar_updates <= budget.maximum_scalar_updates &&
             budget.ledger_bytes <= budget.maximum_ledger_bytes &&
             campaign_bytes <= budget.maximum_campaign_bytes,
         "peer original all-attempt work and external evidence caps respected");
    const std::array<const MasslessFDTransferValue *, 3> values{
        &row.comoving_cdm, &row.spatial_potential, &row.lapse_potential};
    for (std::size_t field = 0; field < values.size(); ++field) {
      const auto &value = *values[field];
      need(value.status == S::ok && value.value && std::isfinite(*value.value),
           "every requested signed native field has an admitted value");
      const W epsilon = 1e-7L + 1e-4L * std::abs(W(*value.value));
      const std::array<W, 5> native_terms{
          value.time_refinement, value.initial_refinement,
          value.hierarchy_refinement, value.background_age_sensitivity,
          value.arithmetic_cast_sensitivity};
      W native_sum = 0, reference_sum = 0;
      for (std::size_t term = 0; term < native_terms.size(); ++term) {
        need(std::isfinite(native_terms[term]) && native_terms[term] >= 0 &&
                 native_terms[term] <= epsilon / 6,
             "native named contribution stays inside its fixed share");
        native_sum += native_terms[term];
        const W contribution = reference.contributions[field][term];
        need(std::isfinite(contribution) && contribution >= 0 &&
                 contribution <= epsilon / 30,
             "peer named contribution stays inside its fixed share");
        reference_sum += contribution;
      }
      need(native_sum <= 5 * epsilon / 6 && reference_sum <= epsilon / 6,
           "native and peer allocations do not widen the requested budget");
      const W native_error = value.absolute_error_estimate;
      need(std::isfinite(native_error) && native_error >= 0 &&
               native_error <= 5 * epsilon / 6,
           "native total numerical diagnostic stays inside fixed allocated "
           "share");
      near(native_error, native_sum, 0,
           8 * std::numeric_limits<double>::epsilon(),
           "native cast total retains exactly the five declared contributions");
      need(reference.error[field] >= reference_sum &&
               reference.error[field] <= epsilon / 6,
           "peer uncertainty retains all contributors inside reference share");
      need(std::isfinite(reference.value[field]) &&
               std::abs(W(*value.value) - reference.value[field]) +
                       reference.error[field] <=
                   epsilon,
           "native and independently angular/trace result agree within fixed "
           "budget");
      std::cout << std::setprecision(std::numeric_limits<W>::max_digits10)
                << "row=" << row_index << " k=" << waves[row_index]
                << " field=" << field << " native=" << *value.value
                << " peer=" << reference.value[field]
                << " peer_error=" << reference.error[field]
                << " native_error=" << value.absolute_error_estimate << '\n';
    }
    need(row.maximum_constraint_residual <= 1e-7 &&
             row.maximum_closure_defect <= 1e-7 &&
             row.maximum_background_identity_defect <= 1e-7,
         "original native Einstein and retained background gates remain");
    for (const auto &witness : row.initial_states) {
      need(witness.status == S::ok && witness.scale_factor > 0 &&
               std::isfinite(witness.projected_delta) &&
               std::isfinite(witness.leading_delta_c) &&
               witness.leading_delta_c > 0,
           "projected and signed gauge-density source witnesses retained");
    }
    // The lowest mode is compared absolutely. No tiny relative coefficient,
    // endpoint zeta renormalization or positive-value replacement is asserted.
    std::cout << "row=" << row_index << " native_rhs=" << row.work.rhs
              << " peer_rhs=" << budget.rhs
              << " peer_bytes=" << budget.ledger_bytes
              << " original_trace=" << reference.diagnostics.maximum_trace
              << " correction=" << reference.diagnostics.maximum_correction
              << '\n';
  }
  return shared_reference;
}

using Background = irred::cosmology::ThermalBackground;
using Model = irred::cosmology::ThermalFlatModel;
using NativeBatch = irred::cosmology::MasslessFDTransferBatch;
using NativePolicy = irred::cosmology::MasslessFDTransferPolicy;

// These additional cases keep their own declared original native work cap.
// They share the campaign's external evidence cap and never retry a refusal.
NativeBatch retained_native_control(const Background &background,
                                    const std::filesystem::path &directory,
                                    std::size_t &campaign_bytes,
                                    const char *name,
                                    std::span<const double> waves,
                                    double target, NativePolicy policy = {}) {
  using namespace irred::cosmology;
  peer::Budget record_budget;
  record_budget.campaign_bytes = &campaign_bytes;
  peer::Ledger ledger((directory / (std::string(name) + ".txt")).string());
  auto record = [&](const std::ostringstream &text) {
    need(ledger.write(text.str(), record_budget),
         "additional native control record retained before admission check");
  };
  auto work = [](std::ostream &out, const MasslessFDTransferWork &count) {
    out << " rhs=" << count.rhs << " scalars=" << count.scalar_updates
        << " background=" << count.background_queries
        << " age_evaluations=" << count.age_quadrature_evaluations
        << " momentum_callbacks=" << count.momentum_callbacks;
  };
  constexpr unsigned outputs = massless_fd_comoving_cdm |
                               massless_fd_spatial_potential |
                               massless_fd_lapse_potential;
  std::ostringstream request;
  const auto &source = background.source();
  request << std::setprecision(std::numeric_limits<W>::max_digits10)
          << "control-request name=" << name
          << " status=started completion=unobserved "
             "incomplete_disposition=interrupted"
          << " actual_counters=unavailable model="
          << massless_fd_transfer_model_id
          << " method=" << massless_fd_transfer_method_id
          << " H0=" << source.h0_km_s_mpc << " gamma=" << source.omega_gamma
          << " massless_nonphoton=" << source.omega_massless_nonphoton
          << " baryon=" << source.omega_b << " cdm=" << source.omega_cdm
          << " initial_a=" << double(1e-14) << " target_a=" << target
          << " output_mask=" << outputs << " waves_count=" << waves.size()
          << " absolute=" << policy.absolute_tolerance
          << " relative=" << policy.relative_tolerance
          << " maximum_log_step=" << policy.maximum_log_step
          << " maximum_phase_step=" << policy.maximum_phase_step
          << " constraint_cap=" << policy.maximum_constraint_residual
          << " points_cap=" << policy.maximum_points
          << " rhs_per_point_cap=" << policy.maximum_rhs_per_point
          << " rhs_total_cap=" << policy.maximum_total_rhs
          << " scalars_per_point_cap="
          << policy.maximum_scalar_updates_per_point
          << " scalars_total_cap=" << policy.maximum_total_scalar_updates
          << " background_cap=" << policy.maximum_background_queries
          << " age_cap=" << policy.maximum_age_quadrature_evaluations
          << " payload_cap=" << policy.maximum_native_bytes
          << " preparation_callbacks=" << background.preparation_callbacks()
          << '\n';
  for (std::size_t j = 0; j < source.species.size(); ++j)
    request << "source-species index=" << j
            << " mass=" << source.species[j].mass_ev
            << " T0=" << source.species[j].temperature_today_ev
            << " weight=" << source.species[j].statistical_weight << '\n';
  for (std::size_t j = 0; j < waves.size(); ++j)
    request << "source-wave index=" << j << " k=" << waves[j] << '\n';
  record(request);
  const auto owner = prepare_massless_fd_transfer(background, 1e-14);
  std::ostringstream prepared;
  prepared << "control-owner-prepare name=" << name
           << " status=" << static_cast<int>(owner.status()) << '\n';
  record(prepared);
  need(owner.status() == S::ok && owner.background(),
       "additional case retains admitted source owner");
  need(owner.background()->source().h0_km_s_mpc == source.h0_km_s_mpc &&
           owner.background()->source().omega_gamma == source.omega_gamma &&
           owner.background()->source().omega_massless_nonphoton ==
               source.omega_massless_nonphoton &&
           owner.background()->source().omega_b == source.omega_b &&
           owner.background()->source().omega_cdm == source.omega_cdm &&
           owner.background()->source().species.size() ==
               source.species.size() &&
           owner.background()->momentum_method() ==
               background.momentum_method(),
       "additional native owner keeps exact supplied source and momentum "
       "method");
  for (std::size_t j = 0; j < source.species.size(); ++j) {
    const auto &retained = owner.background()->source().species[j];
    need(retained.mass_ev == source.species[j].mass_ev &&
             retained.temperature_today_ev ==
                 source.species[j].temperature_today_ev &&
             retained.statistical_weight ==
                 source.species[j].statistical_weight,
         "additional native owner keeps physical species tuple order");
  }
  const auto batch = owner.evaluate(waves, target, outputs, policy);
  std::ostringstream returned;
  returned << "control-return name=" << name
           << " completion=returned status=" << static_cast<int>(batch.status)
           << " rows=" << batch.rows.size();
  work(returned, batch.work);
  returned << '\n';
  record(returned);
  for (std::size_t j = 0; j < batch.rows.size(); ++j) {
    const auto &row = batch.rows[j];
    std::ostringstream summary;
    summary << std::setprecision(std::numeric_limits<W>::max_digits10)
            << "control-row index=" << j << " k=" << row.wavenumber_mpc_inverse
            << " eta=" << row.conformal_age_mpc
            << " eta_error=" << row.conformal_age_error_mpc
            << " constraint=" << row.maximum_constraint_residual
            << " closure=" << row.maximum_closure_defect
            << " background_identity=" << row.maximum_background_identity_defect
            << " attempts=" << row.attempts_recorded;
    work(summary, row.work);
    summary << '\n';
    record(summary);
    const std::array<const MasslessFDTransferValue *, 3> values{
        &row.comoving_cdm, &row.spatial_potential, &row.lapse_potential};
    for (std::size_t field = 0; field < values.size(); ++field) {
      const auto &value = *values[field];
      std::ostringstream line;
      line << std::setprecision(std::numeric_limits<W>::max_digits10)
           << "control-field row=" << j << " field=" << field
           << " status=" << static_cast<int>(value.status)
           << " available=" << bool(value.value);
      if (value.value)
        line << " value=" << *value.value;
      line << " error=" << value.absolute_error_estimate
           << " time=" << value.time_refinement
           << " initial=" << value.initial_refinement
           << " hierarchy=" << value.hierarchy_refinement
           << " background_age=" << value.background_age_sensitivity
           << " arithmetic=" << value.arithmetic_cast_sensitivity << '\n';
      record(line);
    }
    for (unsigned trial = 0;
         trial <
         std::min<std::size_t>(row.attempts_recorded, row.attempts.size());
         ++trial) {
      const auto &a = row.attempts[trial];
      std::ostringstream line;
      line << std::setprecision(std::numeric_limits<W>::max_digits10)
           << "control-attempt row=" << j << " trial=" << trial
           << " status=" << static_cast<int>(a.status)
           << " hierarchy_l=" << a.hierarchy_l
           << " time_divisor=" << a.time_divisor
           << " control_kind=" << a.control_kind
           << " control_sign=" << a.control_sign
           << " endpoint_available=" << a.endpoint_available
           << " lapse_available=" << a.lapse_available
           << " endpoint_scale_factor=" << a.endpoint_scale_factor
           << " metric_epoch_scale_factor=" << a.metric_epoch_scale_factor
           << " maximum_stage_phase_bound=" << a.maximum_stage_phase_bound
           << " maximum_phase_increment_upper="
           << a.maximum_phase_increment_upper;
      if (a.endpoint_available)
        line << " Delta=" << a.endpoint[0] << " phi=" << a.endpoint[1]
             << " endpoint_scaled_shear=" << a.endpoint_scaled_shear;
      else
        line << " Delta=unavailable phi=unavailable "
                "endpoint_scaled_shear=unavailable";
      if (a.endpoint_available && a.lapse_available)
        line << " psi=" << a.endpoint[2];
      else
        line << " psi=unavailable";
      line << " initial_status=" << static_cast<int>(a.initial.status);
      if (a.initial.status == S::invalid_input)
        line << " initial_fields=unavailable";
      else
        line << " initial_fields=available initial_a=" << a.initial.scale_factor
             << " initial_eta=" << a.initial.eta_mpc
             << " initial_phi=" << a.initial.projected_phi
             << " initial_psi=" << a.initial.projected_psi
             << " initial_Vc=" << a.initial.projected_vc
             << " initial_Delta=" << a.initial.projected_delta
             << " initial_Z=" << a.initial.projected_phi_n;
      work(line, a.work);
      line << '\n';
      record(line);
    }
  }
  need(batch.work.rhs <= policy.maximum_total_rhs &&
           batch.work.scalar_updates <= policy.maximum_total_scalar_updates &&
           batch.work.background_queries <= policy.maximum_background_queries &&
           batch.work.age_quadrature_evaluations <=
               policy.maximum_age_quadrature_evaluations,
       "additional controls retain original all-attempt batch caps");
  for (const auto &row : batch.rows)
    need(row.work.rhs <= policy.maximum_rhs_per_point &&
             row.work.scalar_updates <= policy.maximum_scalar_updates_per_point,
         "additional controls retain original all-attempt point caps");
  return batch;
}

bool admitted_control_row(const irred::cosmology::MasslessFDTransferRow &row,
                          bool allow_numerical_refusal,
                          bool allow_boundary_domain_refusal = false) {
  using namespace irred::cosmology;
  const NativePolicy policy;
  const std::array<const MasslessFDTransferValue *, 3> values{
      &row.comoving_cdm, &row.spatial_potential, &row.lapse_potential};
  bool admitted = true;
  for (const auto *value : values) {
    if (value->status != S::ok) {
      need(!value->value && allow_numerical_refusal &&
               (value->status == S::work_limit ||
                value->status == S::conditioning_budget_exceeded ||
                (allow_boundary_domain_refusal &&
                 value->status == S::outside_domain &&
                 row.attempts_recorded > 0)),
           "adversarial numerical refusal retains unqualified status and "
           "withholds field");
      admitted = false;
      continue;
    }
    need(value->value && std::isfinite(*value->value),
         "admitted additional control signed field available");
    const W epsilon = 1e-7L + 1e-4L * std::abs(W(*value->value));
    const std::array<W, 5> terms{
        value->time_refinement, value->initial_refinement,
        value->hierarchy_refinement, value->background_age_sensitivity,
        value->arithmetic_cast_sensitivity};
    W sum = 0;
    for (W term : terms) {
      need(std::isfinite(term) && term >= 0 && term <= epsilon / 6,
           "additional admitted field keeps each original numerical share");
      sum += term;
    }
    need(sum <= 5 * epsilon / 6 &&
             std::isfinite(value->absolute_error_estimate) &&
             value->absolute_error_estimate >= 0 &&
             value->absolute_error_estimate <= 5 * epsilon / 6,
         "additional admitted field keeps original total allocation");
    near(value->absolute_error_estimate, sum, 0,
         8 * std::numeric_limits<double>::epsilon(),
         "additional admitted error retains all named contributions");
  }
  if (!admitted) {
    std::cout << "UNQUALIFIED additional native control k="
              << row.wavenumber_mpc_inverse
              << " status=" << static_cast<int>(row.comoving_cdm.status)
              << '\n';
    return false;
  }
  need(row.attempts_recorded == row.attempts.size() &&
           row.maximum_constraint_residual <= 1e-7 &&
           row.maximum_closure_defect <= 1e-7 &&
           row.maximum_background_identity_defect <= 1e-7,
       "additional admitted case retains original trials and "
       "Einstein/background gates");
  for (const auto &a : row.attempts)
    need(a.status == S::ok && a.endpoint_available && a.lapse_available &&
             a.endpoint_scale_factor == a.metric_epoch_scale_factor &&
             std::isfinite(a.maximum_stage_phase_bound) &&
             a.maximum_stage_phase_bound >= 0 &&
             a.maximum_stage_phase_bound <= W(policy.maximum_phase_step) &&
             std::isfinite(a.maximum_phase_increment_upper) &&
             a.maximum_phase_increment_upper >= 0 &&
             a.maximum_phase_increment_upper <= W(policy.maximum_phase_step),
         "additional admitted trial phase and metric epoch witnesses remain "
         "bound");
  return true;
}

void compare_reference_point(const irred::cosmology::MasslessFDTransferRow &row,
                             const ReferencePoint &reference) {
  const std::array<const irred::cosmology::MasslessFDTransferValue *, 3> values{
      &row.comoving_cdm, &row.spatial_potential, &row.lapse_potential};
  for (std::size_t field = 0; field < values.size(); ++field) {
    const W value = *values[field]->value;
    const W epsilon = 1e-7L + 1e-4L * std::abs(value);
    need(std::isfinite(reference.value[field]) &&
             std::isfinite(reference.error[field]) &&
             reference.error[field] >= 0 &&
             reference.error[field] <= epsilon / 6 &&
             std::abs(value - reference.value[field]) +
                     reference.error[field] <=
                 epsilon,
         "additional native case meets unchanged independent reference "
         "allocation");
  }
}

void species_and_epoch_controls(const Background &baseline,
                                const ReferencePoint &reference,
                                const std::filesystem::path &directory,
                                std::size_t &campaign_bytes) {
  using namespace irred::cosmology;
  const std::array<double, 1> wave{.01};
  peer::Budget identity_budget;
  identity_budget.campaign_bytes = &campaign_bytes;
  peer::Ledger identity_record(
      (directory / "species-identity-setup.txt").string());
  auto retain_identity = [&](const std::ostringstream &line) {
    need(identity_record.write(line.str(), identity_budget),
         "species identity inputs/setup counters retained before checks");
  };
  need(identity_budget.charge(identity_budget.background_queries, 1,
                              identity_budget.maximum_background_queries),
       "baseline species-identity background query charged");
  const auto base_radiation = baseline.scaled_expansion(0);
  identity_budget.momentum_callbacks += base_radiation.callbacks;
  std::ostringstream base_record;
  base_record << std::setprecision(std::numeric_limits<W>::max_digits10)
              << "baseline-radiation status="
              << static_cast<int>(base_radiation.status)
              << " P0=" << base_radiation.a4_e2
              << " error=" << base_radiation.error_estimate
              << " background=" << identity_budget.background_queries
              << " momentum_callbacks=" << identity_budget.momentum_callbacks
              << '\n';
  retain_identity(base_record);
  need(base_radiation.status == S::ok,
       "baseline retained radiation for species identity");
  auto radiation_weight = [&](const Model &m) {
    need(identity_budget.scalars(1 + 4 * m.species.size()),
         "independent species-weight construction charged");
    W sum = 0;
    for (const auto &species : m.species) {
      const W t = species.temperature_today_ev;
      sum += W(species.statistical_weight) * t * t * t * t;
    }
    return sum;
  };
  const Model split{70, 0, 0, 0, .3, {{0, .0002, 1}, {0, .0002, 1}}};
  const Model degenerate{70, 0, 0, 0, .3, {{0, .0002, 1}, {0, .0001, 16}}};
  Model reversed = degenerate;
  std::reverse(reversed.species.begin(), reversed.species.end());
  const std::array<const Model *, 3> models{&split, &degenerate, &reversed};
  const std::array<const char *, 3> names{
      "species-split", "species-temperature-weight", "species-reversed"};
  for (std::size_t variation = 0; variation < models.size(); ++variation) {
    const auto &model = *models[variation];
    std::ostringstream input;
    input << std::setprecision(std::numeric_limits<W>::max_digits10)
          << "species-case name=" << names[variation]
          << " status=started completion=unobserved "
             "incomplete_disposition=interrupted"
          << " H0=" << model.h0_km_s_mpc << " cdm=" << model.omega_cdm
          << " initial_a=" << double(1e-14) << " target_a=" << double(1e-4)
          << " k=" << wave[0] << '\n';
    for (std::size_t j = 0; j < model.species.size(); ++j)
      input << "species index=" << j << " mass=" << model.species[j].mass_ev
            << " T0=" << model.species[j].temperature_today_ev
            << " weight=" << model.species[j].statistical_weight << '\n';
    retain_identity(input);
    near(radiation_weight(model), radiation_weight(baseline.source()), 0,
         arithmetic, "independent massless g*T0^4 degeneracy");
    const auto background = prepare_thermal_background(model);
    identity_budget.momentum_callbacks += background.preparation_callbacks();
    std::ostringstream prepared;
    prepared << "species-prepared name=" << names[variation]
             << " status=" << static_cast<int>(background.status())
             << " preparation_callbacks=" << background.preparation_callbacks()
             << '\n';
    retain_identity(prepared);
    need(background.status() == S::ok &&
             background.source().species.size() == 2,
         "explicit split/temperature-weight source admitted without source "
         "collapse");
    for (std::size_t j = 0; j < model.species.size(); ++j)
      need(background.source().species[j].mass_ev == model.species[j].mass_ev &&
               background.source().species[j].temperature_today_ev ==
                   model.species[j].temperature_today_ev &&
               background.source().species[j].statistical_weight ==
                   model.species[j].statistical_weight,
           "thermal owner preserves exact distinct physical tuple order");
    need(identity_budget.charge(identity_budget.background_queries, 1,
                                identity_budget.maximum_background_queries),
         "species identity retained-radiation query charged");
    const auto radiation = background.scaled_expansion(0);
    identity_budget.momentum_callbacks += radiation.callbacks;
    std::ostringstream radiation_record;
    radiation_record << std::setprecision(std::numeric_limits<W>::max_digits10)
                     << "species-radiation name=" << names[variation]
                     << " status=" << static_cast<int>(radiation.status)
                     << " P0=" << radiation.a4_e2
                     << " error=" << radiation.error_estimate
                     << " background=" << identity_budget.background_queries
                     << " scalars=" << identity_budget.scalar_updates
                     << " momentum_callbacks="
                     << identity_budget.momentum_callbacks << '\n';
    retain_identity(radiation_record);
    need(radiation.status == S::ok &&
             std::abs(radiation.a4_e2 - base_radiation.a4_e2) <=
                 radiation.error_estimate + base_radiation.error_estimate,
         "retained thermal radiation matches mathematical massless degeneracy");
    const auto got = retained_native_control(
        background, directory, campaign_bytes, names[variation], wave, 1e-4);
    need(got.status == S::ok && got.rows.size() == 1 &&
             got.rows[0].wavenumber_mpc_inverse == wave[0],
         "one-k species control preserves ordered required request");
    need(admitted_control_row(got.rows[0], false),
         "required species transfer comparison earns numerical admission");
    // The source algebra fixes the same total radiation. Reuse the existing
    // baseline angular/TRACE reference rather than repeat its full campaign.
    compare_reference_point(got.rows[0], reference);
  }
  std::ostringstream near_input;
  near_input << "near-radiation-prepare status=started completion=unobserved"
             << " incomplete_disposition=interrupted H0=70 gamma=0 "
                "massless_nonphoton=0"
             << " baryon=0 cdm=1e-6 species0_mass=0 species0_T0=.0002 "
                "species0_weight=2\n";
  retain_identity(near_input);
  const auto near_radiation =
      prepare_thermal_background({70, 0, 0, 0, 1e-6, {{0, .0002, 2}}});
  identity_budget.momentum_callbacks += near_radiation.preparation_callbacks();
  std::ostringstream near_prepared;
  near_prepared << "near-radiation-prepared status="
                << static_cast<int>(near_radiation.status())
                << " preparation_callbacks="
                << near_radiation.preparation_callbacks() << '\n';
  retain_identity(near_prepared);
  need(near_radiation.status() == S::ok,
       "positive-CDM near-radiation physical source admitted");
  const auto got = retained_native_control(
      near_radiation, directory, campaign_bytes, "near-radiation", wave, 1e-4);
  need(got.status == S::ok && got.rows.size() == 1 &&
           admitted_control_row(got.rows[0], false),
       "required positive-CDM near-radiation transfer earns native admission");
  peer::Budget budget;
  budget.campaign_bytes = &campaign_bytes;
  const auto independent =
      peer::campaign(near_radiation, wave[0], 1e-14, 1e-4,
                     (directory / "near-radiation-peer.txt").string(), budget);
  // Preserve a first failed reference before asserting its numerical status.
  peer::Budget record_budget;
  record_budget.campaign_bytes = &campaign_bytes;
  peer::Ledger record((directory / "near-radiation-peer-return.txt").string());
  std::ostringstream summary;
  summary << "near-radiation-peer status="
          << static_cast<int>(independent.status)
          << " diagnostics_partial=" << (independent.status == S::ok ? 0 : 1)
          << " rhs=" << budget.rhs << " scalars=" << budget.scalar_updates
          << " background=" << budget.background_queries
          << " age_evaluations=" << budget.age_quadrature_evaluations
          << " ledger_bytes=" << budget.ledger_bytes
          << " logging_refusals=" << budget.logging_refusals;
  for (std::size_t field = 0; field < independent.value.size(); ++field) {
    if (independent.status == S::ok)
      summary << " value" << field << '=' << independent.value[field];
    summary << " diagnostic_error" << field << '=' << independent.error[field];
    for (std::size_t term = 0; term < independent.contributions[field].size();
         ++term)
      summary << " diagnostic_contribution" << field << '_' << term << '='
              << independent.contributions[field][term];
  }
  summary << '\n';
  need(record.write(summary.str(), record_budget),
       "near-radiation peer refusal or result retained");
  need(independent.status == S::ok && budget.status == S::ok,
       "required near-radiation angular TRACE peer earns admission within "
       "original caps");
  const ReferencePoint near_reference{independent.value, independent.error};
  for (std::size_t field = 0; field < independent.value.size(); ++field) {
    const auto &row = got.rows[0];
    const std::array<const MasslessFDTransferValue *, 3> values{
        &row.comoving_cdm, &row.spatial_potential, &row.lapse_potential};
    const W epsilon = 1e-7L + 1e-4L * std::abs(W(*values[field]->value));
    for (W term : independent.contributions[field])
      need(std::isfinite(term) && term >= 0 && term <= epsilon / 30,
           "near-radiation peer keeps each original numerical share");
  }
  compare_reference_point(got.rows[0], near_reference);
}

void late_and_phase_controls(const Background &background,
                             const std::filesystem::path &directory,
                             std::size_t &campaign_bytes) {
  using namespace irred::cosmology;
  // One batch earns an a=1 trial and acquires the native Big-Bang-age witness
  // on an outside-phase row before any hierarchy trial for that row.
  const std::array<double, 2> late_waves{1e-7, .01};
  const auto late = retained_native_control(
      background, directory, campaign_bytes, "late-inclusive-and-phase-outside",
      late_waves, 1);
  need(late.status == S::ok && late.rows.size() == late_waves.size(),
       "a1 inclusive request retains both original rows");
  admitted_control_row(late.rows[0], true);
  const auto &outside = late.rows[1];
  need(outside.comoving_cdm.status == S::outside_domain &&
           outside.spatial_potential.status == S::outside_domain &&
           outside.lapse_potential.status == S::outside_domain &&
           !outside.comoving_cdm.value && !outside.spatial_potential.value &&
           !outside.lapse_potential.value && outside.attempts_recorded == 0 &&
           outside.work.rhs == 0 &&
           W(late_waves[1]) * (outside.conformal_age_mpc +
                               outside.conformal_age_error_mpc) >
               20,
       "actual outside-phase Big-Bang-age gate refuses before hierarchy "
       "attempts");
  const W eta = outside.conformal_age_mpc,
          error = outside.conformal_age_error_mpc;
  need(std::isfinite(eta) && std::isfinite(error) && eta > error && error > 0,
       "actual native phase witness has finite positive age uncertainty");
  peer::Budget age_budget;
  age_budget.campaign_bytes = &campaign_bytes;
  peer::Ledger phase_record(
      (directory / "phase-construction-and-independent-age.txt").string());
  std::ostringstream age_request;
  age_request << "independent-age-request status=started completion=unobserved"
              << " incomplete_disposition=interrupted "
                 "actual_counters=unavailable target_a=1"
              << " background_cap=" << age_budget.maximum_background_queries
              << " age_cap=" << age_budget.maximum_age_quadrature_evaluations
              << " scalars_cap=" << age_budget.maximum_scalar_updates << '\n';
  need(phase_record.write(age_request.str(), age_budget),
       "independent age pending request retained before setup and evaluation");
  need(age_budget.charge(age_budget.background_queries, 1,
                         age_budget.maximum_background_queries),
       "independent late-age radiation query charged");
  const auto radiation = background.scaled_expansion(0);
  age_budget.momentum_callbacks += radiation.callbacks;
  need(radiation.status == S::ok,
       "independent late-age retained radiation acquired");
  const auto independent_age = peer::endpoint_age(
      background, 1, radiation.a4_e2, radiation.error_estimate, age_budget);
  std::ostringstream age_record;
  age_record << std::setprecision(std::numeric_limits<W>::max_digits10)
             << "independent-age status="
             << static_cast<int>(independent_age.status)
             << " background=" << age_budget.background_queries
             << " age_evaluations=" << age_budget.age_quadrature_evaluations
             << " scalars=" << age_budget.scalar_updates
             << " momentum_callbacks=" << age_budget.momentum_callbacks;
  if (independent_age.status == S::ok)
    age_record << " eta=" << independent_age.value
               << " eta_error=" << independent_age.error;
  else
    age_record << " eta=unavailable eta_error=unavailable";
  age_record << '\n';
  need(phase_record.write(age_record.str(), age_budget),
       "independent age result/refusal and its actual counters retained");
  need(independent_age.status == S::ok &&
           std::abs(independent_age.value - eta) <=
               independent_age.error + error,
       "different endpoint-age algorithm agrees within retained age "
       "diagnostics");
  need(age_budget.scalars(12),
       "actual derived phase input construction charged");
  const double threshold = static_cast<double>(20 / (eta + error));
  const double inside = static_cast<double>(19 / (eta + error));
  const double ambiguous_inside = std::nextafter(threshold, 0.0);
  const double beyond =
      std::nextafter(threshold, std::numeric_limits<double>::infinity());
  const double uncertainty_crossing =
      std::nextafter(static_cast<double>(20 / eta), 0.0);
  std::ostringstream phase_inputs;
  phase_inputs
      << std::setprecision(std::numeric_limits<W>::max_digits10)
      << "phase-inputs native_eta=" << eta << " native_eta_error=" << error
      << " target_a=1 original_phase_limit=20\n"
      << "phase-mode index=0 role=physical-interior-19 k=" << inside << '\n'
      << "phase-mode index=1 role=central-upper-neighbor-inside-unresolved k="
      << ambiguous_inside << '\n'
      << "phase-mode index=2 role=central-upper-neighbor-outside k=" << beyond
      << '\n'
      << "phase-mode index=3 role=central-inside-uncertainty-outside k="
      << uncertainty_crossing << '\n';
  need(phase_record.write(phase_inputs.str(), age_budget),
       "each actual phase input and its selection reason retained before "
       "request");
  need(inside >= 1e-7 && beyond <= .01 && uncertainty_crossing <= .01 &&
           W(inside) * (eta + error) <= 20 &&
           W(ambiguous_inside) * (eta + error) <= 20 &&
           W(beyond) * (eta + error) > 20 && W(beyond) * error <= 1e-6L &&
           W(uncertainty_crossing) * eta < 20 &&
           W(uncertainty_crossing) * (eta + error) > 20,
       "binary64 adversarial phase controls straddle actual age upper bound");
  const std::array<double, 4> phases{inside, ambiguous_inside, beyond,
                                     uncertainty_crossing};
  const auto boundary =
      retained_native_control(background, directory, campaign_bytes,
                              "phase-boundary-and-age-uncertainty", phases, 1);
  need(boundary.status == S::ok && boundary.rows.size() == phases.size(),
       "phase controls retain actual original source axis order");
  admitted_control_row(boundary.rows[0], true);
  // The neighboring central inside mode can cross a matched input trial's
  // own uncertainty envelope. Keep that refusal explicitly unqualified.
  admitted_control_row(boundary.rows[1], true, true);
  for (std::size_t j : {std::size_t(2), std::size_t(3)}) {
    const auto &row = boundary.rows[j];
    need(row.wavenumber_mpc_inverse == phases[j] &&
             row.attempts_recorded == 0 && row.work.rhs == 0 &&
             row.comoving_cdm.status == S::outside_domain &&
             row.spatial_potential.status == S::outside_domain &&
             row.lapse_potential.status == S::outside_domain &&
             !row.comoving_cdm.value && !row.spatial_potential.value &&
             !row.lapse_potential.value,
         "outside/uncertainty-crossing phase controls refuse with no "
         "fabricated hierarchy value");
  }
  const std::array<double, 1> low_wave{1e-7};
  const auto too_late = retained_native_control(
      background, directory, campaign_bytes, "endpoint-above-inclusive-a1",
      low_wave, std::nextafter(1.0, std::numeric_limits<double>::infinity()));
  need(too_late.status == S::outside_domain && too_late.rows.empty() &&
           too_late.work.rhs == 0,
       "next binary64 endpoint above a1 refuses original physical domain");
}

} // namespace

int main() {
  try {
    const auto evidence = evidence_directory();
    std::cout << "massless FD peer evidence=" << evidence.string() << '\n';
    const auto background = irred::cosmology::prepare_thermal_background(
        {70, 0, 0, 0, .3, {{0, .0002, 2}}});
    need(background.status() == S::ok,
         "single acquired exact request background");
    std::size_t campaign_bytes = 0;
    characteristic_controls();
    nodal_transport_controls();
    translation_controls();
    initializer_controls();
    constraint_controls([](W g, W x, W kappa, W energy, W momentum,
                           W energy_force, W momentum_force) {
      const auto derivative = peer::constraint_rhs(
          g, x, kappa, energy, momentum, energy_force, momentum_force);
      return std::pair<W, W>{derivative[0], derivative[1]};
    });
    resource_controls(evidence);
    eds_equation_controls();
    full_trace_controls(background, evidence, campaign_bytes);
    const auto reference =
        native_comparison(background, evidence, campaign_bytes);
    species_and_epoch_controls(background, reference, evidence, campaign_bytes);
    late_and_phase_controls(background, evidence, campaign_bytes);
    std::cout << "PASS " << checks
              << " massless FD characteristic/nodal/trace controls\n";
  } catch (const std::exception &error) {
    std::cerr << "FAIL " << error.what() << '\n';
    return 1;
  }
}
