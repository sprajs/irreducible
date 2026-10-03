#pragma once
#include "../src/ideal_acoustic_bridge.hpp"
#include "../src/ideal_acoustic_initial_bounds.hpp"
#include "../src/ideal_acoustic_positive_radius.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace ideal_acoustic_transport_test {
namespace native = irred::cosmology;
namespace detail = irred::cosmology::detail;
using W = long double;
using S = irred::numerics::Status;

inline void require(bool good, const char *message) {
  if (!good) {
    std::cerr << "FAIL ideal acoustic transport: " << message << '\n';
    std::exit(1);
  }
}
inline void relative(W actual, W expected, const char *message) {
  const W allowance = 512 * std::numeric_limits<W>::epsilon() *
                      (std::abs(actual) + std::abs(expected));
  require(std::isfinite(actual) && std::isfinite(expected) &&
              std::abs(actual - expected) <= allowance,
          message);
}
inline void enclosed(W before, W after, W bound, const char *message) {
  // Only the two direct scalar-control evaluations/subtraction are covered by
  // this normal-Wide allowance. It is not a whole-trajectory output budget.
  const W allowance = 512 * std::numeric_limits<W>::epsilon() *
                      (std::abs(before) + std::abs(after));
  require(std::isfinite(before) && std::isfinite(after) &&
              std::isfinite(bound) && bound >= 0 &&
              std::abs(after - before) <= bound + allowance,
          message);
}
struct Counters {
  std::size_t writes = 0, diagnostics = 0;
  static bool write(void *opaque, std::size_t n) noexcept {
    auto &self = *static_cast<Counters *>(opaque);
    if (n > 4096 || self.writes > 4096 - n)
      return false;
    self.writes += n;
    return true;
  }
  static bool diagnostic(void *opaque) noexcept {
    auto &self = *static_cast<Counters *>(opaque);
    if (self.diagnostics == 4096)
      return false;
    ++self.diagnostics;
    return true;
  }
  detail::IdealAcousticTransportAccounting accounting() noexcept {
    return {this, write, diagnostic};
  }
};

inline void signed_rhs_fusion_controls(
    const detail::IdealAcousticTransportFrame &frame,
    const detail::IdealAcousticSourceUncertainty &u) {
  namespace arithmetic = detail::ideal_acoustic_transport_internal;
  // This is an affected-consumer association/ownership check through the
  // existing vector interfaces, not an independent physical algorithm.
  struct RhsCounter {
    std::size_t cap, writes=0, diagnostics=0;
    static bool write(void *opaque,std::size_t n) noexcept {
      auto &self=*static_cast<RhsCounter *>(opaque);
      if (n>self.cap || self.writes>self.cap-n) return false;
      self.writes+=n;
      return true;
    }
    static bool diagnostic(void *opaque) noexcept {
      ++static_cast<RhsCounter *>(opaque)->diagnostics;
      return true;
    }
    detail::IdealAcousticTransportAccounting accounting() noexcept {
      return {this,write,diagnostic};
    }
  };
  const std::array<W,6> y{.25L,-.125L,.0625L,-.5L,.75L,1};
  std::array<W,24> z;
  for (unsigned j=0;j<4;++j)
    for (unsigned i=0;i<6;++i)
      z[6*j+i]=(j%2 ? -1 : 1)*W((i+1)*(j+1))/32;
  for (unsigned i=0;i<6;++i) z[i]=y[i];
  for (unsigned fixture=0;fixture<3;++fixture) {
    auto selected=frame;
    if (fixture) {
      // Synthetic opposite operator produces exact signed cancellation in
      // S and dV, and in Delta when the separately retained defect is zero.
      // It is not a new physical source or a source-accuracy certificate.
      const auto &e=selected.actual;
      selected.directions.gradient[0]={-e.x2,-e.F,-e.B,-e.L,
          -e.loading_over_one_plus_loading,-e.sound_speed_squared};
      selected.actual.acceleration_defect=fixture==1 ? 0 : .125L;
    }
    std::array<W,24> expected,dz;
    std::array<W,6> tau;
    arithmetic::Arithmetic reference;
    arithmetic::Bands bands;
    require(arithmetic::bands(selected,u,bands,reference)==S::ok,
            "fusion reference retains admitted diagnostic bands");
    for (unsigned j=0;j<4;++j) {
      std::array<W,5> partial;
      detail::ideal_acoustic_derivative(
          std::span<const W,5>(z.data()+6*j,5),selected.actual,partial);
      const auto force=detail::ideal_acoustic_source_force(
          std::span<const W,5>(y.data(),5),selected.directions.gradient[j]);
      for (unsigned i=0;i<5;++i)
        expected[6*j+i]=reference.add(partial[i],force[i]);
      expected[6*j+5]=-reference.div(reference.mul(bands.inverse_h,bands.pj[j]),
                                    reference.mul(2,selected.center.shadow_p));
    }
    require(reference.status==S::ok,"materialized reference uses normal arithmetic");
    RhsCounter exact{54};
    require(detail::ideal_acoustic_transport_rhs(y,z,selected,u,dz,tau,
                exact.accounting())==S::ok && exact.writes==54 &&
                exact.diagnostics==1,
            "fused signed RHS stores exactly its 54 actual destinations");
    for (unsigned i=0;i<24;++i)
      require(std::isfinite(dz[i]) && dz[i]==expected[i] &&
                  std::signbit(dz[i])==std::signbit(expected[i]),
              "fused final signed derivative preserves vector-route Wide result");
    for (W v:tau)
      require(arithmetic::radius(v),"original literal RHS allowances remain available");
    if (fixture)
      require(dz[1]==0 && dz[2]==0 &&
                  (fixture==1 ? dz[0]==0 : dz[0]!=0),
              "cancellation retains the separately supplied acceleration defect");
  }
  std::array<W,24> dz;
  std::array<W,6> tau;
  dz.fill(7); tau.fill(7);
  const auto untouched_dz=dz;
  const auto untouched_tau=tau;
  RhsCounter before_first_store{5};
  require(detail::ideal_acoustic_transport_rhs(y,z,frame,u,dz,tau,
              before_first_store.accounting())==S::work_limit &&
              before_first_store.writes==0 && dz==untouched_dz &&
              tau==untouched_tau,
          "initial six-slot reservation refuses before mutation");
  RhsCounter one_short{53};
  require(detail::ideal_acoustic_transport_rhs(y,z,frame,u,dz,tau,
              one_short.accounting())==S::work_limit && one_short.writes==48,
          "one-short RHS cap retains actual prefix and refuses final tau update");
  dz=untouched_dz; tau=untouched_tau;
  auto invalid=y;
  invalid[2]=std::numeric_limits<W>::quiet_NaN();
  RhsCounter invalid_input{54};
  require(detail::ideal_acoustic_transport_rhs(invalid,z,frame,u,dz,tau,
              invalid_input.accounting())==S::outside_domain &&
              invalid_input.writes==0 && invalid_input.diagnostics==0 &&
              dz==untouched_dz && tau==untouched_tau,
          "nonnormal state refuses before work or output mutation");
}

inline void signed_radius_assembly_controls() {
  namespace arithmetic = detail::ideal_acoustic_transport_internal;
  // Powers of two make these exact real sums independent of RK/source
  // ancestry. Exercise the same completion operation used by native combine.
  const W r = .125L;
  const W infinity = std::numeric_limits<W>::infinity();
  arithmetic::Arithmetic adequate;
  const W positive = adequate.assemble_signed_radius(-r, 2*r);
  require(adequate.status == S::ok && arithmetic::radius(positive) &&
              positive >= r && positive <= std::nextafter(r, infinity),
          "signed provisional plus adequate allowance encloses its exact sum");
  arithmetic::Arithmetic inadequate;
  const W negative = inadequate.assemble_signed_radius(-r, r/2);
  require(inadequate.status == S::ok && negative < 0 &&
              !arithmetic::radius(negative),
          "completed negative radius remains refused rather than clipped");
  arithmetic::Arithmetic cancellation;
  require(cancellation.assemble_signed_radius(-r, r) == 0 &&
              cancellation.status == S::ok,
          "exact signed cancellation retains zero");
  const std::array<std::array<W,2>,3> positive_pairs{{{r,r}, {r,0}, {0,r}}};
  for (const auto &operands : positive_pairs) {
    arithmetic::Arithmetic nonnegative;
    const W bound = nonnegative.assemble_signed_radius(operands[0], operands[1]);
    require(nonnegative.status == S::ok && arithmetic::radius(bound) &&
                bound >= operands[0]+operands[1],
            "positive and zero-operand sums stay outward and admissible");
  }
  arithmetic::Arithmetic zero;
  require(zero.assemble_signed_radius(0,0) == 0 && zero.status == S::ok,
          "zero allowance and provisional stay exactly zero");
  for (W invalid : {std::numeric_limits<W>::quiet_NaN(),
                    std::numeric_limits<W>::infinity()}) {
    arithmetic::Arithmetic nonnormal;
    (void)nonnormal.assemble_signed_radius(invalid, r);
    require(nonnormal.status == S::outside_domain,
            "nonnormal signed addition cannot become admitted through zero");
  }
  arithmetic::Arithmetic invalid_allowance;
  (void)invalid_allowance.assemble_signed_radius(r,-r);
  require(invalid_allowance.status == S::outside_domain,
          "assembly allowance itself must be nonnegative");
  if (std::numeric_limits<W>::denorm_min() > 0) {
    const W minimum = std::numeric_limits<W>::min();
    arithmetic::Arithmetic subnormal;
    (void)subnormal.assemble_signed_radius(
        -std::nextafter(minimum, infinity), minimum);
    require(subnormal.status == S::outside_domain,
            "produced subnormal signed sum retains domain refusal");
  }
}
inline void radius_failure_record_controls() {
  namespace arithmetic = detail::ideal_acoustic_transport_internal;
  const W provisional=-.125L, allowance=.0625L, completed=-.0625L;
  const auto negative=arithmetic::record_radius_failure(
      2,native::IdealAcousticRadiusChannel::arithmetic,provisional,allowance,
      completed,S::ok);
  require(negative.coordinate==native::IdealAcousticStateCoordinate::scaled_common_minus_cdm_velocity &&
              negative.channel==native::IdealAcousticRadiusChannel::arithmetic &&
              negative.provisional_radius==provisional &&
              negative.local_assembly_allowance==allowance &&
              negative.completed_radius==completed && negative.arithmetic_status==S::ok &&
              negative.completed_finite && negative.completed_normal_or_zero &&
              !negative.completed_nonnegative,
          "failed signed radius and allowance remain literal without clipping");
  const auto source=arithmetic::record_radius_failure(
      5,native::IdealAcousticRadiusChannel::source,0,std::nullopt,0,S::outside_domain);
  require(source.coordinate==native::IdealAcousticStateCoordinate::conformal_age &&
              source.channel==native::IdealAcousticRadiusChannel::source &&
              source.provisional_radius==0 && !source.local_assembly_allowance &&
              source.completed_radius==0 && source.arithmetic_status==S::outside_domain &&
              source.completed_finite && source.completed_normal_or_zero &&
              source.completed_nonnegative,
          "source has no donated assembly allowance and preserves earlier arithmetic status");
  for (W invalid : {std::numeric_limits<W>::quiet_NaN(),
                    std::numeric_limits<W>::infinity()}) {
    const auto nonfinite=arithmetic::record_radius_failure(
        4,native::IdealAcousticRadiusChannel::arithmetic,invalid,invalid,invalid,S::outside_domain);
    require(!nonfinite.provisional_radius && !nonfinite.local_assembly_allowance &&
                !nonfinite.completed_radius && !nonfinite.completed_finite &&
                !nonfinite.completed_normal_or_zero &&
                nonfinite.arithmetic_status==S::outside_domain,
            "nonfinite diagnostic scalars remain absent for portable records");
  }
  if (std::numeric_limits<W>::denorm_min()>0) {
    const W tiny=std::numeric_limits<W>::denorm_min();
    const auto subnormal=arithmetic::record_radius_failure(
        1,native::IdealAcousticRadiusChannel::source,tiny,std::nullopt,tiny,S::outside_domain);
    require(subnormal.provisional_radius==tiny && subnormal.completed_radius==tiny &&
                subnormal.completed_finite && !subnormal.completed_normal_or_zero &&
                subnormal.completed_nonnegative,
            "finite subnormal witness remains visible without admission");
  }
}

// Direct fixed-source scalar equations for test comparisons only. They do not
// replace the retained production H or use its momentum perturbation law.
struct Primitive {
  W photon, baryon, cdm, vacuum;
};
struct RationalControl {
  std::array<W, 6> coefficient; // x2,F,B,L,T,cs2
  W inverse_h;
};
inline RationalControl rational(Primitive v, W a, W k, W alpha) {
  const W a2 = a * a, a4 = a2 * a2, matter = (v.baryon + v.cdm) * a;
  const W P = v.photon + matter + v.vacuum * a4;
  const W D = 4 * v.photon + 3 * v.baryon * a;
  require(a > 0 && P > 0 && D > 0 && v.photon > 0 && v.baryon > 0 &&
              v.cdm > 0 && v.vacuum > 0,
          "direct positive primitive family");
  const W ka = k * a / alpha;
  return {{{ka * ka / P, (matter + 4 * v.photon / 3) / P,
            (v.baryon * a + 4 * v.photon / 3) / P,
            (matter + 4 * v.vacuum * a4) / (2 * P) - 1, 3 * v.baryon * a / D,
            4 * v.photon / (3 * D)}},
          a / (alpha * std::sqrt(P))};
}

inline void radiation_source_limit() {
  // Formal endpoint, not an admitted baryon-free Thomson background. The
  // physical Gamma direction retains its correlated -dGamma vacuum shift.
  const W a = 1e-8L, gamma = 1e-4L, a4 = a * a * a * a, x2 = 1e-7L;
  const detail::IdealAcousticSourceCenter center{
      a, gamma, 0, 0, gamma, 0, x2, 4.L / 3, 4.L / 3, -1, 0, 1.L / 3};
  const auto directions = detail::ideal_acoustic_source_directions(center);
  require(directions.status == S::ok, "formal radiation directions");
  relative(directions.gradient[0].F, (4.L / 3) * a4 / gamma,
           "Gamma F cancellation leaves vacuum support");
  relative(directions.gradient[0].B, (4.L / 3) * a4 / gamma,
           "Gamma B cancellation leaves vacuum support");
  relative(directions.gradient[0].L, -2 * a4 / gamma,
           "Gamma L cancellation leaves vacuum support");
  const std::array<W, 5> y{x2 / 4, 0, x2 / 36, -1.L / 3, -2.L / 3};
  const auto force =
      detail::ideal_acoustic_source_force(y, directions.gradient[0]);
  relative(force[3], 2 * a4 / (3 * gamma),
           "regular radiation velocity source has no arbitrary growth force");
  require(std::abs(force[3]) < a4 / gamma,
          "radiation source retains a4 rather than an independent L maximum");
}

inline void initial_controls(const detail::IdealAcousticTransportFrame &frame,
                             const detail::IdealAcousticSourceUncertainty &u,
                             Counters &counter) {
  const W p0 = -2.L / 3, x2 = frame.actual.x2;
  const W m = (frame.center.baryon + frame.center.cdm) * frame.center.a /
              frame.center.photon;
  const W phi = p0 * (1 - m / 16 + 7 * m * m / 160 - x2 / 30);
  const W vc = p0 * (.5L + m / 16 - 7 * m * m / 160 - x2 / 120);
  const W entropy = -p0 * x2 * x2 / 96, dv = -p0 * x2 / 24;
  const W delta =
      (-2 * x2 * phi / 3 + frame.actual.B * entropy - 3 * frame.actual.B * dv) /
      frame.actual.F;
  const std::array<W, 6> y{delta, entropy, dv, vc, phi, 1};
  const auto original = y;
  detail::IdealAcousticTransportInitialBounds seed;
  seed.status = S::ok;
  // These are explicitly supplied test enclosures, not a claim that the
  // production initial/early-age owner has earned these six coordinates.
  seed.source_family_radius.fill(1e-5L);
  seed.arithmetic_radius.fill(1e-12L);
  detail::IdealAcousticResponseState response{};
  require(detail::ideal_acoustic_transport_initial(
              y, frame, u, seed, response, counter.accounting()) == S::ok,
          "explicit full initial seed accepted");
  require(y == original,
          "response seed does not reproject or renormalize mode");
  for (unsigned j = 0; j < 4; ++j) {
    const W xj = frame.directions.gradient[j].x2;
    relative(response[6 * j + 1], -p0 * x2 * xj / 48,
             "initial entropy derivative retains x4 support");
    relative(response[6 * j + 2], -p0 * xj / 24,
             "initial relative velocity derivative retains x2 support");
    const W wave = -(2 * (xj * phi + x2 * response[6 * j + 4])) / 3;
    const W density = frame.directions.gradient[j].B * entropy +
                      frame.actual.B * response[6 * j + 1];
    const W momentum = -3 * (frame.directions.gradient[j].B * dv +
                             frame.actual.B * response[6 * j + 2]);
    const W denominator = -frame.directions.gradient[j].F * delta;
    const W projected_derivative = wave + density + momentum + denominator;
    const W lhs = frame.actual.F * response[6 * j];
    // Vacuum support cancels at leading radiation order. The independently
    // grouped scalar residual uses term magnitudes, not its small result.
    const W allowance = 512 * std::numeric_limits<W>::epsilon() *
                        (std::abs(lhs) + std::abs(wave) + std::abs(density) +
                         std::abs(momentum) + std::abs(denominator));
    require(std::isfinite(lhs) && std::isfinite(projected_derivative) &&
                std::abs(lhs - projected_derivative) <= allowance,
            "Delta-only projection has its full quotient derivative");
    require(response[6 * j + 5] == 0,
            "initial eta source is in its supplied full remainder");
  }
  for (unsigned i = 0; i < 6; ++i)
    require(response[24 + i] >= seed.source_family_radius[i] &&
                response[30 + i] == seed.arithmetic_radius[i],
            "all initial family/arithmetic coordinates retained");
  seed.status = S::invalid_input;
  require(detail::ideal_acoustic_transport_initial(
              y, frame, u, seed, response, counter.accounting()) != S::ok,
          "missing full initial seed refuses");
  seed.status = S::ok;
  seed.arithmetic_radius[4] = 0;
  require(detail::ideal_acoustic_transport_initial(
              y, frame, u, seed, response, counter.accounting()) != S::ok,
          "nonzero initial potential cannot receive absent arithmetic bound");
}

// Independently factored scalar initial series and projection. This does not
// call the bound helper's Map/variation/Graph and is not a trajectory oracle.
inline std::array<W, 5> direct_initial(Primitive p, W a,
                                       const std::array<W, 6> &coefficient) {
  const W unit = -2.L / 3;
  const W m = ((p.baryon + p.cdm) / p.photon) * a;
  const W q = coefficient[0];
  const W phi = unit * (1 + m * (-1.L / 16 + 7 * m / 160) - q / 30);
  const W vc = unit * (.5L + m * (1.L / 16 - 7 * m / 160) - q / 120);
  const W entropy = (-unit / 96) * q * q;
  const W dv = (-unit / 24) * q;
  const W numerator =
      (-2.L / 3) * (q * phi) + coefficient[2] * (entropy - 3 * dv);
  return {numerator / coefficient[1], entropy, dv, vc, phi};
}

inline void finite_initial_family_controls(
    Primitive nominal, W k, W alpha,
    const detail::IdealAcousticTransportFrame &fine,
    const detail::IdealAcousticTransportFrame &broad,
    const detail::IdealAcousticSourceUncertainty &family, W fine_clock_radius,
    W broad_clock_radius, Counters &counter) {
  const W a = fine.center.a, p0 = -2.L / 3, q = fine.actual.x2;
  const W m = (nominal.baryon + nominal.cdm) * a / nominal.photon;
  std::array<W, 6> actual{-3 * p0 * q / 8,
                          -p0 * q * q / 96,
                          -p0 * q / 24,
                          p0 * (.5L + m / 16 - 7 * m * m / 160 - q / 120),
                          p0 * (1 - m / 16 + 7 * m * m / 160 - q / 30),
                          1};
  actual[0] = (-(2.L / 3) * q * actual[4] + fine.actual.B * actual[1] -
               3 * fine.actual.B * actual[2]) /
              fine.actual.F;
  const auto original = actual;
  // Positive test-only age estimates: no integral-age accuracy is asserted by
  // this scalar map comparison. The production caller owns the actual age.
  constexpr W age_source = 1e-4L, age_arithmetic = 1e-10L;
  detail::IdealAcousticTransportInitialBounds fine_seed, broad_seed;
  const auto writes_before = counter.writes;
  const auto diagnostics_before = counter.diagnostics;
  require(detail::ideal_acoustic_initial_bounds(
              actual, fine, family, nominal.vacuum, age_source, age_arithmetic,
              fine_seed, counter.accounting()) == S::ok &&
              detail::ideal_acoustic_initial_bounds(
                  actual, broad, family, nominal.vacuum, age_source,
                  age_arithmetic, broad_seed, counter.accounting()) == S::ok,
          "actual borrowed initial map has finite source and clock bounds");
  require(actual == original && counter.writes - writes_before == 24 &&
              counter.diagnostics - diagnostics_before == 2,
          "initial bounds preserve fields and charge twelve outputs per call");
  require(fine_seed.source_family_radius == broad_seed.source_family_radius,
          "initial SOURCE family does not acquire a second clock owner");
  for (unsigned i = 0; i < 6; ++i)
    require(std::isfinite(fine_seed.source_family_radius[i]) &&
                fine_seed.source_family_radius[i] > 0 &&
                std::isfinite(fine_seed.arithmetic_radius[i]) &&
                fine_seed.arithmetic_radius[i] > 0 &&
                broad_seed.arithmetic_radius[i] >=
                    fine_seed.arithmetic_radius[i],
            "all six initial source/arithmetic coordinates remain available");
  require(fine_seed.source_family_radius[5] == age_source &&
              broad_seed.source_family_radius[5] == age_source &&
              fine_seed.arithmetic_radius[5] >= age_arithmetic &&
              broad_seed.arithmetic_radius[5] > fine_seed.arithmetic_radius[5],
          "eta keeps supplied source age and positive full-family clock "
          "arithmetic");

  const auto nominal_rational = rational(nominal, a, k, alpha);
  const auto nominal_initial =
      direct_initial(nominal, a, nominal_rational.coefficient);
  for (unsigned i = 0; i < 5; ++i)
    enclosed(actual[i], nominal_initial[i], fine_seed.arithmetic_radius[i],
             "actual rounded map minus direct smooth nominal is arithmetic");
  for (unsigned corner = 0; corner < 16; ++corner) {
    std::array<W, 4> delta;
    for (unsigned j = 0; j < 4; ++j)
      delta[j] =
          (corner & (1u << j) ? 1 : -1) *
          (std::abs(family.signed_shift[j]) + family.operation_radius[j]);
    const Primitive shifted{nominal.photon + delta[0],
                            nominal.baryon + delta[1], nominal.cdm + delta[2],
                            nominal.vacuum + delta[3] - delta[0] - delta[1] -
                                delta[2]};
    const auto source_rational = rational(shifted, a, k, alpha);
    const auto source_initial =
        direct_initial(shifted, a, source_rational.coefficient);
    for (unsigned i = 0; i < 5; ++i)
      enclosed(
          nominal_initial[i], source_initial[i],
          fine_seed.source_family_radius[i],
          "all projected initial fields enclose correlated source corners");
    for (unsigned clock = 0; clock < 2; ++clock) {
      const auto &seed = clock ? broad_seed : fine_seed;
      const W radius = clock ? broad_clock_radius : fine_clock_radius;
      require(radius > 0 && radius < .01L,
              "initial direct clock controls retain finite log support");
      const W rho = .99L * radius;
      for (W sign : {-1.L, 1.L}) {
        const W shifted_a = a * (1 + sign * rho);
        const auto clock_rational = rational(shifted, shifted_a, k, alpha);
        const auto clock_initial =
            direct_initial(shifted, shifted_a, clock_rational.coefficient);
        for (unsigned i = 0; i < 5; ++i) {
          enclosed(
              source_initial[i], clock_initial[i], seed.arithmetic_radius[i],
              "initial arithmetic encloses both clocks across source family");
          enclosed(
              actual[i], clock_initial[i],
              seed.source_family_radius[i] + seed.arithmetic_radius[i],
              "initial SOURCE plus ARITHMETIC encloses direct joint change");
        }
      }
    }
  }

  detail::IdealAcousticResponseState response{};
  require(detail::ideal_acoustic_transport_initial(
              actual, fine, family, fine_seed, response,
              counter.accounting()) == S::ok,
          "finite initial helper seeds the actual response consumer");
  for (unsigned j = 0; j < 4; ++j)
    require(response[6 * j + 5] == 0,
            "initial eta has no invented signed sensitivity");
  require(
      response[29] >= age_source &&
          response[35] == fine_seed.arithmetic_radius[5],
      "initial eta age and arithmetic stay in their declared remainder slots");

  detail::IdealAcousticTransportInitialBounds refused;
  require(detail::ideal_acoustic_initial_bounds(
              actual, fine, family, nominal.vacuum, 0, age_arithmetic, refused,
              counter.accounting()) != S::ok &&
              detail::ideal_acoustic_initial_bounds(
                  actual, fine, family, nominal.vacuum, age_source, 0, refused,
                  counter.accounting()) != S::ok,
          "each missing positive initial-age owner refuses");
  auto singular = family;
  singular.operation_radius[0] = nominal.photon;
  require(detail::ideal_acoustic_initial_bounds(
              actual, fine, singular, nominal.vacuum, age_source,
              age_arithmetic, refused, counter.accounting()) != S::ok,
          "initial singular primitive family refuses");
  auto no_clock = fine;
  no_clock.arithmetic.coefficient_radius[0] =
      no_clock.arithmetic.coefficient_at_actual_a_radius[0];
  require(detail::ideal_acoustic_initial_bounds(
              actual, no_clock, family, nominal.vacuum, age_source,
              age_arithmetic, refused, counter.accounting()) != S::ok,
          "initial unearned wave clock cannot become an exact zero bound");
  no_clock = fine;
  no_clock.arithmetic.conditional_wide_clock_profile = false;
  require(detail::ideal_acoustic_initial_bounds(
              actual, no_clock, family, nominal.vacuum, age_source,
              age_arithmetic, refused, counter.accounting()) != S::ok,
          "initial missing conditional profile refuses");
  refused.source_family_radius.fill(7);
  refused.arithmetic_radius.fill(11);
  const auto untouched_source = refused.source_family_radius;
  const auto untouched_arithmetic = refused.arithmetic_radius;
  Counters denied;
  denied.writes = 4096;
  require(
      detail::ideal_acoustic_initial_bounds(
          actual, fine, family, nominal.vacuum, age_source, age_arithmetic,
          refused, denied.accounting()) == S::work_limit &&
          refused.status == S::work_limit && denied.writes == 4096 &&
          denied.diagnostics == 1 &&
          refused.source_family_radius == untouched_source &&
          refused.arithmetic_radius == untouched_arithmetic,
      "initial result-write refusal precedes all twelve output assignments");
  denied = {};
  denied.diagnostics = 4096;
  require(detail::ideal_acoustic_initial_bounds(
              actual, fine, family, nominal.vacuum, age_source, age_arithmetic,
              refused, denied.accounting()) == S::work_limit &&
              denied.writes == 0 && denied.diagnostics == 4096 &&
              refused.source_family_radius == untouched_source &&
              refused.arithmetic_radius == untouched_arithmetic,
          "initial diagnostic refusal precedes graph and result mutation");
}

inline void endpoint_controls(const detail::IdealAcousticTransportFrame &frame,
                              const detail::IdealAcousticSourceUncertainty &u,
                              Counters &counter) {
  // Independent constant-coefficient exponential series is an algorithm
  // control for the finite-flow bound, not another cosmological trajectory.
  const W d=1e-4L;
  const std::array<W,6> y{.2L,.03L,-.01L,-.3L,-.6L,1};
  detail::IdealAcousticResponseState response{};
  const auto before_writes=counter.writes,before_diagnostics=counter.diagnostics;
  require(detail::ideal_acoustic_transport_fixed_a(y,response,frame,u,d,
              counter.accounting())==S::ok,
          "finite whole-family endpoint flow admitted");
  require(counter.writes-before_writes==6 &&
              counter.diagnostics-before_diagnostics==1,
          "endpoint correction owns six radii and one diagnostic");
  const auto c=frame.actual;
  const auto derivative=[&](const std::array<W,5> &v) {
    return std::array<W,5>{-c.x2*v[3]+4.5L*c.B*v[2]+
                              3*c.acceleration_defect*v[3],
        c.x2*v[2],(c.L-c.loading_over_one_plus_loading)*v[2]+
            c.sound_speed_squared*(v[0]-v[1]),
        (c.L-1)*v[3]+v[4],-v[4]+1.5L*c.F*v[3]+1.5L*c.B*v[2]};
  };
  std::array<W,5> term{y[0],y[1],y[2],y[3],y[4]},sum=term;
  for (unsigned n=1;n<=6;++n) {
    term=derivative(term);
    for (unsigned i=0;i<5;++i) { term[i]*=d/n; sum[i]+=term[i]; }
  }
  // Nominal constant matrix norm <=8 in this actual initial radiation slice.
  // A geometric remainder after six terms bounds its exponential tail.
  require(c.x2+4.5L*c.B+3*std::abs(c.acceleration_defect)<8 && c.x2<8 &&
              std::abs(c.L)+c.loading_over_one_plus_loading+
                  2*c.sound_speed_squared<8 &&
              std::abs(c.L)+2<8 && 1+1.5L*(c.F+c.B)<8,
          "constant-matrix analytic control norm");
  W tail=.6L;
  for (unsigned n=0;n<7;++n) tail*=8*d;
  tail/=1-8*d;
  for (unsigned i=0;i<5;++i)
    require(std::abs(sum[i]-y[i])+tail<=response[30+i],
            "independent exponential control and its tail fit endpoint radius");
  require(d/frame.hcal<=response[35] && response[35]>0,
          "inhomogeneous eta endpoint has its separate unit allowance");
  for (unsigned i=0;i<30;++i)
    require(response[i]==0,"endpoint flow does not rewrite source or mode");
  auto refused=response;
  require(detail::ideal_acoustic_transport_fixed_a(y,refused,frame,u,1,
              counter.accounting())!=S::ok && refused==response,
          "singular finite endpoint support refuses without radius mutation");
  require(detail::ideal_acoustic_transport_fixed_a(y,refused,frame,u,0,
              counter.accounting())!=S::ok && refused==response,
          "absent endpoint support cannot become zero error");
}

inline void positive_transport_controls(
    const detail::IdealAcousticTransportFrame &frame,
    const detail::IdealAcousticSourceUncertainty &u,Counters &counter) {
  const std::array<W,6> y{.125L,.0625L,.03125L,.25L,.5L,1};
  detail::IdealAcousticResponseState response{};
  for (unsigned i=0;i<6;++i) { response[24+i]=1e-6L; response[30+i]=1e-12L; }
  detail::IdealAcousticRadiusDiagonal diagonal;
  require(detail::ideal_acoustic_radius_diagonal(frame.actual,diagonal,
              counter.accounting())==S::ok,"actual borrowed radius diagonal");
  detail::IdealAcousticRadiusVector p0,p1,bridge,action;
  require(detail::ideal_acoustic_transport_forcing(y,response,frame,u,
              diagonal.positive,p0,counter.accounting())==S::ok,
          "complete positive forcing retains inherited and affine owners");
  require(p0[5]>0 && p0[11]>0,"global source and arithmetic eta forces remain positive");
  constexpr W bend=.125L,h=.01L;
  std::array<W,6> q{};
  q[3]=bend;
  auto final_y=y;
  final_y[3]+=bend; // Exact binary addition: only FINAL central assembly changes.
  require(detail::ideal_acoustic_transport_forcing(final_y,response,frame,u,
              diagonal.positive,p1,counter.accounting())==S::ok,
          "matching endpoint forcing observes final central pulse");
  require(detail::ideal_acoustic_transport_endpoint_bridge(q,frame,u,h,
              diagonal.mu_lower,bridge,counter.accounting())==S::ok,
          "finite PC endpoint source bridge has its own actual denominator");
  // phi diagonal is exactly -1. The independent endpoint partial difference
  // is scaled by h/2/(1+h/2); no full-step continuous-bound claim is made.
  const W direct=(h/2)*std::abs(p1[4]-p0[4])/(1+h/2);
  require(direct>0 && bridge[4]>0,
          "final bend changes endpoint SOURCE; omitted bridge returns zero");
  enclosed(0,direct,bridge[4],"PC endpoint source partial is enclosed by owned bridge");
  detail::IdealAcousticRadiusVector pure_q{};
  pure_q[9]=bend;
  require(detail::ideal_acoustic_transport_local_action(pure_q,frame,u,action,
              counter.accounting())==S::ok && action[4]>0 && action[10]>0,
          "SOURCE cross sees an intermediate fresh arithmetic pulse");
  require(action[5]==0 && action[11]==0 && bridge[5]==0 && bridge[11]==0,
          "fixed-frame eta has zero local gain without removing affine forcing");
  std::array<W,6> zero{};
  require(detail::ideal_acoustic_transport_endpoint_bridge(zero,frame,u,h,
              diagonal.mu_lower,bridge,counter.accounting())==S::ok,
          "exact zero pulse has no endpoint bridge");
  for (W v:bridge) require(v==0,"zero bridge stays exact zero");
  detail::IdealAcousticRadiusVector output{};
  output.fill(7);
  auto untouched=output;
  q[3]=std::numeric_limits<W>::quiet_NaN();
  require(detail::ideal_acoustic_transport_endpoint_bridge(q,frame,u,h,
              diagonal.mu_lower,output,counter.accounting())==S::outside_domain &&
              output==untouched,"invalid pulse refuses before output mutation");
  q[3]=bend;
  Counters denied;
  denied.writes=4096;
  require(detail::ideal_acoustic_transport_endpoint_bridge(q,frame,u,h,
              diagonal.mu_lower,output,denied.accounting())==S::work_limit &&
              output==untouched && denied.writes==4096,
          "bridge output reservation refuses without hidden writes");
}

inline void controls() {
  signed_radius_assembly_controls();
  radius_failure_record_controls();
  radiation_source_limit();
  const auto mapping =
      native::map_thermal_physical_model({70, .02, .10, 2.7, 0, {}});
  require(mapping.status == S::ok && mapping.model && mapping.scalar_witnesses,
          "one actual no-species physical map");
  const auto background = native::prepare_thermal_background(*mapping.model);
  const auto retained =
      detail::ThermalRetainedCoefficientAccess::capture(background);
  require(background.status() == S::ok && retained.has_value(),
          "one actual background preparation and retained capture");
  const auto loading = detail::thermal_baryon_loading(background);
  require(loading.status == S::ok, "same retained baryon loading owner");
  const auto &other = (*mapping.scalar_witnesses)[3];
  require(other.wide_value == 0 && other.emitted_value == 0 &&
              other.measured_absolute_cast_loss == 0 &&
              other.wide_operation_estimate == 0,
          "original fourth mapper coordinate is exactly zero");
  const Primitive nominal{
      W(mapping.model->omega_gamma), W(mapping.model->omega_b),
      W(mapping.model->omega_cdm), retained->lambda_retained};
  detail::IdealAcousticSourceUncertainty family;
  family.status = S::ok;
  family.original_other_zero_witness = true;
  const std::array<W, 4> magnitude{nominal.photon, nominal.baryon, nominal.cdm,
                                   nominal.vacuum};
  for (unsigned j = 0; j < 4; ++j) {
    // Deliberately synthetic finite primitive box for the algebra comparison;
    // these do not replace actual deterministic mapper uncertainties.
    family.signed_shift[j] = (j % 2 ? -1 : 1) * magnitude[j] * .01L;
    family.operation_radius[j] = magnitude[j] * .04L;
  }
  const W alpha =
      W(background.source().h0_km_s_mpc) / detail::thermal_conformal_c_km_s;
  Counters counter;
  constexpr W k = .01L;
  unsigned stage_index = 0;
  for (W a : {1e-10L, 1e-3L, W(.01)}) {
    auto epoch = detail::thermal_conformal_epoch(
        background, a, k, nominal.photon, native::ThermalPolicy{});
    require(detail::thermal_conformal_no_species_diagnostics(
                background, *retained, *mapping.scalar_witnesses, a, epoch) ==
                S::ok,
            "borrowed same-query forward and shadow diagnostics");
    const auto original_p = epoch.p, original_h = epoch.hcal;
    const auto callbacks = epoch.raw_scaled_query->callbacks;
    const auto sound = detail::thermal_baryon_sound(loading, a);
    require(sound.status == S::ok, "same retained loading/sound epoch");
    const detail::IdealAcousticCoefficients actual{
        epoch.x2,
        epoch.g,
        epoch.enthalpy_fraction,
        epoch.fb + 4 * epoch.fg / 3,
        sound.loading / sound.denominator,
        sound.sound_speed_squared_over_c_squared,
        epoch.acceleration_defect};
    const detail::IdealAcousticSourceCenter center{
        a,
        nominal.photon,
        nominal.baryon,
        nominal.cdm,
        epoch.shadow->p_shadow,
        epoch.shadow->p_shadow_n,
        actual.x2,
        actual.F,
        actual.B,
        actual.L,
        actual.loading_over_one_plus_loading,
        actual.sound_speed_squared};
    const auto directions = detail::ideal_acoustic_source_directions(center);
    const W N = std::log(a);
    const detail::IdealAcousticBridgeClock fine_clock{
        S::ok, N, 512 * std::numeric_limits<W>::epsilon() * (1 + std::abs(N)),
        true};
    auto broad_clock = fine_clock;
    broad_clock.accumulated_N_absolute_radius += .001L;
    require(broad_clock.accumulated_N_absolute_radius < .01L,
            "near-boundary synthetic clock shifts retain their log support");
    detail::IdealAcousticTransportFrame fine, broad;
    require(detail::ideal_acoustic_transport_frame(
                epoch, loading, sound, actual, center, directions,
                nominal.vacuum, family, fine_clock, fine,
                counter.accounting()) == S::ok &&
                detail::ideal_acoustic_transport_frame(
                    epoch, loading, sound, actual, center, directions,
                    nominal.vacuum, family, broad_clock, broad,
                    counter.accounting()) == S::ok,
            "positive physical family admits conditional clock diagnostic");
    require(
        epoch.p == original_p && epoch.hcal == original_h &&
            epoch.raw_scaled_query->callbacks == callbacks && fine.x == 0 &&
            broad.x == 0,
        "bridge borrows unchanged H/query and leaves sqrt witness unavailable");
    require(fine.arithmetic.coefficient_at_actual_a_radius ==
                    broad.arithmetic.coefficient_at_actual_a_radius &&
                fine.arithmetic.gradient_radius ==
                    broad.arithmetic.gradient_radius &&
                fine.arithmetic.eta_gradient_radius ==
                    broad.arithmetic.eta_gradient_radius,
            "actual-a base and gradients do not repeat family clock support");
    for (unsigned i = 0; i < 6; ++i)
      require(broad.arithmetic.coefficient_radius[i] >=
                      fine.arithmetic.coefficient_radius[i] &&
                  broad.arithmetic.coefficient_radius[i] >=
                      broad.arithmetic.coefficient_at_actual_a_radius[i],
              "full clock band retains base coefficient radius");
    require(broad.arithmetic.hcal_radius > fine.arithmetic.hcal_radius &&
                broad.arithmetic.coefficient_radius[0] >
                    fine.arithmetic.coefficient_radius[0],
            "finite clock grows full H and wave bounds");
    detail::ideal_acoustic_transport_internal::Arithmetic arithmetic;
    detail::ideal_acoustic_transport_internal::Bands bands;
    require(detail::ideal_acoustic_transport_internal::bands(
                broad, family, bands, arithmetic) == S::ok,
            "full positive source quotient and inverse-H bands");
    detail::ideal_acoustic_transport_internal::Arithmetic fine_arithmetic;
    detail::ideal_acoustic_transport_internal::Bands fine_bands;
    require(
        detail::ideal_acoustic_transport_internal::bands(
            fine, family, fine_bands, fine_arithmetic) == S::ok &&
            fine_bands.change == bands.change &&
            fine_bands.remainder == bands.remainder,
        "fixed-actual-a source quotient bands do not acquire clock forcing");
    for (unsigned corner = 0; corner < 16; ++corner) {
      std::array<W, 4> delta;
      for (unsigned j = 0; j < 4; ++j)
        delta[j] =
            (corner & (1u << j) ? 1 : -1) *
            (std::abs(family.signed_shift[j]) + family.operation_radius[j]);
      const Primitive shifted{nominal.photon + delta[0],
                              nominal.baryon + delta[1], nominal.cdm + delta[2],
                              nominal.vacuum + delta[3] - delta[0] - delta[1] -
                                  delta[2]};
      const auto before = rational(shifted, a, k, alpha);
      for (W sign : {-1.L, 1.L}) {
        // |log(1+-rho)| <= rho/(1-rho) < clock radius for rho=.99*radius
        // and radius<.01. The five-percent primitive box and near-boundary
        // shift exercise family-clock cross terms rather than a loose margin.
        const W rho = .99L * broad_clock.accumulated_N_absolute_radius;
        const auto after = rational(shifted, a * (1 + sign * rho), k, alpha);
        for (unsigned i = 0; i < 6; ++i)
          enclosed(
              before.coefficient[i], after.coefficient[i],
              broad.arithmetic.coefficient_radius[i],
              "full-family rational clock jump is enclosed for both signs");
        enclosed(
            before.inverse_h, after.inverse_h, bands.inverse_h_arithmetic,
            "effective H radius encloses full-family inverse-H clock jump");
      }
    }
    auto refused_clock = fine_clock;
    refused_clock.conditional_profile = false;
    detail::IdealAcousticTransportFrame refused;
    require(detail::ideal_acoustic_transport_frame(
                epoch, loading, sound, actual, center, directions,
                nominal.vacuum, family, refused_clock, refused,
                counter.accounting()) != S::ok &&
                refused.status != S::ok,
            "missing conditional clock profile cannot become zero error");
    auto missing_other = family;
    missing_other.original_other_zero_witness = false;
    require(detail::ideal_acoustic_transport_frame(
                epoch, loading, sound, actual, center, directions,
                nominal.vacuum, missing_other, fine_clock, refused,
                counter.accounting()) != S::ok,
            "missing original zero-other witness refuses");
    missing_other = family;
    missing_other.original_other_shift = 1e-8L;
    require(detail::ideal_acoustic_transport_frame(
                epoch, loading, sound, actual, center, directions,
                nominal.vacuum, missing_other, fine_clock, refused,
                counter.accounting()) != S::ok,
            "nonzero original other direction refuses changed physics");
    auto singular = family;
    singular.operation_radius[0] = nominal.photon;
    require(detail::ideal_acoustic_transport_frame(
                epoch, loading, sound, actual, center, directions,
                nominal.vacuum, singular, fine_clock, refused,
                counter.accounting()) != S::ok,
            "singular positive physical family refuses");
    auto absent = epoch;
    absent.shadow.reset();
    require(detail::ideal_acoustic_transport_frame(
                absent, loading, sound, actual, center, directions,
                nominal.vacuum, family, fine_clock, refused,
                counter.accounting()) != S::ok,
            "missing same-source shadow refuses");
    if (stage_index++ == 0) {
      initial_controls(fine, family, counter);
      finite_initial_family_controls(nominal, k, alpha, fine, broad, family,
                                     fine_clock.accumulated_N_absolute_radius,
                                     broad_clock.accumulated_N_absolute_radius,
                                     counter);
      endpoint_controls(broad,family,counter);
      positive_transport_controls(broad,family,counter);
      signed_rhs_fusion_controls(broad,family);
    }
  }
  require(counter.diagnostics > 0 && counter.writes > 0,
          "test derivative/diagnostic owners were actually charged");
}
} // namespace ideal_acoustic_transport_test
