// Independent synthetic grey-magnitude contract. Three-dimensional exact
// adjugate covariance, analytic de Sitter geometry; shared libm disclosed.
// Frozen assembled absolute budget 1e-10, magnitude budget 2e-14. No inference.
#include "irred/supernova.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
int checks = 0;
long double worst = 0;
void require(bool v, const char *why) {
  ++checks;
  if (!v)
    throw std::runtime_error(why);
}
void close(long double a, long double b, long double budget = 1e-10L) {
  worst = std::max(worst, std::abs(a - b));
  require(std::abs(a - b) <= budget, "independent analytic budget");
}
auto bits(double x) { return std::bit_cast<std::uint64_t>(x); }
supernova::Policy policy() {
  supernova::Policy p;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.maximum_forward_sensitivity = 1e-10;
  p.maximum_models = 8;
  p.maximum_native_bytes = 1048576;
  p.requested = 63;
  p.background = {{numerics::IntegrationPolicy{1e-12, 1e-12, 100000, 30}},
                  16,
                  10000,
                  1000,
                  1048576};
  return p;
}
constexpr std::array<double, 3> z{.125, .5, 1}, observer{.3, .75, 1.25};
constexpr std::array<long double, 3> noise{.02L, -.03L, .01L};
long double shape(std::size_t i, bool constant) {
  return 5 *
         std::log10((1 + (long double)observer[i]) * (constant ? .5L : z[i]));
}
long double basis(std::size_t i, bool constant) {
  return std::log1p(constant ? .5L : (long double)z[i]) / std::log(2.L);
}
supernova::Consumer
make(bool constant = false, double intercept = 17, bool fifth_bin = false,
     observations::Role role = observations::Role::synthetic_control) {
  observations::Input in{};
  in.profile = observations::Profile::typed_magnitude_covariance;
  in.role = role;
  in.unit = observations::Unit::magnitude;
  in.calibration = observations::Calibration::not_applicable;
  in.uncertainty = observations::Uncertainty::covariance;
  in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  in.component = observations::Component::total;
  in.ordering_provenance = "synthetic independent adjugate order";
  in.table_sha256 = std::string(64, 'c');
  in.uncertainty_sha256 = std::string(64, 'd');
  in.measurement_ids = {"short", "peer-long-identity-aaaaaaaaaaaaaaaa",
                        "third"};
  in.event_ids = in.measurement_ids;
  in.uncertainty_axis_ids = in.measurement_ids;
  in.uncertainty_matrix = {4, 1, 0, 1, 3, 1, 0, 1, 2};
  in.missing = {0, 0, 0};
  in.source_selection = {1, 1, 1};
  in.quality = {0, 0, 0};
  std::vector<supernova::MagnitudeCoordinate> qs;
  const std::array<double, 3> high_z{1.25, 2., 2.5},
      high_observer{1.5, 1.7, 2.1};
  for (std::size_t i = 0; i < 3; ++i) {
    const auto zi = fifth_bin ? high_z[i] : constant ? .5 : z[i];
    const auto oi = fifth_bin ? high_observer[i] : observer[i];
    qs.push_back({zi, {oi, cosmology::Convention::released_zhd_zhel}});
    const auto geometric = fifth_bin
                               ? 5 * std::log10((1 + (long double)oi) * zi)
                               : shape(i, constant);
    const auto effect = fifth_bin ? std::log1p((long double)zi) / std::log(2.L)
                                  : basis(i, constant);
    in.values.push_back(
        static_cast<double>(geometric + intercept + .13L * effect + noise[i]));
  }
  auto ids = in.measurement_ids;
  auto source = observations::prepare(std::move(in), {3, 9, 4096});
  require(source.status() == observations::Status::ok, "typed synthetic");
  auto shared =
      std::make_shared<const observations::Prepared>(std::move(source));
  supernova::SelectedMagnitudeSource selected{shared, {0, 1, 2}, ids, qs};
  auto c = supernova::prepare(
      std::move(selected),
      {numerics::Arithmetic::longdouble_cpu_v1, 3, 9, 1e-10, 1000000});
  require(c.selected_source().source.get() == shared.get(),
          "same immutable source retained without full source copy");
  require(c.status() == supernova::Status::ok, "retained peer consumer");
  return c;
}
std::pair<long double, long double>
exact_profile(std::array<long double, 3> r) {
  // C^-1 = [[5,-2,1],[-2,8,-4],[1,-4,11]]/18;
  // W1=[4,2,8]/18 and 1'W1=14/18. No production solves.
  auto h = 4 * r[0] + 2 * r[1] + 8 * r[2];
  auto q = (5 * r[0] * r[0] + 8 * r[1] * r[1] + 11 * r[2] * r[2] -
            4 * r[0] * r[1] + 2 * r[0] * r[2] - 8 * r[1] * r[2] - h * h / 14) /
           18;
  return {q, h / 14};
}
// Independent split-z Simpson; production uses segment antiderivatives.
long double piecewise_integral(double z, const std::array<double, 5> &q,
                               unsigned panels) {
  auto E = [&](long double x) {
    long double value = 1;
    for (size_t k = 0; k < 5; ++k) {
      auto lo = (long double)cosmology::piecewise_q_edges[k],
           hi = std::min(x, (long double)cosmology::piecewise_q_edges[k + 1]);
      if (hi <= lo)
        break;
      value *= std::pow((1 + hi) / (1 + lo), 1 + (long double)q[k]);
    }
    return value;
  };
  long double total = 0;
  for (size_t k = 0; k < 5; ++k) {
    auto lo = (long double)cosmology::piecewise_q_edges[k],
         hi = std::min((long double)z,
                       (long double)cosmology::piecewise_q_edges[k + 1]);
    if (hi <= lo)
      break;
    auto h = (hi - lo) / panels;
    long double sum = 0;
    for (unsigned j = 0; j <= panels; ++j)
      sum += (j == 0 || j == panels ? 1 : j % 2 ? 4 : 2) / E(lo + j * h);
    total += h * sum / 3;
  }
  return total;
}
void absent(const supernova::Slot &s) {
  require(s.status != supernova::Status::ok, "failed slot");
  require(!s.score && s.magnitude_shifts.empty() && s.geometric_shape.empty() &&
              s.corrected_residuals.empty() && s.profiled_residuals.empty(),
          "failed prediction carries no finite arrays");
}
} // namespace
int main() {
  try {
    auto measured =
        make(false, 17, false, observations::Role::observed_measurement);
    auto c = make();

    auto shared_source = c.selected_source().source;
    auto prep_bound = supernova::preparation_payload_bound(
        *shared_source, 3, numerics::Arithmetic::longdouble_cpu_v1);
    require(prep_bound.has_value(), "authoritative SN preparation bound");
    require(
        !supernova::preparation_payload_bound(
            *shared_source, SIZE_MAX, numerics::Arithmetic::longdouble_cpu_v1),
        "SN preparation square overflow rejected");
    supernova::PreparationPolicy prep{numerics::Arithmetic::longdouble_cpu_v1,
                                      3, 9, 1e-10, *prep_bound};
    auto descriptor = c.selected_source();
    require(descriptor.source_indices.capacity() == 3 &&
                descriptor.coordinates.capacity() == 3 &&
                descriptor.ordered_ids.capacity() == 3,
            "exact-bound fixture has compact descriptor storage");
    auto bounded = supernova::prepare(descriptor, prep);
    require(bounded.status() == supernova::Status::ok,
            "exact conservative preparation bound admitted");
    --prep.maximum_native_bytes;
    auto capped = supernova::prepare(descriptor, prep);
    require(capped.status() == supernova::Status::work_limit,
            "preparation bound minus one rejects before factor");
    prep.maximum_native_bytes = 1000000;
    descriptor.source_indices.reserve(200000);
    auto excess = supernova::prepare(std::move(descriptor), prep);
    require(excess.status() == supernova::Status::work_limit,
            "transferred selected-index reserve charged");
    auto p = policy();
    const auto expansion = cosmology::ExpansionSpec(cosmology::ConstantQ(-1));
    std::array<supernova::ModelPoint, 1> none{
        {{expansion, cosmology::FlatFLRW{}, supernova::NoMagnitudeEffect{}}}};
    auto old = c.evaluate_batch(none, p);
    auto fitted =
        make(false, 17, false, observations::Role::released_fitted_summary);
    require(measured.probability_metadata().source_semantics !=
                    fitted.probability_metadata().source_semantics &&
                measured.probability_metadata().source_semantics !=
                    c.probability_metadata().source_semantics,
            "probability metadata preserves distinct actual source roles");
    auto posterior_input = measured.selected_source().source->source();
    posterior_input.role = observations::Role::posterior_summary;
    auto posterior =
        observations::prepare(std::move(posterior_input), {3, 9, 4096});
    require(posterior.status() == observations::Status::incompatible_semantics,
            "posterior is not a generic independent magnitude measurement");
    require(measured.selected_source().source_indices ==
                std::vector<std::size_t>({0, 1, 2}),
            "generic measured admission applies no Pantheon cut");
    auto measured_result = measured.evaluate_batch(none, p);
    auto fitted_result = fitted.evaluate_batch(none, p);
    require(measured.selected_source().source->source().role ==
                observations::Role::observed_measurement,
            "measured role preserved without release upgrade");
    require(fitted.selected_source().source->source().role ==
                observations::Role::released_fitted_summary,
            "fitted role preserved");
    require(bits(measured_result.slots[0].score->relative_profile_score) ==
                    bits(old.slots[0].score->relative_profile_score) &&
                bits(fitted_result.slots[0].score->relative_profile_score) ==
                    bits(old.slots[0].score->relative_profile_score),
            "same numeric values role parity");

    for (double eps : {-.5, -.2, -0., 0., .13, .2, .5}) {
      std::array<supernova::ModelPoint, 1> points{
          {{expansion, cosmology::FlatFLRW{},
            supernova::GreyLog1pMagnitude(eps)}}};
      auto out = c.evaluate_batch(points, p);
      require(out.status == supernova::Status::ok && out.slots.size() == 1,
              "owned batch");
      auto &s = out.slots[0];
      require(s.status == supernova::Status::ok, "finite grey");
      require(
          bits(std::get<supernova::GreyLog1pMagnitude>(s.source.source_effect)
                   .epsilon_mag) == bits(eps),
          "raw epsilon bits");
      std::array<long double, 3> r;
      for (std::size_t i = 0; i < 3; ++i) {
        r[i] = 17 + (.13L - eps) * basis(i, false) + noise[i];
        close(s.magnitude_shifts[i], eps * basis(i, false), 2e-14L);
        close(s.geometric_shape[i], shape(i, false), 2e-14L);
        close(s.corrected_residuals[i], r[i], 2e-14L);
        if (eps != 0)
          require(std::abs(s.magnitude_shifts[i] -
                           eps * std::log1p(observer[i]) / std::log(2.)) > 1e-3,
                  "zhd distinct from observer");
      }
      auto [q, a] = exact_profile(r);
      close(s.score->quadratic, q);
      close(s.score->relative_profile_score, -q / 2);
      close(s.score->offset_coefficient, a);
      if (eps == 0) {
        require(bits(s.score->quadratic) == bits(old.slots[0].score->quadratic),
                "zero q exact");
        require(s.corrected_residuals == old.slots[0].corrected_residuals,
                "zero residual exact");
        require(s.profiled_residuals == old.slots[0].profiled_residuals,
                "zero canonical exact");
        require(s.work.callbacks == old.slots[0].work.callbacks,
                "zero work exact");
      }
    }
    auto same = make(true), shifted = make(true, 117);
    for (double eps : {-.5, 0., .5}) {
      std::array<supernova::ModelPoint, 1> points{
          {{expansion, cosmology::FlatFLRW{},
            supernova::GreyLog1pMagnitude(eps)}}};
      auto a = same.evaluate_batch(points, p),
           b = shifted.evaluate_batch(points, p);
      auto expected = exact_profile(noise).first;
      close(a.slots[0].score->quadratic, expected);
      close(b.slots[0].score->quadratic, expected);
      close(b.slots[0].score->offset_coefficient -
                a.slots[0].score->offset_coefficient,
            100);
    }
    for (double eps : {std::nextafter(.5, 1.), std::nextafter(-.5, -1.),
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::denorm_min()}) {
      std::array<supernova::ModelPoint, 1> points{
          {{expansion, cosmology::FlatFLRW{},
            supernova::GreyLog1pMagnitude(eps)}}};
      auto out = c.evaluate_batch(points, p);
      require(out.slots.size() == 1, "attempt materialized");
      if (eps == std::numeric_limits<double>::denorm_min()) {
        const auto &s = out.slots[0];
        require(s.geometry.availability == cosmology::Availability::available &&
                    s.geometric_shape.size() == 3,
                "independent geometry survives failed effect");
        require(!s.score && s.magnitude_shifts.empty() &&
                    s.corrected_residuals.empty() &&
                    s.profiled_residuals.empty(),
                "failed dependent outputs carry no payload");
        require(s.effect.numerical_status == numerics::Status::outside_domain &&
                    s.corrected.numerical_status ==
                        numerics::Status::outside_domain &&
                    s.profile.numerical_status ==
                        numerics::Status::outside_domain,
                "dependent requested states retain effect failure cause");
      } else {
        absent(out.slots[0]);
      }
      require(bits(std::get<supernova::GreyLog1pMagnitude>(
                       out.slots[0].source.source_effect)
                       .epsilon_mag) == bits(eps),
              "failed source bits");
      if (eps == std::numeric_limits<double>::denorm_min())
        require(out.slots[0].numerical_status ==
                    numerics::Status::outside_domain,
                "tiny shift normal-output policy");
    }
    std::array<supernova::ModelPoint, 2> points{
        {{expansion, cosmology::FlatFLRW{}, supernova::GreyLog1pMagnitude(.2)},
         {expansion, cosmology::FlatFLRW{},
          supernova::GreyLog1pMagnitude(-.2)}}};
    auto no = p;
    no.maximum_native_bytes = 0;
    auto out = c.evaluate_batch(points, no);
    require(out.status == supernova::Status::work_limit && out.slots.empty(),
            "preallocation byte cap");
    // Current requested storage differs from obsolete five-array wrapper. Test
    // the declared host-layout boundary without preserving an old layout
    // golden.
    size_t lo = 0, hi = p.maximum_native_bytes;
    while (lo + 1 < hi) {
      auto mid = (lo + hi) / 2;
      no = p;
      no.maximum_native_bytes = mid;
      auto v = c.evaluate_batch(points, no);
      if (v.status == supernova::Status::ok)
        hi = mid;
      else
        lo = mid;
    }
    no = p;
    no.maximum_native_bytes = lo;
    require(c.evaluate_batch(points, no).status ==
                supernova::Status::work_limit,
            "storage cap minus one");
    no.maximum_native_bytes = hi;
    out = c.evaluate_batch(points, no);
    require(out.slots.size() == 2 &&
                out.slots[0].status == supernova::Status::ok,
            "storage exactly admitted");
    no = p;
    no.maximum_models = 0;
    out = c.evaluate_batch(points, no);
    require(out.status == supernova::Status::work_limit && out.slots.empty(),
            "model cap");
    no = p;
    no.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    out = c.evaluate_batch(points, no);
    require(out.status == supernova::Status::incompatible_metadata &&
                out.slots.empty(),
            "precision mismatch");
    auto first = c.evaluate_batch(std::span(points).first(1), p);
    auto work = first.slots[0].work.callbacks;
    require(work > 0, "actual callbacks");
    no = p;
    no.background.maximum_callbacks = work;
    out = c.evaluate_batch(points, no);
    require(out.slots.size() == 2 &&
                out.slots[0].status == supernova::Status::ok,
            "first consumes global callback budget");
    require(!out.slots[1].score && out.slots[1].geometric_shape.empty() &&
                out.slots[1].corrected_residuals.empty() &&
                out.slots[1].profiled_residuals.empty(),
            "failed geometry dependencies carry no payload");
    require(out.slots[1].status == supernova::Status::work_limit &&
                out.slots[1].numerical_status == numerics::Status::work_limit,
            "truthful global work failure");
    // New composition: fixed-five-bin expansion and the same magnitude effect,
    // with independent radial integration and the unchanged exact adjugate.
    const std::array<double, 5> qbins{-1, 0, -.5, .2, 1};
    for (double eps : {-.2, 0., .2}) {
      std::array<supernova::ModelPoint, 1> composed{
          {{cosmology::FixedFiveBinQ(qbins), cosmology::FlatFLRW{},
            supernova::GreyLog1pMagnitude(eps)}}};
      auto result = c.evaluate_batch(composed, p);
      require(result.slots.size() == 1 &&
                  result.slots[0].status == supernova::Status::ok,
              "piecewise plus grey composition");
      const auto &slot = result.slots[0];
      require(result.work.callbacks == 0 && result.work.segment_visits > 0,
              "analytic composition accounts segments not callbacks");
      std::array<long double, 3> residual;
      for (size_t i = 0; i < 3; ++i) {
        auto a = piecewise_integral(z[i], qbins, 4096),
             b = piecewise_integral(z[i], qbins, 8192);
        require(std::abs(a - b) <= 2e-15L,
                "composition independent reference refinement");
        auto geometric = 5 * std::log10((1 + (long double)observer[i]) * b);
        close(slot.geometric_shape[i], geometric, 2e-14L);
        close(slot.magnitude_shifts[i], eps * basis(i, false), 2e-14L);
        residual[i] =
            (long double)c.selected_source().source->source().values[i] -
            geometric - eps * basis(i, false);
        close(slot.corrected_residuals[i], residual[i], 2e-14L);
      }
      auto [q, a] = exact_profile(residual);
      close(slot.score->quadratic, q);
      close(slot.score->relative_profile_score, -q / 2);
      close(slot.score->offset_coefficient, a);
    }
    // Companion composition reaches the entire fifth bin and its upper
    // endpoint.
    auto high = make(false, 17, true);
    const std::array<double, 3> high_z{1.25, 2., 2.5},
        high_observer{1.5, 1.7, 2.1};
    std::array<supernova::ModelPoint, 1> high_none{
        {{cosmology::FixedFiveBinQ(qbins), cosmology::FlatFLRW{},
          supernova::NoMagnitudeEffect{}}}};
    auto high_baseline = high.evaluate_batch(high_none, p);
    for (double eps : {-.2, -0., 0., .2}) {
      std::array<supernova::ModelPoint, 1> model{
          {{cosmology::FixedFiveBinQ(qbins), cosmology::FlatFLRW{},
            supernova::GreyLog1pMagnitude(eps)}}};
      auto prediction = high.evaluate_batch(model, p);
      const auto &slot = prediction.slots[0];
      require(slot.status == supernova::Status::ok,
              "fifth-bin composition finite");
      std::array<long double, 3> residual;
      for (size_t i = 0; i < 3; ++i) {
        auto coarse = piecewise_integral(high_z[i], qbins, 4096);
        auto fine = piecewise_integral(high_z[i], qbins, 8192);
        close(coarse, fine, 2e-15L);
        const auto geometric =
            5 * std::log10((1 + (long double)high_observer[i]) * fine);
        const auto B = eps * std::log1p((long double)high_z[i]) / std::log(2.L);
        close(slot.geometric_shape[i], geometric, 2e-14L);
        close(slot.magnitude_shifts[i], B, 2e-14L);
        residual[i] =
            high.selected_source().source->source().values[i] - geometric - B;
        close(slot.corrected_residuals[i], residual[i], 2e-14L);
      }
      auto [quadratic, coefficient] = exact_profile(residual);
      close(slot.score->quadratic, quadratic);
      close(slot.score->offset_coefficient, coefficient);
      close(slot.score->relative_profile_score, -quadratic / 2);
      if (eps == 0) {
        require(bits(slot.score->relative_profile_score) ==
                    bits(high_baseline.slots[0].score->relative_profile_score),
                "fifth-bin none/zero score bits");
        require(slot.corrected_residuals ==
                        high_baseline.slots[0].corrected_residuals &&
                    slot.profiled_residuals ==
                        high_baseline.slots[0].profiled_residuals,
                "fifth-bin none/zero residual parity");
      }
    }
    // Requested-only dependencies: geometric output cannot be erased by an
    // unrequested tiny effect; effect-only does no expansion integration.
    std::array<supernova::ModelPoint, 1> tiny{
        {{expansion, cosmology::FlatFLRW{},
          supernova::GreyLog1pMagnitude(
              std::numeric_limits<double>::denorm_min())}}};
    no = p;
    no.requested = (uint32_t)supernova::Output::geometric_shape;
    out = c.evaluate_batch(tiny, no);
    require(out.slots[0].status == supernova::Status::ok &&
                out.slots[0].geometric_shape.size() == 3 &&
                out.slots[0].magnitude_shifts.empty(),
            "unrequested effect underflow isolated");
    no = p;
    no.requested = (uint32_t)supernova::Output::magnitude_effect;
    out = c.evaluate_batch(std::span(points).first(1), no);
    require(out.slots[0].status == supernova::Status::ok &&
                out.slots[0].magnitude_shifts.size() == 3 &&
                out.work.callbacks == 0 && out.work.segment_visits == 0 &&
                out.slots[0].geometric_shape.empty(),
            "effect-only zero expansion work");
    no = p;
    no.background.maximum_callbacks = 0;
    out = c.evaluate_batch(tiny, no);
    require(out.slots[0].geometry.numerical_status ==
                numerics::Status::work_limit,
            "geometry retains its own work-limit cause");
    require(out.slots[0].effect.numerical_status ==
                numerics::Status::outside_domain,
            "effect retains its own underflow cause");
    std::printf("grey magnitude independent PASS %d maxerror %.18Lg\n", checks,
                worst);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "grey independent FAIL %s after%d\n", e.what(),
                 checks);
    return 1;
  }
}
