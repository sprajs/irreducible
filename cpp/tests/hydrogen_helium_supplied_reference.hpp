#pragma once
#include "hydrogen_helium_history_reference.hpp"
#include <cstdint>
#include <string>
// This reference owns no engine dependency. Source doubles are coefficients,
// never sampled native H, temperatures or histories. Rate/constant/libm ancestry
// is explicitly shared with the old independent reference, whose original
// initializer and driver remain available under their original identity.
namespace hydrogen_helium_supplied_reference {
namespace original = hydrogen_helium_history_reference;
using original::W;
using original::State;
using original::Model;
using original::Row;
using original::Stats;
inline constexpr const char *method =
    "direct-SI-emitted-source-log-expansion-RadauIIA3-FD6x6-bounded-halving/v1";
struct EmittedSource {
  double H0 = 0, photon = 0, other = 0, baryon = 0, cdm = 0;
  // Exact ordered emitted mass/eV-temperature/statistical-weight triples.
  std::vector<std::array<double, 3>> species;
};
struct Boundary { double hydrogen = 0, helium = 0, temperature = 0; };
struct Policy {
  unsigned refinement = 1;
  W root_tolerance = 2e-16L;
  std::size_t maximum_attempts = 1u << 20, maximum_rhs = 1u << 27,
              maximum_bytes = 32u * 1024 * 1024;
};
struct Physics {
  original::Physics original_physics;
  EmittedSource source;
  W H0, critical_energy, radiation, matter, lambda;
  bool original_input_lane = false;
  mutable std::optional<std::array<W, 2>> attempted_temperature;
  mutable std::optional<W> maximum_log_heiii;
  Physics(Model model, EmittedSource emitted)
      : original_physics(std::move(model)), source(std::move(emitted)) {
    for (const double value : {source.H0, source.photon, source.other, source.baryon, source.cdm})
      if (!std::isfinite(value)) throw std::runtime_error("reference nonfinite emitted source");
    constexpr W hbar = original::h / (2 * std::numbers::pi_v<W>);
    H0 = W(source.H0) * 1000 / original::mpc();
    critical_energy = 3 * H0 * H0 * original::c * original::c /
                      (8 * std::numbers::pi_v<W> * original::G);
    radiation = W(source.photon) + source.other;
    matter = W(source.baryon) + source.cdm;
    for (const auto &item : source.species) {
      if (item[0] != 0 || !(item[1] > 0) || !(item[2] > 0))
        throw std::runtime_error("reference massless source domain");
      const W density = 7 * W(item[2]) * std::numbers::pi_v<W> *
          std::numbers::pi_v<W> / 240 * std::pow(W(item[1]) * original::ev, 4) /
          std::pow(hbar * original::c, 3);
      if (!(density > 0) || !std::isfinite(density))
        throw std::runtime_error("reference species representation");
      radiation += density / critical_energy;
    }
    lambda = 1 - radiation - matter;
    if (!(H0 > 0 && critical_energy > 0 && radiation > 0 && matter > 0 && source.photon > 0 &&
          source.other >= 0 && source.baryon >= 0 && source.cdm >= 0 && lambda >= 0) ||
          !std::isfinite(lambda) || !std::isfinite(H0) || !std::isfinite(critical_energy))
      throw std::runtime_error("reference emitted flat closure");
  }
  W hubble(W u) const {
    if (original_input_lane) return original_physics.hubble(u);
    return H0 * std::sqrt(radiation * u * u * u * u + matter * u * u * u + lambda);
  }
  State operator()(W s, const State &state) const {
    const auto &m = original_physics.model;
    const W u = (1 + W(m.initial)) * std::exp(-s), T = state[2] * W(m.T0) * u;
    if (!(T > 0) || !std::isfinite(T))
      throw std::runtime_error("reference temperature representation");
    if (!attempted_temperature) attempted_temperature = std::array<W, 2>{T, T};
    else {
      (*attempted_temperature)[0] = std::min((*attempted_temperature)[0], T);
      (*attempted_temperature)[1] = std::max((*attempted_temperature)[1], T);
    }
    const W ne = u * u * u * (W(m.nH0) * state[0] + W(m.nHe0) * state[1]),
        quantum = std::pow(2 * std::numbers::pi_v<W> * original::me * original::kb * T /
                          (original::h * original::h), 1.5L),
        activity = std::log(quantum) - 54.4177655282L * original::ev /
                   (original::kb * T) - std::log(ne);
    if (!maximum_log_heiii || activity > *maximum_log_heiii) maximum_log_heiii = activity;
    const W photon_energy = W(source.photon) * critical_energy * u * u * u * u;
    return original_physics.rhs(s, state, hubble(u),
        original_input_lane ? std::optional<W>{} : std::optional<W>{photon_energy});
  }
  W opacity(W z, const State &state) const {
    const auto &m = original_physics.model;
    const W u = 1 + z;
    return original::c * original::sigma * u * u *
           (W(m.nH0) * state[0] + W(m.nHe0) * state[1]) / hubble(u);
  }
};
struct Result {
  std::vector<Row> rows;
  Stats stats;
  bool complete = false;
  std::string failure;
  std::optional<std::array<W, 2>> attempted_temperature;
  std::optional<W> maximum_log_heiii;
};
// A failing trial leaves the last accepted state and every work counter intact.
// Binary subdivision covers the exact nominal interval; it cannot skip a stop.
inline void advance(const Physics &physics, W begin, W end, State &state,
                    Stats &stats, const Policy &policy, unsigned depth = 0) {
  if (stats.attempts >= policy.maximum_attempts || stats.rhs >= policy.maximum_rhs)
    throw std::runtime_error("reference work cap");
  if (!(end > begin)) throw std::runtime_error("reference nonadvancing step");
  ++stats.attempts;
  try {
    const auto next = original::step(physics, begin, end - begin, state, stats,
        {policy.root_tolerance, policy.maximum_rhs, true});
    state = next;
    return;
  } catch (const std::runtime_error &) {
    ++stats.rejected;
    if (depth == 40 || stats.rhs >= policy.maximum_rhs || policy.maximum_rhs - stats.rhs < 2 ||
        stats.attempts >= policy.maximum_attempts) throw;
  }
  const W midpoint = begin + (end - begin) / 2;
  if (!(midpoint > begin && midpoint < end))
    throw std::runtime_error("reference subdivision representation");
  advance(physics, begin, midpoint, state, stats, policy, depth + 1);
  advance(physics, midpoint, end, state, stats, policy, depth + 1);
}
inline Result integrate(Model model, EmittedSource source, Boundary boundary,
                        const std::vector<double> &requested, Policy policy = {},
                        std::optional<State> wide_initial = {},
                        bool original_input_lane = false) {
  Result result;
  std::optional<Physics> physics;
  try {
    if ((policy.refinement != 1 && policy.refinement != 2 && policy.refinement != 4) ||
        !(policy.root_tolerance > 0) || !std::isfinite(policy.root_tolerance) ||
        requested.size() > 4096 || source.species.size() > 16)
      throw std::runtime_error("reference request domain");
    // Three owned species copies at peak, stops, scheduled rows, original-order
    // output and fixed-depth stack/strings. This is requested payload, not RSS.
    const std::size_t bound = 65536 + 3 * 16 * sizeof(std::array<double, 3>) +
        (requested.size() + 3) * (sizeof(W) + 2 * sizeof(Row));
    if (bound > policy.maximum_bytes) throw std::runtime_error("reference byte cap");
    auto payload = [&](std::size_t stop_capacity, std::size_t scheduled_capacity,
                       std::size_t output_capacity) {
      std::size_t total = 65536;
      auto add = [&](std::size_t count, std::size_t width) {
        if (count > (policy.maximum_bytes - std::min(total, policy.maximum_bytes)) / width)
          throw std::runtime_error("reference actual-capacity byte cap");
        total += count * width;
      };
      std::size_t species_capacity = std::max(model.massless_species.capacity(), source.species.capacity());
      if (physics) species_capacity = std::max({species_capacity, physics->source.species.capacity(),
          physics->original_physics.model.massless_species.capacity()});
      // Conservative simultaneous constructor source/model acquisition copies.
      add(species_capacity, 4 * sizeof(std::array<double, 3>));
      add(stop_capacity, sizeof(W)); add(scheduled_capacity, sizeof(Row)); add(output_capacity, sizeof(Row));
      if (total > policy.maximum_bytes) throw std::runtime_error("reference actual byte cap");
    };
    payload(0, 0, 0);
    if (!wide_initial && !(boundary.hydrogen > 0 && boundary.hydrogen < 1 &&
          boundary.helium > 0 && boundary.helium < 1 && boundary.temperature >= 4000 &&
          W(boundary.temperature) <= W(model.T0) * (1 + W(model.initial))))
      throw std::runtime_error("reference supplied boundary domain");
    for (const double query : requested)
      if (!std::isfinite(query) || query < model.late || query > model.initial)
        throw std::runtime_error("reference query domain");
    physics.emplace(model, std::move(source));
    payload(0, 0, 0);
    physics->original_input_lane = original_input_lane;
    // Optional wide p,q,Tm is a separately labelled reference input, never a
    // native expected value. Preserve its Tm rather than a theta round trip.
    const W initial_temperature = wide_initial ? (*wide_initial)[2] : W(boundary.temperature),
        initial_radiation_temperature = W(model.T0) * (1 + W(model.initial));
    State state{wide_initial ? (*wide_initial)[0] : W(boundary.hydrogen),
                wide_initial ? (*wide_initial)[1] : W(boundary.helium),
                initial_temperature / initial_radiation_temperature};
    if (!(state[0] > 0 && state[0] < 1 && state[1] > 0 && state[1] < 1 &&
          state[2] > 0 && state[2] <= 1 && initial_temperature >= 4000 &&
          initial_temperature <= initial_radiation_temperature))
      throw std::runtime_error("reference wide initial domain");
    std::vector<W> stops;
    stops.reserve(requested.size() + 3);
    payload(stops.capacity(), 0, 0);
    for (const double query : requested) stops.push_back(W(query));
    stops.push_back(W(model.initial)); stops.push_back(W(model.late));
    stops.push_back(W(model.initial) - .1L);
    std::sort(stops.begin(), stops.end(), std::greater<W>());
    stops.erase(std::unique(stops.begin(), stops.end()), stops.end());
    std::vector<Row> scheduled;
    scheduled.reserve(stops.size());
    payload(stops.capacity(), scheduled.capacity(), 0);
    auto append = [&](W z) {
      const W u = 1 + z;
      const Row row{z, state[0], state[1], z == W(model.initial)
          ? (wide_initial ? initial_temperature : W(boundary.temperature))
          : state[2] * W(model.T0) * u,
          u * u * u * (W(model.nH0) * state[0] + W(model.nHe0) * state[1]),
          physics->opacity(z, state)};
      for (const W x : {row.hydrogen, row.helium, row.temperature, row.electrons, row.opacity})
        if (!(x > 0) || !std::isfinite(x)) throw std::runtime_error("reference output representation");
      scheduled.push_back(row);
    };
    append(stops.front());
    W s = 0;
    for (std::size_t j = 1; j < stops.size(); ++j) {
      const W end = std::log((1 + W(model.initial)) / (1 + stops[j])), length = end - s;
      std::size_t count = std::max<std::size_t>(1, std::ceil(length / 1e-4L));
      if (stops[j] >= W(model.initial) - .1L) count = std::max<std::size_t>(128, count);
      if (count > policy.maximum_attempts / policy.refinement)
        throw std::runtime_error("reference nominal work cap");
      count *= policy.refinement;
      const W begin = s;
      for (std::size_t k = 0; k < count; ++k) {
        const W next = k + 1 == count ? end : begin + (k + 1) * length / count;
        advance(*physics, s, next, state, result.stats, policy);
        s = next;
      }
      append(stops[j]);
    }
    result.rows.reserve(requested.size());
    payload(stops.capacity(), scheduled.capacity(), result.rows.capacity());
    for (const double query : requested) {
      auto found = std::find_if(scheduled.begin(), scheduled.end(),
          [&](const Row &row) { return row.z == W(query); });
      if (found == scheduled.end()) throw std::runtime_error("reference source order");
      result.rows.push_back(*found);
    }
    result.complete = true;
  } catch (const std::exception &error) { result.failure = error.what(); }
  if (physics) {
    result.attempted_temperature = physics->attempted_temperature;
    result.maximum_log_heiii = physics->maximum_log_heiii;
  }
  return result;
}
} // namespace hydrogen_helium_supplied_reference
