// Independent named controls: Hogg astro-ph/9905116 eq.23 and the photon
// energy law. Frequency-coordinate Simpson is independent of the production
// rectangular-band algebra. Constants are SI definitions, shared ancestry,
// not a separately measured physical oracle. No copied upstream code.
#include "irred/photometry.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
unsigned checks = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
constexpr long double c = 299792458.0L, h = 6.62607015e-34L;
photometry::Policy policy(unsigned mask = 7) { return {mask, 64, 1u << 20}; }
photometry::Input ordinary() {
  return {1e15, 4e-7, 8e-7, 3e-7, 1.5e-6, 1e10, 0.5, 2, 0.75, 120};
}
const photometry::Outcome &out(const photometry::Row &r, unsigned bit) {
  if (bit == 1)
    return r.flux_watt_per_square_metre;
  if (bit == 2)
    return r.energy_joule;
  return r.expected_photons;
}
void available(const photometry::Outcome &o) {
  check(o.availability == photometry::Availability::available, "available tag");
  check(o.numerical_status == numerics::Status::ok,
        "available numerical cause");
  check(o.value && std::isfinite(*o.value), "finite payload");
}
void near(double v, long double ref) {
  check(std::fabs(static_cast<long double>(v) - ref) <=
            2e-12L * std::fabs(ref) + 1e-300L,
        "frozen observable allocation");
}
photometry::Row eval(photometry::Input x, unsigned mask = 7) {
  auto b = photometry::evaluate(std::span(&x, 1), policy(mask));
  check(b.status == numerics::Status::ok && b.rows.size() == 1,
        "materialized batch");
  return b.rows[0];
}
std::array<long double, 3> frequency(const photometry::Input &x, unsigned n) {
  const long double u = 1 + static_cast<long double>(x.redshift),
                    d = x.luminosity_distance_metre;
  const long double a = std::max(
                        static_cast<long double>(x.observed_lower_metre),
                        u * x.rest_lower_metre),
                    b = std::min(
                        static_cast<long double>(x.observed_upper_metre),
                        u * x.rest_upper_metre);
  if (b <= a)
    return {};
  const long double lo = c / b, hi = c / a, step = (hi - lo) / n;
  long double energy = 0, photons = 0;
  const long double scale =
      x.luminosity_watt_per_metre * c / (4 * std::acos(-1.0L) * d * d * u);
  for (unsigned i = 0; i <= n; ++i) {
    const long double nu = lo + i * step,
                      weight = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
    const long double f = scale / (nu * nu);
    energy += weight * f;
    photons += weight * f / (h * nu);
  }
  const long double flux = energy * step / 3,
                    at = x.collecting_area_square_metre *
                         x.optical_transmission * x.observer_exposure_second;
  return {flux, at * flux, at * photons * step / 3};
}
void positive_control(photometry::Input x) {
  auto r = eval(x);
  check(r.admission_status == numerics::Status::ok, "valid admission");
  const auto coarse = frequency(x, 8192), fine = frequency(x, 16384);
  for (unsigned k = 0; k < 3; ++k) {
    available(out(r, 1u << k));
    check(fine[k] > 0, "positive independent prediction");
    check(std::fabs(fine[k] - coarse[k]) <= 2e-13L * std::fabs(fine[k]),
          "independent reference refinement");
    near(*out(r, 1u << k).value, fine[k]);
  }
}
void binary_edge(photometry::Input x, long double width) {
  auto r = eval(x);
  check(r.admission_status == numerics::Status::ok, "edge input admitted");
  const long double f = x.luminosity_watt_per_metre /
                        (4 * std::acos(-1.0L) *
                         static_cast<long double>(x.luminosity_distance_metre) *
                         x.luminosity_distance_metre *
                         (1 + static_cast<long double>(x.redshift)));
  const long double energy = x.collecting_area_square_metre *
                             x.optical_transmission *
                             x.observer_exposure_second * f * width;
  for (unsigned bit : {1u, 2u, 4u}) {
    const auto &o = out(r, bit);
    if (o.availability == photometry::Availability::failed) {
      check(o.numerical_status ==
                numerics::Status::conditioning_budget_exceeded,
            "unresolved boundary precise conditioning cause");
      check(!o.value, "conditioning no payload");
    } else {
      available(o);
      check(*o.value > 0, "exact positive overlap never zero");
      near(*o.value, bit == 1 ? f * width
                              : (bit == 2 ? energy
                                          : energy *
                                                (static_cast<long double>(
                                                     x.observed_lower_metre) +
                                                 width / 2) /
                                                (h * c)));
    }
  }
}
} // namespace
int main() {
  try {
    auto x = ordinary();
    positive_control(x);
    for (double z : {0.0, 1.0, 3.0, 10.0}) {
      x = ordinary();
      x.redshift = z;
      x.observed_lower_metre = (1 + z) * 5e-7;
      x.observed_upper_metre = (1 + z) * 7e-7;
      positive_control(x);
    }
    x = ordinary();
    x.observed_lower_metre = 7e-7;
    x.observed_upper_metre = 9e-7;
    positive_control(x);
    x = ordinary();
    auto base = eval(x);
    for (unsigned field = 0; field < 4; ++field) {
      auto y = x;
      if (field == 0)
        y.luminosity_watt_per_metre *= 2;
      if (field == 1)
        y.collecting_area_square_metre *= 2;
      if (field == 2)
        y.observer_exposure_second *= 2;
      if (field == 3)
        y.luminosity_distance_metre *= 2;
      auto r = eval(y);
      for (unsigned bit : {1u, 2u, 4u}) {
        available(out(r, bit));
        const long double factor =
            field == 3 ? 0.25L
                       : ((bit == 1 && (field == 1 || field == 2)) ? 1 : 2);
        near(*out(r, bit).value, factor * *out(base, bit).value);
      }
    }
    for (unsigned field = 0; field < 4; ++field) {
      auto y = ordinary();
      if (field == 0)
        y.luminosity_watt_per_metre = 0;
      if (field == 1)
        y.collecting_area_square_metre = 0;
      if (field == 2)
        y.optical_transmission = 0;
      if (field == 3)
        y.observer_exposure_second = 0;
      auto r = eval(y);
      for (unsigned bit : {1u, 2u, 4u}) {
        available(out(r, bit));
        if (field == 0 || bit != 1)
          check(*out(r, bit).value == 0, "mathematical zero");
        else
          check(*out(r, bit).value > 0,
                "detector zero preserves incident flux");
      }
    }
    x = ordinary();
    x.redshift = 1;
    x.observed_lower_metre = 2 * x.rest_upper_metre;
    x.observed_upper_metre = 3 * x.rest_upper_metre;
    for (unsigned bit : {1u, 2u, 4u}) {
      auto r = eval(x);
      available(out(r, bit));
      check(*out(r, bit).value == 0, "exact tangent zero");
    }
    const double s = std::ldexp(1.0, -20), e = std::ldexp(1.0, -52);
    x = ordinary();
    x.redshift = e;
    x.rest_lower_metre = s;
    x.rest_upper_metre = s * (1 + e);
    x.observed_lower_metre = s * (1 + 2 * e);
    x.observed_upper_metre = std::nextafter(x.observed_lower_metre, INFINITY);
    binary_edge(x, std::ldexp(1.0L, -124));
    x = ordinary();
    x.redshift = std::ldexp(1.0, -65);
    x.rest_lower_metre = std::nextafter(s, 0.0);
    x.rest_upper_metre = s;
    x.observed_lower_metre = s;
    x.observed_upper_metre = std::nextafter(s, INFINITY);
    binary_edge(x, std::ldexp(1.0L, -85));
    x = ordinary();
    x.redshift = 0;
    x.observed_lower_metre = 6e-7;
    x.observed_upper_metre = std::nextafter(x.observed_lower_metre, INFINITY);
    auto thin = eval(x);
    for (unsigned bit : {1u, 2u, 4u}) {
      available(out(thin, bit));
      check(*out(thin, bit).value > 0, "thin interior positive");
    }
    x = ordinary();
    x.optical_transmission = 0;
    x.luminosity_distance_metre = std::numeric_limits<double>::quiet_NaN();
    auto invalid = eval(x);
    check(invalid.admission_status != numerics::Status::ok,
          "full admission before zero");
    for (unsigned bit : {1u, 2u, 4u})
      check(!out(invalid, bit).value, "invalid no finite payload");
    x = ordinary();
    x.luminosity_watt_per_metre = 1e205;
    x.rest_lower_metre = 1e100;
    x.rest_upper_metre = 2e100;
    x.observed_lower_metre = 1.2e100;
    x.observed_upper_metre = std::nextafter(x.observed_lower_metre, INFINITY);
    x.luminosity_distance_metre = 1e3;
    x.redshift = 0;
    x.collecting_area_square_metre = 1;
    x.optical_transmission = 1;
    x.observer_exposure_second = 1;
    auto separated = eval(x);
    available(separated.energy_joule);
    check(separated.expected_photons.availability ==
                  photometry::Availability::failed &&
              !separated.expected_photons.value,
          "photon overflow does not erase energy");
    auto energy_only = eval(x, 2);
    available(energy_only.energy_joule);
    near(*energy_only.energy_joule.value, *separated.energy_joule.value);
    x = ordinary();
    x.luminosity_watt_per_metre = std::numeric_limits<double>::denorm_min();
    auto under = eval(x);
    check(under.flux_watt_per_square_metre.availability ==
                  photometry::Availability::failed &&
              !under.flux_watt_per_square_metre.value,
          "positive underflow is not zero");
    x = ordinary();
    auto one = eval(x, 2);
    check(one.flux_watt_per_square_metre.availability ==
                  photometry::Availability::omitted &&
              !one.flux_watt_per_square_metre.value,
          "unrequested flux absent");
    check(one.expected_photons.availability ==
                  photometry::Availability::omitted &&
              !one.expected_photons.value,
          "unrequested photons absent");
    available(one.energy_joule);
    auto bound = photometry::output_payload_bound(1);
    check(bound.has_value(), "payload bound");
    auto p = policy();
    p.maximum_output_bytes = *bound;
    auto admitted = photometry::evaluate(std::span(&x, 1), p);
    check(admitted.status == numerics::Status::ok, "exact cap admitted");
    p.maximum_output_bytes = *bound - 1;
    auto capped = photometry::evaluate(std::span(&x, 1), p);
    check(capped.status == numerics::Status::work_limit && capped.rows.empty(),
          "cap minus one owned failure");
    check(!photometry::output_payload_bound(
              std::numeric_limits<std::size_t>::max()),
          "count bound overflow");
    auto source = x;
    auto owned = photometry::evaluate(std::span(&x, 1), policy());
    x.redshift = 9;
    check(std::bit_cast<std::uint64_t>(owned.rows[0].source.redshift) ==
              std::bit_cast<std::uint64_t>(source.redshift),
          "result source owns bits");
    auto bad = source;
    bad.redshift = -1;
    std::array mixed{source, bad, source};
    auto batch = photometry::evaluate(mixed, policy());
    check(batch.rows.size() == 3, "mixed order retained");
    available(batch.rows[0].energy_joule);
    check(batch.rows[1].admission_status != numerics::Status::ok,
          "mixed failure retained");
    available(batch.rows[2].energy_joule);
    std::cout << "PASS " << checks << " independent photometry controls\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL check " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
