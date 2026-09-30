// Closed rectangular-band controls derived in docs/photometry.md.
// Shared SI constants/pi are declared ancestry, not independent measurements.
#include "irred/photometry.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
namespace {
unsigned checks = 0;
void check(bool good) { ++checks; if (!good) throw std::runtime_error("photometry owner control"); }
void near(const irred::photometry::Outcome &o, long double expected) {
  check(o.availability == irred::photometry::Availability::available);
  check(o.numerical_status == irred::numerics::Status::ok && o.value.has_value());
  check(std::abs(static_cast<long double>(*o.value)-expected) <= 2e-12L*std::abs(expected)+1e-300L);
}
}
int main() {
  try {
    using namespace irred;
    photometry::Input x{static_cast<double>(4*std::numbers::pi_v<long double>),
                       1e-7, 3e-6, 1e-6, 2e-6, 1, 0, 2, 0.5, 3};
    // Lλ/(4πDL²)=1 W/m³; collect through area*transmission*time=3 m²s.
    auto batch = photometry::evaluate(std::span(&x,1),{7,16,1u<<20});
    check(batch.status == numerics::Status::ok && batch.rows.size()==1);
    near(batch.rows[0].flux_watt_per_square_metre,1e-6L);
    near(batch.rows[0].energy_joule,3e-6L);
    near(batch.rows[0].expected_photons,4.5e-12L/(6.62607015e-34L*299792458.0L));
    // Observer exposure is not divided by (1+z) a second time.
    x.redshift=1;
    batch=photometry::evaluate(std::span(&x,1),{7,16,1u<<20});
    near(batch.rows[0].energy_joule,1.5e-6L);
    near(batch.rows[0].expected_photons,2.25e-12L/(6.62607015e-34L*299792458.0L));
    // Scaling of supplied luminosity distance is inverse square.
    x.luminosity_distance_metre=2;
    batch=photometry::evaluate(std::span(&x,1),{7,16,1u<<20});
    near(batch.rows[0].energy_joule,0.375e-6L);
    x.collecting_area_square_metre=0;
    batch=photometry::evaluate(std::span(&x,1),{7,16,1u<<20});
    near(batch.rows[0].energy_joule,0);
    check(*batch.rows[0].flux_watt_per_square_metre.value>0);
    std::cout<<"PASS "<<checks<<" photometry owner controls\n";
  } catch(const std::exception &e) { std::cerr<<e.what()<<'\n';return 1; }
}
