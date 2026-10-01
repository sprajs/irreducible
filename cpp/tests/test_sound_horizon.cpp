#include "irred/sound_horizon.hpp"
#include "irred/quantities.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred;
int main() {
  unsigned checks = 0;
  auto require = [&](bool yes, const char *what) {
    ++checks;
    if (!yes) throw std::runtime_error(what);
  };
  cosmology::SoundHorizonPolicy policy{1e-9, 2e-11, 1000000, 40, 64,
                                      2000000, 1024 * 1024};
  // Radiation-only and Lambda=0 no-baryon-loading analytic controls. These
  // identities are derived directly by integrating 1/sqrt(r+m*a).
  for (double redshift : {0., 999., 1e6}) {
    for (double matter : {0., .75}) {
      cosmology::SoundHorizonRequest request{{70., matter, 1.-matter, 0., .1},
                                            redshift, "analytic control"};
      const auto batch = cosmology::evaluate_sound_horizon({&request, 1}, policy);
      require(batch.status == numerics::Status::ok && batch.rows.size()==1,
              "owned analytic row");
      const auto &row = batch.rows.front();
      require(row.status == numerics::Status::ok && row.sound_horizon_mpc,
              "analytic payload");
      const long double a = 1.L/(1.L+redshift), radiation=1.L-matter;
      const long double integral = 2*a/(std::sqrt(radiation+matter*a)+std::sqrt(radiation));
      const long double expected = (long double)speed_of_light_m_per_s/
          1000/70/std::sqrt(3.L)*integral;
      require(std::isfinite(*row.sound_horizon_mpc) &&
                  std::abs((long double)*row.sound_horizon_mpc-expected) <=
                      1e-9L+2e-11L*std::abs(expected), "analytic fixed budget");
      require(row.callbacks == batch.callbacks && row.callbacks > 0,
              "actual callback accounting");
    }
  }
  for (bool reverse : {false, true}) {
    const double tiny=std::numeric_limits<double>::denorm_min();
    cosmology::SoundHorizonRequest request{{70., reverse?tiny:1., reverse?1.:tiny,
                                          0., reverse?1.:tiny}, 1000., "closure adversary"};
    const auto batch=cosmology::evaluate_sound_horizon({&request,1},policy);
    require(batch.rows.front().status==numerics::Status::outside_domain &&
            !batch.rows.front().sound_horizon_mpc && batch.callbacks==0,
            "exact positive excess closure rejected");
  }
  std::cout << "sound_horizon owner checks=" << checks << '\n';
}
