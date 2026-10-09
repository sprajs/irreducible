#include <irred/plummer_population.hpp>
#include <array>
#include <cmath>
#include <cstdlib>
#include <utility>
int main() {
  using namespace irred::gravity;
  using S=irred::numerics::Status;
  PlummerSphere source{2e41,3e19,6.67430e-11};
  std::array<double,3> radii{0,3e19,3e20};
  auto population=prepare_plummer_population(source,radii);
  if(population.status()!=S::ok || population.native_evaluations()!=6) return EXIT_FAILURE;
  std::array<PlummerVelocityRequest,3> velocities{{{0,0},{1,population.radii()[1].escape_speed_m_s.value*.5},
      {2,population.radii()[2].escape_speed_m_s.value*1.001}}};
  auto batch=population.evaluate(velocities);
  if(batch.status!=S::ok || batch.rows.size()!=3 || batch.velocity_evaluations!=3) return EXIT_FAILURE;
  if(batch.rows[0].speed_density_s_m.value!=0 || batch.rows[0].distribution_function_kg_s3_m6.value<=0) return EXIT_FAILURE;
  if(batch.rows[2].support!=PlummerEnergySupport::outside || batch.rows[2].speed_density_s_m.value!=0) return EXIT_FAILURE;
  if(std::abs(population.radii()[1].one_axis_variance_m2_s2.value/(population.radii()[1].relative_potential_m2_s2.value/6)-1)>1e-14) return EXIT_FAILURE;
  auto owned=std::move(population);
  if(population.status()!=S::invalid_input || owned.evaluate(velocities).status!=S::ok) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
