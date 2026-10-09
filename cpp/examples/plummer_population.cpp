#include "irred/plummer_population.hpp"
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
int main() {
  using namespace irred::gravity;
  using S=irred::numerics::Status;
  PlummerSphere source{2e41,3e19,6.67430e-11};
  std::array<double,4> radii{0,3e19,6e19,3e20};
  auto population=prepare_plummer_population(source,radii);
  if(population.status()!=S::ok) return EXIT_FAILURE;
  std::vector<PlummerVelocityRequest> requests;
  for(std::size_t i=0;i<radii.size();++i)
    for(double fraction:{0.0,.25,.5,.75,.9,1.0,1.001})
      requests.push_back({i,population.radii()[i].escape_speed_m_s.value*fraction});
  auto batch=population.evaluate(requests);
  if(batch.rows.size()!=requests.size()) return EXIT_FAILURE;
  std::cout<<std::setprecision(17)
           <<"# synthetic self-gravitating isotropic collisionless mass-traces-tracer SI Plummer population\n"
           <<"# model_id="<<plummer_population_id<<'\n'
           <<"# M_kg="<<source.total_mass_kg<<" b_m="<<source.scale_radius_metres
           <<" G_m3_kg_s2="<<source.gravitational_coupling_m3_kg_s2<<'\n'
           <<"radius_m\tspeed_m_s\tsupport\trho_kg_m3\tsigma_axis2_m2_s2\tsigma_LOS2_m2_s2\tf_kg_s3_m6\tvector_PDF_s3_m3\tspeed_PDF_s_m\n";
  auto field=[](const irred::numerics::ScalarResult& value) {
    if(value.status==S::ok) std::cout<<value.value;else std::cout<<"NA";
  };
  for(const auto& row:batch.rows) {
    const auto& state=population.radii()[row.request.radius_index];
    const char* support=row.support==PlummerEnergySupport::bound?"bound":
                        row.support==PlummerEnergySupport::outside?"outside":
                        row.support==PlummerEnergySupport::ambiguous?"ambiguous":"invalid";
    if(row.status!=S::ok && !(row.support==PlummerEnergySupport::ambiguous && row.status==S::conditioning_budget_exceeded))
      return EXIT_FAILURE;
    std::cout<<row.radius_metres<<'\t'<<row.request.speed_m_s<<'\t'<<support<<'\t';
    field(state.density_kg_m3);std::cout<<'\t';field(state.one_axis_variance_m2_s2);std::cout<<'\t';
    field(state.projected_los_variance_m2_s2);std::cout<<'\t';field(row.distribution_function_kg_s3_m6);std::cout<<'\t';
    field(row.vector_velocity_density_s3_m3);std::cout<<'\t';field(row.speed_density_s_m);std::cout<<'\n';
  }
  return EXIT_SUCCESS;
}
