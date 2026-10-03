// Fresh installed SDK consumer: retain one physical scene and one Gaussian
// factor, then use the full ordered vector. No CLI, image assets or RNG.
#include "irred/two_deflector_forward.hpp"
#include "irred/statistics.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
using namespace irred::lensing;
using S = irred::numerics::Status;
using W = long double;
namespace {
void need(bool ok,const char *why) {
  if(!ok) { std::cerr<<"FAIL installed "<<why<<'\n'; std::exit(1); }
}
TwoDeflectorScene scene() {
  TwoDeflectorScene s;
  s.deflectors={SoftenedPotentialComponent{{-.3,.1},.6,.3,.8,.125},
                SoftenedPotentialComponent{{.4,-.15},.4,.2,.9,-.25}};
  s.source={{.0625,.125},1.5,1,.25,1e13}; s.shear_1=.04; s.shear_2=.03;
  s.theta_scale_radians=1e-5; s.exposure_seconds=600;
  s.uniform_sky_electrons_per_second_per_radian_squared=1e9;
  s.origin="installed-synthetic-two-potential-source/affine-electron-response/no-measured-ACS";
  return s;
}
std::array<ForwardPixelRectangle,16> pixels() {
  std::array<ForwardPixelRectangle,16> out{};
  constexpr double centers[]{-48,-16,16,48};
  for(unsigned iy=0;iy<4;++iy) for(unsigned ix=0;ix<4;++ix) {
    const unsigned i=4*iy+ix; const auto x=centers[ix],y=centers[iy];
    out[i]={x-.05,x+.05,y-.05,y+.05,i};
  }
  return out;
}
bool complete(const ForwardPixelMeans &means,std::size_t n) {
  if(means.status!=S::ok||means.rows().size()!=n) return false;
  for(const auto &row:means.rows()) if(row.status!=S::ok||!row.electrons||!row.errors) return false;
  return true;
}
void print_work(const ForwardWork &w) {
  std::cout<<'{';
#define IRRED_PRINT_FORWARD_WORK(field) <<"\"" #field "\":"<<w.field
  std::cout IRRED_PRINT_FORWARD_WORK(preparation_started)
      <<',' IRRED_PRINT_FORWARD_WORK(pixels_started)
      <<',' IRRED_PRINT_FORWARD_WORK(index_comparisons_started)
      <<',' IRRED_PRINT_FORWARD_WORK(psf_preflight_started)
      <<',' IRRED_PRINT_FORWARD_WORK(field_samples_started)
      <<',' IRRED_PRINT_FORWARD_WORK(deflector_evaluations_started)
      <<',' IRRED_PRINT_FORWARD_WORK(outer_integrations_started)
      <<',' IRRED_PRINT_FORWARD_WORK(inner_integrations_started)
      <<',' IRRED_PRINT_FORWARD_WORK(outer_callbacks_started)
      <<',' IRRED_PRINT_FORWARD_WORK(inner_callbacks_started)
      <<',' IRRED_PRINT_FORWARD_WORK(failed_starts);
#undef IRRED_PRINT_FORWARD_WORK
  std::cout<<'}';
}
} // namespace
int main() {
  static_assert(!std::is_copy_constructible_v<PreparedTwoDeflectorForward>);
  static_assert(!std::is_copy_constructible_v<ForwardPixelMeans>);
  const AffineDetectorCutout cutout{{0,0},{.02,.004,.002,.015}};
  const std::array<ShiftPSFComponent,3> psf{{{{-.02,0},.25},{{0,0},.5},{{.02,0},.25}}};
  const auto original=pixels();
  auto acquired=scene();
  auto lens=prepare_two_deflector_forward(std::move(acquired),cutout,psf);
  need(lens.status()==S::ok&&lens.scene()&&lens.preparation_work().preparation_started==1,
       "independent installed archive creates retained scene");
  const auto retained=lens.retained_payload_bound();
  need(retained.has_value(),"retained payload identity");
  auto means=lens.means(original);
  need(complete(means,16),"complete required means before joint score");

  // Synthetic observation law: C=2I+J/4 electrons^2; d_i=1+i/8 electrons.
  // The chosen data are a fixed vector, not a claimed independently emitted
  // draw or measured asset. Every original pixel is retained, no selection.
  std::vector<double> covariance(256,.25),data(16),residual(16);
  irred::statistics::Metadata metadata;
  metadata.measure="ordered electron product measure";
  metadata.table_identity="synthetic-fixed-original-16-pixel-vector/v1";
  metadata.uncertainty_identity="supplied-fixed-SPD-C=2I+J/4-electrons^2/v1";
  metadata.ordering_provenance="y-outer/x-inner [-48,-16,16,48]; original indices0..15";
  metadata.calibration_provenance="fixed supplied electron response; no marginalized latent prior";
  metadata.dependence_provenance="complete supplied cross-pixel C; no scalar-product independence";
  metadata.source_semantics="synthetic fixed vector; all parent pixels retained; no conditioned selection";
  metadata.input_matrix_convention="covariance in electrons squared";
  for(std::size_t i=0;i<16;++i) {
    covariance[16*i+i]=2.25; data[i]=1+double(i)/8;
    metadata.ordered_ids.push_back("synthetic-original-pixel/"+std::to_string(i));
    need(means.rows()[i].original_index==original[i].original_index,"original source pixel lineage");
    residual[i]=data[i]-*means.rows()[i].electrons;
  }
  const auto ids=metadata.ordered_ids;
  auto gaussian=irred::statistics::prepare_gaussian(covariance,irred::statistics::MatrixKind::covariance,
      std::move(metadata),256,1e-12,irred::numerics::Arithmetic::longdouble_cpu_v1);
  need(gaussian.status()==irred::statistics::DensityStatus::finite,"one retained ordered Gaussian factor");
  const auto score=gaussian.evaluate(residual,ids,1e-12);
  need(score.density.status==irred::statistics::DensityStatus::finite,"complete joint normalized density");
  W sum=0,squares=0;
  for(double r:residual) { sum+=r; squares+=W(r)*r; }
  const W quadratic=squares/2-sum*sum/48;
  const W logdet=15*std::log(2.L)+std::log(6.L);
  const W normalized=-(quadratic+logdet+16*std::log(2*std::numbers::pi_v<W>))/2;
  const W reference_error=1024*std::numeric_limits<W>::epsilon()*
      (1+std::abs(quadratic)+std::abs(logdet)+std::abs(normalized));
  const W total=2e-8L*(1+std::abs(normalized));
  need(reference_error<=.05L*total&&std::abs(W(score.density.log_value)-normalized)+reference_error<=total,
       "independent rank-one inverse/determinant within fixed total log budget");

  ForwardPixelPolicy cap; cap.maximum_field_samples=0;
  const auto unavailable=lens.means(original,cap);
  need(unavailable.status==S::work_limit&&!complete(unavailable,16),"missing required mean causally withholds density");
  // This caller never invokes Gaussian::evaluate on a refused forward batch.
  bool dependent_score_available=false;
  if(complete(unavailable,16)) dependent_score_available=true;
  need(!dependent_score_available,"no joint result after a required mean refusal");

  auto zero=scene(); zero.source.peak_electrons_per_second_per_radian_squared=0;
  zero.uniform_sky_electrons_per_second_per_radian_squared=0;
  const auto exact=prepare_two_deflector_forward(std::move(zero),cutout,psf).means(original,cap);
  need(complete(exact,16)&&exact.work.field_samples_started==0,"structural-zero mean law");
  for(std::size_t i=0;i<16;++i) { need(*exact.rows()[i].electrons==0,"exact zero retained"); residual[i]=data[i]; }
  const auto joint=gaussian.evaluate(residual,ids,1e-12);
  need(joint.density.status==irred::statistics::DensityStatus::finite,"same retained Gaussian factor for second mean law");
  const W joint_quad=38.L/3,product_quad=523.L/18;
  const W joint_log=-(joint_quad+logdet+16*std::log(2*std::numbers::pi_v<W>))/2;
  const W product_log=-(product_quad+16*std::log(2.25L)+16*std::log(2*std::numbers::pi_v<W>))/2;
  const W difference=(16*std::log(2.25L)-logdet+295.L/18)/2;
  need(std::abs(W(joint.quadratic)-joint_quad)<1e-12L&&
       std::abs(W(joint.density.log_value)-joint_log)<2e-8L*(1+std::abs(joint_log))&&
       std::abs((joint_log-product_log)-difference)<1e-14L&&std::abs(difference)>1,
       "joint covariance witness differs meaningfully from scalar marginal product");
  auto reordered=ids; std::swap(reordered[0],reordered[1]);
  need(gaussian.evaluate(residual,reordered,1e-12).density.status==irred::statistics::DensityStatus::incompatible_metadata,
       "ordered observation IDs cannot silently permute");

  std::cout<<std::setprecision(17)
      <<"{\"consumer\":\"installed-fixed-two-deflector-joint-Gaussian/v1\","
      <<"\"model_id\":\""<<two_deflector_forward_id<<"\","
      <<"\"method_id\":\""<<two_deflector_forward_method_id<<"\","
      <<"\"arithmetic_id\":\""<<two_deflector_forward_arithmetic_id<<"\","
      <<"\"mean_unit\":\"electron\",\"covariance_unit\":\"electron^2\","
      <<"\"measure\":\"ordered electron product measure\","
      <<"\"selection\":\"all fixed parent pixels retained\","
      <<"\"status\":\"finite\",\"qualification\":\"synthetic directed control only\","
      <<"\"joint_log_density\":"<<score.density.log_value
      <<",\"joint_quadratic\":"<<score.quadratic
      <<",\"joint_log_determinant\":"<<score.log_determinant
      <<",\"joint_normalization\":"<<score.normalization
      <<",\"joint_backward_residual\":"<<score.backward_residual
      <<",\"joint_forward_sensitivity\":"<<score.estimated_forward_sensitivity
      <<",\"zero_mean_joint_minus_marginal_log\":"<<double(joint_log-product_log)
      <<",\"retained_native_payload_bytes\":"<<*retained
      <<",\"batch_admitted_payload_bytes\":"<<*means.admitted_payload_bytes
      <<",\"preparation_work\":";
  print_work(lens.preparation_work()); std::cout<<",\"batch_work\":";
  print_work(means.work); std::cout<<",\"refused_batch_work\":";
  print_work(unavailable.work); std::cout<<",\"rows\":[";
  for(std::size_t i=0;i<16;++i) {
    if(i) std::cout<<',';
    const auto &r=means.rows()[i]; const auto &e=*r.errors;
    std::cout<<"{\"original_index\":"<<r.original_index<<",\"data_electrons\":"<<data[i]
        <<",\"mean_electrons\":"<<*r.electrons<<",\"status\":\"ok\",\"errors_electrons\":["
        <<e.inner_quadrature_electrons<<','<<e.outer_quadrature_electrons<<','
        <<e.field_psf_arithmetic_electrons<<','<<e.projection_area_electrons<<','<<e.total_electrons
        <<"],\"work\":";
    print_work(r.work); std::cout<<'}';
  }
  std::cout<<"]}\n";
  return 0;
}
