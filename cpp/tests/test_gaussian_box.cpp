// Original synthetic controls. Narrow intervals concern the reported numerical
// Gaussian completion; independent source accuracy is a separate gate.
#include "irred/gaussian_box.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::statistics;
namespace {
int checks=0;
void check(bool b,const char *message) {++checks;if(!b)throw std::runtime_error(message);}
void near(double actual,double expected,const char *message) {
  check(std::abs(actual-expected)<=2e-12*(1+std::abs(expected)),message);
}
bool contains(BoxInterval a,double x) {return a.lower<=x&&x<=a.upper;}
Metadata metadata() {
  Metadata m;m.ordered_ids={"r0","r1"};m.measure="product d(residual)";
  m.ordering_provenance="synthetic declared order";m.calibration_provenance="none synthetic";
  m.dependence_provenance="supplied synthetic full covariance";return m;
}
DesignMetadata design_metadata() {
  return {{"original0","original2"},{"u0","u2"},{},"u",
          "synthetic active0,2;fixed1=0","supplied synthetic full covariance"};
}
DesignProfile profile(std::vector<double> c={1,0,0,1},
                      std::vector<double> x={1,0,0,1}) {
  auto m=metadata();
  auto g=prepare_gaussian(c,MatrixKind::covariance,m,4,1e-10,
                          irred::numerics::Arithmetic::longdouble_cpu_v1);
  check(g.status()==DensityStatus::finite,"synthetic Gaussian prepared");
  return DesignProfile::prepare(std::move(g),x,m.ordered_ids,design_metadata());
}
BoxSupport support(double lo0=-8,double hi0=8,double lo2=-10,double hi2=10) {
  return {{"original0","original2"},{lo0,lo2},{hi0,hi2},
          "product d(original0)d(original2)","normalized uniform finite box",
          "original1 literal0 point mass outside active2 measure"};
}
GaussianBoxResult evaluate(GaussianBox &b,BoxMarginalRequest r={0,.5},BoxPolicy p={}) {
  const std::vector<double> y={0,0};return b.evaluate(y,metadata().ordered_ids,r,p);
}
} // namespace
int main() {try {
  auto d=profile();auto b=GaussianBox::prepare(std::move(d),support());
  check(b.status()==DensityStatus::finite&&d.status()==DensityStatus::invalid_input,
        "successful box consumes retained design");
  auto r=evaluate(b);
  check(r.status==DensityStatus::finite,"broad separable box median admitted");
  check(contains(r.requested_quantile,0),"analytic symmetric marginal median");
  check(r.requested_quantile.upper-r.requested_quantile.lower<=1e-8,
        "original coordinate quantile width");
  near(r.minimum_quadratic,0,"exact full square q minimum");
  near(r.log_design_precision_determinant,0,"exact identity precision logdet");
  near(r.source_covariance_log_determinant,0,"exact identity covariance logdet");
  near(r.unboxed_variance[0],1,"first analytic marginal variance");
  near(r.unboxed_variance[1],1,"second analytic marginal variance");
  // Q(8)=6.2209605742717841235e-16, Q(10)=7.619853024160526066e-24.
  // Exact separable mass is (1-2Q8)(1-2Q10), below1 and above union lower.
  const double exact_mass=1-2*6.2209605742717841235e-16;
  check(contains(r.box_probability,exact_mass),"analytic product mass enclosed");
  check(r.excluded_mass_upper>0&&r.log_box_probability.lower<0,
        "positive omitted tails never silently zero");
  const double exact_log_integral=1.83787706640934548356-2*6.2209605742717841235e-16;
  check(contains(r.log_relative_box_integral,exact_log_integral),
        "analytic relative normalization and explicit missing mass");
  check(contains(r.log_prior_volume,std::log(320.)),"two active widths only; fixed point excluded");
  check(r.lower_endpoint_cdf.lower<=.5&&r.upper_endpoint_cdf.upper>=.5,
        "median endpoint CDF witnesses");
  auto lower=evaluate(b,{0,.025}),upper=evaluate(b,{0,.975});
  check(lower.status==DensityStatus::finite&&upper.status==DensityStatus::finite,
        "two tail quantiles admitted");
  // Independently tabulated standard Gaussian quantile. Removed8sigma mass
  // alters these quantiles by less3e-14; expected interval width is1e-8.
  check(contains(lower.requested_quantile,-1.9599639845400542355),
        "analytic lower Gaussian quantile limit");
  check(contains(upper.requested_quantile,1.9599639845400542355),
        "analytic upper Gaussian quantile limit");
  auto correlated_profile=profile({1,.6,.6,1});
  auto correlated=GaussianBox::prepare(std::move(correlated_profile),support(-9,9,-9,9));
  auto cr=evaluate(correlated);
  check(cr.status==DensityStatus::finite&&contains(cr.requested_quantile,0),
        "correlated centrally symmetric median and correlation valid tails");
  near(cr.unboxed_variance[0],1,"correlated first diagonal variance");
  near(cr.unboxed_variance[1],1,"correlated second diagonal variance");
  near(cr.log_design_precision_determinant,-std::log(.64),"independent2x2 cofactor determinant");
  // Original-coordinate Jacobian: X'=X diag(2,.5), beta'=diag(.5,2)beta.
  auto scaled_profile=profile({1,0,0,1},{2,0,0,.5});
  auto scaled=GaussianBox::prepare(std::move(scaled_profile),support(-4,4,-20,20));
  auto sr=evaluate(scaled);
  check(sr.status==DensityStatus::finite&&contains(sr.requested_quantile,0),"unit scaled median");
  near(sr.unboxed_variance[0],.25,"first unit Jacobian variance");
  near(sr.unboxed_variance[1],4,"second unit Jacobian variance");
  near(sr.log_prior_volume.lower,r.log_prior_volume.lower,"unit product Jacobian unity");
  auto moved=std::move(b);
  check(b.status()==DensityStatus::invalid_input&&evaluate(moved).status==DensityStatus::finite,
        "move invalidates old owner");moved=std::move(moved);
  check(moved.status()==DensityStatus::finite,"self move retains owner");
  auto invalid_profile=profile();auto invalid=support();invalid.lower[0]=invalid.upper[0];
  check(GaussianBox::prepare(std::move(invalid_profile),invalid).status()==DensityStatus::invalid_input&&
        invalid_profile.status()==DensityStatus::finite,"zero width refuses without consuming design");
  invalid=support();std::swap(invalid.ordered_parameter_ids[0],invalid.ordered_parameter_ids[1]);
  check(GaussianBox::prepare(std::move(invalid_profile),invalid).status()==DensityStatus::incompatible_metadata,
        "exact active source axis order required");
  auto cut_profile=profile();auto cut=GaussianBox::prepare(std::move(cut_profile),support(0,8,-10,10));
  auto refused=evaluate(cut);
  check(refused.status==DensityStatus::numerical_failure&&refused.unboxed_mean.empty(),
        "boundary cutting proper box refused at frozen tail allocation");
  auto negative_profile=profile();auto negative=GaussianBox::prepare(std::move(negative_profile),support(1,8,-10,10));
  check(evaluate(negative).status==DensityStatus::numerical_failure,
        "negative signed distance never squared into small tail");
  auto almost_profile=profile();auto almost=GaussianBox::prepare(std::move(almost_profile),support(-6,8,-10,10));
  check(evaluate(almost).status==DensityStatus::numerical_failure,
        "informative mass cannot evade stricter log tail width");
  BoxPolicy quota;quota.maximum_cdf_nodes=2;
  check(evaluate(moved,{0,.975},quota).numerical_status==irred::numerics::Status::work_limit,
        "actual Simpson endpoints charged to node quota");
  quota={};quota.maximum_cdf_evaluations=1;
  auto fq=evaluate(moved,{0,.5},quota);
  check(fq.numerical_status==irred::numerics::Status::work_limit&&fq.cdf_evaluations>0&&fq.unboxed_mean.empty(),
        "failed CDF attempts charged and payload withheld");
  quota={};quota.maximum_bisections=1;
  check(evaluate(moved,{0,.975},quota).numerical_status==irred::numerics::Status::work_limit,
        "insufficient inversion refuses");
  quota={};quota.maximum_quantile_width=1e-30;
  check(evaluate(moved,{0,.5},quota).status==DensityStatus::numerical_failure,
        "unrepresentable requested quantile width refuses");
  const auto old=std::fegetround();std::fesetround(FE_UPWARD);
  check(evaluate(moved).status==DensityStatus::unsupported_domain,"rounding environment contract");
  std::fesetround(old);
  std::cout<<checks<<" Gaussian box checks passed\n";return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
