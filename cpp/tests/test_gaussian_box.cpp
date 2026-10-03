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
  // Every observation remains in the likelihood even when its design row is
  // exactly zero. The small isolated variance affects logdet(C), not H.
  auto isolated_metadata=metadata();
  isolated_metadata.ordered_ids.push_back("isolated-r2");
  const std::vector<double> isolated_covariance={1,0,0,0,1,0,0,0,0x1p-40};
  const std::vector<double> isolated_design={1,0,0,1,0,0};
  auto isolated_gaussian=prepare_gaussian(isolated_covariance,MatrixKind::covariance,
      isolated_metadata,9,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
  check(isolated_gaussian.status()==DensityStatus::finite,"dyadic isolated covariance admitted");
  const std::vector<double> isolated_zero(3,0);
  const auto unnecessary_solve=isolated_gaussian.evaluate(isolated_zero,
      isolated_metadata.ordered_ids,1e-10);
  // The preserved refusal is specific to the qualified 64-bit mantissa
  // profile. A wider mantissa can admit this separate inverse request.
  if constexpr(std::numeric_limits<long double>::digits==64)
    check(unnecessary_solve.density.status==DensityStatus::numerical_failure&&
        unnecessary_solve.density.numerical_status==irred::numerics::Status::conditioning_budget_exceeded,
        "full covariance inverse sensitivity is a separate unchanged solve gate");
  auto isolated_profile=DesignProfile::prepare(std::move(isolated_gaussian),
      isolated_design,isolated_metadata.ordered_ids,design_metadata());
  auto isolated_box=GaussianBox::prepare(std::move(isolated_profile),
      support(.125-8,.125+8,-.25-8,-.25+8));
  const std::vector<double> isolated_y={.125,-.25,0};
  const auto isolated=isolated_box.evaluate(isolated_y,isolated_metadata.ordered_ids,{0,.5});
  check(isolated.status==DensityStatus::finite&&isolated.stage==BoxStage::complete&&
      isolated.completion_step==BoxCompletionStep::complete&&!isolated.completion_parameter_index,
      "retained determinant completes without an artificial inverse solve");
  near(isolated.unboxed_mean[0],.125,"dyadic first fitted mean");
  near(isolated.unboxed_mean[1],-.25,"dyadic second fitted mean");
  near(isolated.unboxed_variance[0],1,"isolated first variance");
  near(isolated.unboxed_variance[1],1,"isolated second variance");
  near(isolated.minimum_quadratic,0,"isolated zero residual contributes exact zero q");
  near(isolated.log_design_precision_determinant,0,"isolated row contributes no design information");
  near(isolated.source_covariance_log_determinant,-40*std::log(2.),
      "isolated row retains its exact dyadic determinant contribution");
  check(contains(isolated.requested_quantile,.125),"isolated symmetric conditional median");
  const auto midpoint=[](BoxInterval v) {return (v.lower+v.upper)/2;};
  // Independent scalar facts from the refined Decimal control: log(2pi) and
  // Q(8). log[(1-2Q8)^2] differs from -4Q8 by less than 2e-30.
  const double isolated_log_evidence=12*std::log(2.)-.5*1.83787706640934548356
      -4*6.2209605742717841235e-16;
  near(midpoint(isolated.log_observation_normalized_evidence),isolated_log_evidence,
      "all three observation dimensions and retained tiny determinant normalize evidence");
  const std::vector<double> isolated_nonzero_y={.125,-.25,0x1p-20};
  const auto isolated_nonzero=isolated_box.evaluate(isolated_nonzero_y,
      isolated_metadata.ordered_ids,{0,.5});
  check(isolated_nonzero.status==DensityStatus::finite,"isolated nonzero residual admitted");
  near(isolated_nonzero.minimum_quadratic,1,"retained isolated residual adds exactly one to q minimum");
  near(midpoint(isolated_nonzero.log_observation_normalized_evidence),isolated_log_evidence-.5,
      "isolated nonzero residual reduces normalized evidence by one half");
  BoxPolicy tight_profile;tight_profile.design.maximum_forward_sensitivity=1e-18;
  const auto inherited_refusal=isolated_box.evaluate(isolated_y,
      isolated_metadata.ordered_ids,{0,.5},tight_profile);
  check(inherited_refusal.status==DensityStatus::numerical_failure&&
      inherited_refusal.numerical_status==irred::numerics::Status::conditioning_budget_exceeded&&
      inherited_refusal.completion_step==BoxCompletionStep::profile_evaluation&&
      !inherited_refusal.completion_parameter_index&&inherited_refusal.stage==BoxStage::unassessed&&
      !inherited_refusal.gaussian_completion_available&&!inherited_refusal.endpoint_margins_available&&
      !inherited_refusal.rectangle_enclosure_available&&!inherited_refusal.normalization_enclosures_available&&
      !inherited_refusal.quantile_enclosure_available&&!inherited_refusal.endpoint_cdf_enclosures_available&&
      inherited_refusal.unboxed_mean.empty()&&inherited_refusal.unboxed_variance.empty(),
      "genuine inherited profile sensitivity refusal preserves attempted gate without outputs");
  auto lower=evaluate(b,{0,.025}),upper=evaluate(b,{0,.975});
  check(lower.status==DensityStatus::finite&&upper.status==DensityStatus::finite,
        "two tail quantiles admitted");
  // Independent90/120-digit Decimal Machin-pi/integrated Gaussian power
  // series/bisection reference, refined difference<3e-70 for these facts.
  // These are the finite[-8,8] marginal quantiles, not unboxed replacements.
  check(contains(lower.requested_quantile,-1.9599639845400441236),
        "analytic lower finite box quantile");
  check(contains(upper.requested_quantile,1.9599639845400441236),
        "analytic upper finite box quantile");
  auto asym_profile=profile();
  auto asym=GaussianBox::prepare(std::move(asym_profile),support(-7,9,-10,11));
  const auto ar=evaluate(asym);
  check(ar.status==DensityStatus::finite&&contains(ar.requested_quantile,1.6040070129182506649e-12),
        "independent asymmetric singly truncated marginal median");
  check(contains(ar.box_probability,.9999999999987201873432553244),
        "independent separable asymmetric finite mass");
  auto correlated_profile=profile({1,.6,.6,1});
  auto correlated=GaussianBox::prepare(std::move(correlated_profile),support(-9,9,-9,9));
  auto cr=evaluate(correlated);
  check(cr.status==DensityStatus::finite&&contains(cr.requested_quantile,0),
        "correlated centrally symmetric median and correlation valid tails");
  near(cr.unboxed_variance[0],1,"correlated first diagonal variance");
  near(cr.unboxed_variance[1],1,"correlated second diagonal variance");
  near(cr.log_design_precision_determinant,-std::log(.64),"independent2x2 cofactor determinant");
  near((cr.log_relative_box_integral.lower+cr.log_relative_box_integral.upper)/2,
       1.61473351509513572779,"independent high precision cofactor normalization; source arithmetic comparison");
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
        "move invalidates old owner");auto *self=&moved;moved=std::move(*self);
  check(moved.status()==DensityStatus::finite,"self move retains owner");
  auto invalid_profile=profile();auto invalid=support();invalid.lower[0]=invalid.upper[0];
  check(GaussianBox::prepare(std::move(invalid_profile),invalid).status()==DensityStatus::invalid_input&&
        invalid_profile.status()==DensityStatus::finite,"zero width refuses without consuming design");
  invalid=support();std::swap(invalid.ordered_parameter_ids[0],invalid.ordered_parameter_ids[1]);
  check(GaussianBox::prepare(std::move(invalid_profile),invalid).status()==DensityStatus::incompatible_metadata,
        "exact active source axis order required");
  auto cut_profile=profile();auto cut=GaussianBox::prepare(std::move(cut_profile),support(0,8,-10,10));
  auto refused=evaluate(cut);
  check(refused.status==DensityStatus::numerical_failure&&refused.gaussian_completion_available&&
        refused.endpoint_margins_available&&!refused.rectangle_enclosure_available&&
        !refused.normalization_enclosures_available&&!refused.quantile_enclosure_available&&
        refused.stage==BoxStage::endpoint_margins&&refused.unboxed_mean.size()==2&&
        refused.completion_step==BoxCompletionStep::complete&&!refused.completion_parameter_index,
        "boundary refusal retains earned completion and every tail margin");
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
  check(fq.numerical_status==irred::numerics::Status::work_limit&&fq.cdf_evaluations>0&&
        fq.normalization_enclosures_available&&!fq.quantile_enclosure_available&&
        fq.gaussian_completion_available&&fq.stage==BoxStage::normalization_enclosures,
        "late inversion refusal retains earned normalization with explicit availability");
  const std::vector<double> invalid_y={std::numeric_limits<double>::quiet_NaN(),0};
  const auto no_completion=moved.evaluate(invalid_y,metadata().ordered_ids,{0,.5});
  check(no_completion.status==DensityStatus::invalid_input&&no_completion.stage==BoxStage::unassessed&&
        !no_completion.gaussian_completion_available&&no_completion.unboxed_mean.empty()&&
        no_completion.completion_step==BoxCompletionStep::profile_evaluation&&
        !no_completion.completion_parameter_index,
        "invalid input has no fabricated zero completion");
  quota={};quota.maximum_bisections=1;
  check(evaluate(moved,{0,.975},quota).numerical_status==irred::numerics::Status::work_limit,
        "insufficient inversion refuses");
  quota={};quota.maximum_quantile_width=1e-30;
  check(evaluate(moved,{0,.5},quota).status==DensityStatus::numerical_failure,
        "unrepresentable requested quantile width refuses");
  auto large_profile=profile();
  auto large=GaussianBox::prepare(std::move(large_profile),support(1e9-10,1e9+10,-10,10));
  const std::vector<double> large_y={1e9,0};
  check(large.evaluate(large_y,metadata().ordered_ids,{0,.5}).status==DensityStatus::numerical_failure,
        "binary64 returned coordinate width independently admitted");
  const auto old=std::fegetround();std::fesetround(FE_UPWARD);
  check(evaluate(moved).status==DensityStatus::unsupported_domain,"rounding environment contract");
  std::fesetround(old);
  std::cout<<checks<<" Gaussian box checks passed\n";return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
