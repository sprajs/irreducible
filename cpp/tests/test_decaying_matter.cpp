// Independent reference integrates physical densities and a in proper time,
// not the production log-scale/comoving-daughter coordinates. No external data.
#include "irred/decaying_matter.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
using namespace irred::cosmology;
using S = irred::numerics::Status;
using W = long double;
namespace {
unsigned checks=0;
void need(bool b,const char *s) { ++checks; if(!b) throw std::runtime_error(s); }
void near(W a,W b,W tol,const char *s) { need(std::abs(a-b)<=tol*std::max(1.0L,std::abs(b)),s); }
using Y=std::array<W,4>; // a/a_i, stable, parent, radiation; rho/rho_ci
Y reference(const DecayingMatterModel &m,W tau,unsigned n) {
  Y y{1,m.stable_matter_fraction,m.parent_fraction,m.daughter_radiation_fraction};
  auto rhs=[&](Y z) {
    const W e=std::sqrt(z[1]+z[2]+z[3]+m.lambda_fraction);
    const W q=m.decay_rate_over_initial_hubble*z[2];
    return Y{z[0]*e,-3*e*z[1],-3*e*z[2]-q,-4*e*z[3]+q};
  };
  const W h=tau/n;
  for(unsigned i=0;i<n;++i) {
    auto add=[&](Y k,W s) { Y z{}; for(unsigned j=0;j<4;++j)z[j]=y[j]+s*k[j];return z; };
    auto k1=rhs(y),k2=rhs(add(k1,h/2)),k3=rhs(add(k2,h/2)),k4=rhs(add(k3,h));
    for(unsigned j=0;j<4;++j)y[j]+=h*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;
  }
  return y;
}
}
int main() {
 try {
  DecayingMatterModel m; m.initial_scale_factor=.2; m.initial_hubble_per_second=2e-18;
  // Pure nondecaying matter/radiation/Lambda analytic E and proper-time limits.
  for(unsigned species=0;species<3;++species) {
    m.parent_fraction=0;m.stable_matter_fraction=species==0?1:0;
    m.daughter_radiation_fraction=species==1?1:0;m.lambda_fraction=species==2?1:0;
    const double scales[]{.2,.4,.8}; auto out=evolve_decaying_matter(m,scales);
    need(out.status==S::ok,"analytic species accepted");
    for(const auto &r:out.rows) {
      W b=W(r.scale_factor)/m.initial_scale_factor;
      W e=species==0?std::pow(b,-1.5L):species==1?1/(b*b):1;
      W t=species==0?2*(std::pow(b,1.5L)-1)/3:species==1?(b*b-1)/2:std::log(b);
      near(r.expansion_over_initial_hubble,e,2e-15L,"analytic expansion");
      near(r.elapsed_initial_hubble_time,t,3e-9L,"analytic elapsed time");
      near(r.deceleration,species==0?.5L:species==1?1:-1,2e-15L,"analytic q");
    }
  }
  m.stable_matter_fraction=.2;m.parent_fraction=.3;m.daughter_radiation_fraction=.1;m.lambda_fraction=.4;
  m.decay_rate_over_initial_hubble=0;
  const double scales[]{.2,.3,.6,1.2};
  auto zero=evolve_decaying_matter(m,scales); need(zero.status==S::ok,"Gamma=0 mixture");
  for(const auto &r:zero.rows) {
    W b=W(r.scale_factor)/m.initial_scale_factor;
    near(r.total_density_over_initial_critical,.5L/std::pow(b,3)+.1L/std::pow(b,4)+.4L,2e-15L,"Gamma zero analytic densities");
    need(r.transfer_over_hubble_total_density==0,"zero transfer exact");
  }
  for(double g:{.01,1.,10.}) {
    m.decay_rate_over_initial_hubble=g;
    auto out=evolve_decaying_matter(m,scales); need(out.status==S::ok,"decay mixture accepted");
    DecayingMatterPolicy tight;tight.absolute_tolerance=1e-13;tight.relative_tolerance=1e-11;
    auto refined=evolve_decaying_matter(m,scales,tight); need(refined.status==S::ok,"decay refinement accepted");
    W prev_comoving=.1;
    for(std::size_t i=0;i<out.rows.size();++i) {
      const auto &r=out.rows[i]; W b=W(r.scale_factor)/m.initial_scale_factor;
      near(r.parent_density_over_initial_critical*std::pow(b,3)*
           std::exp(g*r.elapsed_initial_hubble_time),m.parent_fraction,2e-15L,"proper-time particle survival");
      W comoving=r.daughter_density_over_initial_critical*std::pow(b,4);
      need(comoving>=prev_comoving && r.parent_fraction>0 && r.daughter_radiation_fraction>0,"positive survival and monotone daughter energy");
      prev_comoving=comoving;
      near(r.stable_matter_fraction+r.parent_fraction+r.daughter_radiation_fraction+r.lambda_fraction,1,2e-15L,"closure");
      need(std::abs(r.continuity_residual)<2e-15L,"continuous total conservation");
      near(r.expansion_over_initial_hubble,refined.rows[i].expansion_over_initial_hubble,3e-9L,"native refinement E");
      near(r.elapsed_initial_hubble_time,refined.rows[i].elapsed_initial_hubble_time,3e-9L,"native refinement time");
      // Named independent-reference budget 2e-8 absolute/relative; proper-time
      // doubling difference must occupy less than 5% of that allocation.
      auto ref=reference(m,r.elapsed_initial_hubble_time,16384);
      auto ref2=reference(m,r.elapsed_initial_hubble_time,32768);
      for(unsigned j=0;j<4;++j) near(ref[j],ref2[j],1e-9L,"independent proper-time refinement");
      near(ref2[0],b,2e-8L,"independent proper-time scale");
      near(ref2[1],r.stable_density_over_initial_critical,2e-8L,"independent stable density");
      near(ref2[2],r.parent_density_over_initial_critical,2e-8L,"independent parent density");
      near(ref2[3],r.daughter_density_over_initial_critical,2e-8L,"independent daughter density");
    }
    need(out.rhs_evaluations>0 && out.accepted_steps>0,"accounted adaptive work");
  }
  // Pure initial parent: radiation emerges, parent depletes, late radiation q.
  m.stable_matter_fraction=0;m.parent_fraction=1;m.daughter_radiation_fraction=0;m.lambda_fraction=0;
  m.decay_rate_over_initial_hubble=3;
  const double pure_scales[]{.2,.22,.4,1.}; auto pure=evolve_decaying_matter(m,pure_scales);
  need(pure.status==S::ok,"pure decay accepted");
  need(pure.rows.front().daughter_radiation_fraction==0 && pure.rows[1].daughter_radiation_fraction>0,"daughter initial boundary and emergence");
  need(pure.rows.back().parent_fraction<1e-6L && pure.rows.back().deceleration>.999L,"late radiation asymptote");
  for(const auto &r:pure.rows) {
    auto ref=reference(m,r.elapsed_initial_hubble_time,32768);
    near(ref[0],W(r.scale_factor)/m.initial_scale_factor,2e-8L,"pure decay proper-time reference");
  }
  auto copy=pure; pure.rows.clear(); need(!copy.rows.empty(),"retained batch independent lifetime");
  auto moved=std::move(copy); need(!moved.rows.empty(),"retained move lifetime");
  const double duplicate[]{.2,.2,.4,.4}; auto duplicates=evolve_decaying_matter(m,duplicate);
  need(duplicates.status==S::ok && duplicates.rows[2].elapsed_seconds==duplicates.rows[3].elapsed_seconds,"ordered duplicates retained");
  DecayingMatterPolicy cap;cap.maximum_rhs_evaluations=1;
  auto fail=evolve_decaying_matter(m,scales,cap);
  need(fail.status==S::work_limit && fail.rows.empty() && fail.rhs_evaluations==1,"callback cap and atomic refusal");
  cap={};cap.maximum_native_bytes=1;need(evolve_decaying_matter(m,scales,cap).status==S::work_limit,"payload cap");
  cap={};cap.maximum_rows=2;need(evolve_decaying_matter(m,scales,cap).status==S::work_limit,"row cap");
  const double backwards[]{.3,.2};need(evolve_decaying_matter(m,backwards).status==S::outside_domain,"backward trajectory refusal");
  const double past[]{.1};need(evolve_decaying_matter(m,past).status==S::outside_domain,"past anchor refusal");
  const double distant[]{1e8};need(evolve_decaying_matter(m,distant).status==S::outside_domain,"scale domain cap");
  m.parent_fraction=.9;need(evolve_decaying_matter(m,scales).status==S::outside_domain,"no closure renormalization");m.parent_fraction=1;
  m.decay_rate_over_initial_hubble=-1;need(evolve_decaying_matter(m,scales).status==S::outside_domain,"negative decay refusal");
  m.decay_rate_over_initial_hubble=1e5;need(evolve_decaying_matter(m,scales).status==S::outside_domain,"rate cap");
  m.decay_rate_over_initial_hubble=1e4; const double over_lifetime[]{1.};
  need(evolve_decaying_matter(m,over_lifetime).status==S::outside_domain,"proper lifetime domain cap");
  m.decay_rate_over_initial_hubble=0;m.initial_hubble_per_second=0;
  need(evolve_decaying_matter(m,scales).status==S::outside_domain,"physical H required");
  m.initial_hubble_per_second=std::numeric_limits<double>::quiet_NaN();
  need(evolve_decaying_matter(m,scales).status==S::nonfinite_input,"nonfinite physical state");
  m.initial_hubble_per_second=2e-18;const double nan[]{std::numeric_limits<double>::quiet_NaN()};
  need(evolve_decaying_matter(m,nan).status==S::nonfinite_input,"nonfinite requested scale");
  cap={};cap.absolute_tolerance=0;need(evolve_decaying_matter(m,scales,cap).status==S::invalid_input,"positive absolute control required");
  std::cout<<checks<<" decaying-matter scientific/resource checks passed\n";
 } catch(const std::exception &e) { std::cerr<<e.what()<<'\n';return 1; }
}
