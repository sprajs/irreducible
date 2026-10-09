#include "irred/decaying_matter.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using W = long double;
using S = numerics::Status;
using State = std::array<W, 2>; // H_i*(t-t_i), rho_r/rho_ci * exp(4x)
struct Solver {
  const DecayingMatterModel &m;
  const DecayingMatterPolicy &p;
  DecayingMatterTrajectory &out;
  S failure = S::ok;
  bool densities(W x, State y, std::array<W, 4> &d, W &e) {
    if (!std::isfinite(y[0]) || !std::isfinite(y[1]) || y[0] < 0 ||
        y[1] < 0 || y[0] > 1e6L || (m.parent_fraction > 0 && m.decay_rate_over_initial_hubble * y[0] > 600)) {
      failure = S::outside_domain; return false;
    }
    d = {m.stable_matter_fraction * std::exp(-3*x),
         m.parent_fraction * std::exp(-3*x-m.decay_rate_over_initial_hubble*y[0]),
         y[1] * std::exp(-4*x), m.lambda_fraction};
    const W sum = d[0]+d[1]+d[2]+d[3];
    e = std::sqrt(sum);
    if (!(e > 0) || !std::isfinite(e) ||
        (m.parent_fraction > 0 && !(d[1] > 0)) ||
        (m.stable_matter_fraction > 0 && !(d[0] > 0)) ||
        (y[1] > 0 && !(d[2] > 0))) {
      failure = S::outside_domain; return false;
    }
    return true;
  }
  bool rhs(W x, State y, State &v) {
    if (out.rhs_evaluations >= p.maximum_rhs_evaluations) {
      failure = S::work_limit; return false;
    }
    ++out.rhs_evaluations;
    std::array<W, 4> d{}; W e = 0;
    if (!densities(x,y,d,e)) return false;
    v = {1/e, m.decay_rate_over_initial_hubble*m.parent_fraction*
                  std::exp(x-m.decay_rate_over_initial_hubble*y[0])/e};
    if (!std::isfinite(v[0]) || !std::isfinite(v[1]) ||
        (m.parent_fraction > 0 && m.decay_rate_over_initial_hubble > 0 && !(v[1] > 0))) {
      failure = S::outside_domain; return false;
    }
    return true;
  }
  bool step(W x, State y, W h, State &z) {
    State a{},b{},c{},d{};
    auto add = [y](State v, W k) { return State{y[0]+k*v[0],y[1]+k*v[1]}; };
    if (!rhs(x,y,a) || !rhs(x+h/2,add(a,h/2),b) ||
        !rhs(x+h/2,add(b,h/2),c) || !rhs(x+h,add(c,h),d)) return false;
    for (unsigned j=0;j<2;++j) z[j]=y[j]+h*(a[j]+2*b[j]+2*c[j]+d[j])/6;
    return true;
  }
  bool row(double scale, W x, State y, State error) {
    std::array<W,4> d{}; W e=0;
    if (!densities(x,y,d,e)) return false;
    const W total=e*e;
    DecayingMatterRow r;
    r.scale_factor=scale; r.elapsed_initial_hubble_time=y[0];
    r.elapsed_seconds=y[0]/m.initial_hubble_per_second;
    r.expansion_over_initial_hubble=e; r.hubble_per_second=e*m.initial_hubble_per_second;
    if (!std::isfinite(r.elapsed_seconds) || !std::isfinite(r.hubble_per_second) ||
        !std::isnormal(r.hubble_per_second) ||
        (y[0] > 0 && !std::isnormal(r.elapsed_seconds))) { failure=S::outside_domain; return false; }
    r.stable_density_over_initial_critical=d[0]; r.parent_density_over_initial_critical=d[1];
    r.daughter_density_over_initial_critical=d[2]; r.lambda_density_over_initial_critical=d[3];
    r.total_density_over_initial_critical=total;
    r.stable_matter_fraction=d[0]/total; r.parent_fraction=d[1]/total;
    r.daughter_radiation_fraction=d[2]/total; r.lambda_fraction=d[3]/total;
    r.total_equation_of_state=r.daughter_radiation_fraction/3-r.lambda_fraction;
    r.deceleration=(1+3*r.total_equation_of_state)/2;
    r.transfer_over_hubble_total_density=m.decay_rate_over_initial_hubble*r.parent_fraction/e;
    if ((d[0] > 0 && !(r.stable_matter_fraction > 0)) ||
        (d[1] > 0 && !(r.parent_fraction > 0)) ||
        (d[2] > 0 && !(r.daughter_radiation_fraction > 0)) ||
        (d[3] > 0 && !(r.lambda_fraction > 0)) ||
        (m.decay_rate_over_initial_hubble > 0 && d[1] > 0 &&
         !(r.transfer_over_hubble_total_density > 0))) {
      failure=S::outside_domain; return false;
    }
    const W transfer=m.decay_rate_over_initial_hubble*d[1]/e;
    const W derivative=-3*d[0]+(-3*d[1]-transfer)+(-4*d[2]+transfer);
    r.continuity_residual=(derivative+3*(d[0]+d[1])+4*d[2])/total;
    r.elapsed_initial_hubble_time_estimate=error[0]; r.comoving_daughter_estimate=error[1];
    r.expansion_relative_estimate=(m.decay_rate_over_initial_hubble*d[1]*error[0]+
        std::exp(-4*x)*error[1])/(2*total)+32*std::numeric_limits<W>::epsilon();
    out.rows.push_back(r); return true;
  }
};
} // namespace
DecayingMatterTrajectory evolve_decaying_matter(const DecayingMatterModel &m,
    std::span<const double> scales, DecayingMatterPolicy p) {
  DecayingMatterTrajectory out; out.initial_state=m;
  const double inputs[]{m.initial_scale_factor,m.initial_hubble_per_second,
    m.stable_matter_fraction,m.parent_fraction,m.daughter_radiation_fraction,
    m.lambda_fraction,m.decay_rate_over_initial_hubble,p.absolute_tolerance,p.relative_tolerance};
  for (double v:inputs) if (!std::isfinite(v)) { out.status=S::nonfinite_input; return out; }
  if (std::fegetround()!=FE_TONEAREST || !(p.absolute_tolerance>0) || p.relative_tolerance<0 ||
      p.absolute_tolerance<128*std::numeric_limits<W>::epsilon() ||
      (p.relative_tolerance>0 && p.relative_tolerance<128*std::numeric_limits<W>::epsilon())) return out;
  if (!(m.initial_scale_factor>0) || !(m.initial_hubble_per_second>0) ||
      m.stable_matter_fraction<0 || m.parent_fraction<0 ||
      m.daughter_radiation_fraction<0 || m.lambda_fraction<0 ||
      m.decay_rate_over_initial_hubble<0 || m.decay_rate_over_initial_hubble>1e4) {
    out.status=S::outside_domain; return out;
  }
  const W closure=W(m.stable_matter_fraction)+m.parent_fraction+m.daughter_radiation_fraction+m.lambda_fraction;
  if (std::abs(closure-1)>8*std::numeric_limits<double>::epsilon()) {
    out.status=S::outside_domain; return out;
  }
  // Requested retained allocation plus operation stack allowance. No caches,
  // endpoint-sort copies or hidden borrowed resources; allocator/RSS excluded.
  constexpr std::size_t overhead=sizeof(DecayingMatterTrajectory)+8192;
  if (scales.size()>4096 || scales.size()>p.maximum_rows ||
      p.maximum_native_bytes<overhead ||
      scales.size()>(p.maximum_native_bytes-overhead)/sizeof(DecayingMatterRow)) {
    out.status=S::work_limit; return out;
  }
  W previous=m.initial_scale_factor;
  for (double a:scales) {
    if (!std::isfinite(a)) { out.status=S::nonfinite_input; return out; }
    if (!(a>0) || a<previous || std::log(W(a)/m.initial_scale_factor)>16) {
      out.status=S::outside_domain; return out;
    }
    previous=a;
  }
  out.rows.reserve(scales.size());
  Solver solver{m,p,out};
  State y{0,m.daughter_radiation_fraction}, error{0,0};
  W x=0,h=0.02L;
  const W span=scales.empty()?0:std::log(W(scales.back())/m.initial_scale_factor);
  for (double a:scales) {
    const W target=std::log(W(a)/m.initial_scale_factor);
    while (x<target) {
      h=std::min(h,target-x);
      if (!(x+h>x)) { solver.failure=S::conditioning_budget_exceeded; break; }
      State coarse{},half{},fine{};
      if (!solver.step(x,y,h,coarse) || !solver.step(x,y,h/2,half) ||
          !solver.step(x+h/2,half,h/2,fine)) break;
      W ratio=0; State local{};
      for (unsigned j=0;j<2;++j) {
        local[j]=std::abs(fine[j]-coarse[j])/15;
        const W budget=(p.absolute_tolerance+p.relative_tolerance*std::abs(fine[j]))*
            h/std::max(span,1.0L);
        ratio=std::max(ratio,local[j]/budget);
      }
      if (ratio<=1) {
        y=fine; x=(h==target-x)?target:x+h;
        for(unsigned j=0;j<2;++j) error[j]+=local[j];
        ++out.accepted_steps;
      } else ++out.rejected_steps;
      h*=ratio==0?2:std::clamp(0.8L*std::pow(ratio,-0.2L),0.1L,2.0L);
    }
    if (solver.failure!=S::ok || !solver.row(a,target,y,error)) {
      out.rows.clear(); out.status=solver.failure; return out;
    }
  }
  out.status=S::ok; return out;
}
} // namespace irred::cosmology
