#include "irred/bianchi_i.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using W=long double;
using S=numerics::Status;
constexpr W eps=std::numeric_limits<double>::epsilon();
struct Context {
  const BianchiIModel &m;
  W sigma, smax;
  mutable S failure=S::ok;
  W scaled_density(W a) const {
    const W a2=a*a,a3=a2*a;
    return m.matter_fraction*a3+m.radiation_fraction*a2+m.lambda_fraction*a3*a3+sigma;
  }
};
double integrand(double x,const void *ptr) {
  const auto &c=*static_cast<const Context *>(ptr);
  const W a=std::exp(W(x)),value=c.smax/std::sqrt(c.scaled_density(a));
  const double result=static_cast<double>(value);
  if (!(value>0) || !std::isnormal(result)) {
    c.failure=S::outside_domain;
    return std::numeric_limits<double>::quiet_NaN();
  }
  return result;
}
bool finite_model(const BianchiIModel &m) {
  return std::isfinite(m.anchor_hubble_per_second) && std::isfinite(m.matter_fraction) &&
      std::isfinite(m.radiation_fraction) && std::isfinite(m.lambda_fraction) &&
      std::isfinite(m.shear_x_over_anchor_hubble) && std::isfinite(m.shear_y_over_anchor_hubble);
}
} // namespace
BianchiITrajectory evolve_bianchi_i(const BianchiIModel &m,
    std::span<const BianchiIQuery> queries,BianchiIPolicy p) {
  BianchiITrajectory out; out.initial_state=m;
  if (!finite_model(m) || !std::isfinite(p.absolute_beta_tolerance) || !std::isfinite(p.relative_tolerance)) {
    out.status=S::nonfinite_input; return out;
  }
  if (std::fegetround()!=FE_TONEAREST || !(p.absolute_beta_tolerance>=128*eps) ||
      p.relative_tolerance<0 || (p.relative_tolerance>0 && p.relative_tolerance<128*eps) ||
      p.maximum_depth>60) return out;
  if (!(m.anchor_hubble_per_second>0) || m.matter_fraction<0 ||
      m.radiation_fraction<0 || m.lambda_fraction<0) {
    out.status=S::outside_domain; return out;
  }
  out.anchor_shear_over_hubble={m.shear_x_over_anchor_hubble,m.shear_y_over_anchor_hubble,
      -(W(m.shear_x_over_anchor_hubble)+m.shear_y_over_anchor_hubble)};
  W sigma=0,smax=0;
  for (W s:out.anchor_shear_over_hubble) { sigma+=s*s/6;smax=std::max(smax,std::abs(s)); }
  out.anchor_shear_fraction=sigma;
  if (!std::isfinite(sigma) || (smax>0 && !std::isnormal(sigma)) ||
      std::abs(W(m.matter_fraction)+m.radiation_fraction+m.lambda_fraction+sigma-1)>8*eps) {
    out.status=S::outside_domain; return out;
  }
  // Retained vector payload and fixed operation/recursive integration allowance.
  // Excludes borrowed queries, caller-held copies, allocator metadata and RSS.
  constexpr std::size_t overhead=sizeof(BianchiITrajectory)+32768;
  if (queries.size()>4096 || queries.size()>p.maximum_rows || p.maximum_native_bytes<overhead ||
      queries.size()>(p.maximum_native_bytes-overhead)/sizeof(BianchiIRow)) {
    out.status=S::work_limit; return out;
  }
  double previous=1;
  for (const auto &q:queries) {
    if (!std::isfinite(q.scale_factor)) { out.status=S::nonfinite_input;return out; }
    if (q.scale_factor<1e-4 || q.scale_factor>previous) { out.status=S::outside_domain;return out; }
    previous=q.scale_factor;
    W norm=0;
    for(double n:q.observed_direction) {
      if (!std::isfinite(n)) { out.status=S::nonfinite_input;return out; }
      norm+=W(n)*n;
    }
    if (std::abs(norm-1)>8*eps) { out.status=S::outside_domain;return out; }
  }
  out.rows.reserve(queries.size());
  Context context{m,sigma,smax};
  const bool pure_shear=m.matter_fraction==0 && m.radiation_fraction==0 && m.lambda_fraction==0;
  W scaled_integral=0,error=0;
  double prior_log=0;
  const double full_span=queries.empty()?0:-std::log(queries.back().scale_factor);
  for (const auto &q:queries) {
    const double x=std::log(q.scale_factor);
    if (smax>0 && x<prior_log) {
      if (pure_shear) {
        scaled_integral=-smax*W(x)/std::sqrt(sigma);
        error=32*eps*std::abs(scaled_integral);
      } else {
        if (p.maximum_callbacks-out.callbacks<3) { out.status=S::work_limit;out.rows.clear();return out; }
        numerics::IntegrationPolicy policy{p.absolute_beta_tolerance*(prior_log-x)/std::max(1.0,full_span),
          p.relative_tolerance,p.maximum_callbacks-out.callbacks,p.maximum_depth};
        const auto part=numerics::integrate(integrand,&context,x,prior_log,policy);
        out.callbacks+=part.evaluations;
        if (context.failure!=S::ok || part.status!=S::ok) {
          out.status=context.failure!=S::ok?context.failure:part.status;
          out.rows.clear();return out;
        }
        if (!(part.value>0)) { out.status=S::outside_domain;out.rows.clear();return out; }
        scaled_integral+=part.value;
        error+=part.error_estimate+32*eps*std::abs(part.value);
      }
    }
    prior_log=x;
    const W a=q.scale_factor,a3=a*a*a,root=std::sqrt(context.scaled_density(a));
    BianchiIRow row;row.query=q;
    row.expansion_over_anchor_hubble=root/a3;
    row.mean_hubble_per_second=m.anchor_hubble_per_second*row.expansion_over_anchor_hubble;
    row.shear_fraction=sigma/(root*root);
    row.maximum_beta_estimate=error;
    W directional_energy2=0;
    bool valid=std::isnormal(row.mean_hubble_per_second) && row.mean_hubble_per_second>0;
    for(unsigned i=0;i<3;++i) {
      const W s=out.anchor_shear_over_hubble[i];
      row.beta[i]=smax==0?0:-(s/smax)*scaled_integral;
      row.beta_trace_residual+=row.beta[i];
      if (!std::isfinite(row.beta[i]) || std::abs(row.beta[i])>64) { valid=false;break; }
      row.directional_scale_factors[i]=a*std::exp(row.beta[i]);
      const W scaled_axis_h=root+s;
      row.axis_hubble_per_second[i]=m.anchor_hubble_per_second*scaled_axis_h/a3;
      if (!std::isnormal(row.directional_scale_factors[i]) ||
          (scaled_axis_h!=0 && !std::isnormal(row.axis_hubble_per_second[i]))) { valid=false;break; }
      const W n=q.observed_direction[i];
      directional_energy2+=n*n*std::exp(-2*row.beta[i]);
    }
    row.one_plus_directional_redshift=std::sqrt(directional_energy2)/a;
    row.directional_redshift=row.one_plus_directional_redshift-1;
    row.redshift_relative_estimate=std::expm1(error)+64*eps;
    if (!valid || !std::isnormal(row.one_plus_directional_redshift) ||
        !std::isfinite(row.redshift_relative_estimate)) {
      out.status=S::outside_domain;out.rows.clear();return out;
    }
    // The same absolute+relative allocation used for the scaled shear integral
    // must also accommodate arithmetic diagnostics of its retained cumulative sum.
    if (error>p.absolute_beta_tolerance+p.relative_tolerance*std::abs(scaled_integral)) {
      out.status=S::conditioning_budget_exceeded;out.rows.clear();return out;
    }
    out.rows.push_back(row);
  }
  out.status=S::ok;return out;
}
} // namespace irred::cosmology
