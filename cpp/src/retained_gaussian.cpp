#include "retained_gaussian.hpp"
#include "box_enclosure.hpp"
#include <utility>
namespace irred::detail {
using namespace statistics;
double RetainedQrAccess::preparation_sensitivity(const DesignProfile&d)noexcept{return d.preparation_sensitivity_;}
std::pair<std::size_t,std::size_t>RetainedQrAccess::preparation_whitenings(const DesignProfile&d)noexcept {
  return {d.preparation_whitening_attempts_,d.preparation_whitenings_completed_};
}
RetainedQrView RetainedQrAccess::view(const DesignProfile &d) noexcept {
  return {d.gaussian_.metadata().ordered_ids.size(),d.scales_.size(),d.qr_,d.scales_,
      d.pivot_,d.triangular_condition_,d.transpose_triangular_condition_};
}
const numerics::Factorization &RetainedQrAccess::factor(const DesignProfile &d) noexcept {
  return d.gaussian_.factor_;
}
const numerics::Factorization &RetainedQrAccess::factor(const Gaussian &g) noexcept {return g.factor_;}
std::span<const double> RetainedQrAccess::design(const DesignProfile &d) noexcept {return d.x_;}
std::span<const long double> RetainedQrAccess::normalized_design(const DesignProfile &d) noexcept {
  return d.whitened_design_;
}
std::span<const long double> RetainedQrAccess::reflectors(const DesignProfile &d) noexcept {return d.tau_;}
QrConditions qr_conditions(std::span<const long double> r,std::size_t p) {
  QrConditions out;
  if(!p||p>SIZE_MAX/p||r.size()<p*p) return out;
  long double norm_r=0,norm_transpose_r=0,norm_inverse_transpose_r=0;
  std::vector<long double> inverse_rows(p,0),inverse_column(p);
  for(std::size_t i=0;i<p;++i) {
    long double row=0;
    for(std::size_t j=i;j<p;++j) row+=std::abs(r[i*p+j]);
    norm_r=std::max(norm_r,row);
  }
  for(std::size_t j=0;j<p;++j) {
    long double column_norm=0;
    for(std::size_t i=0;i<=j;++i) column_norm+=std::abs(r[i*p+j]);
    norm_transpose_r=std::max(norm_transpose_r,column_norm);
  }
  for(std::size_t j=0;j<p;++j) {
    std::fill(inverse_column.begin(),inverse_column.end(),0);
    for(std::size_t k=p;k>0;--k) {
      const auto i=k-1;long double x=i==j ? 1 : 0;
      for(std::size_t t=i+1;t<p;++t) x-=r[i*p+t]*inverse_column[t];
      inverse_column[i]=x/r[i*p+i];inverse_rows[i]+=std::abs(inverse_column[i]);
    }
    long double column_norm=0;
    for(const auto x:inverse_column) column_norm+=std::abs(x);
    norm_inverse_transpose_r=std::max(norm_inverse_transpose_r,column_norm);
  }
  out.triangular=std::max(1.L,norm_r* *std::max_element(inverse_rows.begin(),inverse_rows.end()));
  out.transpose=std::max(1.L,norm_transpose_r*norm_inverse_transpose_r);
  out.status=std::isfinite(out.triangular)&&std::isfinite(out.transpose) ? numerics::Status::ok : numerics::Status::overflow;
  return out;
}
void qr_transform(const RetainedQrView &v,std::span<const long double> tau,std::span<long double> w) {
  for(std::size_t k=0;k<v.p;++k) {
    long double dot=w[k];
    for(std::size_t i=k+1;i<v.n;++i) dot+=v.qr[i*v.p+k]*w[i];
    dot*=tau[k];w[k]-=dot;
    for(std::size_t i=k+1;i<v.n;++i) w[i]-=v.qr[i*v.p+k]*dot;
  }
}
long double qr_tail(std::span<const long double> w,std::size_t p) {
  long double q=0;for(std::size_t i=p;i<w.size();++i) q+=w[i]*w[i];return q;
}
long double qr_precision_logdet(const RetainedQrView &v) {
  long double value=0;
  for(std::size_t j=0;j<v.p;++j)
    value+=2*(std::log(v.scales[j])+std::log(std::abs(v.qr[j*v.p+j])));
  return value;
}
EstimatorVarianceResult qr_variance(const RetainedQrView &v,std::span<const double> w,double budget) {
  EstimatorVarianceResult out;
  auto fail=[](numerics::Status s) {EstimatorVarianceResult r;r.status=DensityStatus::numerical_failure;
    r.numerical_status=s;return r;};
  const auto p=v.p,n=v.n;
  std::vector<long double> u(p),z(p);
  for(std::size_t i=0;i<p;++i) {
    u[i]=static_cast<long double>(w[v.pivot[i]])/v.scales[v.pivot[i]];
    if(!std::isfinite(u[i])||std::fpclassify(u[i])==FP_SUBNORMAL||(w[v.pivot[i]]!=0&&u[i]==0))
      return fail(numerics::Status::outside_domain);
    long double x=u[i];
    for(std::size_t j=0;j<i;++j) x-=v.qr[j*p+i]*z[j];
    z[i]=x/v.qr[i*p+i];
    if(!std::isfinite(z[i])||std::fpclassify(z[i])==FP_SUBNORMAL) return fail(numerics::Status::outside_domain);
  }
  long double variance=0,norm_z=0,norm_u=0,norm_rt=0,defect=0;
  for(std::size_t i=0;i<p;++i) {
    variance+=z[i]*z[i];norm_z=std::max(norm_z,std::abs(z[i]));norm_u=std::max(norm_u,std::abs(u[i]));
    long double row_norm=0,residual=-u[i];
    for(std::size_t j=0;j<=i;++j) {row_norm+=std::abs(v.qr[j*p+i]);residual+=v.qr[j*p+i]*z[j];}
    norm_rt=std::max(norm_rt,row_norm);defect=std::max(defect,std::abs(residual));
  }
  const auto denominator=norm_rt*norm_z+norm_u;
  const auto backward=denominator==0 ? 0 : defect/denominator;
  if(!std::isfinite(variance)||(variance==0&&norm_u!=0)||!qr_output(variance,out.variance))
    return fail(numerics::Status::outside_domain);
  const auto rounding=variance==0 ? 0 : std::abs(variance-static_cast<long double>(out.variance))/variance;
  const auto eta=v.transpose_condition*(backward+(static_cast<long double>(n)+p)*std::numeric_limits<long double>::epsilon());
  const auto sensitivity=2*eta+eta*eta+p*std::numeric_limits<long double>::epsilon()+rounding;
  if(!std::isfinite(sensitivity)||sensitivity>budget) return fail(numerics::Status::conditioning_budget_exceeded);
  out.triangular_backward_residual=static_cast<double>(backward);
  out.estimated_forward_sensitivity=static_cast<double>(sensitivity);
  out.output_rounding_error_relative=static_cast<double>(rounding);
  out.status=DensityStatus::finite;out.numerical_status=numerics::Status::ok;return out;
}
void qr_complete(GaussianBoxResult &out,const RetainedQrView &v,DesignResult &&fit,
    long double q,double logc,double budget) {
  out.completion_step=BoxCompletionStep::qr_completion;
  out.minimum_quadratic=box_scalar(q);out.log_design_precision_determinant=box_scalar(qr_precision_logdet(v));
  out.reported_profile_quadratic=fit.quadratic;out.profile_stationarity=fit.normalized_normal_equation_residual;
  out.unboxed_mean=std::move(fit.coefficients);out.unboxed_variance.resize(v.p);
  std::vector<double> weights(v.p,0);
  out.completion_step=BoxCompletionStep::marginal_variances;
  for(std::size_t j=0;j<v.p;++j) {
    out.completion_parameter_index=j;weights[j]=1;
    const auto variance=qr_variance(v,weights,budget);weights[j]=0;
    if(variance.status!=DensityStatus::finite) throw BoxRefusal{variance.numerical_status,{}};
    out.unboxed_variance[j]=variance.variance;
    out.maximum_variance_sensitivity=std::max(out.maximum_variance_sensitivity,variance.estimated_forward_sensitivity);
  }
  out.completion_parameter_index.reset();out.completion_step=BoxCompletionStep::covariance_determinant;
  out.source_covariance_log_determinant=box_scalar(logc);
  out.completion_step=BoxCompletionStep::complete;out.gaussian_completion_available=true;
  out.stage=BoxStage::gaussian_completion;
}
} // namespace irred::detail
