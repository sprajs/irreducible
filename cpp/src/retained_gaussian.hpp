#pragma once
#include "irred/gaussian_box.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace irred::detail {
// Private factor-free view. n is the original observation dimension; only the
// first p triangular rows are required for coefficients/variance/completion.
struct RetainedQrView {
  std::size_t n=0,p=0;
  std::span<const long double> qr,scales;
  std::span<const std::size_t> pivot;
  double triangular_condition=0,transpose_condition=0;
};
struct RetainedQrAccess {
  static double preparation_sensitivity(const statistics::DesignProfile &) noexcept;
  static std::pair<std::size_t,std::size_t> preparation_whitenings(const statistics::DesignProfile &) noexcept;
  static RetainedQrView view(const statistics::DesignProfile &) noexcept;
  static const numerics::Factorization &factor(const statistics::DesignProfile &) noexcept;
  static const numerics::Factorization &factor(const statistics::Gaussian &) noexcept;
  static std::span<const double> design(const statistics::DesignProfile &) noexcept;
  static std::span<const long double> normalized_design(const statistics::DesignProfile &) noexcept;
  static std::span<const long double> reflectors(const statistics::DesignProfile &) noexcept;
};
struct QrConditions {
  numerics::Status status=numerics::Status::invalid_input;
  long double triangular=0,transpose=0;
};
QrConditions qr_conditions(std::span<const long double>,std::size_t p);
void qr_transform(const RetainedQrView &,std::span<const long double> tau,
                  std::span<long double>);
long double qr_tail(std::span<const long double>,std::size_t p);
long double qr_precision_logdet(const RetainedQrView &);
statistics::EstimatorVarianceResult qr_variance(const RetainedQrView &,
    std::span<const double>,double maximum_forward_sensitivity);
void qr_complete(statistics::GaussianBoxResult &,const RetainedQrView &,
    statistics::DesignResult &&,long double minimum_quadratic,double covariance_logdet,
    double maximum_forward_sensitivity);

inline bool qr_output(long double v,double &d) noexcept {
  if(!std::isfinite(v)||std::abs(v)>std::numeric_limits<double>::max()) return false;
  d=static_cast<double>(v);
  return v==0 ? d==0 : std::fpclassify(d)==FP_NORMAL;
}
// One compiled equation/diagnostic owner used by ordinary and appended-row
// consumers. The callbacks are bounded compiled matrix/whitening operations,
// never input expressions or runtime physical models.
template<class Input,class RawDesign,class NormalizedDesign,class Whitening>
statistics::DesignResult qr_finish_fit(const RetainedQrView &v,
    std::span<const Input> residual,std::span<const long double> transformed_head,
    std::span<const long double> original_rhs_column_scales,
    RawDesign raw_design,NormalizedDesign normalized_design,Whitening whiten_adjusted,
    double maximum_forward_sensitivity) {
  using namespace statistics;
  auto fail=[](numerics::Status s) { DesignResult r;r.status=DensityStatus::numerical_failure;
    r.numerical_status=s;return r; };
  if(residual.size()!=v.n||transformed_head.size()<v.p||
      original_rhs_column_scales.size()!=v.p||v.scales.size()!=v.p||v.pivot.size()!=v.p)
    return fail(numerics::Status::invalid_input);
  DesignResult out;
  std::vector<long double> beta(v.p);
  for(std::size_t k=v.p;k>0;--k) {
    const auto i=k-1;long double x=transformed_head[i];
    for(std::size_t j=i+1;j<v.p;++j) x-=v.qr[i*v.p+j]*beta[j];
    beta[i]=x/v.qr[i*v.p+i];
  }
  out.coefficients.resize(v.p);out.adjusted_residuals.resize(v.n);
  for(std::size_t j=0;j<v.p;++j)
    if(!qr_output(beta[j]/v.scales[v.pivot[j]],out.coefficients[v.pivot[j]]))
      return fail(numerics::Status::outside_domain);
  for(std::size_t i=0;i<v.n;++i) {
    long double x=residual[i];
    for(std::size_t j=0;j<v.p;++j) x-=static_cast<long double>(raw_design(i,j))*out.coefficients[j];
    if(!qr_output(x,out.adjusted_residuals[i])) return fail(numerics::Status::outside_domain);
  }
  const auto we=whiten_adjusted(std::span<const double>(out.adjusted_residuals));
  if(we.status!=numerics::Status::ok) return fail(we.status);
  long double q=0,defect=0,denom=0;
  for(std::size_t i=0;i<v.n;++i) q+=we.value[i]*we.value[i];
  for(std::size_t j=0;j<v.p;++j) {
    long double d=0;
    for(std::size_t i=0;i<v.n;++i) d+=normalized_design(i,j)*we.value[i];
    defect=std::max(defect,std::abs(d));
    denom=std::max(denom,original_rhs_column_scales[j]);
  }
  const long double normal=denom==0 ? defect : defect/denom;
  if(!std::isfinite(normal)||normal>maximum_forward_sensitivity)
    return fail(numerics::Status::conditioning_budget_exceeded);
  if(q<0||!qr_output(q,out.quadratic)||!qr_output(-q/2,out.relative_log_score))
    return fail(numerics::Status::outside_domain);
  out.normalized_normal_equation_residual=static_cast<double>(normal);
  out.covariance_whitening_backward_residual=we.backward_residual;
  out.covariance_whitening_rounding_estimate=we.arithmetic_rounding_estimate;
  out.status=DensityStatus::finite;out.numerical_status=numerics::Status::ok;
  return out;
}
} // namespace irred::detail
