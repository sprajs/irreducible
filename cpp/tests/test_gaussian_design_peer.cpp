// Synthetic exact anchor frozen before owner outputs. Independent augmented KKT
// elimination [C X; X^T 0][w;a]=[r;0], not the production QR route.
// C tridiagonal(2,1), X=[1,i], r=[1,0,0]: w=[1/4,-1/2,1/4],
// a=[1,-1/2], adjusted=[0,-1/2,0], Q=1/4. No external copied code.
// Reference absolute allocation 1e-12; product coefficients/residuals 2e-12
// times max(1,reference), Q absolute1e-10. Synthetic numerical qualification
// only; relative profile score is not a normalized density or H0 inference.
#include "irred/gaussian_design.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cfenv>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>
namespace {
using W=long double;using namespace irred::statistics;
unsigned checks=0;
void check(bool b,const char* n){++checks;if(!b)throw std::runtime_error(n);}
void near(double a,W b,const char* n){check(std::abs(W(a)-b)<=2e-12L*std::max(1.L,std::abs(b)),n);}
Metadata metadata(size_t n){Metadata m;for(size_t i=0;i<n;++i)m.ordered_ids.push_back("synthetic-"+std::to_string(i));m.measure="synthetic residual measure";m.ordering_provenance="explicit synthetic order";return m;}
DesignMetadata design(){return {{"intercept","slope"},{"u","u/index"},{"intercept"},"u","synthetic-two-column","synthetic shared covariance"};}
struct Oracle{std::vector<W> a,adjusted,w;W q=0;};
Oracle kkt(std::span<const double> c,std::span<const double> x,std::span<const double> r,size_t p){
 size_t n=r.size(),d=n+p;std::vector<W> a(d*(d+1));
 for(size_t i=0;i<n;++i){for(size_t j=0;j<n;++j)a[i*(d+1)+j]=c[i*n+j];for(size_t j=0;j<p;++j){a[i*(d+1)+n+j]=x[i*p+j];a[(n+j)*(d+1)+i]=x[i*p+j];}a[i*(d+1)+d]=r[i];}
 for(size_t k=0;k<d;++k){size_t pivot=k;for(size_t j=k+1;j<d;++j)if(std::abs(a[j*(d+1)+k])>std::abs(a[pivot*(d+1)+k]))pivot=j;
 check(a[pivot*(d+1)+k]!=0,"independent KKT nonsingular");for(size_t j=k;j<=d;++j)std::swap(a[k*(d+1)+j],a[pivot*(d+1)+j]);
 const W v=a[k*(d+1)+k];for(size_t j=k;j<=d;++j)a[k*(d+1)+j]/=v;
 for(size_t i=0;i<d;++i)if(i!=k){W factor=a[i*(d+1)+k];for(size_t j=k;j<=d;++j)a[i*(d+1)+j]-=factor*a[k*(d+1)+j];}}
 Oracle o;o.a.resize(p);o.w.resize(n);o.adjusted.resize(n);
 for(size_t j=0;j<p;++j)o.a[j]=a[(n+j)*(d+1)+d];for(size_t i=0;i<n;++i){o.w[i]=a[i*(d+1)+d];o.adjusted[i]=r[i];for(size_t j=0;j<p;++j)o.adjusted[i]-=x[i*p+j]*o.a[j];o.q+=o.adjusted[i]*o.w[i];}return o;
}
void compare(std::span<const double> c,std::span<const double> x,std::span<const double> r){
 auto m=metadata(r.size());auto ref=kkt(c,x,r,2);
 auto g=prepare_gaussian(c,MatrixKind::covariance,m,c.size(),1e-10);
 auto retained=DesignProfile::prepare(std::move(g),x,m.ordered_ids,design());
 check(retained.status()==DensityStatus::finite,"design prepared");check(g.status()==DensityStatus::invalid_input,"Gaussian consumed");
 auto got=retained.evaluate(r,m.ordered_ids);check(got.status==DensityStatus::finite,"design evaluated");
 check(got.coefficients.size()==2&&got.adjusted_residuals.size()==r.size(),"result dimensions");
 for(size_t j=0;j<2;++j)near(got.coefficients[j],ref.a[j],"KKT coefficients");
 for(size_t i=0;i<r.size();++i){near(got.adjusted_residuals[i],ref.adjusted[i],"KKT adjusted residual");W actual=r[i];for(size_t j=0;j<2;++j)actual-=W(x[i*2+j])*got.coefficients[j];check(got.adjusted_residuals[i]==static_cast<double>(actual),"reported coefficients own adjustment");}
 check(std::abs(W(got.quadratic)-ref.q)<=1e-10L,"KKT objective");near(got.relative_log_score,-W(got.quadratic)/2,"relative score");
 auto independent_gaussian=prepare_gaussian(c,MatrixKind::covariance,m,c.size(),1e-10);
 auto scored=independent_gaussian.evaluate(got.adjusted_residuals,m.ordered_ids,1e-10);
 check(scored.density.status==DensityStatus::finite,"actual rounded residual solve");near(got.quadratic,scored.quadratic,"objective on rounded residuals");
 check(got.normalized_normal_equation_residual>=0&&got.normalized_normal_equation_residual<1e-10,"normal equation diagnostic");
 auto bad=m.ordered_ids;std::swap(bad[0],bad[1]);check(retained.evaluate(r,bad).status==DensityStatus::incompatible_metadata,"wrong row order");
}
// Frozen exact rational control, also solved with an independent 90-decimal
// direct SVD before native evaluation. delta=2^-14, t=(-3,-1,1,3),
// X=[1,1+delta*t], r=(2,-1,4,0). beta=(820.45,-819.2),
// e=(.6,-2.3,2.8,-1.1), q=14.7. Values in X are exactly binary64.
// KKT normal elimination would square the condition here, so this anchor
// deliberately compares the analytic rational result, not its rounded KKT.
void dyadic_stress(int scale0,int scale1,bool reverse){
 std::array<double,16> c{};for(size_t i=0;i<4;++i)c[i*4+i]=1;
 std::array<int,4> t{-3,-1,1,3};std::array<double,8> x{};
 const double s0=std::ldexp(1.,scale0),s1=std::ldexp(1.,scale1),delta=std::ldexp(1.,-14);
 for(size_t i=0;i<4;++i){x[i*2]=s0;x[i*2+1]=s1*(1+delta*t[i]);if(reverse)std::swap(x[i*2],x[i*2+1]);}
 std::array<double,4> r{2,-1,4,0};std::array<W,2> beta{820.45L/s0,-819.2L/s1};
 std::array<W,4> expected{.6L,-2.3L,2.8L,-1.1L};if(reverse)std::swap(beta[0],beta[1]);
 auto m=metadata(4);auto d=design();d.parameter_units={"u/scaled-column0","u/scaled-column1"};
 if(reverse){std::swap(d.ordered_parameter_ids[0],d.ordered_parameter_ids[1]);std::swap(d.parameter_units[0],d.parameter_units[1]);}
 d.design_identity="synthetic-dyadic-scaled-two-column";
 auto g=prepare_gaussian(c,MatrixKind::covariance,m,16,1e-10);
 auto retained=DesignProfile::prepare(std::move(g),x,m.ordered_ids,d);
 check(retained.status()==DensityStatus::finite,"dyadic stress admitted without squared condition");
 check(retained.rank()==DesignRank::full_within_conditioning_contract,"dyadic rank scoped conditioning");
 check(g.status()==DensityStatus::invalid_input,"stress consumes source only on success");
 check(std::string_view(retained.method_id())=="retained-whitened-pivoted-householder-qr/v1","QR method identity");
 check(std::string_view(retained.qr_arithmetic_id())=="longdouble-cpu/v1","QR arithmetic identity");
 const auto condition=retained.equilibrated_triangular_condition_inf();
 check(std::isfinite(condition)&&condition>=1&&condition<1e6,"triangular condition finite scoped estimate");
 auto moved=std::move(retained);check(retained.status()!=DensityStatus::finite,"moved profile invalidated");
 // The retained operator owns X and metadata, not the caller's mutable spans.
 const auto original=x;std::fill(x.begin(),x.end(),std::numeric_limits<double>::quiet_NaN());
 for(int repeat=0;repeat<2;++repeat){auto got=moved.evaluate(r,m.ordered_ids);
  check(got.status==DensityStatus::finite,"retained repeated stress evaluation");
  for(size_t j=0;j<2;++j)near(got.coefficients[j],beta[j],"exact rational stress coefficients");
  W q=0,numerator=0,denominator=0;
  for(size_t i=0;i<4;++i){near(got.adjusted_residuals[i],expected[i],"exact rational stress residual");
   W actual=r[i];for(size_t j=0;j<2;++j)actual-=W(original[i*2+j])*got.coefficients[j];
   check(got.adjusted_residuals[i]==static_cast<double>(actual),"stress reported beta owns binary64 residual");q+=W(got.adjusted_residuals[i])*got.adjusted_residuals[i];}
  for(size_t j=0;j<2;++j){W norm=0,dot=0,den=0;for(size_t i=0;i<4;++i)norm+=W(original[i*2+j])*original[i*2+j];norm=std::sqrt(norm);
   for(size_t i=0;i<4;++i){W a=W(original[i*2+j])/norm;dot+=a*got.adjusted_residuals[i];den+=std::abs(a*r[i]);}numerator=std::max(numerator,std::abs(dot));denominator=std::max(denominator,den);}
  check(std::abs(W(got.quadratic)-14.7L)<=1e-10L,"exact rational stress objective");
  check(std::abs(W(got.quadratic)-q)<=1e-12L,"stress objective on actual rounded residual");
  const W defect=denominator?numerator/denominator:numerator;
  check(std::abs(W(got.normalized_normal_equation_residual)-defect)<=1e-15L,"independent actual residual stationarity metric");
  check(defect<=1e-10L,"actual residual stationarity admitted");
  check(got.covariance_whitening_backward_residual>=0&&std::isfinite(got.covariance_whitening_rounding_estimate),"whitening diagnostics finite");
 }
}
}
int main(){try{
 std::array<double,9> c{2,1,0,1,2,1,0,1,2};std::array<double,6> x{1,0,1,1,1,2};std::array<double,3> r{1,0,0};
 auto exact=kkt(c,x,r,2);check(std::abs(exact.a[0]-1)<=1e-12L&&std::abs(exact.a[1]+.5L)<=1e-12L&&std::abs(exact.q-.25L)<=1e-12L,"frozen exact KKT anchor");compare(c,x,r);
 std::array<size_t,3> order{2,0,1};std::array<double,9> cp{};std::array<double,6> xp{};std::array<double,3> rp{};
 for(size_t i=0;i<3;++i){rp[i]=r[order[i]];for(size_t j=0;j<3;++j)cp[i*3+j]=c[order[i]*3+order[j]];for(size_t j=0;j<2;++j)xp[i*2+j]=x[order[i]*2+j];}compare(cp,xp,rp);
 for(size_t i=0;i<3;++i){xp[i*2]=x[i*2+1];xp[i*2+1]=x[i*2];}compare(c,xp,r);
 xp=x;for(size_t i=0;i<3;++i)xp[i*2+1]*=8;compare(c,xp,r);
 std::array<double,16> correlated{4,1,1,1,1,3,1,1,1,1,5,1,1,1,1,2};std::array<double,8> multi{1,-2,1,-1,1,1,1,3};std::array<double,4> rr{2,-3,1,4};compare(correlated,multi,rr);
 dyadic_stress(0,0,false);dyadic_stress(20,-20,false);dyadic_stress(-20,20,true);
 auto m=metadata(3);auto g=prepare_gaussian(c,MatrixKind::covariance,m,9,1e-10);
 std::array<double,6> zero{1,0,1,0,1,0};auto failed=DesignProfile::prepare(std::move(g),zero,m.ordered_ids,design());check(failed.status()!=DensityStatus::finite&&failed.rank()==DesignRank::deficient,"proven zero column rank");check(g.status()==DensityStatus::finite,"failed preparation preserves source");
 std::array<double,6> duplicate{1,1,1,1,1,1};failed=DesignProfile::prepare(std::move(g),duplicate,m.ordered_ids,design());check(failed.status()!=DensityStatus::finite&&failed.rank()==DesignRank::unresolved,"dependent columns unresolved not jittered");
 std::array<double,6> close{1,1,1,1+1e-10,1,1-1e-10};failed=DesignProfile::prepare(std::move(g),close,m.ordered_ids,design());check(failed.status()!=DensityStatus::finite&&failed.rank()==DesignRank::unresolved,"near-collinear refused");
 check(g.status()==DensityStatus::finite,"near-collinear failure preserves source");
 auto nan=x;nan[1]=std::numeric_limits<double>::quiet_NaN();failed=DesignProfile::prepare(std::move(g),nan,m.ordered_ids,design());check(failed.status()!=DensityStatus::finite,"nonfinite design");
 auto prior=g.proper_offset(std::array<double,3>{1,1,1},m.ordered_ids,0,1,"proper-shared",true,9,1e-10);failed=DesignProfile::prepare(std::move(prior),x,m.ordered_ids,design());check(failed.status()==DensityStatus::incompatible_metadata,"proper prior separate from profile");
 std::array<double,3> tiny{1e-170,0,0};auto op=DesignProfile::prepare(std::move(g),x,m.ordered_ids,design());check(op.status()==DensityStatus::finite,"underflow operator");auto low=op.evaluate(tiny,m.ordered_ids);check(low.status==DensityStatus::numerical_failure,"positive Q underflow not zero");
 // Long-double arithmetic contract applies even with legacy binary64 C factor.
 auto rounding_source=prepare_gaussian(c,MatrixKind::covariance,m,9,1e-10);
 auto saved_rounding=std::fegetround();check(std::fesetround(FE_DOWNWARD)==0,"set hostile rounding");
 auto rounding_evaluation=op.evaluate(r,m.ordered_ids);
 auto rounding_preparation=DesignProfile::prepare(std::move(rounding_source),x,m.ordered_ids,design());
 check(std::fesetround(saved_rounding)==0,"restore rounding");
 check(rounding_evaluation.status!=DensityStatus::finite&&rounding_evaluation.numerical_status==irred::numerics::Status::outside_domain,"evaluation rejects unsupported rounding");
 check(rounding_preparation.status()!=DensityStatus::finite&&rounding_preparation.numerical_status()==irred::numerics::Status::outside_domain&&rounding_source.status()==DensityStatus::finite,"preparation rounding rejects preserves source");
 std::cout<<"PASS "<<checks<<" Gaussian design peer controls\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
