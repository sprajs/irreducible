// Optional standalone high-precision reproduction of the permanent fixture.
// Derived from the independently reviewed f4ab4d99 caller and retained c6ed18f
// output; no production Bessel/Simpson calls occur in the reference.
// Requires explicit MPFR/GMP linking and run; excluded from default CI.
#include "../fixtures/continuous_cmb_class_pl_604.hpp"
#include <mpfr.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace irred::projection;
using S = irred::numerics::Status;
mpfr_prec_t precision = 192;
std::uint64_t operations = 0, cell_attempts = 0, native_nodes = 0;
unsigned native_calls=0;
void tick() {
  if (operations >= 600000000) throw std::runtime_error("reference MPFR call cap");
  ++operations;
}
#define MP(name, ...) (tick(), mpfr_##name(__VA_ARGS__))
void need(bool ok, const char *why) { if (!ok) throw std::runtime_error(why); }
struct R {
  mpfr_t v;
  R() { mpfr_init2(v, precision); mpfr_set_zero(v, 0); }
  explicit R(long x) : R() { MP(set_si, v, x, MPFR_RNDN); }
  explicit R(double x) : R() { MP(set_d, v, x, MPFR_RNDN); }
  R(const R &o) { mpfr_init2(v, mpfr_get_prec(o.v)); MP(set, v, o.v, MPFR_RNDN); }
  R(R &&o) noexcept { mpfr_init2(v, mpfr_get_prec(o.v)); mpfr_swap(v, o.v); }
  R &operator=(const R &o) { if (this != &o) { mpfr_set_prec(v, mpfr_get_prec(o.v)); MP(set, v, o.v, MPFR_RNDN); } return *this; }
  R &operator=(R &&o) noexcept { if (this != &o) mpfr_swap(v, o.v); return *this; }
  ~R() { mpfr_clear(v); }
};
void sane(const R &x) {
  need(mpfr_number_p(x.v), "nonfinite reference arithmetic");
  if (!mpfr_zero_p(x.v)) need(std::abs(mpfr_get_exp(x.v)) < 50000, "reference exponent fence");
}
R abs_r(const R &a) { R r; MP(abs, r.v, a.v, MPFR_RNDU); return r; }
R up_add(const R &a, const R &b) { R r; MP(add, r.v, a.v, b.v, MPFR_RNDU); sane(r); return r; }
R up_mul(const R &a, const R &b) { R r; MP(mul, r.v, a.v, b.v, MPFR_RNDU); sane(r); need(mpfr_zero_p(a.v) || mpfr_zero_p(b.v) || !mpfr_zero_p(r.v), "reference positive product underflow"); return r; }
R up_div(const R &a, const R &b) { R r; need(mpfr_sgn(b.v) > 0, "reference bound denominator"); MP(div, r.v, a.v, b.v, MPFR_RNDU); sane(r); return r; }
R ulp(const R &a) {
  R r;
  if (!mpfr_zero_p(a.v)) MP(set_ui_2exp, r.v, 1, mpfr_get_exp(a.v)-mpfr_get_prec(a.v), MPFR_RNDU);
  return r;
}
struct B {
  R c, r;
  B() = default;
  explicit B(long x) : c(x) {}
  explicit B(double x) : c(x) {}
  explicit B(const R &x) : c(x) {} // exact emitted dyadic, not a root certificate
};
B neg(const B &a) { B o; MP(neg, o.c.v, a.c.v, MPFR_RNDN); o.r=a.r; return o; }
B operator+(const B &a, const B &b) {
  B o; MP(add, o.c.v, a.c.v, b.c.v, MPFR_RNDN); sane(o.c);
  o.r=up_add(up_add(a.r,b.r),ulp(o.c)); return o;
}
B operator-(const B &a, const B &b) { return a+neg(b); }
B operator*(const B &a, const B &b) {
  B o; MP(mul, o.c.v, a.c.v, b.c.v, MPFR_RNDN); sane(o.c);
  need(mpfr_zero_p(a.c.v)||mpfr_zero_p(b.c.v)||!mpfr_zero_p(o.c.v), "reference central product underflow");
  o.r=up_add(up_add(up_mul(abs_r(a.c),b.r),up_mul(abs_r(b.c),a.r)),up_add(up_mul(a.r,b.r),ulp(o.c))); return o;
}
B operator/(const B &a, const B &b) {
  B o; R den, ab=abs_r(b.c);
  MP(sub,den.v,ab.v,b.r.v,MPFR_RNDD);sane(den);need(mpfr_sgn(den.v)>0,"reference interval divides across zero");
  MP(div,o.c.v,a.c.v,b.c.v,MPFR_RNDN); sane(o.c);
  need(mpfr_zero_p(a.c.v)||!mpfr_zero_p(o.c.v),"reference central quotient underflow");
  const R u=ulp(o.c);
  o.r=up_add(up_div(up_add(a.r,up_mul(up_add(abs_r(o.c),u),b.r)),den),u); return o;
}
R upper_abs(const B &a) { return up_add(abs_r(a.c),a.r); }
R upper(const B &a) { return up_add(a.c,a.r); }
R lower(const B &a) { R r; MP(sub,r.v,a.c.v,a.r.v,MPFR_RNDD);sane(r);return r; }
B sin_b(const B &a) {
  B o; MP(sin,o.c.v,a.c.v,MPFR_RNDN); sane(o.c);
  need(mpfr_zero_p(a.c.v)||!mpfr_zero_p(o.c.v),"reference sine rounded to zero");
  o.r=up_add(a.r,ulp(o.c)); return o;
}
B cos_b(const B &a) { B o; MP(cos,o.c.v,a.c.v,MPFR_RNDN); sane(o.c); need(!mpfr_zero_p(o.c.v),"reference cosine rounded to zero"); o.r=up_add(a.r,ulp(o.c)); return o; }
B sqrt_b(const B &a) {
  B o; R lower, root;
  MP(sub,lower.v,a.c.v,a.r.v,MPFR_RNDD);sane(lower);need(mpfr_sgn(lower.v)>0,"reference sqrt domain");
  MP(sqrt,root.v,lower.v,MPFR_RNDD); MP(sqrt,o.c.v,a.c.v,MPFR_RNDN);sane(root);sane(o.c);
  need(mpfr_sgn(root.v)>0&&mpfr_sgn(o.c.v)>0,"reference sqrt rounded to zero");
  o.r=up_add(up_div(a.r,root),ulp(o.c)); return o;
}
bool contains_zero(const B &a) { R v=abs_r(a.c); return mpfr_cmp(v.v,a.r.v)<=0; }
std::string str(const R &x) { char b[256]; int n=mpfr_snprintf(b,sizeof b,"%Ra",x.v); need(n>0&&std::size_t(n)<sizeof b,"reference encoding fence"); return b; }
std::array<B,3> legendre(unsigned n,const B &x) {
  B p(1L),d(0L),dd(0L),oldp(0L),oldd(0L),olddd(0L);
  for(unsigned j=1;j<=n;++j) {
    const B a{long(2*j-1)}, b{long(j-1)}, den{long(j)};
    B np=(a*x*p-b*oldp)/den, nd=(a*(p+x*d)-b*oldd)/den;
    B ndd=(a*(B(2L)*d+x*dd)-b*olddd)/den;
    oldp=p;oldd=d;olddd=dd;p=np;d=nd;dd=ndd;
  }
  return {p,d,dd};
}
// Root generation is point arithmetic only. Its emitted dyadics are certified
// as a finite rule by moments, not assumed to be exact Gauss roots/weights.
std::pair<R,R> point_legendre(unsigned n,const R &x) {
  R p(1L),old(0L),next,tmp,a,b,den;
  for(unsigned j=1;j<=n;++j) {
    MP(mul,tmp.v,x.v,p.v,MPFR_RNDN); MP(mul_ui,a.v,tmp.v,2*j-1,MPFR_RNDN);
    MP(mul_ui,b.v,old.v,j-1,MPFR_RNDN); MP(sub,next.v,a.v,b.v,MPFR_RNDN);
    MP(div_ui,next.v,next.v,j,MPFR_RNDN); old=p;p=next;
  }
  MP(mul,a.v,x.v,p.v,MPFR_RNDN); MP(sub,a.v,a.v,old.v,MPFR_RNDN);
  MP(mul_ui,a.v,a.v,n,MPFR_RNDN); MP(mul,den.v,x.v,x.v,MPFR_RNDN);
  MP(sub_ui,den.v,den.v,1,MPFR_RNDN);need(mpfr_sgn(den.v)<0,"reference Legendre denominator");
  MP(div,b.v,a.v,den.v,MPFR_RNDN);sane(p);sane(b);
  return {p,b};
}
struct Rule { unsigned n; std::vector<std::pair<R,R>> points; R defect,weight_sum; };
Rule rule(unsigned n) {
  need(n>=32&&n<=256&&n%2==0,"reference rule order");
  Rule out{n,{},R(),R()}; out.points.reserve(n);
  R pi,seed,x,delta,tmp,a,b,w,tol;
  MP(const_pi,pi.v,MPFR_RNDN); MP(set_ui_2exp,tol.v,1,-precision+8,MPFR_RNDN);
  for(unsigned i=0;i<n/2;++i) {
    MP(mul_ui,seed.v,pi.v,4*i+3,MPFR_RNDN); MP(div_ui,seed.v,seed.v,4*n+2,MPFR_RNDN);
    MP(cos,x.v,seed.v,MPFR_RNDN); bool converged=false;
    for(unsigned it=0;it<12;++it) {
      const auto pd=point_legendre(n,x);need(!mpfr_zero_p(pd.second.v),"reference Newton derivative");
      MP(div,delta.v,pd.first.v,pd.second.v,MPFR_RNDN);sane(delta);
      MP(sub,x.v,x.v,delta.v,MPFR_RNDN); MP(abs,tmp.v,delta.v,MPFR_RNDN);
      sane(x);
      if(mpfr_cmp(tmp.v,tol.v)<=0) {converged=true;break;}
    }
    need(converged,"reference Newton cap"); const auto pd=point_legendre(n,x);
    MP(mul,a.v,x.v,x.v,MPFR_RNDN); MP(ui_sub,a.v,1,a.v,MPFR_RNDN);
    MP(mul,b.v,pd.second.v,pd.second.v,MPFR_RNDN); MP(mul,a.v,a.v,b.v,MPFR_RNDN);
    MP(ui_div,w.v,2,a.v,MPFR_RNDN); sane(x);sane(w);
    need(mpfr_sgn(x.v)>0&&mpfr_cmp_ui(x.v,1)<0&&mpfr_sgn(w.v)>0,"reference node/weight domain");
    out.points.emplace_back(x,w); R nx;MP(neg,nx.v,x.v,MPFR_RNDN);out.points.emplace_back(nx,w);
  }
  std::sort(out.points.begin(),out.points.end(),[](const auto &a,const auto &b){return mpfr_cmp(a.first.v,b.first.v)<0;});
  for(std::size_t j=1;j<out.points.size();++j)need(mpfr_cmp(out.points[j-1].first.v,out.points[j].first.v)<0,"reference duplicate node");
  std::vector<B> moments(2*n);
  for(const auto &point:out.points) {
    out.weight_sum=up_add(out.weight_sum,point.second); B power(1L);
    for(unsigned j=0;j<2*n;++j) {moments[j]=moments[j]+B(point.second)*power;power=power*B(point.first);}
  }
  for(unsigned j=0;j<2*n;++j) {
    const B exact=j%2?B(0L):B(2L)/B(long(j+1)); const R e=upper_abs(moments[j]-exact);
    if(mpfr_cmp(e.v,out.defect.v)>0)out.defect=e;
  }
  return out;
}
struct C { B re,im; };
using Values=std::array<C,4>;
Values reference(const ContinuousCmbSource &s,unsigned ell,const Rule &q) {
  need(s.k_mpc_inverse.size()==1&&s.eta_mpc.size()>=2&&ell>=2&&ell<=19,"reference fixed shape");
  Values total;
  const B k(s.k_mpc_inverse[0]),obs(s.observer_eta_mpc);
  const std::array<const std::vector<double>*,4> channels{&s.t0,&s.t1,&s.t2,&s.polarization};
  for(const auto &point:q.points) {
    const B mu(point.first), weight(point.second); const auto p=legendre(ell,mu);
    const B p2=(B(3L)*mu*mu-B(1L))/B(2L), one_minus=B(1L)-mu*mu;
    Values f;
    for(std::size_t j=0;j+1<s.eta_mpc.size();++j) {
      need(cell_attempts<1000000,"reference source-cell attempt cap");++cell_attempts;
      const B left(s.eta_mpc[j]),right(s.eta_mpc[j+1]);
      const B h=(right-left)/B(2L), mid=(right+left)/B(2L);
      const B v=k*mu*h, y=k*mu*(obs-mid);
      const B sv=sin_b(v), cv=cos_b(v), sy=sin_b(y), cy=cos_b(y);
      const B sinc=sv/v, derivative=(v*cv-sv)/(v*v);
      for(unsigned c=0;c<4;++c) {
        const B a((*channels[c])[j]),b((*channels[c])[j+1]);
        const B i0=h*(a+b)*sinc, i1=h*(b-a)*derivative;
        f[c].re=f[c].re+cy*i0-sy*i1;f[c].im=f[c].im+sy*i0+cy*i1;
      }
    }
    const std::array<B,4> factor{weight*p[0],weight*mu*p[0],neg(weight*p2*p[0]),weight*one_minus*one_minus*p[2]};
    for(unsigned c=0;c<4;++c) {
      const B re=c==1?neg(f[c].im):f[c].re, im=c==1?f[c].re:f[c].im;
      total[c].re=total[c].re+factor[c]*re;total[c].im=total[c].im+factor[c]*im;
    }
  }
  const B product(long((ell-1)*ell*(ell+1)*(ell+2)));
  const B spin=sqrt_b(B(3L)*product/B(8L));
  for(unsigned c=0;c<4;++c) {
    C v=total[c];
    switch(ell%4) {
      case 0:break;
      case 1:std::swap(v.re,v.im);v.im=neg(v.im);break;
      case 2:v.re=neg(v.re);v.im=neg(v.im);break;
      case 3:std::swap(v.re,v.im);v.re=neg(v.re);break;
    }
    B scale(0.5);if(c==3)scale=neg(scale*spin/product);
    total[c]={v.re*scale,v.im*scale};
  }
  // Explicit actual-rule polynomial moment defect + real-phase Taylor tail.
  B l0(1L),l1(1L);
  for(unsigned j=2;j<=ell;++j){const B next=(B(long(2*j-1))*l1+B(long(j-1))*l0)/B(long(j));l0=l1;l1=next;}
  const B L=ell?l1:l0;
  const std::array<B,4> coefficient{L/B(2L),L/B(2L),L,B(long(2*ell*(ell-1)))*L*spin/product};
  const R X=upper(k*(obs-B(s.eta_mpc.front()))); R exponential;MP(exp,exponential.v,X.v,MPFR_RNDU);sane(exponential);
  const unsigned d=ell+2;need(2*q.n>d,"reference Taylor shape");const unsigned power=2*q.n-d;
  R tail(1L);for(unsigned j=1;j<=power;++j){tail=up_mul(tail,X);tail=up_div(tail,R(long(j)));}
  const R kernel_error=up_add(up_mul(q.defect,exponential),up_mul(up_add(R(2L),q.weight_sum),tail));
  for(unsigned c=0;c<4;++c) {
    B mass;
    for(std::size_t j=0;j+1<s.eta_mpc.size();++j)
      mass=mass+(B(s.eta_mpc[j+1])-B(s.eta_mpc[j]))*B(up_add(abs_r(R((*channels[c])[j])),abs_r(R((*channels[c])[j+1]))))/B(2L);
    const R e=up_mul(up_mul(upper(coefficient[c]),upper(mass)),kernel_error);
    total[c].re.r=up_add(total[c].re.r,e);total[c].im.r=up_add(total[c].im.r,e);
    std::cout<<"ref_detail ell="<<ell<<" channel="<<c<<" N="<<q.n<<" bits="<<precision
             <<" center="<<str(total[c].re.c)<<" radius="<<str(total[c].re.r)
             <<" imaginary_center="<<str(total[c].im.c)<<" imaginary_radius="<<str(total[c].im.r)
             <<" moment_defect="<<str(q.defect)<<" taylor_tail="<<str(tail)<<" angular_bound="<<str(e)<<'\n';
    need(contains_zero(total[c].im),"reference angular imaginary enclosure excludes zero");
  }
  return total;
}
std::array<B,2> projected(const Values &v) {return {v[0].re+v[1].re+v[2].re,v[3].re};}
bool overlap(const B &a,const B &b) {
  const B d=a-b;
  const R distance=up_add(abs_r(d.c),ulp(d.c));
  R radius;MP(add,radius.v,a.r.v,b.r.v,MPFR_RNDD);sane(radius);
  return mpfr_cmp(distance.v,radius.v)<=0;
}
ContinuousCmbResult native_call(const ContinuousCmbProjection &owner,std::span<const unsigned> ell,unsigned mask,ContinuousCmbPolicy policy) {
  need(native_calls<4,"four scoped native calls");++native_calls;
  auto result=project_continuous_cmb(owner,ell,mask,policy);
  need(result.kernel_evaluations<=2003072-native_nodes,"whole native node cap");native_nodes+=result.kernel_evaluations;
  std::cout<<"native_call="<<native_calls<<" status="<<int(result.status)<<" nodes="<<result.kernel_evaluations<<" payload="<<result.payload_bound<<'\n';
  return result;
}
void controls() {
  precision=256;const Rule q=rule(32);ContinuousCmbPolicy policy;policy.maximum_kernel_evaluations=1024;
  for(double amplitude:{1.,-3.}) {
    ContinuousCmbSource s;s.k_mpc_inverse={.5};s.eta_mpc={1,3};s.observer_eta_mpc=12;
    s.t0={0,0};s.t1={amplitude,amplitude};s.t2={0,0};s.polarization={0,0};
    s.producer_id="synthetic-signed-constant-Doppler";s.signed_mode_id="declared-signed-amplitude";s.normalization_id="dimensionless";
    auto owner=prepare_continuous_cmb_projection(std::move(s));need(owner.status()==S::ok,"Doppler control acquisition");
    const auto ref=projected(reference(*owner.source(),2,q));const std::array<unsigned,1> ell{2};
    const auto result=native_call(owner,ell,continuous_temperature,policy);need(result.status==S::ok&&result.rows[0].temperature,"Doppler native control");
    auto j2=[](double x){const B b(x);return (B(3L)/(b*b*b)-B(1L)/b)*sin_b(b)-B(3L)/(b*b)*cos_b(b);};
    const B exact=B(amplitude)*(j2(5.5)-j2(4.5))/B(.5);
    need(overlap(ref[0],exact),"independent angular/trigonometric finite-endpoint control");
    const R difference=upper_abs(B(*result.rows[0].temperature)-exact);
    const B allowed=B(policy.absolute_tolerance)+B(policy.relative_tolerance)*B(abs_r(R(*result.rows[0].temperature)));
    need(mpfr_cmp(difference.v,lower(allowed).v)<=0,"signed native finite-endpoint gate");
    std::cout<<"control=signed-Doppler amplitude="<<amplitude<<" native="<<*result.rows[0].temperature<<" exact="<<str(exact.c)<<" radius="<<str(exact.r)<<'\n';
  }
  ContinuousCmbSource s;s.k_mpc_inverse={1};s.eta_mpc={0,1e-12};s.observer_eta_mpc=1e-12;
  s.t0={0,0};s.t1={0,0};s.t2={1e12,1e12};s.polarization={1e12,1e12};
  s.producer_id="synthetic-tiny-observer-endpoint";s.signed_mode_id="positive-regular-limit";s.normalization_id="explicit-1e12-amplitude";
  auto owner=prepare_continuous_cmb_projection(std::move(s));need(owner.status()==S::ok,"boundary control acquisition");
  const auto ref=projected(reference(*owner.source(),2,q));const std::array<unsigned,1> ell{2};
  const auto result=native_call(owner,ell,continuous_temperature|continuous_e_mode,policy);
  need(result.status==S::ok&&result.rows[0].temperature&&result.rows[0].e_mode,"boundary native control");
  const B leading=B(1e-12)*B(1e12)/B(5L), limit_allowance=B(1e-23);
  for(unsigned c=0;c<2;++c) {
    const R difference=upper_abs(ref[c]-leading);need(mpfr_cmp(difference.v,lower(limit_allowance).v)<=0,"regular one-fifth observer endpoint");
    const double value=c?*result.rows[0].e_mode:*result.rows[0].temperature;
    const R error=upper_abs(B(value)-ref[c]);need(mpfr_cmp(error.v,lower(B(policy.absolute_tolerance)+B(policy.relative_tolerance)*B(std::abs(value))).v)<=0,"boundary geometry gate");
    std::cout<<"control=observer-endpoint output="<<c<<" native="<<value<<" reference="<<str(ref[c].c)<<" radius="<<str(ref[c].r)<<'\n';
  }
}
} // namespace
int main() {
  try {
    need(mpfr_get_emin()<=-100000&&mpfr_get_emax()>=100000,"MPFR exponent profile");
    ContinuousCmbPreparationPolicy pp;pp.maximum_eta=20000; // historical explicit ceiling, default16384
    auto owner=prepare_continuous_cmb_projection(irred::test_fixtures::continuous_cmb_class_pl::make_source(),pp);need(owner.status()==S::ok,"original native source acquisition");
    const auto *identity=owner.source();std::cout<<std::setprecision(21);
    std::array<std::array<std::array<B,2>,2>,3> refs;
    for(unsigned attempt=0;attempt<3;++attempt) {
      precision=attempt==2?256:192;const Rule q=rule(attempt==0?128:256);
      for(unsigned li=0;li<2;++li) {
        refs[attempt][li]=projected(reference(*identity,li?19:2,q));
        for(unsigned c=0;c<2;++c)
          std::cout<<"projected_reference attempt="<<attempt<<" ell="<<(li?19:2)<<" output="<<(c?"E":"T")
                   <<" N="<<q.n<<" bits="<<precision<<" center="<<str(refs[attempt][li][c].c)<<" radius="<<str(refs[attempt][li][c].r)<<'\n';
      }
    }
    precision=256;controls();
    const ContinuousCmbPolicy policy;
    need(policy.absolute_tolerance==1e-9&&policy.relative_tolerance==1e-5&&policy.maximum_kernel_evaluations==2000000,"frozen native policy");
    const std::array<unsigned,2> ell{2,19};const auto native=native_call(owner,ell,continuous_temperature|continuous_e_mode,policy);
    need(native.source_owner.get()==identity&&native.status==S::ok&&native.rows.size()==2,"actual native/owner/status gate");
    bool agreement=true,earned=true;
    for(unsigned li=0;li<2;++li)for(unsigned c=0;c<2;++c) {
      const auto &row=native.rows[li];need(row.completed_source_cells==603&&row.temperature&&row.e_mode,"actual whole source support");
      const double value=c?*row.e_mode:*row.temperature;
      const double qe=c?row.e_quadrature_estimate:row.temperature_quadrature_estimate;
      const double ae=c?row.e_arithmetic_estimate:row.temperature_arithmetic_estimate;
      const B target=B(policy.absolute_tolerance)+B(policy.relative_tolerance)*B(std::abs(value));
      const B &r=refs[2][li][c];const B discrepancy=B(value)-r;
      const R difference=abs_r(discrepancy.c);
      R refinement=abs_r((refs[0][li][c]-refs[1][li][c]).c),paxis=abs_r((refs[1][li][c]-r).c);
      if(mpfr_cmp(paxis.v,refinement.v)>0)refinement=paxis;
      const R empirical=up_add(up_add(difference,refinement),upper(B(qe)+B(ae)));
      const R combined=up_add(upper_abs(discrepancy),upper(B(qe)+B(ae)));
      const bool agree=mpfr_cmp(empirical.v,lower(target).v)<=0;
      const bool complete=mpfr_cmp(combined.v,lower(target).v)<=0&&overlap(refs[0][li][c],refs[1][li][c])&&overlap(refs[1][li][c],r)&&overlap(refs[0][li][c],r);
      agreement&=agree;earned&=complete;
      std::cout<<"actual ell="<<ell[li]<<" output="<<(c?"E":"T")<<" native="<<value<<" reference="<<str(r.c)<<" reference_radius="<<str(r.r)
               <<" difference="<<str(difference)<<" refinement="<<str(refinement)<<" combined="<<str(combined)<<" allowed="<<str(lower(target))
               <<" native_time="<<qe<<" native_arithmetic="<<ae<<" native_radial=absent source_grid=absent agreement_control="<<agree<<" complete_reference_budget="<<complete<<'\n';
    }
    std::cout<<"agreement_control="<<agreement<<" complete_reference_budget="<<earned<<" MPFR_calls="<<operations<<" reference_cell_attempts="<<cell_attempts
             <<" native_nodes="<<native.kernel_evaluations<<" whole_native_nodes="<<native_nodes<<" native_payload="<<native.payload_bound<<" full_CLASS_CMB_qualification=false\n";
    return agreement&&earned?0:1;
  } catch(const std::exception &e) {
    std::cerr<<"REFUSED "<<e.what()<<" MPFR_calls="<<operations<<" reference_cell_attempts="<<cell_attempts<<" native_calls="<<native_calls<<" native_nodes="<<native_nodes<<'\n';return 1;
  }
}
