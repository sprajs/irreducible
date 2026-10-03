// Optional test-only independent finite-IVP reference. Production stays free
// of GMP/MPFR. No private finite-opacity equation/Schur/Radau header is used.
// SOURCE admission: e1e7c3d / independent review0264c4f8 (ignored evidence).
// Native actual theta/k core is compared to independently evolved F/theta.
// Arithmetic assumes the pinned MPFR directed-rounding/custom-storage contract;
// controls challenge that assumption, rather than proving a library/compiler.
#include "irred/finite_opacity_source.hpp"
#include "irred/quantities.hpp"
#include <mpfr.h>
#include <gmp.h>
#include <algorithm>
#include <array>
#include <charconv>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error "The conditional full reference requires strict floating arithmetic"
#endif

namespace reference {
constexpr unsigned L=192, D=2*L+8, degree=5;
constexpr unsigned B=2*(L+1), TB=B+1, C=B+2, TC=B+3, PH=B+4, SCALE=B+5;
constexpr unsigned G(unsigned l){return L+1+l;}
using W=long double;
using S=irred::numerics::Status;
struct Refusal:std::exception {
  const char *reason;
  explicit Refusal(const char *value):reason(value){}
  const char *what() const noexcept override{return reason;}
};
struct Ledger {
  std::size_t primitives=0,endpoints=0,guards=0,copies=0,status_calls=0;
  std::size_t steps=0,tube_checks=0,denied=0,live_intervals=0,peak_intervals=0;
  std::size_t byte_owners=0,peak_bytes=0,serializations=0;
  std::size_t encoded_copy_bytes=0,completed_steps=0;
  std::size_t maximum_step_primitives=0,maximum_step_endpoints=0;
  std::size_t nonstep_primitives=0,nonstep_endpoints=0;
  bool step_active=false;
  std::size_t primitive_limit=100000000,endpoint_limit=800000000;
  unsigned precision=192;
  void charge(std::size_t p,std::size_t e) {
    if(p>primitive_limit || primitives>primitive_limit-p ||
       e>endpoint_limit || endpoints>endpoint_limit-e) {
      ++denied;throw Refusal("reference work prefix limit");
    }
    if(!step_active && (p>1000000||nonstep_primitives>1000000-p)){
      ++denied;throw Refusal("one-time/setup/output/control prefix allowance");}
    primitives+=p;endpoints+=e;
    if(!step_active){nonstep_primitives+=p;nonstep_endpoints+=e;}
  }
  void guard(std::size_t e=0){charge(1,e);++guards;}
  void status_call(){charge(1,0);++status_calls;}
  void copy(std::size_t e=2){charge(1,e);++copies;}
  void acquire_interval() {
    charge(1,0);
    if(live_intervals==8192){++denied;throw Refusal("live interval limit");}
    ++live_intervals;peak_intervals=std::max(peak_intervals,live_intervals);
  }
  void bytes(std::size_t n) {
    guard();
    if(n>64*1024*1024 || byte_owners>64*1024*1024-n) {
      ++denied;throw Refusal("requested owned byte limit");
    }
    byte_owners+=n;peak_bytes=std::max(peak_bytes,byte_owners);
  }
  void step(){guard();if(steps==1536){++denied;throw Refusal("begun step limit");}++steps;}
};
Ledger *active=nullptr;
Ledger &work(){if(!active)throw Refusal("missing arithmetic context");return *active;}
// Requested byte owners are distinct from provider-internal MPFR/GMP scratch.
// Reserve local byte buffers before their initialization and hold the complete
// native owner's declared payload while any reference is alive.
struct Bytes {
  std::size_t amount;
  explicit Bytes(std::size_t n):amount(n){work().bytes(n);}
  ~Bytes(){work().byte_owners-=amount;}
  Bytes(const Bytes &)=delete;Bytes &operator=(const Bytes &)=delete;
};
struct StepReservation {
  std::size_t previous,start_primitives=0,start_endpoints=0;
  StepReservation():previous(work().primitive_limit) {
    // Reserve before beginning. A denied reservation computes no next step.
    work().guard();
    if(work().step_active||work().nonstep_primitives>980000||work().primitives>previous || previous-work().primitives<65536){
      ++work().denied;throw Refusal("step reservation/failure-publication headroom denied");}
    work().primitive_limit=work().primitives+65536;
    start_primitives=work().primitives;start_endpoints=work().endpoints;
    work().step_active=true;
    try{work().step();}catch(...){work().primitive_limit=previous;work().step_active=false;throw;}
  }
  ~StepReservation(){work().maximum_step_primitives=std::max(work().maximum_step_primitives,work().primitives-start_primitives);
    work().maximum_step_endpoints=std::max(work().maximum_step_endpoints,work().endpoints-start_endpoints);
    work().primitive_limit=previous;work().step_active=false;}
};
// Caller-owned mantissas: no mpfr_clear/set_prec/reallocation on these values.
// MPFR/GMP internal temporary allocations are not these requested owners;
// the separately required whole-process900s/2GiB supervisor covers them.
struct Interval {
  mpfr_t lo,hi;
  std::array<mp_limb_t,4> low_words,high_words;
  Interval() {
    auto &w=work();w.acquire_interval();
    try{w.bytes(sizeof(*this));}catch(...){--w.live_intervals;throw;}
    low_words.fill(0);high_words.fill(0);
    mpfr_custom_init_set(lo,MPFR_ZERO_KIND,0,w.precision,low_words.data());
    mpfr_custom_init_set(hi,MPFR_ZERO_KIND,0,w.precision,high_words.data());
  }
  Interval(const Interval &x):Interval(){set(x);}
  Interval(Interval &&x):Interval(){set(x);} // complete no-elide copy ownership
  Interval &operator=(const Interval &x){if(this!=&x)set(x);return *this;}
  Interval &operator=(Interval &&x){if(this!=&x)set(x);return *this;}
  ~Interval(){--work().live_intervals;work().byte_owners-=sizeof(*this);}
  void set(const Interval &x){work().copy();mpfr_set(lo,x.lo,MPFR_RNDD);mpfr_set(hi,x.hi,MPFR_RNDU);}
};
using Jet=std::array<Interval,degree+1>;
using State=std::array<Interval,D>;
using Table=std::array<Jet,D>;

struct Arithmetic {
  std::array<Interval,4> temporary;
  void check(const Interval &x) {
    work().guard(3);
    const bool finite_low=mpfr_number_p(x.lo),finite_high=mpfr_number_p(x.hi);
    const int order=mpfr_cmp(x.lo,x.hi);
    if(!finite_low||!finite_high||order>0)
      throw Refusal("invalid interval");
    const auto bad=MPFR_FLAGS_UNDERFLOW|MPFR_FLAGS_OVERFLOW|MPFR_FLAGS_NAN|
                   MPFR_FLAGS_ERANGE|MPFR_FLAGS_DIVBY0;
    work().status_call();if(mpfr_flags_test(bad))throw Refusal("MPFR exceptional arithmetic flag");
  }
  void integer(Interval &r,long value) {
    work().charge(1,2);mpfr_set_si(r.lo,value,MPFR_RNDD);mpfr_set_si(r.hi,value,MPFR_RNDU);
  }
  void stored(Interval &r,W value) {
    work().charge(1,2);
    if(!std::isfinite(value))throw Refusal("nonfinite stored scalar");
    const int a=mpfr_set_ld(r.lo,value,MPFR_RNDD),b=mpfr_set_ld(r.hi,value,MPFR_RNDU);
    if(a||b)throw Refusal("stored scalar import was not exact");
  }
  void add(Interval &r,const Interval &a,const Interval &b) {
    work().charge(1,2);
    mpfr_add(r.lo,a.lo,b.lo,MPFR_RNDD);mpfr_add(r.hi,a.hi,b.hi,MPFR_RNDU);
  }
  void sub(Interval &r,const Interval &a,const Interval &b) {
    work().charge(1,2);
    if(&r==&b){mpfr_sub(temporary[0].lo,a.lo,b.hi,MPFR_RNDD);
      mpfr_sub(r.hi,a.hi,b.lo,MPFR_RNDU);work().copy(1);mpfr_set(r.lo,temporary[0].lo,MPFR_RNDD);}
    else{mpfr_sub(r.lo,a.lo,b.hi,MPFR_RNDD);mpfr_sub(r.hi,a.hi,b.lo,MPFR_RNDU);}
  }
  void scale(Interval &r,const Interval &a,long n) {
    work().charge(1,2);
    if(n<0 && &r==&a){mpfr_mul_si(temporary[0].lo,a.hi,n,MPFR_RNDD);
      mpfr_mul_si(r.hi,a.lo,n,MPFR_RNDU);work().copy(1);mpfr_set(r.lo,temporary[0].lo,MPFR_RNDD);}
    else{mpfr_mul_si(r.lo,n>=0?a.lo:a.hi,n,MPFR_RNDD);
      mpfr_mul_si(r.hi,n>=0?a.hi:a.lo,n,MPFR_RNDU);}
  }
  void divide_integer(Interval &r,const Interval &a,unsigned n) {
    if(!n)throw Refusal("zero rational denominator");
    work().charge(1,2);mpfr_div_ui(r.lo,a.lo,n,MPFR_RNDD);mpfr_div_ui(r.hi,a.hi,n,MPFR_RNDU);
  }
  void multiply(Interval &r,const Interval &a,const Interval &b) {
    // Four downward and four upward candidates. Comparisons are separately
    // charged guards; all actual numeric MPFR calls enter endpoint work.
    work().charge(1,8);
    mpfr_srcptr av[2]{a.lo,a.hi},bv[2]{b.lo,b.hi};
    for(unsigned i=0;i<4;++i) {
      mpfr_mul(temporary[i].lo,av[i/2],bv[i%2],MPFR_RNDD);
      mpfr_mul(temporary[i].hi,av[i/2],bv[i%2],MPFR_RNDU);
    }
    unsigned low=0,high=0;work().guard(6);
    for(unsigned i=1;i<4;++i){if(mpfr_cmp(temporary[i].lo,temporary[low].lo)<0)low=i;
                             if(mpfr_cmp(temporary[i].hi,temporary[high].hi)>0)high=i;}
    work().copy();mpfr_set(r.lo,temporary[low].lo,MPFR_RNDD);mpfr_set(r.hi,temporary[high].hi,MPFR_RNDU);
  }
  void reciprocal(Interval &r,const Interval &a) {
    work().guard(2);
    const int low=mpfr_cmp_si(a.lo,0),high=mpfr_cmp_si(a.hi,0);
    if(low<=0 && high>=0)
      throw Refusal("interval denominator contains zero");
    work().charge(1,2);
    if(&r==&a){mpfr_ui_div(temporary[0].lo,1,a.hi,MPFR_RNDD);
      mpfr_ui_div(r.hi,1,a.lo,MPFR_RNDU);work().copy(1);mpfr_set(r.lo,temporary[0].lo,MPFR_RNDD);}
    else{mpfr_ui_div(r.lo,1,a.hi,MPFR_RNDD);mpfr_ui_div(r.hi,1,a.lo,MPFR_RNDU);}
  }
  void square_root(Interval &r,const Interval &a) {
    work().guard(1);if(mpfr_cmp_si(a.lo,0)<0)throw Refusal("negative sqrt support");
    work().charge(1,2);mpfr_sqrt(r.lo,a.lo,MPFR_RNDD);mpfr_sqrt(r.hi,a.hi,MPFR_RNDU);
  }
  void exponential(Interval &r,const Interval &a) {
    work().charge(1,2);mpfr_exp(r.lo,a.lo,MPFR_RNDD);mpfr_exp(r.hi,a.hi,MPFR_RNDU);
  }
  // Signed linear tail coefficients have their fixed sign checked at creation.
  // Select the two extremal products, then accumulate with correctly rounded
  // MPFR FMA. Each actual FMA is an endpoint call, each sign pair a guard.
  // r is disjoint from coefficient/value; its old endpoints are the addend.
  void linear_term(Interval &r,const Interval &coefficient,const Interval &value,bool negative) {
    work().guard(2);
    const bool low_negative=mpfr_cmp_si(value.lo,0)<0;
    const bool high_negative=mpfr_cmp_si(value.hi,0)<0;
    mpfr_srcptr cl=negative?(high_negative?coefficient.hi:coefficient.lo):
                            (low_negative?coefficient.hi:coefficient.lo);
    mpfr_srcptr ch=negative?(low_negative?coefficient.lo:coefficient.hi):
                            (high_negative?coefficient.lo:coefficient.hi);
    work().charge(1,2);
    mpfr_fma(r.lo,cl,negative?value.hi:value.lo,r.lo,MPFR_RNDD);
    mpfr_fma(r.hi,ch,negative?value.lo:value.hi,r.hi,MPFR_RNDU);
  }
  void positive_point_product(Interval &r,const Interval &value,const Interval &point) {
    // Caller checks the positive point once, before this fixed graph. Result
    // may alias value; the two endpoints are independent destinations.
    work().charge(1,2);mpfr_mul(r.lo,value.lo,point.lo,MPFR_RNDD);
    mpfr_mul(r.hi,value.hi,point.lo,MPFR_RNDU);
  }
  void nonnegative_product(Interval &r,const Interval &factor,const Interval &value) {
    work().guard(2);
    const bool l=mpfr_cmp_si(value.lo,0)<0,h=mpfr_cmp_si(value.hi,0)<0;
    work().charge(1,2);mpfr_mul(r.lo,l?factor.hi:factor.lo,value.lo,MPFR_RNDD);
    mpfr_mul(r.hi,h?factor.lo:factor.hi,value.hi,MPFR_RNDU);
  }
  void expanded_box(Interval &r,const Interval &value,const Interval &radius) {
    work().charge(1,2);mpfr_sub(r.lo,value.lo,radius.hi,MPFR_RNDD);
    mpfr_add(r.hi,value.hi,radius.hi,MPFR_RNDU);
  }
  void absolute(Interval &r,const Interval &a) {
    auto &t=temporary[0];work().guard(2);
    const int low=mpfr_cmp_si(a.lo,0),high=mpfr_cmp_si(a.hi,0);
    if(low>=0){r.set(a);return;}
    if(high<=0){scale(r,a,-1);return;}
    work().charge(1,2);mpfr_set_zero(t.lo,1);mpfr_neg(t.hi,a.lo,MPFR_RNDU);
    work().guard(1);if(mpfr_cmp(t.hi,a.hi)<0){work().copy(1);mpfr_set(t.hi,a.hi,MPFR_RNDU);}r.set(t);
  }
  bool contains(const Interval &outer,const Interval &inner,bool strict=false) {
    work().guard(2);const int l=mpfr_cmp(outer.lo,inner.lo),h=mpfr_cmp(outer.hi,inner.hi);
    return strict?(l<0&&h>0):(l<=0&&h>=0);
  }
  bool positive(const Interval &x){work().guard(1);return mpfr_cmp_si(x.lo,0)>0;}
  void hull(Interval &r,const Interval &x,const Interval &y) {
    work().guard(2);work().copy();
    mpfr_set(r.lo,mpfr_cmp(x.lo,y.lo)<=0?x.lo:y.lo,MPFR_RNDD);
    mpfr_set(r.hi,mpfr_cmp(x.hi,y.hi)>=0?x.hi:y.hi,MPFR_RNDU);
  }
  void upper(Interval &r,const Interval &x) {
    work().copy();mpfr_set(r.lo,x.hi,MPFR_RNDD);mpfr_set(r.hi,x.hi,MPFR_RNDU);
  }
  void maximum(Interval &r,const Interval &x,const Interval &y) {
    work().guard(2);work().copy();
    mpfr_set(r.lo,mpfr_cmp(x.lo,y.lo)>=0?x.lo:y.lo,MPFR_RNDD);
    mpfr_set(r.hi,mpfr_cmp(x.hi,y.hi)>=0?x.hi:y.hi,MPFR_RNDU);
  }
  void require_bound(const Interval &x,W low,W high) {
    stored(temporary[1],low);stored(temporary[2],high);work().guard(2);
    const int below=mpfr_cmp(x.lo,temporary[1].lo),above=mpfr_cmp(x.hi,temporary[2].hi);
    if(below<0||above>0)
      throw Refusal("declared finite source domain guard");
  }
};

struct Input {
  W photon=0,baryon=0,cdm=0,lambda=0,h0=0,initial_a=0;
  std::array<double,5> time{};
  std::array<double,2> opacity{},k{};
  std::array<std::array<double,11>,2> emitted_seed{};
};

// This finite-degree recurrence is written only for the selected equations.
// It is not a runtime AD system or an alternative production H owner.
struct Equations {
  Arithmetic arithmetic;
  std::array<Interval,20> t;
  Jet p,a2,a4,s,ia,ip,h,h2,fg,fb,fc,ratio,ieta,kjet,psi,pp,momentum;
  std::array<Interval,L-2> stream_previous,stream_next;
  std::array<Interval,2> negative_opacity;
  Interval photon,baryon,cdm,lambda,beta,k,k2,initial_a,eta_left,eta_right;
  Interval opacity_left,slope,one;
  explicit Equations(const Input &in,unsigned ki) {
    auto &a=arithmetic;
    a.stored(photon,in.photon);a.stored(baryon,in.baryon);a.stored(cdm,in.cdm);
    a.stored(lambda,in.lambda);a.stored(k,in.k[ki]);a.multiply(k2,k,k);
    a.stored(initial_a,in.initial_a);a.stored(eta_left,in.time.front());a.stored(eta_right,in.time.back());
    a.stored(opacity_left,in.opacity[0]);a.stored(t[0],in.opacity[1]);
    a.sub(t[0],t[0],opacity_left);a.sub(t[1],eta_right,eta_left);
    a.reciprocal(t[1],t[1]);a.multiply(slope,t[0],t[1]);
    a.integer(t[0],299792458);a.divide_integer(t[0],t[0],1000);
    a.reciprocal(t[0],t[0]);a.stored(t[1],in.h0);a.multiply(beta,t[1],t[0]);
    a.integer(one,1);
    if(!a.positive(k))throw Refusal("fixed positive k domain");
    a.require_bound(slope,0,1);
    for(unsigned l=3;l<L;++l){a.scale(stream_previous[l-3],k,l);a.divide_integer(stream_previous[l-3],stream_previous[l-3],2*l+1);
      a.scale(stream_next[l-3],k,-static_cast<long>(l+1));a.divide_integer(stream_next[l-3],stream_next[l-3],2*l+1);}
  }
  void product(Jet &out,const Jet &x,const Jet &y,unsigned n) {
    auto &a=arithmetic;
    for(unsigned j=0;j<=n;++j){a.integer(out[j],0);
      for(unsigned r=0;r<=j;++r){a.multiply(t[0],x[r],y[j-r]);a.add(out[j],out[j],t[0]);}}
  }
  void inverse(Jet &out,const Jet &x,unsigned n) {
    auto &a=arithmetic;a.reciprocal(out[0],x[0]);
    for(unsigned j=1;j<=n;++j){a.integer(t[2],0);
      for(unsigned r=1;r<=j;++r){a.multiply(t[0],x[r],out[j-r]);a.add(t[2],t[2],t[0]);}
      a.multiply(out[j],out[0],t[2]);a.scale(out[j],out[j],-1);}
  }
  void square_root_jet(Jet &out,const Jet &x,unsigned n) {
    auto &a=arithmetic;a.square_root(out[0],x[0]);a.scale(t[3],out[0],2);a.reciprocal(t[3],t[3]);
    for(unsigned j=1;j<=n;++j){a.integer(t[2],0);
      for(unsigned r=1;r<j;++r){a.multiply(t[0],out[r],out[j-r]);a.add(t[2],t[2],t[0]);}
      a.sub(t[2],x[j],t[2]);a.multiply(out[j],t[2],t[3]);}
  }
  void coefficients(Table &y,const Interval &eta,unsigned n) {
    auto &a=arithmetic;
    work().guard(3);const int positive=mpfr_cmp_si(eta.lo,0),left=mpfr_cmp(eta.lo,eta_left.lo),right=mpfr_cmp(eta.hi,eta_right.hi);
    if(positive<=0||left<0||right>0)
      throw Refusal("original positive opacity-cell time support");
    product(a2,y[SCALE],y[SCALE],n);product(a4,a2,a2,n);
    a.add(t[3],baryon,cdm);
    for(unsigned j=0;j<=n;++j){a.multiply(p[j],t[3],y[SCALE][j]);
      a.multiply(t[0],lambda,a4[j]);a.add(p[j],p[j],t[0]);
      if(!j)a.add(p[j],p[j],photon);}
    if(!a.positive(p[0])||!a.positive(y[SCALE][0]))throw Refusal("nonpositive clock support");
    square_root_jet(s,p,n);
    inverse(ia,y[SCALE],n);inverse(ip,p,n);product(h,s,ia,n);
    for(unsigned j=0;j<=n;++j){a.multiply(h[j],h[j],beta);a.multiply(fg[j],ip[j],photon);}
    product(fb,y[SCALE],ip,n);product(fc,y[SCALE],ip,n);
    for(unsigned j=0;j<=n;++j){a.multiply(fb[j],fb[j],baryon);a.multiply(fc[j],fc[j],cdm);}
    product(h2,h,h,n);
    a.reciprocal(t[3],baryon);a.multiply(t[3],t[3],photon);a.scale(t[3],t[3],4);a.divide_integer(t[3],t[3],3);
    for(unsigned j=0;j<=n;++j)a.multiply(ratio[j],ia[j],t[3]);
    // Absolute eta, not eta-eta_i, belongs to the original terminal closure.
    a.reciprocal(ieta[0],eta);
    for(unsigned j=1;j<=n;++j){a.multiply(ieta[j],ieta[j-1],ieta[0]);a.scale(ieta[j],ieta[j],-1);}
    a.sub(t[0],eta,eta_left);a.multiply(kjet[0],slope,t[0]);a.add(kjet[0],kjet[0],opacity_left);
    if(n)kjet[1].set(slope);
    for(unsigned j=2;j<=n;++j)a.integer(kjet[j],0);
    a.require_bound(kjet[0],0,1);
    a.scale(negative_opacity[0],kjet[0],-1);a.scale(negative_opacity[1],slope,-1);
  }
  // coefficient of a convolution, with separately owned destination.
  void convolution(Interval &out,const Jet &x,const Jet &y,unsigned j) {
    auto &a=arithmetic;a.integer(out,0);
    for(unsigned r=0;r<=j;++r){a.multiply(t[0],x[r],y[j-r]);a.add(out,out,t[0]);}
  }
  void affine_convolution(Interval &out,const Jet &x,const Jet &y,unsigned j) {
    auto &a=arithmetic;a.integer(out,0);
    for(unsigned r=0;r<=std::min(1u,j);++r){a.multiply(t[0],x[r],y[j-r]);a.add(out,out,t[0]);}
  }
  void metric(Table &y,unsigned n) {
    auto &a=arithmetic;
    // Direct metric construction in F2/theta coordinates.
    product(psi,h2,fg,n);product(pp,psi,y[2],n);
    a.reciprocal(t[3],k2);
    for(unsigned j=0;j<=n;++j){a.multiply(t[0],pp[j],t[3]);a.scale(t[0],t[0],3);a.sub(psi[j],y[PH][j],t[0]);}
    product(momentum,fc,y[TC],n);product(pp,fb,y[TB],n);
    for(unsigned j=0;j<=n;++j)a.add(momentum[j],momentum[j],pp[j]);
    product(pp,fg,y[1],n);
    for(unsigned j=0;j<=n;++j){a.multiply(t[0],pp[j],k);a.add(momentum[j],momentum[j],t[0]);}
    product(pp,h2,momentum,n);product(momentum,h,psi,n);
    for(unsigned j=0;j<=n;++j){a.multiply(pp[j],pp[j],t[3]);a.scale(pp[j],pp[j],3);a.divide_integer(pp[j],pp[j],2);
      a.sub(pp[j],pp[j],momentum[j]);}
  }
  void fill(Table &y,const State &start,const Interval &eta,unsigned order) {
    auto &a=arithmetic;
    for(unsigned i=0;i<D;++i)y[i][0].set(start[i]);
    for(unsigned j=0;j<order;++j) {
      coefficients(y,eta,j);metric(y,j);
      a.multiply(y[SCALE][j+1],beta,s[j]);a.divide_integer(y[SCALE][j+1],y[SCALE][j+1],j+1);
      // Temporaries t[4..] never alias t[0..3] used by convolution.
      a.multiply(t[4],k,y[1][j]);a.scale(t[4],t[4],-1);a.scale(t[5],pp[j],4);a.add(t[4],t[4],t[5]);
      a.divide_integer(y[0][j+1],t[4],j+1);
      a.scale(t[4],y[2][j],2);a.sub(t[4],y[0][j],t[4]);a.multiply(t[4],t[4],k);a.divide_integer(t[4],t[4],3);
      a.multiply(t[5],psi[j],k);a.scale(t[5],t[5],4);a.divide_integer(t[5],t[5],3);a.add(t[4],t[4],t[5]);
      for(unsigned r=0;r<=std::min(1u,j);++r){a.reciprocal(t[6],k);a.multiply(t[5],y[TB][j-r],t[6]);a.scale(t[5],t[5],4);a.divide_integer(t[5],t[5],3);
        a.sub(t[5],t[5],y[1][j-r]);a.multiply(t[5],t[5],kjet[r]);a.add(t[4],t[4],t[5]);}
      a.divide_integer(y[1][j+1],t[4],j+1);
      for(unsigned i:{2u,G(0),G(1),G(2)}) {
        if(i==2){a.scale(t[4],y[1][j],2);a.scale(t[5],y[3][j],3);a.sub(t[4],t[4],t[5]);a.divide_integer(t[4],t[4],5);}
        else if(i==G(0)){a.scale(t[4],y[G(1)][j],-1);}
        else {const unsigned l=i-G(0);a.scale(t[4],y[G(l-1)][j],l);a.scale(t[5],y[G(l+1)][j],l+1);
              a.sub(t[4],t[4],t[5]);a.divide_integer(t[4],t[4],2*l+1);}
        a.multiply(t[4],t[4],k);
        for(unsigned r=0;r<=std::min(1u,j);++r){
          a.add(t[5],y[2][j-r],y[G(0)][j-r]);a.add(t[5],t[5],y[G(2)][j-r]);
          if(i==G(1))a.integer(t[5],0);
          else a.divide_integer(t[5],t[5],i==G(0)?2:10);
          a.sub(t[5],y[i][j-r],t[5]);a.multiply(t[5],kjet[r],t[5]);a.sub(t[4],t[4],t[5]);}
        a.divide_integer(y[i][j+1],t[4],j+1);
      }
      for(unsigned base:{0u,G(0)})for(unsigned l=3;l<=L;++l) {
        const unsigned i=base+l;
        if(l<L){a.integer(t[4],0);a.linear_term(t[4],stream_previous[l-3],y[i-1][j],false);
          a.linear_term(t[4],stream_next[l-3],y[i+1][j],true);
          a.linear_term(t[4],negative_opacity[0],y[i][j],true);
          if(j)a.linear_term(t[4],negative_opacity[1],y[i][j-1],true);}
        else {a.multiply(t[4],y[i-1][j],k);affine_convolution(t[5],kjet,y[i],j);a.sub(t[4],t[4],t[5]);
          convolution(t[5],ieta,y[i],j);a.scale(t[5],t[5],L+1);a.sub(t[4],t[4],t[5]);}
        a.divide_integer(y[i][j+1],t[4],j+1);
      }
      for(unsigned density:{B,C}) {
        a.scale(t[4],y[density+1][j],-1);a.scale(t[5],pp[j],3);a.add(t[4],t[4],t[5]);
        a.divide_integer(y[density][j+1],t[4],j+1);
        convolution(t[4],h,y[density+1],j);a.scale(t[4],t[4],-1);a.multiply(t[5],k2,psi[j]);a.add(t[4],t[4],t[5]);
        if(density==B)for(unsigned r=0;r<=j;++r)for(unsigned v=0;v<=std::min(1u,j-r);++v){
          a.multiply(t[5],k,y[1][j-r-v]);a.scale(t[5],t[5],3);a.divide_integer(t[5],t[5],4);
          a.sub(t[5],t[5],y[TB][j-r-v]);a.multiply(t[5],t[5],ratio[r]);a.multiply(t[5],t[5],kjet[v]);a.add(t[4],t[4],t[5]);}
        a.divide_integer(y[density+1][j+1],t[4],j+1);
      }
      a.divide_integer(y[PH][j+1],pp[j],j+1);
    }
  }
};

struct Stepper {
  Equations e;
  Table start_jets,tube_jets;
  State tube,image,next,radius;
  std::array<Interval,24> t;
  Interval time_box,time_span,next_eta;
  explicit Stepper(const Input &in,unsigned ki):e(in,ki){}
  void scale_unit(Interval &out,unsigned i) {
    if(i==SCALE)out.set(e.initial_a);
    else if(i==TB||i==TC)out.set(e.k);
    else e.arithmetic.integer(out,1);
  }
  void initial(State &y,const Input &in,unsigned ki) {
    auto &a=e.arithmetic;
    for(auto &v:y)a.integer(v,0);
    const auto &c=in.emitted_seed[ki];
    a.stored(y[0],c[0]);a.stored(y[1],c[1]);a.scale(y[1],y[1],4);a.divide_integer(y[1],y[1],3);
    a.stored(y[2],c[2]);a.scale(y[2],y[2],2);
    for(unsigned l=0;l<3;++l)a.stored(y[G(l)],c[3+l]);
    a.stored(y[B],c[6]);a.stored(y[TB],c[7]);a.multiply(y[TB],y[TB],e.k);
    a.stored(y[C],c[8]);a.stored(y[TC],c[9]);a.multiply(y[TC],y[TC],e.k);
    a.stored(y[PH],c[10]);y[SCALE].set(e.initial_a);
  }
  // Independently derived row-sum bounds in F and theta/k,a/a_i coordinates.
  // Include the CLOCK column of the complete392-dimensional Jacobian.
  void jacobian_guard(const State &y,const Interval &eta,const Interval &hstep) {
    auto &a=e.arithmetic;
    a.require_bound(y[SCALE],.000999L,.001011L);
    a.require_bound(e.photon,5e-5L,6e-5L);a.require_bound(e.baryon,.048L,.051L);
    a.add(t[0],e.baryon,e.cdm);a.require_bound(t[0],.30L,.32L);
    a.require_bound(e.lambda,0,1);a.require_bound(e.h[0],.0041L,.0045L);
    // P_a=b+c+4lambda*a^3; H_a=H*(P_a/(2P)-1/a).
    a.multiply(t[1],e.a2[0],y[SCALE]);a.multiply(t[1],t[1],e.lambda);a.scale(t[1],t[1],4);a.add(t[1],t[1],t[0]);
    a.multiply(t[2],t[1],e.ip[0]);a.divide_integer(t[2],t[2],2);a.sub(t[2],t[2],e.ia[0]);a.multiply(t[2],t[2],e.h[0]);
    // derivatives of fg,fb,fc in t3,t4,t5 respectively
    a.multiply(t[3],e.ip[0],t[1]);a.multiply(t[3],t[3],e.fg[0]);a.scale(t[3],t[3],-1);
    for(unsigned i=0;i<2;++i){
      const auto &fraction=i?e.fc[0]:e.fb[0];const auto &coefficient=i?e.cdm:e.baryon;
      a.multiply(t[4+i],coefficient,e.ip[0]);a.multiply(t[6],fraction,e.ip[0]);a.multiply(t[6],t[6],t[1]);a.sub(t[4+i],t[4+i],t[6]);}
    a.multiply(t[6],e.h[0],t[2]);a.scale(t[6],t[6],2);a.multiply(t[6],t[6],e.fg[0]);
    a.multiply(t[7],e.h2[0],t[3]);a.add(t[6],t[6],t[7]);a.multiply(t[6],t[6],y[2]);
    a.reciprocal(t[7],e.k2);a.multiply(t[6],t[6],t[7]);a.scale(t[6],t[6],-3); // psi_a
    a.multiply(t[8],t[2],e.psi[0]);a.scale(t[8],t[8],-1);
    a.multiply(t[9],e.h[0],t[6]);a.sub(t[8],t[8],t[9]);
    a.multiply(t[9],t[4],y[TB]);a.multiply(t[10],t[5],y[TC]);a.add(t[9],t[9],t[10]);
    a.multiply(t[10],t[3],e.k);a.multiply(t[10],t[10],y[1]);a.add(t[9],t[9],t[10]);
    a.multiply(t[9],t[9],e.h2[0]);
    a.multiply(t[10],e.fb[0],y[TB]);a.multiply(t[11],e.fc[0],y[TC]);a.add(t[10],t[10],t[11]);
    a.multiply(t[11],e.fg[0],e.k);a.multiply(t[11],t[11],y[1]);a.add(t[10],t[10],t[11]);
    a.multiply(t[10],t[10],e.h[0]);a.multiply(t[10],t[10],t[2]);a.scale(t[10],t[10],2);
    a.add(t[9],t[9],t[10]);a.multiply(t[9],t[9],t[7]);a.scale(t[9],t[9],3);a.divide_integer(t[9],t[9],2);a.add(t[8],t[8],t[9]); // phi'_a
    a.absolute(t[8],t[8]);a.multiply(t[8],t[8],e.initial_a);
    a.absolute(t[6],t[6]);a.multiply(t[6],t[6],e.initial_a);
    // pp row bound, psi row bound, then explicit species/hierarchy maxima
    a.multiply(t[9],e.h2[0],e.fg[0]);a.multiply(t[9],t[9],t[7]);a.scale(t[9],t[9],3);a.add(t[9],t[9],e.one); // psi norm
    a.multiply(t[10],e.h[0],t[9]);a.add(t[11],e.fb[0],e.fc[0]);a.add(t[11],t[11],e.fg[0]);
    a.multiply(t[11],t[11],e.h2[0]);a.reciprocal(t[12],e.k);a.multiply(t[11],t[11],t[12]);a.scale(t[11],t[11],3);a.divide_integer(t[11],t[11],2);
    a.add(t[10],t[10],t[11]);a.add(t[13],t[10],t[8]); // phi full row
    a.scale(t[14],t[13],4);a.add(t[14],t[14],e.k);a.maximum(t[13],t[13],t[14]); // F0
    a.multiply(t[14],e.k,t[9]);a.multiply(t[15],e.k,t[6]);a.add(t[14],t[14],t[15]);a.scale(t[14],t[14],4);a.divide_integer(t[14],t[14],3);
    a.scale(t[15],e.kjet[0],7);a.divide_integer(t[15],t[15],3);a.add(t[14],t[14],t[15]);a.add(t[14],t[14],e.k);a.maximum(t[13],t[13],t[14]);
    a.scale(t[14],e.kjet[0],5);a.divide_integer(t[14],t[14],2);a.add(t[14],t[14],e.k);a.maximum(t[13],t[13],t[14]); // covers G0 and F2/G2
    for(unsigned velocity:{TB,TC}) {
      a.reciprocal(t[14],e.k);a.multiply(t[14],t[14],y[velocity]);a.absolute(t[14],t[14]);a.absolute(t[15],t[2]);a.multiply(t[14],t[14],t[15]);
      a.multiply(t[14],t[14],e.initial_a);a.multiply(t[15],t[9],e.k);a.add(t[14],t[14],t[15]);a.multiply(t[15],t[6],e.k);a.add(t[14],t[14],t[15]);a.add(t[14],t[14],e.h[0]);
      if(velocity==TB){a.scale(t[15],e.ratio[0],7);a.divide_integer(t[15],t[15],4);a.multiply(t[15],t[15],e.kjet[0]);a.add(t[14],t[14],t[15]);
        a.reciprocal(t[15],e.k);a.multiply(t[15],t[15],y[TB]);a.scale(t[16],y[1],3);a.divide_integer(t[16],t[16],4);a.sub(t[15],t[16],t[15]);a.absolute(t[15],t[15]);
        a.multiply(t[15],t[15],e.ratio[0]);a.multiply(t[15],t[15],e.ia[0]);a.multiply(t[15],t[15],e.kjet[0]);a.multiply(t[15],t[15],e.initial_a);a.add(t[14],t[14],t[15]);}
      a.maximum(t[13],t[13],t[14]);}
    a.add(t[14],t[10],t[8]);a.scale(t[14],t[14],3);a.add(t[14],t[14],e.k);a.maximum(t[13],t[13],t[14]);
    a.reciprocal(t[14],eta);a.scale(t[14],t[14],L+1);a.add(t[14],t[14],e.kjet[0]);a.add(t[14],t[14],e.k);a.maximum(t[13],t[13],t[14]);
    a.multiply(t[14],t[1],e.beta);a.reciprocal(t[15],e.s[0]);a.multiply(t[14],t[14],t[15]);a.divide_integer(t[14],t[14],2);a.maximum(t[13],t[13],t[14]);
    a.require_bound(t[13],0,2);a.multiply(t[13],t[13],hstep);a.require_bound(t[13],0,.125L);
  }
  void step(State &y,Interval &eta,const Interval &hstep) {
    auto &a=e.arithmetic;StepReservation reservation;
    work().guard(2);const int point=mpfr_cmp(hstep.lo,hstep.hi),positive=mpfr_cmp_si(hstep.lo,0);
    if(point||positive<=0)
      throw Refusal("step must be an exact positive dyadic point");
    e.fill(start_jets,y,eta,4);
    a.integer(t[0],0);a.hull(time_span,t[0],hstep);a.add(time_box,eta,time_span);
    a.integer(t[3],0);
    for(unsigned i=0;i<D;++i){scale_unit(t[4],i);a.reciprocal(t[4],t[4]);a.multiply(t[5],y[i],t[4]);a.absolute(t[5],t[5]);a.maximum(t[3],t[3],t[5]);}
    a.multiply(t[3],t[3],e.k2);a.positive_point_product(t[3],t[3],hstep);
    a.stored(t[4],1e-30L);a.add(t[3],t[3],t[4]);
    for(unsigned i=0;i<D;++i){a.absolute(radius[i],start_jets[i][1]);scale_unit(t[4],i);
      a.multiply(t[4],t[4],t[3]);a.add(radius[i],radius[i],t[4]);a.positive_point_product(radius[i],radius[i],hstep);a.scale(radius[i],radius[i],2);}
    bool included=false;
    // One checked tube in this first fixed implementation. The source proposal
    // allows at most four; its fourth-success graph exceeded65536 and is kept
    // as a source finding. No retry, width discard or resource-cap increase.
    for(unsigned attempt=0;attempt<1;++attempt){++work().tube_checks;
      for(unsigned i=0;i<D;++i)a.expanded_box(tube[i],y[i],radius[i]);
      e.fill(tube_jets,tube,time_box,1);jacobian_guard(tube,time_box,hstep);
      included=true;
      for(unsigned i=0;i<D;++i){a.nonnegative_product(image[i],time_span,tube_jets[i][1]);a.add(image[i],image[i],y[i]);
        if(!a.contains(tube[i],image[i],true))included=false;}
      if(included)break;
    }
    if(!included)throw Refusal("strict Picard tube inclusion failed");
    // This is J5 over EVERY point of the tube and time interval. It is not the
    // point-start J5 and does not discard previous endpoint wrapping widths.
    e.fill(tube_jets,tube,time_box,5);
    t[6].set(hstep);for(unsigned j=1;j<5;++j)a.positive_point_product(t[6],t[6],hstep);
    for(unsigned i=0;i<D;++i){next[i].set(start_jets[i][4]);
      for(unsigned j=4;j-->0;){a.positive_point_product(next[i],next[i],hstep);a.add(next[i],next[i],start_jets[i][j]);}
      a.nonnegative_product(t[7],t[6],tube_jets[i][5]);a.add(next[i],next[i],t[7]);a.check(next[i]);}
    a.add(next_eta,eta,hstep);work().guard();
    if(work().primitives>work().primitive_limit || work().primitive_limit-work().primitives<D+1 || work().endpoint_limit<2*(D+1) ||
       work().endpoints>work().endpoint_limit-2*(D+1))throw Refusal("atomic full-state/time commit reservation");
    for(unsigned i=0;i<D;++i)y[i].set(next[i]);
    eta.set(next_eta);
    ++work().completed_steps;
  }
};

// The completed runs own encoded bytes, NOT thousands of live MPFR intervals.
// Endpoint hex digits and binary exponent decode EXACTLY at their emitted
// precision. Oversized or inexact export/import refuses instead of truncating.
struct EncodedEndpoint {std::array<char,128> bytes{};};
struct WireBox {EncodedEndpoint low,high;};
struct WireRecord {unsigned kind=0,ki=0,ti=0,coordinate=0;WireBox box;};
struct Records {
  std::unique_ptr<WireRecord[]> owner;
  std::size_t used=0,reserved=0;
  std::size_t size() const noexcept{return used;}
  std::size_t capacity() const noexcept{return reserved;}
  WireRecord &operator[](std::size_t i){return owner[i];}
  const WireRecord &operator[](std::size_t i) const{return owner[i];}
  const WireRecord *begin() const noexcept{return owner.get();}
  const WireRecord *end() const noexcept{return used?owner.get()+used:owner.get();}
};
struct ClosedRun {
  unsigned subdivisions=0,precision=0;
  bool complete=false;
  Records records;
  std::size_t accounted_bytes=0;
  ~ClosedRun(){if(active)work().byte_owners-=accounted_bytes;}
  ClosedRun()=default;
  ClosedRun(const ClosedRun &)=delete;
  ClosedRun &operator=(const ClosedRun &)=delete;
  void reserve(std::size_t count) {
    work().guard();if(count>1000||accounted_bytes)throw Refusal("closed record count/reservation ceiling");
    const auto request=count*sizeof(WireRecord);work().bytes(request);
    // Exactly count requested record objects, admitted before array allocation
    // and value initialization. Allocator book-keeping is provider-internal,
    // separately bounded by the ordinary whole-process supervisor.
    try{work().copy(0);work().encoded_copy_bytes+=request;
      records.owner.reset(new WireRecord[count]);}
    catch(...){work().byte_owners-=request;throw;}
    records.reserved=count;accounted_bytes=request;
  }
  void append(WireRecord &&r) {
    work().guard();if(records.size()==records.capacity())throw Refusal("unreserved record publication");
    work().copy(0);work().encoded_copy_bytes+=sizeof(r);records.owner[records.used]=std::move(r);++records.used;
  }
};
void encode_endpoint(EncodedEndpoint &out,mpfr_srcptr value) {
  Bytes buffers(128+32+sizeof(mpfr_exp_t));work().copy(0);work().encoded_copy_bytes+=160;
  std::array<char,128> digits{};mpfr_exp_t exponent=0;
  work().charge(1,1);++work().serializations;
  if(!mpfr_get_str(digits.data(),&exponent,16,work().precision/4+1,value,MPFR_RNDN))
    throw Refusal("endpoint export refused");
  const auto length=std::strlen(digits.data());const auto sign=digits[0]=='-'?1u:0u;
  const auto e2=4*(exponent-static_cast<mpfr_exp_t>(length-sign));
  std::array<char,32> power{};
  const auto result=std::to_chars(power.data(),power.data()+power.size(),e2);
  const auto power_size=static_cast<std::size_t>(result.ptr-power.data());
  work().guard();if(result.ec!=std::errc{}||length+1+power_size+1>out.bytes.size())throw Refusal("endpoint wire128 limit");
  work().copy(0);work().encoded_copy_bytes+=length+power_size+2;
  std::copy_n(digits.data(),length,out.bytes.data());out.bytes[length]='@';
  std::copy_n(power.data(),power_size,out.bytes.data()+length+1);out.bytes[length+1+power_size]='\0';
}
void decode_endpoint(mpfr_ptr out,const EncodedEndpoint &wire) {
  Bytes buffer(128);
  work().copy(0);work().encoded_copy_bytes+=128;std::array<char,128> digits=wire.bytes;
  const auto end=std::find(digits.begin(),digits.end(),'\0');
  const auto mark=std::find(digits.begin(),end,'@');
  if(mark==end||mark==digits.begin()||end==digits.end())throw Refusal("missing dyadic exponent/terminator");
  const auto position=static_cast<std::size_t>(mark-digits.begin());digits[position]='\0';
  long exponent=0;const auto parsed=std::from_chars(digits.data()+position+1,digits.data()+(end-digits.begin()),exponent);
  work().guard();if(parsed.ec!=std::errc{}||parsed.ptr!=digits.data()+(end-digits.begin())||exponent< -16384||exponent>16384)
    throw Refusal("endpoint exponent/profile refusal");
  work().charge(1,1);if(mpfr_set_str(out,digits.data(),16,MPFR_RNDN)!=0)throw Refusal("inexact endpoint integer decode");
  work().charge(0,1);if(mpfr_mul_2si(out,out,exponent,MPFR_RNDN)!=0)throw Refusal("inexact endpoint exponent decode");
}
void encode(WireBox &out,const Interval &value,Arithmetic &a,Interval &roundtrip) {
  a.check(value);encode_endpoint(out.low,value.lo);encode_endpoint(out.high,value.hi);
  decode_endpoint(roundtrip.lo,out.low);decode_endpoint(roundtrip.hi,out.high);
  work().guard(2);const int low=mpfr_cmp(value.lo,roundtrip.lo),high=mpfr_cmp(value.hi,roundtrip.hi);
  if(low||high)throw Refusal("endpoint changed in byte roundtrip");
}
void decode(Interval &out,const WireBox &wire,Arithmetic &a){decode_endpoint(out.lo,wire.low);decode_endpoint(out.hi,wire.high);a.check(out);}

void canonical_core(Interval &out,unsigned coordinate,const State &y,Equations &e) {
  auto &a=e.arithmetic;
  if(coordinate==0)out.set(y[0]);
  else if(coordinate==1){a.scale(out,y[1],3);a.divide_integer(out,out,4);}
  else if(coordinate==2)a.divide_integer(out,y[2],2);
  else if(coordinate<6)out.set(y[G(coordinate-3)]);
  else if(coordinate==6)out.set(y[B]);
  else if(coordinate==7){a.reciprocal(e.t[4],e.k);a.multiply(out,y[TB],e.t[4]);}
  else if(coordinate==8)out.set(y[C]);
  else if(coordinate==9){a.reciprocal(e.t[4],e.k);a.multiply(out,y[TC],e.t[4]);}
  else out.set(y[PH]);
}
// Exact one-cell opacity integral, raw sources and Hamiltonian are independent
// of native helpers. The actual fixture has two original knots. Other inputs
// causally refuse, rather than crossing or inventing an opacity cell.
void raw_outputs(std::array<Interval,4> &channels,Interval &constraint,
                 Interval &survival,Interval &tau,const State &y,const Interval &eta,
                 Table &jets,Equations &e) {
  auto &a=e.arithmetic;e.fill(jets,y,eta,1);
  a.sub(e.t[4],e.eta_right,eta);a.sub(e.t[5],e.eta_right,e.eta_left);
  a.multiply(e.t[5],e.t[5],e.slope);a.add(e.t[5],e.t[5],e.opacity_left);
  a.add(e.t[5],e.t[5],e.kjet[0]);a.multiply(tau,e.t[4],e.t[5]);a.divide_integer(tau,tau,2);
  a.scale(e.t[4],tau,-1);a.exponential(survival,e.t[4]);
  if(!a.positive(survival))throw Refusal("positive finite survival unavailable");
  a.multiply(e.t[4],survival,e.kjet[0]); // g
  a.divide_integer(e.t[5],y[0],4);a.multiply(e.t[5],e.t[5],e.t[4]);
  a.multiply(channels[0],survival,e.pp[0]);a.add(channels[0],channels[0],e.t[5]);
  a.multiply(channels[1],survival,e.psi[0]);a.multiply(channels[1],channels[1],e.k);
  a.reciprocal(e.t[5],e.k);a.multiply(e.t[5],e.t[5],y[TB]);a.multiply(e.t[5],e.t[5],e.t[4]);a.add(channels[1],channels[1],e.t[5]);
  a.add(e.t[5],y[2],y[G(0)]);a.add(e.t[5],e.t[5],y[G(2)]);a.multiply(channels[2],e.t[5],e.t[4]);a.divide_integer(channels[2],channels[2],8);
  a.integer(e.t[5],6);a.square_root(e.t[5],e.t[5]);a.multiply(channels[3],e.t[5],channels[2]);
  a.multiply(constraint,e.k2,y[PH]);
  a.multiply(e.t[5],e.h[0],e.psi[0]);a.add(e.t[5],e.pp[0],e.t[5]);a.multiply(e.t[5],e.t[5],e.h[0]);a.scale(e.t[5],e.t[5],3);a.add(constraint,constraint,e.t[5]);
  a.multiply(e.t[5],e.fc[0],y[C]);a.multiply(e.t[6],e.fb[0],y[B]);a.add(e.t[5],e.t[5],e.t[6]);
  a.multiply(e.t[6],e.fg[0],y[0]);a.add(e.t[5],e.t[5],e.t[6]);a.multiply(e.t[5],e.t[5],e.h2[0]);a.scale(e.t[5],e.t[5],3);a.divide_integer(e.t[5],e.t[5],2);a.add(constraint,constraint,e.t[5]);
}
void require(bool value,const char *message){work().guard();if(!value)throw Refusal(message);}
template<class F> void refusal_control(F &&f,const char *name) {
  bool refused=false;try{f();}catch(const Refusal &){refused=true;}
  require(refused,name);work().status_call();mpfr_clear_flags();
}
void profile() {
  require(GMP_NUMB_BITS==64 && GMP_NAIL_BITS==0 && sizeof(mp_limb_t)==8,"fixed custom limb profile");
  require(std::numeric_limits<W>::radix==2 && std::numeric_limits<W>::digits<=192 &&
          std::numeric_limits<double>::is_iec559 && std::fegetround()==FE_TONEAREST,"stored IEEE input profile");
  work().status_call();const auto *version=mpfr_get_version();
  require(std::strcmp(version,"4.2.2")==0 && std::strcmp(gmp_version,"6.3.0")==0,"pinned backend version");
  work().status_call();const auto emin=mpfr_get_emin();work().status_call();const auto emax=mpfr_get_emax();
  require(emin<=-16384 && emax>=16384,"pinned exponent support");
  require(sizeof(Interval)==128&&sizeof(WireRecord)==272&&sizeof(Stepper)==872960,"declared custom/frame layout");
  work().status_call();mpfr_clear_flags();
  std::cout<<"PROFILE MPFR="<<version<<" GMP="<<gmp_version<<" Wdigits="<<std::numeric_limits<W>::digits
           <<" directed-basic/fma/sqrt/exp/custom-storage=conditional assumptions;controls are not compiler proof\n";
  std::cout<<"LAYOUT interval="<<sizeof(Interval)<<" record="<<sizeof(WireRecord)<<" stepper="<<sizeof(Stepper)
           <<" interval_live_source_upper=7240 no-elide-owned-copies=true;provider-internal scratch excluded\n";
}
void append(ClosedRun &run,unsigned kind,unsigned ki,unsigned ti,unsigned coordinate,
            const Interval &value,Arithmetic &a,Interval &roundtrip) {
  Bytes local(sizeof(WireRecord));work().copy(0);work().encoded_copy_bytes+=sizeof(WireRecord);
  WireRecord record;record.kind=kind;record.ki=ki;record.ti=ti;record.coordinate=coordinate;
  encode(record.box,value,a,roundtrip);run.append(std::move(record));
  const auto &published=run.records[run.records.size()-1];
  // Publish as soon as this exact record is complete. A later refinement,
  // allocation or profile refusal cannot destroy earlier completed bytes.
  std::cout<<"REFERENCE_RECORD n="<<run.subdivisions<<" precision="<<run.precision<<" kind="<<kind
    <<" k="<<ki<<" ti="<<ti<<" coordinate="<<coordinate<<" partial_run=true interval="
    <<published.box.low.bytes.data()<<','<<published.box.high.bytes.data()<<'\n';
  require(std::cout.good(),"exact reference record output refusal");
}
void print_endpoint(mpfr_srcptr x) {
  Bytes local(sizeof(EncodedEndpoint));work().copy(0);work().encoded_copy_bytes+=sizeof(EncodedEndpoint);
  EncodedEndpoint value;encode_endpoint(value,x);std::cout<<value.bytes.data();
}
void print_box(const Interval &x){print_endpoint(x.lo);std::cout<<',';print_endpoint(x.hi);}
void center(Interval &out,const Interval &x) {
  work().charge(1,2);mpfr_add(out.lo,x.lo,x.hi,MPFR_RNDN);mpfr_div_2ui(out.lo,out.lo,1,MPFR_RNDN);
  work().copy(1);mpfr_set(out.hi,out.lo,MPFR_RNDN);
}
bool within(const Interval &error,const Interval &allowance) {
  work().guard(1);return mpfr_cmp(error.hi,allowance.lo)<=0;
}
bool overlaps(const Interval &x,const Interval &y) {
  work().guard(2);const int left=mpfr_cmp(x.lo,y.hi),right=mpfr_cmp(y.lo,x.hi);return left<=0 && right<=0;
}
void rational_control(Arithmetic &a,Interval &x,unsigned numerator,unsigned denominator) {
  // This exact Fraction oracle is GMP rational, separate from the interval
  // convolution. Reserve the struct and a conservative128-byte numerator/
  // denominator owner; GMP internal instructions are not claimed as counted.
  Bytes fraction(sizeof(mpq_t)+128);mpq_t value;work().status_call();mpq_init(value);
  work().charge(1,0);mpq_set_ui(value,numerator,denominator);
  work().guard(2);const int left=mpfr_cmp_q(x.lo,value),right=mpfr_cmp_q(x.hi,value);const bool ok=left<=0 && right>=0;
  work().status_call();mpq_clear(value);require(ok,"independent exact rational oracle");a.check(x);
}
// Meaningful bounded operator/interval controls. They do not select a new
// physical source, substitute its seed, or establish a global libm certificate.
void arithmetic_and_operator_controls() {
  Arithmetic a;std::array<Interval,24> x;
  a.integer(x[0],1);a.divide_integer(x[1],x[0],3);rational_control(a,x[1],1,3);
  a.scale(x[2],x[1],-3);a.integer(x[3],-1);require(a.contains(x[2],x[3]),"signed rational interval");
  a.stored(x[4],.125L);a.square_root(x[5],x[4]);a.multiply(x[6],x[5],x[5]);
  require(a.contains(x[6],x[4]),"sqrt residual enclosure");
  a.integer(x[4],4);a.square_root(x[5],x[4]);a.integer(x[6],2);
  require(a.contains(x[5],x[6])&&a.contains(x[6],x[5]),"exact dyadic sqrt");
  a.integer(x[7],0);a.exponential(x[8],x[7]);require(a.contains(x[8],x[0]),"exact exp0");
  a.integer(x[9],-1);a.integer(x[10],2);a.hull(x[11],x[9],x[10]);
  a.integer(x[12],0);a.linear_term(x[12],x[1],x[11],false);a.multiply(x[13],x[1],x[11]);
  require(a.contains(x[13],x[12]),"signed FMA endpoint corner control");
  a.scale(x[14],x[1],-1);a.integer(x[12],0);a.linear_term(x[12],x[14],x[11],true);
  a.multiply(x[13],x[14],x[11]);require(a.contains(x[13],x[12]),"negative FMA coefficient");
  refusal_control([&]{a.reciprocal(x[12],x[11]);},"zero-support inverse refusal");
  refusal_control([&]{a.square_root(x[12],x[9]);},"negative sqrt refusal");
  refusal_control([&]{a.stored(x[12],std::numeric_limits<W>::quiet_NaN());},"nonfinite import refusal");
  work().charge(1,1);mpfr_set_nan(x[12].lo);refusal_control([&]{a.check(x[12]);},"nonfinite result refusal before aggregate");
  a.stored(x[12],0x1p-1000L);a.reciprocal(x[13],x[12]);a.multiply(x[14],x[12],x[13]);
  require(a.contains(x[14],x[0]),"near-positive denominator retained without jitter");
  refusal_control([&]{a.require_bound(x[11],0,1);},"too-wide declared support refusal");
  a.integer(x[12],0);require(!a.contains(x[12],x[1]),"deliberately wrong narrow rational caught");
  // Polynomial normalized jets and reciprocal/sqrt recurrences have exact
  // independently chosen coefficients. Horner is checked against Fraction.
  a.stored(x[12],.125L);a.integer(x[13],1);a.integer(x[14],2);
  a.multiply(x[15],x[12],x[12]);a.add(x[15],x[15],x[14]);a.multiply(x[15],x[15],x[12]);
  a.add(x[15],x[15],x[13]);rational_control(a,x[15],641,512);
  // Taylor4 for exp over[0,1/8] needs tube-wide J5=exp(B)/120.
  a.integer(x[15],1);for(unsigned j=4;j>0;--j){a.divide_integer(x[16],x[0],j);
    a.multiply(x[15],x[15],x[16]);a.multiply(x[15],x[15],x[12]);a.add(x[15],x[15],x[0]);}
  // Above nested coefficients are 1+h(1+h/2(1+h/3(1+h/4))).
  a.exponential(x[16],x[12]);a.multiply(x[17],x[12],x[12]);a.multiply(x[18],x[17],x[17]);a.multiply(x[18],x[18],x[12]);
  a.multiply(x[19],x[18],x[16]);a.divide_integer(x[19],x[19],120);
  a.integer(x[20],0);a.hull(x[19],x[20],x[19]);a.add(x[19],x[19],x[15]);
  require(a.contains(x[19],x[16]),"whole-tube fifth derivative Taylor remainder");
  a.divide_integer(x[18],x[18],120);a.add(x[18],x[18],x[15]);
  require(!a.contains(x[18],x[16]),"point-start fifth derivative is inadequate");
  // Thomson exchange: fg=3/8,fb=1/4 gives 4fg/(3fb)=2. Slip
  // has eigenvalue -3K and enthalpy-weighted momentum is conserved.
  a.integer(x[1],2);a.integer(x[2],3);a.stored(x[3],.25L);a.stored(x[4],.375L);
  a.integer(x[5],1);a.integer(x[6],-1);a.sub(x[7],x[5],x[6]);a.scale(x[8],x[7],-1);a.multiply(x[8],x[8],x[3]);
  a.scale(x[9],x[8],-2);a.scale(x[10],x[8],2);a.add(x[10],x[10],x[9]);a.integer(x[11],0);
  require(a.contains(x[10],x[11]),"Thomson weighted exchange exact");
  a.sub(x[10],x[8],x[9]);a.multiply(x[12],x[7],x[3]);a.scale(x[12],x[12],-3);
  require(overlaps(x[10],x[12]),"Thomson slip rate has negative sign and baryon factor");
  a.scale(x[13],x[12],-1);require(!overlaps(x[10],x[13]),"reversed collision sign caught");
  // Collision-only Pi'=-(3/10)K Pi; vector(1,5,1) is its slow
  // mode, while (1,-1,0) and(1,0,-1) have Pi0 and eigenvalue-K.
  for(unsigned mode=0;mode<3;++mode){a.integer(x[5],1);a.integer(x[6],mode==0?5:mode==1?-1:0);
    a.integer(x[7],mode==0?1:mode==2?-1:0);a.add(x[8],x[5],x[6]);a.add(x[8],x[8],x[7]);
    for(unsigned i=0;i<3;++i){a.divide_integer(x[9],x[8],i==1?2:10);a.sub(x[9],x[9],x[5+i]);
      a.multiply(x[9],x[9],x[3]);a.multiply(x[10],x[5+i],x[3]);a.scale(x[10],x[10],-1);
      if(mode==0){a.scale(x[10],x[10],3);a.divide_integer(x[10],x[10],10);}
      require(overlaps(x[9],x[10]),"three polarized Thomson eigenmodes");}}
  // Stationary forced moments sigma=16ku/(45K),G0=5sigma/2,G2=sigma/2.
  a.stored(x[5],.125L);a.stored(x[6],.25L);a.multiply(x[7],x[5],x[6]);a.scale(x[7],x[7],16);
  a.divide_integer(x[7],x[7],45);a.reciprocal(x[8],x[3]);a.multiply(x[7],x[7],x[8]);
  a.scale(x[8],x[7],5);a.divide_integer(x[8],x[8],2);a.divide_integer(x[9],x[7],2);
  a.scale(x[10],x[7],2);a.add(x[10],x[10],x[8]);a.add(x[10],x[10],x[9]);
  a.multiply(x[11],x[5],x[6]);a.scale(x[11],x[11],8);a.divide_integer(x[11],x[11],15);
  a.divide_integer(x[12],x[10],10);a.scale(x[13],x[7],2);a.sub(x[12],x[12],x[13]);
  a.multiply(x[12],x[12],x[3]);a.add(x[12],x[12],x[11]);require(a.contains(x[12],x[20]),"stationary polarized16/45 shear");
  // Positive survival, correctly signed tau derivative and Psrc normalization.
  a.scale(x[5],x[3],-1);a.exponential(x[6],x[5]);require(a.positive(x[6]),"strictly positive surviving boundary");
  work().guard(1);require(mpfr_cmp(x[6].hi,x[0].lo)<0,"positive opacity survival below1");
  a.multiply(x[7],x[6],x[3]);a.integer(x[8],6);a.square_root(x[8],x[8]);a.multiply(x[9],x[7],x[8]);
  a.multiply(x[10],x[9],x[9]);a.multiply(x[11],x[7],x[7]);a.scale(x[11],x[11],6);
  require(overlaps(x[10],x[11])&&!overlaps(x[9],x[7]),"Psrc sqrt6 versus missing factor");
  // Byte-prefix and actual no-elide interval lifetime/codec controls.
  const auto live=work().live_intervals;{Interval copy=x[1];Interval moved=std::move(copy);a.check(moved);}
  require(work().live_intervals==live,"complete no-elide value lifetime");
  {Bytes local(sizeof(WireBox));WireBox wire;encode(wire,x[1],a,x[22]);decode(x[23],wire,a);require(overlaps(x[1],x[23]),"exact endpoint wire roundtrip");
    wire.low.bytes[0]='?';refusal_control([&]{decode(x[23],wire,a);},"malformed wire causal refusal");}
  const auto old_limit=work().primitive_limit;work().primitive_limit=work().primitives;
  bool denied=false;try{a.integer(x[12],0);}catch(const Refusal &){denied=true;}
  work().primitive_limit=old_limit;require(denied,"work prefix denies before assignment");
  std::cout<<"CONTROLS exact-Fraction/interval/Taylor/Thomson/polarization/source-sign/codec/lifetime/refusal passed\n";
}
void passive_bessel_control() {
  Arithmetic a;std::array<Interval,7> value,derivative;std::array<Interval,12> t;
  a.stored(t[0],.125L);a.multiply(t[1],t[0],t[0]);a.reciprocal(t[2],t[0]);a.integer(t[11],0);
  // j_l(z)=z^l/(2l+1)!! sum_r(-z²)^r/[2^r r!
  // (2l+3)(2l+5)...]. At z=1/8 both value and derivative terms
  // strictly decrease; the signed next term bounds the alternating remainder.
  for(unsigned l=0;l<7;++l){a.integer(t[3],1);
    for(unsigned n=0;n<l;++n)a.multiply(t[3],t[3],t[0]);
    for(unsigned n=1;n<=2*l+1;n+=2)a.divide_integer(t[3],t[3],n);
    a.integer(t[4],0);a.integer(t[5],0);
    for(unsigned r=0;r<8;++r){a.add(t[4],t[4],t[3]);a.scale(t[6],t[3],l+2*r);a.multiply(t[6],t[6],t[2]);a.add(t[5],t[5],t[6]);
      a.multiply(t[3],t[3],t[1]);a.scale(t[3],t[3],-1);a.divide_integer(t[3],t[3],2*(r+1)*(2*l+2*r+3));}
    a.add(t[6],t[4],t[3]);a.hull(value[l],t[4],t[6]);
    a.scale(t[6],t[3],l+16);a.multiply(t[6],t[6],t[2]);a.add(t[6],t[5],t[6]);a.hull(derivative[l],t[5],t[6]);
    require(a.positive(value[l]),"small-phase Bessel sign");
  }
  for(unsigned l=1;l<6;++l){a.scale(t[3],value[l],2*l+1);a.multiply(t[3],t[3],t[2]);
    a.sub(t[3],value[l-1],t[3]);a.add(t[3],t[3],value[l+1]);require(a.contains(t[3],t[11]),"signed Bessel streaming recurrence");
    a.scale(t[3],value[l],l+1);a.multiply(t[3],t[3],t[2]);a.sub(t[3],value[l-1],t[3]);
    require(overlaps(t[3],derivative[l]),"absolute-eta passive terminal law");
    a.scale(t[4],value[l],l+1);a.multiply(t[4],t[4],t[2]);a.add(t[4],t[4],value[l-1]);
    require(!overlaps(t[4],derivative[l]),"wrong terminal damping sign caught");}
  std::cout<<"CONTROLS passive small-phase Bessel signed streaming/absolute-eta terminal passed; not finite seed\n";
}

// Inherited synthetic axis construction is an input builder ONLY. Its rounded
// Simpson age and native thermal callback are not an independent clock oracle.
irred::cosmology::FiniteOpacityRequest fixture(std::size_t &calls) {
  using namespace irred::cosmology;FiniteOpacityRequest r;
  r.model={67.4,.0224,.12,2.7255,0,{}};r.initial_scale_factor=.001;r.observer_scale_factor=.00101;
  const auto mapped=map_thermal_physical_model(r.model);require(mapped.model.has_value(),"fixture map");
  const auto bg=prepare_thermal_background(*mapped.model);require(bg.status()==S::ok,"fixture background");
  const auto age=[&](W limit){constexpr unsigned n=16384;W sum=0;
    const W c=W(irred::speed_of_light_m_per_s)/1000;
    for(unsigned i=0;i<=n;++i){++calls;const auto p=bg.scaled_expansion(limit*i/n);
      require(p.status==S::ok&&p.a4_e2>0,"inherited fixture age callback");
      sum+=(i==0||i==n?1:i%2?4:2)*c/bg.source().h0_km_s_mpc/std::sqrt(p.a4_e2);}
    return limit*sum/(3*n);};
  r.opacity.eta_mpc={double(age(r.initial_scale_factor)),double(age(r.observer_scale_factor))};
  r.opacity.differential_opacity_mpc_inverse={.1,.2};r.opacity.origin="synthetic-no-electron-history";
  r.opacity.law_id="emitted-linear-conformal-opacity/v1";
  r.opacity.exact_member_digest="synthetic-fixture-axis-construction-not-content-certification";
  r.k_mpc_inverse={.005,.01};r.source_origin="synthetic-finite-IVP-native-contract";return r;
}
void capture(Input &in,const irred::cosmology::FiniteOpacitySourceResult &native) {
  require(native.status==S::conditioning_budget_exceeded && native.source && native.identity &&
          !native.source_numerically_admitted,"complete but unadmitted native source");
  require(native.attempts.size()==9 && native.source->eta_mpc.size()==5 &&
          native.source->t0.size()==10 && native.source->t1.size()==10 && native.source->t2.size()==10 &&
          native.source->polarization.size()==10 && native.identity->seeds.size()==2 &&
          native.identity->original.k_mpc_inverse.size()==2 && native.identity->original.opacity.eta_mpc.size()==2 &&
          native.identity->original.opacity.differential_opacity_mpc_inverse.size()==2,"full native source shape");
  require(native.source->k_mpc_inverse.size()==2 &&
          native.source->k_mpc_inverse==native.identity->original.k_mpc_inverse &&
          native.source->eta_mpc.front()==native.identity->original.opacity.eta_mpc.front() &&
          native.source->eta_mpc.back()==native.identity->original.opacity.eta_mpc.back(),"same original source axis/order");
  for(const auto &attempt:native.attempts)require(attempt.status==S::ok&&attempt.reached_wavenumbers==2&&
      attempt.nodes.size()==10&&attempt.completed_steps==attempt.attempted_steps,"complete native attempt");
  const auto &last=native.attempts.back();require(last.hierarchy==192&&last.time_refinement==2&&
      last.final_temperature_tail.size()==380&&last.final_polarization_tail.size()==380,"endpoint-only native tails");
  for(unsigned ki=0;ki<2;++ki)for(unsigned ti=0;ti<5;++ti){const auto &node=last.nodes[ki*5+ti];
    require(node.k_index==ki&&node.eta_index==ti&&node.eta_mpc==native.source->eta_mpc[ti],"actual ordered node lineage");}
  for(const auto &d:native.diagnostics)require(!d.common_background_clock_error && !d.arithmetic_linear_error &&
      !d.source_grid_error,"native upstream missing error stays absent");
  const auto &id=*native.identity;const auto &bg=id.background.source();
  work().charge(37,0);work().copies+=37;
  in.photon=bg.omega_gamma;in.baryon=bg.omega_b;in.cdm=bg.omega_cdm;
  in.lambda=id.lambda_retained;in.h0=bg.h0_km_s_mpc;in.initial_a=id.original.initial_scale_factor;
  std::copy_n(native.source->eta_mpc.begin(),5,in.time.begin());
  std::copy_n(id.original.opacity.differential_opacity_mpc_inverse.begin(),2,in.opacity.begin());
  std::copy_n(id.original.k_mpc_inverse.begin(),2,in.k.begin());
  for(unsigned ki=0;ki<2;++ki)in.emitted_seed[ki]=id.seeds[ki].emitted_core;
  require(in.time[0]>=225&&in.time[0]<=640&&in.time[4]-in.time[0]<=2.5&&in.time[4]>in.time[0],"selected short low-phase source axis");
  for(unsigned ti=1;ti<5;++ti)require(in.time[ti]>in.time[ti-1],"actual emitted quarter order");
  require(native.boundary()&&native.boundary()->survival_i>0&&native.boundary()->tau_i>0&&
          native.boundary()->photon_monopole.size()==2&&native.boundary()->photon_dipole_theta_over_k.size()==2&&
          native.boundary()->omitted_temperature_absolute_bound.size()==2,"native positive boundary retained");
  for(unsigned ki=0;ki<2;++ki)require(std::isfinite(native.boundary()->photon_monopole[ki])&&
      std::isfinite(native.boundary()->photon_dipole_theta_over_k[ki])&&
      std::isfinite(native.boundary()->omitted_temperature_absolute_bound[ki])&&
      native.boundary()->omitted_temperature_absolute_bound[ki]>0,"finite positive original boundary state");
}
void recurrence_and_zero_opacity_controls(const Input &original) {
  require(work().live_intervals==0,"isolated operator control lifetime");
  Bytes input_bytes(sizeof(Input));work().charge(37,0);work().copies+=37;Input control=original;
  work().charge(2,0);control.opacity={0,0};Stepper stepper(control,0);auto &e=stepper.e;auto &a=e.arithmetic;
  State y;Interval eta,h,expected;std::array<Interval,4> channels;Interval constraint,survival,tau;
  Jet first,product,inverse;
  for(unsigned j=0;j<=5;++j)a.integer(first[j],j<2?1:0);
  e.product(product,first,first,5);
  for(unsigned j=0;j<=5;++j){a.integer(expected,j==0||j==2?1:j==1?2:0);require(a.contains(product[j],expected),"actual jet affine-square recurrence");}
  e.inverse(inverse,first,5);
  for(unsigned j=0;j<=5;++j){a.integer(expected,j%2?-1:1);require(a.contains(inverse[j],expected),"actual jet reciprocal recurrence");}
  // sqrt((1+x)^2)=1+x, using the selected normalized sqrt recurrence.
  e.square_root_jet(inverse,product,5);
  for(unsigned j=1;j<=5;++j){
    a.integer(expected,j==1?1:0);require(a.contains(inverse[j],expected),"actual normalized sqrt jet exact polynomial");}
  stepper.initial(y,control,0);a.stored(eta,control.time[0]);a.stored(h,control.time[1]);a.sub(h,h,eta);a.divide_integer(h,h,64);
  // A separately labelled shear operator input; no change to the actual seed.
  a.stored(y[2],.002L);stepper.step(y,eta,h);e.fill(stepper.start_jets,y,eta,1);
  a.integer(expected,0);
  for(unsigned l=0;l<=L;++l)require(a.contains(y[G(l)],expected)&&a.contains(expected,y[G(l)]),"K0 keeps all193 polarization moments exactly zero");
  a.sub(e.t[8],y[B],y[C]);require(a.contains(e.t[8],expected),"K0 equal cold densities invariant");
  a.sub(e.t[8],y[TB],y[TC]);require(a.contains(e.t[8],expected),"K0 equal cold velocities invariant");
  work().guard(1);require(mpfr_cmp(e.psi[0].hi,y[PH].lo)<0,"positive photon shear implies psi below phi");
  raw_outputs(channels,constraint,survival,tau,y,eta,stepper.start_jets,e);
  for(unsigned i:{2u,3u})require(a.contains(channels[i],expected)&&a.contains(expected,channels[i]),"K0 raw quadrupole/polarization source exactly zero");
  a.integer(expected,1);require(a.contains(survival,expected)&&a.contains(expected,survival),"K0 survival exactly1");
  // Adversarial time step is rejected without modifying current full state or
  // reached eta; all work and the begun attempt remain in the global ledger.
  a.integer(h,1);const auto completed=work().completed_steps;
  refusal_control([&]{stepper.step(y,eta,h);},"large-h cell/Jacobian preflight refusal");
  require(work().completed_steps==completed,"refused full step has no partial commit");
  std::cout<<"CONTROLS actual jets/zero-opacity full392 operator/shear/large-h refusal passed; actual fixture unchanged\n";
}
void output(ClosedRun &run,Stepper &stepper,const State &y,const Interval &eta,
            unsigned ki,unsigned ti) {
  auto &e=stepper.e;auto &a=e.arithmetic;
  std::array<Interval,4> channels;Interval core,constraint,survival,tau,roundtrip;
  for(unsigned i=0;i<11;++i){canonical_core(core,i,y,e);append(run,0,ki,ti,i,core,a,roundtrip);}
  raw_outputs(channels,constraint,survival,tau,y,eta,stepper.start_jets,e);
  for(unsigned i=0;i<4;++i)append(run,1,ki,ti,i,channels[i],a,roundtrip);
  append(run,2,ki,ti,0,y[SCALE],a,roundtrip);
  // Supplementary nominal Hamiltonian invariant and boundary records have no
  // unearned upstream error allocation. They are separate from acceptance.
  a.multiply(constraint,constraint,y[SCALE]);append(run,4,ki,ti,0,constraint,a,roundtrip);
  if(!ti){append(run,5,ki,ti,0,tau,a,roundtrip);append(run,5,ki,ti,1,survival,a,roundtrip);
    a.divide_integer(core,y[0],4);a.absolute(core,core);a.scale(e.t[7],y[1],3);a.divide_integer(e.t[7],e.t[7],4);
    a.absolute(e.t[7],e.t[7]);a.add(core,core,e.t[7]);a.multiply(core,core,survival);
    require(a.positive(core),"positive omitted boundary absolute bound");append(run,5,ki,ti,2,core,a,roundtrip);}
  if(ti==4)for(unsigned role=0;role<2;++role)for(unsigned l=3;l<=L;++l)
    append(run,3,ki,ti,role*(L-2)+l-3,y[role?G(l):l],a,roundtrip);
}
void campaign(ClosedRun &run,const Input &in,unsigned subdivisions,unsigned precision) {
  require(work().live_intervals==0,"precision switch requires no live interval values");
  require(precision==192||precision==256,"selected interval precision");work().precision=precision;
  work().status_call();mpfr_clear_flags();run.subdivisions=subdivisions;run.precision=precision;run.reserve(1000);
  for(unsigned ki=0;ki<2;++ki){Stepper stepper(in,ki);State y;Interval eta,h,left,increment;
    stepper.initial(y,in,ki);auto &a=stepper.e.arithmetic;a.stored(eta,in.time.front());
    unsigned reached_ti=0,reached_sub=0;output(run,stepper,y,eta,ki,0);
    try{for(unsigned ti=1;ti<5;++ti){a.stored(left,in.time[ti-1]);a.stored(increment,in.time[ti]);
        a.sub(h,increment,left);a.divide_integer(h,h,subdivisions);
        for(unsigned sub=0;sub<subdivisions;++sub){reached_ti=ti;reached_sub=sub;stepper.step(y,eta,h);}
        a.stored(increment,in.time[ti]);work().guard(2);
        const int low=mpfr_cmp(eta.lo,increment.lo),high=mpfr_cmp(eta.hi,increment.hi);
        require(low==0&&high==0,"exact emitted target landing without reset");
        output(run,stepper,y,eta,ki,ti);}
    }catch(const Refusal &error){
      // The campaign reserves one million primitives of the unchanged100M
      // global ceiling for failure serialization. A denied step never consumes
      // it. This is a stricter compute prefix, not a counter reset or new cap.
      work().primitive_limit=100000000;
      std::cout<<"REFERENCE_REFUSAL n="<<subdivisions<<" precision="<<precision<<" k="<<ki
               <<" requested_ti="<<reached_ti<<" completed_sub="<<reached_sub<<" reason="<<error.what()
               <<" earlier_rows="<<run.records.size()<<" complete=false eta=";print_box(eta);std::cout<<'\n';
      // Full392-state current prefix is retained in the log, never re-labelled
      // as an endpoint result. No serialization/truncation failure is hidden.
      for(unsigned i=0;i<D;++i){std::cout<<"PREFIX "<<i<<' ';print_box(y[i]);std::cout<<'\n';}throw;
    }
  }
  run.complete=true;std::cout<<"REFERENCE_RUN n="<<subdivisions<<" precision="<<precision<<" complete=true rows="<<run.records.size()<<'\n';
}

W nominal(const WireRecord &r,const irred::cosmology::FiniteOpacitySourceResult &native) {
  const auto &finest=native.attempts.back();const auto &node=finest.nodes[r.ki*5+r.ti];
  if(r.kind==0)return node.core[r.coordinate];
  if(r.kind==1){const auto index=r.ti*2+r.ki;const auto &s=*native.source;
    return r.coordinate==0?s.t0[index]:r.coordinate==1?s.t1[index]:r.coordinate==2?s.t2[index]:s.polarization[index];}
  if(r.kind==2)return node.scale_factor;
  if(r.kind==3){const unsigned role=r.coordinate/(L-2),l=r.coordinate%(L-2);
    return role?finest.final_polarization_tail[r.ki*(L-2)+l]:finest.final_temperature_tail[r.ki*(L-2)+l];}
  throw Refusal("supplementary record has no output allocation");
}
bool compare(const std::array<ClosedRun,4> &runs,const irred::cosmology::FiniteOpacitySourceResult &native) {
  require(work().live_intervals==0,"comparison precision switch");work().precision=256;
  Arithmetic a;std::array<Interval,22> x;bool accepted=true;std::size_t compared=0;
  for(const auto &run:runs)require(run.complete && run.records.size()==runs.back().records.size(),"all fixed complete refinements");
  for(std::size_t j=0;j<runs.back().records.size();++j){const auto &r=runs.back().records[j];
    for(unsigned run=0;run<4;++run){const auto &old=runs[run].records[j];
      require(old.kind==r.kind&&old.ki==r.ki&&old.ti==r.ti&&old.coordinate==r.coordinate,"original comparison mask/order");
      decode(x[run],old.box,a);}
    bool refinement=true;for(unsigned left=0;left<4;++left)for(unsigned right=left+1;right<4;++right)
      refinement=overlaps(x[left],x[right])&&refinement;
    if(r.kind>=4){
      std::cout<<"SUPPLEMENTARY kind="<<r.kind<<" k="<<r.ki<<" ti="<<r.ti<<" coordinate="<<r.coordinate
               <<" refinement_overlap="<<refinement<<" interval="<<r.box.low.bytes.data()<<','<<r.box.high.bytes.data()
               <<" upstream_error=unavailable\n";
      if(r.kind==5){const auto &boundary=*native.boundary();
        a.stored(x[19],r.coordinate==0?boundary.tau_i:r.coordinate==1?boundary.survival_i:boundary.omitted_temperature_absolute_bound[r.ki]);
        center(x[18],x[3]);a.sub(x[19],x[19],x[18]);a.absolute(x[19],x[19]);
        std::cout<<"BOUNDARY_NOMINAL_DIFFERENCE k="<<r.ki<<" coordinate="<<r.coordinate<<" interval=";
        print_box(x[19]);std::cout<<" allocation=unavailable;positive boundary is retained,not composed\n";
      }
      accepted=accepted&&refinement;
      if(r.kind==4){ // aE must retain the same emitted seed's small E_i.
        const WireRecord *initial=nullptr;for(const auto &candidate:runs.back().records)
          if(candidate.kind==4&&candidate.ki==r.ki&&candidate.ti==0)initial=&candidate;
        require(initial,"initial Hamiltonian record present");decode(x[20],initial->box,a);
        require(overlaps(x[20],x[3]),"nominal aE constraint invariant without reset");
      }continue;
    }
    a.stored(x[4],nominal(r,native));center(x[5],x[3]);center(x[6],x[2]);
    a.sub(x[7],x[3],x[5]);a.absolute(x[7],x[7]);a.upper(x[7],x[7]);
    a.sub(x[8],x[6],x[5]);a.absolute(x[8],x[8]);a.upper(x[8],x[8]);
    a.absolute(x[9],x[4]);
    // Original declared decimal policy constants enter as exact rationals.
    a.integer(x[10],3);a.divide_integer(x[10],x[10],100000);a.multiply(x[9],x[9],x[10]);
    a.integer(x[11],1);a.divide_integer(x[11],x[11],r.kind==1?1000000000:100000000);
    if(r.kind==1)a.divide_integer(x[11],x[11],100); //1e-11
    a.add(x[9],x[9],x[11]);
    a.scale(x[12],x[9],4);a.divide_integer(x[12],x[12],100); //radius4%
    a.divide_integer(x[13],x[9],200); //precision0.5%
    a.divide_integer(x[14],x[9],20); //combined5%
    a.add(x[15],x[7],x[8]);
    a.sub(x[16],x[4],x[5]);a.absolute(x[16],x[16]);a.add(x[17],x[16],x[15]);
    const bool row=refinement&&within(x[7],x[12])&&within(x[8],x[13])&&within(x[15],x[14])&&within(x[17],x[9]);
    // e_export=0 is earned only by every exact byte roundtrip and exact decode;
    // no decimal print value is used in the acceptance arithmetic.
    ++compared;accepted=accepted&&row;
    std::cout<<(row?"COMPARISON_PASS":"COMPARISON_REFUSAL")<<" kind="<<r.kind<<" k="<<r.ki
             <<" ti="<<r.ti<<" coordinate="<<r.coordinate<<" native=";
    print_box(x[4]);std::cout<<" interval="<<r.box.low.bytes.data()<<','<<r.box.high.bytes.data()<<" radius=";
    print_box(x[7]);std::cout<<" precision=";print_box(x[8]);std::cout<<" export=exact-roundtrip allocation=";
    print_box(x[9]);std::cout<<" outer_gap=";print_box(x[17]);std::cout<<'\n';
  }
  require(compared==920,"all11core/4channels/clock at5times and380endpointtails per k");return accepted;
}
void ledger_report(const Ledger &w) {
  std::cout<<"ACTUAL_REFERENCE_LEDGER primitives="<<w.primitives<<" endpoints="<<w.endpoints<<" guards="<<w.guards
    <<" copies="<<w.copies<<" status_calls="<<w.status_calls<<" begun_steps="<<w.steps<<" completed_steps="<<w.completed_steps
    <<" maximum_step_primitives="<<w.maximum_step_primitives<<" maximum_step_endpoints="<<w.maximum_step_endpoints
    <<" nonstep_primitives="<<w.nonstep_primitives<<" nonstep_endpoints="<<w.nonstep_endpoints
    <<" tube_checks="<<w.tube_checks<<" denied="<<w.denied<<" peak_intervals="<<w.peak_intervals
    <<" peak_requested_bytes="<<w.peak_bytes<<" live_intervals="<<w.live_intervals<<" requested_bytes="<<w.byte_owners
    <<" serializations="<<w.serializations<<" encoded_copy_bytes="<<w.encoded_copy_bytes<<'\n';
}
} // namespace reference

int main() {
  using namespace reference;Ledger ledger;active=&ledger;int code=1;
  try{
    // Native8MiB, input builder128KiB, bounded stdio/control metadata64KiB.
    // They are conservative requested-owner reservations, not actual RSS.
    Bytes fixed(8*1024*1024+128*1024+64*1024+sizeof(Input)+4*sizeof(ClosedRun));
    profile();arithmetic_and_operator_controls();passive_bessel_control();std::size_t fixture_calls=0;
    auto request=fixture(fixture_calls);irred::cosmology::FiniteOpacityPolicy policy;
    policy.maximum_attempted_steps=2000;policy.maximum_background_clock_calls=400000;
    policy.maximum_coupled_stage_solves=2000;policy.maximum_destination_writes=100000000;policy.maximum_native_bytes=8*1024*1024;
    auto producer=irred::cosmology::prepare_finite_opacity_source(std::move(request),policy);
    require(producer.status()==S::ok,"original native prepare");auto native=producer.produce(policy);
    require(native.peak_owned_payload_bound&&*native.peak_owned_payload_bound<=8*1024*1024,"native payload bound present");
    work().charge(37,0);work().copies+=37;Input input;capture(input,native);
    recurrence_and_zero_opacity_controls(input);
    std::cout<<"NATIVE fixture_callbacks="<<fixture_calls<<" producer_callbacks="<<native.work.background_clock_calls
      <<" begun_steps="<<native.work.attempted_steps<<" coupled_solves="<<native.work.coupled_stage_solves
      <<" destination_writes="<<native.work.destination_writes<<" payload="<<*native.peak_owned_payload_bound
      <<" admission=false compound_clock/source/physical_error=unavailable\n";
    require(native.work.background_clock_calls<=400000 && fixture_calls<=400000-native.work.background_clock_calls,"whole actual fixture/native callback cap");
    std::array<ClosedRun,4> runs;
    work().primitive_limit=99000000;
    campaign(runs[0],input,16,192);campaign(runs[1],input,32,192);
    campaign(runs[2],input,64,192);campaign(runs[3],input,64,256);
    work().primitive_limit=100000000;const bool passed=compare(runs,native);code=passed?0:1;
    std::cout<<(passed?"conditional finite-IVP full-state interval reference passed":"conditional finite-IVP comparison refused")
             <<"; native admission remains false; no upstream/physical/fullCLASS transfer\n";
  }catch(const Refusal &error){std::cerr<<"FINITE_REFERENCE_REFUSAL "<<error.what()<<'\n';}
  catch(const std::bad_alloc &){std::cerr<<"FINITE_REFERENCE_REFUSAL requested/provider allocation failure\n";}
  catch(const std::exception &error){std::cerr<<"FINITE_REFERENCE_REFUSAL "<<error.what()<<'\n';}
  ledger_report(ledger);if(ledger.live_intervals||ledger.byte_owners)code=1;active=nullptr;return code;
}
