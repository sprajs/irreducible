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
#include <cfenv>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace reference {
constexpr unsigned L=192, D=2*L+8, degree=5;
constexpr unsigned B=2*(L+1), TB=B+1, C=B+2, TC=B+3, PH=B+4, SCALE=B+5;
constexpr unsigned G(unsigned l){return L+1+l;}
using W=long double;
using S=irred::numerics::Status;
struct Refusal:std::runtime_error {using std::runtime_error::runtime_error;};
struct Ledger {
  std::size_t primitives=0,endpoints=0,guards=0,copies=0;
  std::size_t steps=0,tube_checks=0,denied=0,live_intervals=0,peak_intervals=0;
  std::size_t byte_owners=0,peak_bytes=0,serializations=0;
  std::size_t primitive_limit=100000000,endpoint_limit=800000000;
  unsigned precision=192;
  void charge(std::size_t p,std::size_t e) {
    if(p>primitive_limit || primitives>primitive_limit-p ||
       e>endpoint_limit || endpoints>endpoint_limit-e) {
      ++denied;throw Refusal("reference work prefix limit");
    }
    primitives+=p;endpoints+=e;
  }
  void guard(std::size_t e=0){charge(1,e);++guards;}
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
    work().guard(1);
    if(!mpfr_number_p(x.lo)||!mpfr_number_p(x.hi)||mpfr_cmp(x.lo,x.hi)>0)
      throw Refusal("invalid interval");
    const auto bad=MPFR_FLAGS_UNDERFLOW|MPFR_FLAGS_OVERFLOW|MPFR_FLAGS_NAN|
                   MPFR_FLAGS_ERANGE|MPFR_FLAGS_DIVBY0;
    if(mpfr_flags_test(bad))throw Refusal("MPFR exceptional arithmetic flag");
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
    auto &t=temporary[0];work().charge(1,2);
    mpfr_add(t.lo,a.lo,b.lo,MPFR_RNDD);mpfr_add(t.hi,a.hi,b.hi,MPFR_RNDU);r.set(t);
  }
  void sub(Interval &r,const Interval &a,const Interval &b) {
    auto &t=temporary[0];work().charge(1,2);
    mpfr_sub(t.lo,a.lo,b.hi,MPFR_RNDD);mpfr_sub(t.hi,a.hi,b.lo,MPFR_RNDU);r.set(t);
  }
  void scale(Interval &r,const Interval &a,long n) {
    auto &t=temporary[0];work().charge(1,2);
    mpfr_mul_si(t.lo,n>=0?a.lo:a.hi,n,MPFR_RNDD);
    mpfr_mul_si(t.hi,n>=0?a.hi:a.lo,n,MPFR_RNDU);r.set(t);
  }
  void divide_integer(Interval &r,const Interval &a,unsigned n) {
    if(!n)throw Refusal("zero rational denominator");
    auto &t=temporary[0];work().charge(1,2);
    mpfr_div_ui(t.lo,a.lo,n,MPFR_RNDD);mpfr_div_ui(t.hi,a.hi,n,MPFR_RNDU);r.set(t);
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
    if(mpfr_cmp_si(a.lo,0)<=0 && mpfr_cmp_si(a.hi,0)>=0)
      throw Refusal("interval denominator contains zero");
    auto &t=temporary[0];work().charge(1,2);
    mpfr_ui_div(t.lo,1,a.hi,MPFR_RNDD);mpfr_ui_div(t.hi,1,a.lo,MPFR_RNDU);r.set(t);
  }
  void square_root(Interval &r,const Interval &a) {
    work().guard(1);if(mpfr_cmp_si(a.lo,0)<0)throw Refusal("negative sqrt support");
    auto &t=temporary[0];work().charge(1,2);
    mpfr_sqrt(t.lo,a.lo,MPFR_RNDD);mpfr_sqrt(t.hi,a.hi,MPFR_RNDU);r.set(t);
  }
  void exponential(Interval &r,const Interval &a) {
    auto &t=temporary[0];work().charge(1,2);
    mpfr_exp(t.lo,a.lo,MPFR_RNDD);mpfr_exp(t.hi,a.hi,MPFR_RNDU);r.set(t);
  }
  void absolute(Interval &r,const Interval &a) {
    auto &t=temporary[0];work().guard(2);
    if(mpfr_cmp_si(a.lo,0)>=0){r.set(a);return;}
    if(mpfr_cmp_si(a.hi,0)<=0){scale(r,a,-1);return;}
    work().charge(1,2);mpfr_set_zero(t.lo,1);mpfr_neg(t.hi,a.lo,MPFR_RNDU);
    work().guard(1);if(mpfr_cmp(t.hi,a.hi)<0){work().copy(1);mpfr_set(t.hi,a.hi,MPFR_RNDU);}r.set(t);
  }
  bool contains(const Interval &outer,const Interval &inner,bool strict=false) {
    work().guard(2);const int l=mpfr_cmp(outer.lo,inner.lo),h=mpfr_cmp(outer.hi,inner.hi);
    return strict?(l<0&&h>0):(l<=0&&h>=0);
  }
  bool positive(const Interval &x){work().guard(1);return mpfr_cmp_si(x.lo,0)>0;}
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
  void coefficients(Table &y,const Interval &eta,unsigned n) {
    auto &a=arithmetic;
    product(a2,y[SCALE],y[SCALE],n);product(a4,a2,a2,n);
    a.add(t[3],baryon,cdm);
    for(unsigned j=0;j<=n;++j){a.multiply(p[j],t[3],y[SCALE][j]);
      a.multiply(t[0],lambda,a4[j]);a.add(p[j],p[j],t[0]);
      if(!j)a.add(p[j],p[j],photon);}
    if(!a.positive(p[0])||!a.positive(y[SCALE][0]))throw Refusal("nonpositive clock support");
    a.square_root(s[0],p[0]);a.scale(t[3],s[0],2);a.reciprocal(t[3],t[3]);
    for(unsigned j=1;j<=n;++j){a.integer(t[2],0);
      for(unsigned r=1;r<j;++r){a.multiply(t[0],s[r],s[j-r]);a.add(t[2],t[2],t[0]);}
      a.sub(t[2],p[j],t[2]);a.multiply(s[j],t[2],t[3]);}
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
  }
  // coefficient of a convolution, with separately owned destination.
  void convolution(Interval &out,const Jet &x,const Jet &y,unsigned j) {
    auto &a=arithmetic;a.integer(out,0);
    for(unsigned r=0;r<=j;++r){a.multiply(t[0],x[r],y[j-r]);a.add(out,out,t[0]);}
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
      for(unsigned r=0;r<=j;++r){a.reciprocal(t[6],k);a.multiply(t[5],y[TB][j-r],t[6]);a.scale(t[5],t[5],4);a.divide_integer(t[5],t[5],3);
        a.sub(t[5],t[5],y[1][j-r]);a.multiply(t[5],t[5],kjet[r]);a.add(t[4],t[4],t[5]);}
      a.divide_integer(y[1][j+1],t[4],j+1);
      for(unsigned i:{2u,G(0),G(1),G(2)}) {
        if(i==2){a.scale(t[4],y[1][j],2);a.scale(t[5],y[3][j],3);a.sub(t[4],t[4],t[5]);a.divide_integer(t[4],t[4],5);}
        else if(i==G(0)){a.scale(t[4],y[G(1)][j],-1);}
        else {const unsigned l=i-G(0);a.scale(t[4],y[G(l-1)][j],l);a.scale(t[5],y[G(l+1)][j],l+1);
              a.sub(t[4],t[4],t[5]);a.divide_integer(t[4],t[4],2*l+1);}
        a.multiply(t[4],t[4],k);
        for(unsigned r=0;r<=j;++r){
          a.add(t[5],y[2][j-r],y[G(0)][j-r]);a.add(t[5],t[5],y[G(2)][j-r]);
          if(i==G(1))a.integer(t[5],0);
          else a.divide_integer(t[5],t[5],i==G(0)?2:10);
          a.sub(t[5],y[i][j-r],t[5]);a.multiply(t[5],kjet[r],t[5]);a.sub(t[4],t[4],t[5]);}
        a.divide_integer(y[i][j+1],t[4],j+1);
      }
      for(unsigned base:{0u,G(0)})for(unsigned l=3;l<=L;++l) {
        const unsigned i=base+l;
        if(l<L){a.scale(t[4],y[i-1][j],l);a.scale(t[5],y[i+1][j],l+1);
          a.sub(t[4],t[4],t[5]);a.divide_integer(t[4],t[4],2*l+1);}
        else t[4].set(y[i-1][j]);
        a.multiply(t[4],t[4],k);convolution(t[5],kjet,y[i],j);a.sub(t[4],t[4],t[5]);
        if(l==L){convolution(t[5],ieta,y[i],j);a.scale(t[5],t[5],L+1);a.sub(t[4],t[4],t[5]);}
        a.divide_integer(y[i][j+1],t[4],j+1);
      }
      for(unsigned density:{B,C}) {
        a.scale(t[4],y[density+1][j],-1);a.scale(t[5],pp[j],3);a.add(t[4],t[4],t[5]);
        a.divide_integer(y[density][j+1],t[4],j+1);
        convolution(t[4],h,y[density+1],j);a.scale(t[4],t[4],-1);a.multiply(t[5],k2,psi[j]);a.add(t[4],t[4],t[5]);
        if(density==B)for(unsigned r=0;r<=j;++r)for(unsigned v=0;v<=j-r;++v){
          a.multiply(t[5],k,y[1][j-r-v]);a.scale(t[5],t[5],3);a.divide_integer(t[5],t[5],4);
          a.sub(t[5],t[5],y[TB][j-r-v]);a.multiply(t[5],t[5],ratio[r]);a.multiply(t[5],t[5],kjet[v]);a.add(t[4],t[4],t[5]);}
        a.divide_integer(y[density+1][j+1],t[4],j+1);
      }
      a.divide_integer(y[PH][j+1],pp[j],j+1);
    }
  }
};
} // namespace reference

// Remaining tube, output/campaign and control bodies follow in the next source
// checkpoint. This file is deliberately unwired until root's complete review.
