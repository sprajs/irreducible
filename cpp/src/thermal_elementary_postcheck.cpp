#include "thermal_elementary_postcheck.hpp"

#include <bit>
#include <cfenv>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>

#if defined(__GNUC__)
#define IRRED_ELEMENTARY_INLINE __attribute__((always_inline)) inline
#else
#define IRRED_ELEMENTARY_INLINE inline
#endif
namespace irred::cosmology::detail {
namespace {
using W = ElementaryWide;
using Work = ElementaryWork;
using Stage = ElementaryStage;
using Status = ElementaryStatus;
constexpr W zero=0.0L, one=1.0L, two=2.0L, three=3.0L;
constexpr W tail=0x1p-144L, range_lower=0x1p-128L, range_upper=0x1p128L;
constexpr W positive_max=std::numeric_limits<W>::max(), negative_max=-positive_max;
constexpr std::uint64_t fingerprint=0x4945454538300002ULL;
enum Count:unsigned { basic,neighbor,scaling,term,library,write,copy,guards,attempt };
enum Op:unsigned { add,subtract,multiply,divide,signed_add };

// Admission does not call itself. No floating object is constructed here.
enum class Admission {ok,invalid,overflow,limited};
IRRED_ELEMENTARY_INLINE Admission admit(ElementaryBudget &b,ElementaryRefusal &r,ElementaryControlWork &c,
           Stage stage,const Work &delta,std::uint64_t scalar) noexcept {
  ++c.transitions;
  const auto &v=b.ledger.served.counters;
  // Integer consistency is part of EVERY transaction before any commit/action.
  // It is not a recursively admitted floating predicate or an assumed cache.
  const bool valid=v[basic]<=UINT64_MAX-v[neighbor] &&
    v[basic]+v[neighbor]<=UINT64_MAX-v[library] &&
    v[basic]+v[neighbor]+v[library]==b.ledger.scalar_total &&
    v[scaling]<=v[basic] && v[copy]<=v[write];
  if(!valid) {
    r.occurred=true;r.stage=stage;r.requested_increment=delta;
    return Admission::invalid;
  }
  bool overflow=false;
  for(unsigned i=0;i<9;++i)
    overflow |= delta.counters[i]>UINT64_MAX-b.ledger.served.counters[i];
  overflow |= scalar>UINT64_MAX-b.ledger.scalar_total;
  const bool limited=b.ledger.scalar_total>b.maximum_scalar ||
    v[write]>b.maximum_writes || v[copy]>b.maximum_copies || v[guards]>b.maximum_guards ||
    (!overflow && (scalar>b.maximum_scalar-b.ledger.scalar_total ||
      delta.counters[write]>b.maximum_writes-v[write] ||
      delta.counters[copy]>b.maximum_copies-v[copy] ||
      delta.counters[guards]>b.maximum_guards-v[guards]));
  if(overflow || limited) {
    r.occurred=true; r.stage=stage; r.requested_increment=delta;
    return overflow?Admission::overflow:Admission::limited;
  }
  for(unsigned i=0;i<9;++i)b.ledger.served.counters[i]+=delta.counters[i];
  b.ledger.scalar_total+=scalar;
  return Admission::ok;
}
IRRED_ELEMENTARY_INLINE void finish_control(ElementaryControlWork &c) noexcept {
  c.integer_operation_upper=128*c.transitions+128*c.predicates+32*c.iterations+512;
}
IRRED_ELEMENTARY_INLINE bool profile() noexcept {
#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__) && defined(__linux__) && !defined(__FAST_MATH__) && (!defined(__FINITE_MATH_ONLY__) || __FINITE_MATH_ONLY__ == 0)
  if(CHAR_BIT!=8 || sizeof(void *)!=8 || sizeof(double)!=8 ||
     std::numeric_limits<double>::digits!=53 || !std::numeric_limits<double>::is_iec559 ||
     sizeof(W)!=16 || alignof(W)>16 || LDBL_MANT_DIG!=64 ||
     LDBL_MIN_EXP!=-16381 || LDBL_MAX_EXP!=16384 ||
     std::endian::native!=std::endian::little || std::fegetround()!=FE_TONEAREST)
    return false;
  unsigned short x87=0;
  unsigned mxcsr=0;
  __asm__ volatile("fnstcw %0":"=m"(x87));
  __asm__ volatile("stmxcsr %0":"=m"(mxcsr));
  if((x87&0x0f3f)!=0x033f || (mxcsr&0xffc0)!=0x1f80)return false;
  std::uint64_t significand=0; std::uint16_t exponent=0;
  std::memcpy(&significand,&one,8);
  std::memcpy(&exponent,reinterpret_cast<const unsigned char *>(&one)+8,2);
  return significand==0x8000000000000000ULL && exponent==0x3fff;
#else
  return false;
#endif
}
IRRED_ELEMENTARY_INLINE bool normal(const W &x) noexcept { return std::isnormal(x); }
IRRED_ELEMENTARY_INLINE bool overlap(const void *a,std::size_t an,const void *b,std::size_t bn) noexcept {
  const auto av=reinterpret_cast<std::uintptr_t>(a), bv=reinterpret_cast<std::uintptr_t>(b);
  return av<=bv ? bv-av<an : av-bv<bn;
}
struct Engine {
  ElementaryBudget &b;
  ElementaryScratch &s;
  ElementaryRefusal &refusal;
  ElementaryControlWork &control;
  Status &status;
  Stage &stage;
  Work &before; // Snapshot in the final integer receipt, not a second stack payload.

  IRRED_ELEMENTARY_INLINE Engine(ElementaryBudget &bb,ElementaryScratch &ss,ElementaryRefusal &rr,
         ElementaryControlWork &cc,Status &st,Stage &sg,Work &call) noexcept
    :b(bb),s(ss),refusal(rr),control(cc),status(st),stage(sg),before(call) {
    before=bb.ledger.served;
    refusal={}; control={}; status=Status::unavailable; stage=Stage::entry;
    s.validity.fill(0);
  }
  IRRED_ELEMENTARY_INLINE bool charge(const Work &delta,std::uint64_t scalar=0) noexcept {
    const auto admission=admit(b,refusal,control,stage,delta,scalar);
    if(admission==Admission::ok)return true;
    status=admission==Admission::invalid?Status::invalid_input:
      admission==Admission::overflow?Status::counter_overflow:Status::work_limit;
    return false;
  }
  template<class Predicate> IRRED_ELEMENTARY_INLINE bool decide(bool &answer,Predicate predicate) noexcept {
    if(!charge(Work{{0,0,0,0,0,0,0,1,0}}))return false;
    ++control.predicates; answer=predicate(); return true;
  }
  template<class Predicate> IRRED_ELEMENTARY_INLINE bool require(Predicate predicate,Status failure=Status::normal_margin) noexcept {
    bool answer=false;
    if(!decide(answer,predicate))return false;
    if(answer)return true;
    status=failure; refusal.occurred=true; refusal.stage=stage;
    refusal.requested_increment={}; return false;
  }
  IRRED_ELEMENTARY_INLINE bool begin() noexcept {
    if(!charge(Work{{0,0,0,0,0,0,0,0,1}}))return false;
    stage=Stage::profile;
    if(!require([]{return profile();},Status::unsupported_arithmetic))return false;
    return require([]{return normal(one)&&normal(two)&&normal(three)&&normal(tail)&&
      normal(range_lower)&&normal(range_upper)&&normal(positive_max)&&normal(negative_max);},
      Status::unsupported_arithmetic);
  }
  IRRED_ELEMENTARY_INLINE bool ownership(const W &x,const W &y,const void *output,std::size_t n,
                 const void *context=nullptr,std::size_t cn=0) noexcept {
    // Five finite leaves, each <=4 overlap calls (<=8 uintptr destinations).
    // Do not fold17 calls into the old16-write semantic-leaf reservation.
    return require([&]{
      return !overlap(&x,sizeof(W),&s,sizeof(s)) && !overlap(&y,sizeof(W),&s,sizeof(s)) &&
        !overlap(&x,sizeof(W),output,n) && !overlap(&y,sizeof(W),output,n);
    },Status::invalid_ownership) && require([&]{
      return !context || (!overlap(&x,sizeof(W),context,cn)&&!overlap(&y,sizeof(W),context,cn));
    },Status::invalid_ownership) && require([&]{
      return !overlap(&s,sizeof(s),output,n) && (!context || (!overlap(context,cn,output,n)&&
          !overlap(context,cn,&s,sizeof(s))));
    },Status::invalid_ownership) && require([&]{
      return !overlap(&x,sizeof(W),&b,sizeof(b)) && !overlap(&y,sizeof(W),&b,sizeof(b)) &&
        !overlap(&x,sizeof(W),&b.ledger,sizeof(b.ledger)) && !overlap(&y,sizeof(W),&b.ledger,sizeof(b.ledger));
    },Status::invalid_ownership) && require([&]{
      return !overlap(&s,sizeof(s),&b,sizeof(b)) && !overlap(&s,sizeof(s),&b.ledger,sizeof(b.ledger)) &&
        !overlap(output,n,&b,sizeof(b)) && !overlap(output,n,&b.ledger,sizeof(b.ledger));
    },Status::invalid_ownership);
  }
  IRRED_ELEMENTARY_INLINE bool copies(W &destination,const W &source) noexcept {
    if(!charge(Work{{0,0,0,0,0,1,1,0,0}}))return false;
    destination=source; return true;
  }
  IRRED_ELEMENTARY_INLINE bool raw(W &destination,const W &a,const W &d,Op op,bool scale=false) noexcept {
    if(!charge(Work{{1,0,scale?1u:0u,0,0,1,0,0,0}},1))return false;
    switch(op) {
      case add: case signed_add:destination=a+d;break; case subtract:destination=a-d;break;
      case multiply:destination=a*d;break; case divide:destination=a/d;break;
    }
    return require([&]{return normal(destination);});
  }
  IRRED_ELEMENTARY_INLINE bool negate(W &destination,const W &value) noexcept {
    if(!charge(Work{{1,0,0,0,0,1,0,0,0}},1))return false;
    destination=-value;
    return require([&]{return normal(destination);});
  }
  IRRED_ELEMENTARY_INLINE bool adjacent(W &out,const W &value,bool up) noexcept {
    if(!require([&]{return up?value<positive_max:value>negative_max;}))return false;
    if(!charge(Work{{0,1,0,0,0,1,0,0,0}},1))return false;
    out=std::nextafter(value,up?positive_max:negative_max);
    return require([&]{return normal(out);});
  }
  IRRED_ELEMENTARY_INLINE bool rounded(W &out,const W &a,const W &d,Op op,bool up) noexcept {
    bool exact_zero=false;
    if(op==subtract) {
      if(!decide(exact_zero,[&]{return a==d;}))return false;
    } else if(op==signed_add) {
      if(!negate(s.scalars[13],d))return false;
      if(!decide(exact_zero,[&]{return a==s.scalars[13];}))return false;
    }
    if(exact_zero)return copies(out,zero);
    return raw(s.scalars[0],a,d,op) && adjacent(out,s.scalars[0],up);
  }
  IRRED_ELEMENTARY_INLINE bool positive_binary(unsigned out,unsigned a,unsigned d,Op op) noexcept {
    if(op==divide && !require([&]{return s.intervals[d].lower>zero;}))return false;
    return rounded(s.intervals[out].lower,s.intervals[a].lower,
                   op==divide?s.intervals[d].upper:s.intervals[d].lower,op,false) &&
      rounded(s.intervals[out].upper,s.intervals[a].upper,
              op==divide?s.intervals[d].lower:s.intervals[d].upper,op,true);
  }
  IRRED_ELEMENTARY_INLINE bool point(unsigned interval,const W &value) noexcept {
    if(!charge(Work{{0,0,0,0,0,2,2,0,0}}))return false;
    s.intervals[interval].lower=value;s.intervals[interval].upper=value;return true;
  }
  IRRED_ELEMENTARY_INLINE bool integer(W &out,std::int32_t value) noexcept {
    if(!charge(Work{{1,0,0,0,0,1,0,0,0}},1))return false;
    out=static_cast<W>(value);
    return require([&]{return normal(out);});
  }
  IRRED_ELEMENTARY_INLINE bool series() noexcept {
    stage=Stage::series;
    if(!point(4,zero))return false;
    if(!charge(Work{{0,0,0,0,0,2,2,0,0}}))return false;
    s.intervals[2].lower=s.intervals[0].lower;s.intervals[2].upper=s.intervals[0].upper;
    unsigned power=2,next_power=3,sum=4,next_sum=5;
    for(std::int32_t j=0;j<48;++j) {
      ++control.iterations;
      if(!charge(Work{{0,0,0,1,0,0,0,0,0}}) || !integer(s.scalars[4],2*j+1))return false;
      if(!rounded(s.intervals[6].lower,s.intervals[power].lower,s.scalars[4],divide,false) ||
         !rounded(s.intervals[6].upper,s.intervals[power].upper,s.scalars[4],divide,true) ||
         !positive_binary(next_sum,sum,6,add))return false;
      std::swap(sum,next_sum);
      if(j!=47) {
        if(!positive_binary(next_power,power,1,multiply))return false;
        std::swap(power,next_power);
      }
    }
    return raw(s.intervals[7].lower,s.intervals[sum].lower,two,multiply) &&
      raw(s.intervals[7].upper,s.intervals[sum].upper,two,multiply) &&
      rounded(s.intervals[7].upper,s.intervals[7].upper,tail,add,true);
  }
  IRRED_ELEMENTARY_INLINE bool logarithm(const ElementaryOwnedInterval &log2,const W &x) noexcept {
    stage=Stage::input;
    if(!require([&]{return normal(x)&&x>zero&&x>=range_lower&&x<=range_upper;},Status::invalid_input) ||
       !copies(s.scalars[1],x))return false;
    std::int32_t exponent=0;
    stage=Stage::scaling;
    bool proceed=false;
    for(;;) {
      if(!decide(proceed,[&]{return s.scalars[1]<one;}))return false;
      if(!proceed)break;
      ++control.iterations;
      if(!require([&]{return exponent>-128;},Status::invalid_input) ||
         !raw(s.scalars[1],s.scalars[1],two,multiply,true))return false;
      --exponent;
    }
    for(;;) {
      if(!decide(proceed,[&]{return s.scalars[1]>=two;}))return false;
      if(!proceed)break;
      ++control.iterations;
      if(!require([&]{return exponent<128;},Status::invalid_input) ||
         !raw(s.scalars[1],s.scalars[1],two,divide,true))return false;
      ++exponent;
    }
    if(!decide(proceed,[&]{return s.scalars[1]==one;}))return false;
    if(proceed) { if(!point(7,zero))return false; }
    else {
      stage=Stage::rational;
      if(!raw(s.scalars[2],s.scalars[1],one,subtract) ||
         !rounded(s.intervals[13].lower,s.scalars[1],one,add,false) ||
         !rounded(s.intervals[13].upper,s.scalars[1],one,add,true) ||
         !rounded(s.intervals[0].lower,s.scalars[2],s.intervals[13].upper,divide,false) ||
         !rounded(s.intervals[0].upper,s.scalars[2],s.intervals[13].lower,divide,true) ||
         !positive_binary(1,0,0,multiply) || !series())return false;
    }
    if(exponent==0) {
      if(!charge(Work{{0,0,0,0,0,2,2,0,0}}))return false;
      s.intervals[9].lower=s.intervals[7].lower;s.intervals[9].upper=s.intervals[7].upper;
      return true;
    }
    if(!integer(s.scalars[3],exponent))return false;
    if(!rounded(s.intervals[8].lower,exponent>0?log2.lower:log2.upper,s.scalars[3],multiply,false) ||
       !rounded(s.intervals[8].upper,exponent>0?log2.upper:log2.lower,s.scalars[3],multiply,true) ||
       !rounded(s.intervals[9].lower,s.intervals[7].lower,s.intervals[8].lower,signed_add,false) ||
       !rounded(s.intervals[9].upper,s.intervals[7].upper,s.intervals[8].upper,signed_add,true))return false;
    return true;
  }
  IRRED_ELEMENTARY_INLINE void finish(Work &call) noexcept {
    for(unsigned i=0;i<9;++i)call.counters[i]=b.ledger.served.counters[i]-before.counters[i];
    finish_control(control);
  }
};
} // namespace

struct ElementaryImplementation {
  static IRRED_ELEMENTARY_INLINE bool context(Engine &e,const ElementaryLogContext &c) noexcept {
    e.stage=Stage::context;
    return e.require([&]{return c.status==Status::ok && c.log2_ &&
      c.profile_fingerprint_==fingerprint;},Status::missing_log_context);
  }
  static IRRED_ELEMENTARY_INLINE void prepare(ElementaryLogContext &c,ElementaryBudget &b,ElementaryScratch &s) noexcept {
    const bool unattempted=c.status==Status::unavailable && !c.log2_;
    c.log2_.reset();c.profile_fingerprint_=0; // Ends prior endpoints; no Wide assignment.
    Engine e(b,s,c.refusal,c.control_work,c.status,c.stage,c.preparation_work);
    if(e.begin() && e.require([&]{return unattempted;},Status::invalid_input)) {
      e.stage=Stage::rational;
      if(e.raw(s.scalars[0],one,three,divide) &&
         e.adjacent(s.intervals[0].lower,s.scalars[0],false) &&
         e.adjacent(s.intervals[0].upper,s.scalars[0],true) &&
         e.positive_binary(1,0,0,multiply) && e.series() &&
         e.charge(Work{{0,0,0,0,0,2,2,0,0}})) {
        c.log2_.emplace(s.intervals[7].lower,s.intervals[7].upper);
        c.profile_fingerprint_=fingerprint;c.status=Status::ok;c.stage=Stage::completed;
      }
    }
    e.finish(c.preparation_work);
  }
  static IRRED_ELEMENTARY_INLINE void bound(Engine &e,ElementaryResult &r,const W &lower,const W &upper,const W &error) noexcept {
    e.stage=Stage::image;
    if(!e.require([&]{return std::isfinite(lower)&&std::isfinite(upper)&&std::isfinite(error)&&
       lower<=upper&&error>=zero;}))return;
    if(!e.charge(Work{{0,0,0,0,0,3,3,0,0}}))return;
    r.bound.emplace(lower,upper,error);r.status=Status::ok;r.stage=Stage::completed;
  }
};

void prepare_elementary_log(ElementaryLogContext &c,ElementaryBudget &b,ElementaryScratch &s) noexcept {
  ElementaryImplementation::prepare(c,b,s);
}
void postcheck_log(const ElementaryLogContext &c,const W &x,const W &y,ElementaryBudget &b,
                  ElementaryScratch &s,ElementaryResult &r) noexcept {
  r.bound.reset();Engine e(b,s,r.refusal,r.control_work,r.status,r.stage,r.call_work);
  if(e.begin() && e.ownership(x,y,&r,sizeof(r),&c,sizeof(c)) &&
     ElementaryImplementation::context(e,c)) {
    e.stage=Stage::input;
    if(e.require([&]{return normal(y)||y==zero;},Status::invalid_input) && e.logarithm(*c.log2(),x)) {
      e.stage=Stage::discrepancy;
      bool left=false;
      if(e.rounded(s.scalars[5],y,s.intervals[9].lower,subtract,true) &&
         e.rounded(s.scalars[6],s.intervals[9].upper,y,subtract,true) &&
         e.decide(left,[&]{return s.scalars[5]>=s.scalars[6];}))
        ElementaryImplementation::bound(e,r,s.intervals[9].lower,s.intervals[9].upper,
                                       s.scalars[left?5:6]);
    }
  }
  e.finish(r.call_work);
}
void postcheck_sqrt(const W &x,const W &y,ElementaryBudget &b,ElementaryScratch &s,
                   ElementaryResult &r) noexcept {
  r.bound.reset();Engine e(b,s,r.refusal,r.control_work,r.status,r.stage,r.call_work);
  if(e.begin() && e.ownership(x,y,&r,sizeof(r))) {
    e.stage=Stage::input;
    bool exact=false,left=false;
    if(e.require([&]{return normal(x)&&x>zero&&normal(y)&&y>zero;},Status::invalid_input) &&
       e.decide(exact,[&]{return x==one&&y==one;})) {
      if(exact)ElementaryImplementation::bound(e,r,one,one,zero);
      else {
        e.stage=Stage::residual;
        if(e.rounded(s.intervals[10].lower,y,y,multiply,false) &&
           e.rounded(s.intervals[10].upper,y,y,multiply,true) &&
           e.rounded(s.intervals[11].lower,s.intervals[10].lower,x,subtract,false) &&
           e.rounded(s.intervals[11].upper,s.intervals[10].upper,x,subtract,true) &&
           e.negate(s.scalars[12],s.intervals[11].lower) &&
           e.decide(left,[&]{return s.scalars[12]>=s.intervals[11].upper;}) &&
           e.rounded(s.scalars[9],left?s.scalars[12]:s.intervals[11].upper,y,divide,true) &&
           e.require([&]{return s.scalars[9]<y;}) && e.raw(s.scalars[11],two,y,multiply)) {
          e.stage=Stage::denominator;
          if(e.rounded(s.scalars[8],s.scalars[11],s.scalars[9],subtract,false) &&
             e.require([&]{return s.scalars[8]>zero;}) &&
             e.rounded(s.scalars[10],left?s.scalars[12]:s.intervals[11].upper,s.scalars[8],divide,true) &&
             e.rounded(s.intervals[12].lower,y,s.scalars[10],subtract,false) &&
             e.rounded(s.intervals[12].upper,y,s.scalars[10],add,true) &&
             e.require([&]{return s.intervals[12].lower>zero;}))
            ElementaryImplementation::bound(e,r,s.intervals[12].lower,s.intervals[12].upper,s.scalars[10]);
        }
      }
    }
  }
  e.finish(r.call_work);
}
void postcheck_exp(const ElementaryLogContext &c,const W &n,const W &y,ElementaryBudget &b,
                  ElementaryScratch &s,ElementaryResult &r) noexcept {
  r.bound.reset();Engine e(b,s,r.refusal,r.control_work,r.status,r.stage,r.call_work);
  if(e.begin() && e.ownership(n,y,&r,sizeof(r),&c,sizeof(c)) &&
     ElementaryImplementation::context(e,c)) {
    e.stage=Stage::input;
    bool exact=false,left=false;
    if(e.require([&]{return normal(n)||n==zero;},Status::invalid_input) &&
       e.require([&]{return normal(y)&&y>zero&&y>=range_lower&&y<=range_upper;},Status::invalid_input) &&
       e.decide(exact,[&]{return n==zero&&y==one;})) {
      if(exact)ElementaryImplementation::bound(e,r,one,one,zero);
      else if(e.logarithm(*c.log2(),y)) {
        e.stage=Stage::discrepancy;
        if(e.rounded(s.scalars[5],s.intervals[9].upper,n,subtract,true) &&
           e.rounded(s.scalars[6],n,s.intervals[9].lower,subtract,true) &&
           e.decide(left,[&]{return s.scalars[5]>=s.scalars[6];})) {
          const unsigned d=left?5:6;e.stage=Stage::denominator;
          if(e.require([&]{return s.scalars[d]>=zero&&s.scalars[d]<one;}) &&
             e.rounded(s.scalars[8],one,s.scalars[d],subtract,false) &&
             e.require([&]{return s.scalars[8]>zero;}) &&
             e.rounded(s.scalars[9],s.scalars[d],s.scalars[8],divide,true) &&
             e.rounded(s.scalars[10],y,s.scalars[9],multiply,true) &&
             e.rounded(s.intervals[12].lower,y,s.scalars[10],subtract,false) &&
             e.rounded(s.intervals[12].upper,y,s.scalars[10],add,true) &&
             e.require([&]{return s.intervals[12].lower>zero;}))
            ElementaryImplementation::bound(e,r,s.intervals[12].lower,s.intervals[12].upper,s.scalars[10]);
        }
      }
    }
  }
  e.finish(r.call_work);
}
bool admit_elementary_candidate(ElementaryBudget &b,ElementaryRefusal &r,ElementaryControlWork &c) noexcept {
  r={};c={};
  if(admit(b,r,c,Stage::profile,Work{{0,0,0,0,0,0,0,1,0}},0)!=Admission::ok) {
    finish_control(c);return false;
  }
  ++c.predicates;
  if(!profile()) {
    r.occurred=true;r.stage=Stage::profile;finish_control(c);return false;
  }
  const bool accepted=admit(b,r,c,Stage::candidate_call,Work{{0,0,0,0,1,1,0,0,0}},1)==Admission::ok;
  finish_control(c);return accepted;
}
bool import_elementary_coordinate(const double &x,W &out,ElementaryBudget &b,
                                 ElementaryRefusal &r,ElementaryControlWork &c) noexcept {
  r={};c={};
  if(admit(b,r,c,Stage::profile,Work{{0,0,0,0,0,0,0,1,0}},0)!=Admission::ok) {
    finish_control(c);return false;
  }
  ++c.predicates;
  if(!profile()) {r.occurred=true;r.stage=Stage::profile;finish_control(c);return false;}
  if(admit(b,r,c,Stage::input,Work{{0,0,0,0,0,0,0,1,0}},0)!=Admission::ok) {
    finish_control(c);return false;
  }
  ++c.predicates;
  if(!(std::isnormal(x)||x==0.0)) {r.occurred=true;r.stage=Stage::input;finish_control(c);return false;}
  const bool accepted=admit(b,r,c,Stage::input,Work{{1,0,0,0,0,1,0,0,0}},1)==Admission::ok;
  if(accepted)out=static_cast<W>(x);
  finish_control(c);return accepted;
}
bool admit_elementary_output_guard(ElementaryBudget &b,ElementaryRefusal &r,
                                  ElementaryControlWork &c) noexcept {
  r={};c={};const bool accepted=admit(b,r,c,Stage::image,Work{{0,0,0,0,0,0,0,1,0}},0)==Admission::ok;
  if(accepted)++c.predicates;
  finish_control(c);return accepted;
}
} // namespace irred::cosmology::detail

#undef IRRED_ELEMENTARY_INLINE
