#include "thermal_elementary_wire.hpp"
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <type_traits>

namespace {
namespace d=irred::cosmology::detail;
using W=d::ElementaryWide;
using Status=d::ElementaryStatus;
bool allocation_active=false;
std::uint64_t new_calls=0,malloc_calls=0;
constexpr W zero=0.0L,one=1.0L,four=4.0L,sixteen=16.0L;
constexpr double clock0=1e-10,clock1=.01;
unsigned failures=0;
unsigned requests_started=0;
void check(bool yes,const char *what) {
  if(!yes){++failures;std::fprintf(stderr,"thermal elementary: %s\n",what);}
}
struct Request {
  d::ElementaryLogContext context;
  d::ElementaryScratch scratch;
  d::ElementaryResult result;
  d::ElementaryLedger ledger;
  d::ElementaryBudget budget;
  W argument,log_candidate,exp_candidate,root_candidate; // Unassigned until charged.
  explicit Request(std::uint64_t scalar=8192,std::uint64_t writes=16384,
                   std::uint64_t copies=16384,std::uint64_t guards=65536)
    :budget{scalar,writes,copies,guards,ledger} {++requests_started;}
};
static_assert(!std::is_copy_constructible_v<d::ElementaryLogContext>);
static_assert(!std::is_move_constructible_v<d::ElementaryLogContext>);
static_assert(!std::is_copy_constructible_v<d::ElementaryResult>);
static_assert(!std::is_move_constructible_v<d::ElementaryResult>);
static_assert(sizeof(d::ElementaryLogContext)<=240 && sizeof(d::ElementaryScratch)<=832);
static_assert(sizeof(d::ElementaryResult)<=256 && sizeof(d::ElementaryLedger)<=80);
static_assert(sizeof(d::ElementaryBudget)<=40);
// All-active frames and opaque primitive frames are separate compiled/profile gates.
constexpr std::size_t envelope=sizeof(Request)+sizeof(elementary_test::Wire)+256+512+176+384;
static_assert(envelope<=4096);

bool available(Request &q,d::ElementaryRefusal &refusal,d::ElementaryControlWork &control) {
  if(!d::admit_elementary_output_guard(q.budget,refusal,control))return false;
  const bool valid=q.result.status==Status::ok && q.result.bound &&
    std::isfinite(q.result.bound->lower) && std::isfinite(q.result.bound->upper) &&
    std::isfinite(q.result.bound->absolute_error) &&
    q.result.bound->lower<=q.result.bound->upper && q.result.bound->absolute_error>=zero;
  if(!valid){refusal.occurred=true;refusal.stage=d::ElementaryStage::image;refusal.requested_increment={};}
  return valid;
}
bool fixture(elementary_test::Wire &wire) {
  Request q;
  auto &refusal=q.result.refusal;auto &control=q.result.control_work;
  wire.text("[\"ieee80-le\",[\"basic\",\"next\",\"scale\",\"term\",\"lib\",\"write\",\"copy\",\"guard\",\"attempt\"],\"clock150249\",\"");
#ifdef IRRED_ELEMENTARY_SOURCE_SHA256
  wire.text(IRRED_ELEMENTARY_SOURCE_SHA256);
#else
  wire.text("unconfigured");
#endif
  wire.text("\"]");if(!wire.flush())return false;
  allocation_active=true;
  d::prepare_elementary_log(q.context,q.budget,q.scratch);
  wire.add_control(q.context.control_work);
  allocation_active=false;
  wire.text("[\"prep\",");wire.integer(static_cast<unsigned>(q.context.status));wire.character(',');
  wire.work(q.context.preparation_work);wire.character(',');wire.control(q.context.control_work);
  wire.character(']');if(!wire.flush())return false;
  if(q.context.status!=Status::ok)return false;
  for(unsigned row=0;row<2;++row) {
    allocation_active=true;
    const bool imported=d::import_elementary_coordinate(row?clock1:clock0,q.argument,q.budget,refusal,control);
    wire.add_control(control);
    if(!imported){allocation_active=false;
      wire.refused(row,0,row?clock1:clock0,nullptr,nullptr,refusal,control,q.ledger);return false;}
    bool called=d::admit_elementary_candidate(q.budget,refusal,control);
    wire.add_control(control);
    bool admitted=false;
    if(!called){allocation_active=false;
      wire.refused(row,0,row?clock1:clock0,&q.argument,nullptr,refusal,control,q.ledger);return false;}
    if(called) {
      q.log_candidate=std::log(q.argument);
      d::postcheck_log(q.context,q.argument,q.log_candidate,q.budget,q.scratch,q.result);
      wire.add_control(q.result.control_work);
    }
    allocation_active=false;
    if(!wire.result(row,0,q.argument,&q.log_candidate,q.result,q.ledger))return false;
    admitted=available(q,refusal,control);wire.add_control(control);
    if(!wire.gate(row,0,admitted,refusal,control,q.ledger) || !admitted)return false;
    allocation_active=true;
    called=d::admit_elementary_candidate(q.budget,refusal,control);
    wire.add_control(control);
    if(!called){allocation_active=false;
      wire.refused(row,1,row?clock1:clock0,&q.log_candidate,nullptr,refusal,control,q.ledger);return false;}
    if(called) {
      q.exp_candidate=std::exp(q.log_candidate);
      d::postcheck_exp(q.context,q.log_candidate,q.exp_candidate,q.budget,q.scratch,q.result);
      wire.add_control(q.result.control_work);
    }
    allocation_active=false;
    if(!wire.result(row,1,q.log_candidate,&q.exp_candidate,q.result,q.ledger))return false;
    admitted=available(q,refusal,control);wire.add_control(control);
    if(!wire.gate(row,1,admitted,refusal,control,q.ledger) || !admitted)return false;
    allocation_active=true;
    called=d::admit_elementary_candidate(q.budget,refusal,control);
    wire.add_control(control);
    if(!called){allocation_active=false;
      wire.refused(row,2,row?clock1:clock0,row?&sixteen:&four,nullptr,refusal,control,q.ledger);return false;}
    if(called) {
      q.root_candidate=std::sqrt(row?sixteen:four);
      d::postcheck_sqrt(row?sixteen:four,q.root_candidate,q.budget,q.scratch,q.result);
      wire.add_control(q.result.control_work);
    }
    allocation_active=false;
    if(!wire.result(row,2,row?sixteen:four,&q.root_candidate,q.result,q.ledger))return false;
    admitted=available(q,refusal,control);wire.add_control(control);
    if(!wire.gate(row,2,admitted,refusal,control,q.ledger) || !admitted)return false;
  }
  check(q.ledger.served.counters[4]==6 && q.ledger.served.counters[8]==7,"fixed original calls/attempts");
  check(new_calls==0 && malloc_calls==0,"native owning fixture allocations");
  wire.text("[\"summary\",");wire.integer(envelope);wire.character(',');wire.work(q.ledger.served);
  wire.character(',');wire.integer(q.ledger.scalar_total);wire.character(',');wire.integer(new_calls);
  wire.character(',');wire.integer(malloc_calls);wire.character(',');wire.integer(wire.io_calls+1);
  wire.character(',');wire.aggregate_control();
  wire.character(',');wire.integer(sizeof(Request));wire.character(',');wire.integer(sizeof(elementary_test::Wire));
  wire.text("]");return wire.flush(); // Completion makes this last one-call count actual.
}
#ifndef IRRED_ELEMENTARY_FIXTURE_ONLY
enum Expect { accepted,refused,either };
struct Case { unsigned op;W x,y;bool library;Expect expectation; };
constexpr Case cases[]={
  {0,1,0,false,accepted},{1,0,1,false,accepted},{2,1,1,false,accepted},
  {0,1,1,false,accepted},{2,1,2,false,accepted},{1,0,1.125L,false,accepted},
  {1,0,2,false,refused},{2,4,0,true,accepted},{2,16,0,true,accepted},
  {2,0x1p128L,0,true,accepted},
  {0,0x1p-128L,0,true,accepted},{0,0x1p-64L,0,true,accepted},
  {0,.5L,0,true,accepted},{0,1,0,true,accepted},{0,2,0,true,accepted},
  {0,0x1p64L,0,true,accepted},{0,0x1p128L,0,true,accepted},
  {0,1.0L-0x1p-64L,0,true,accepted},{0,1.0L+0x1p-63L,0,true,accepted},
  {2,1.0L+0x1p-62L,1.0L+0x1p-63L,false,accepted},
  {0,0x1p-128L-0x1p-192L,0,false,refused},
  {0,0x1p128L+0x1p65L,0,false,refused},
  {0,0,0,false,refused},{0,-1,0,false,refused},{2,0,1,false,refused},
  {2,1,0,false,refused},{1,0,0,false,refused},{1,-1,0,true,accepted},
  {0,2,-1,false,accepted},{0,2,0,false,accepted},
  {2,std::numeric_limits<W>::min(),1,false,refused},
  {2,1,std::numeric_limits<W>::max(),false,refused},
  {0,std::numeric_limits<W>::denorm_min(),0,false,refused},
  {0,std::numeric_limits<W>::infinity(),0,false,refused},
  {2,1,std::numeric_limits<W>::quiet_NaN(),false,refused},
  {1,-std::numeric_limits<W>::infinity(),1,false,refused},
  {0,1,-0.0L,false,accepted},{1,-0.0L,1,false,accepted},
};
constexpr unsigned case_count=sizeof(cases)/sizeof(cases[0]);
bool run_case(unsigned id,elementary_test::Wire *wire=nullptr) {
  Request q;
  d::prepare_elementary_log(q.context,q.budget,q.scratch);
  if(q.context.status!=Status::ok)return false;
  const Case &c=cases[id];
  const W *candidate=&c.y;
  if(c.library) {
    d::ElementaryRefusal refusal;d::ElementaryControlWork control;
    if(!d::admit_elementary_candidate(q.budget,refusal,control))return false;
    if(c.op==0)q.log_candidate=std::log(c.x);
    else if(c.op==1)q.log_candidate=std::exp(c.x);
    else q.log_candidate=std::sqrt(c.x);
    candidate=&q.log_candidate;
  }
  if(c.op==0)d::postcheck_log(q.context,c.x,*candidate,q.budget,q.scratch,q.result);
  else if(c.op==1)d::postcheck_exp(q.context,c.x,*candidate,q.budget,q.scratch,q.result);
  else d::postcheck_sqrt(c.x,*candidate,q.budget,q.scratch,q.result);
  const bool ok=q.result.status==Status::ok && q.result.bound;
  check(c.expectation==either || ok==(c.expectation==accepted),"named case availability");
  if(id<=2 || id==36 || id==37)check(ok && q.result.bound->absolute_error==zero,"exact identity radius");
  if(id==3 || id==4 || id==19)check(ok && q.result.bound->absolute_error>zero,"wrong candidate/rounded-square nonzero radius");
  if(!ok)check(!q.result.bound && q.result.refusal.occurred,"refused case absence/causal receipt");
  return !wire || wire->result(id,c.op,c.x,candidate,q.result,q.ledger);
}
void prefixes() {
  Request original;
  d::prepare_elementary_log(original.context,original.budget,original.scratch);
  const auto served=original.ledger.served;
  const auto scalar=original.ledger.scalar_total;
  for(unsigned dimension=0;dimension<4;++dimension) {
    for(unsigned boundary=0;boundary<3;++boundary) {
      const std::uint64_t required=dimension==0?scalar:served.counters[dimension==1?5:dimension==2?6:7];
      const std::uint64_t cap=boundary==0?0:boundary==1?required-1:required;
      Request q(dimension==0?cap:8192,dimension==1?cap:16384,
                dimension==2?cap:16384,dimension==3?cap:65536);
      d::prepare_elementary_log(q.context,q.budget,q.scratch);
      check((q.context.status==Status::ok)==(boundary==2),"exact context cap boundary");
      if(boundary!=2) {
        check(q.context.status==Status::work_limit && q.context.refusal.occurred,"cap causal refusal");
        check(q.ledger.scalar_total<=q.budget.maximum_scalar &&
          q.ledger.served.counters[5]<=q.budget.maximum_writes &&
          q.ledger.served.counters[6]<=q.budget.maximum_copies &&
          q.ledger.served.counters[7]<=q.budget.maximum_guards,"served cap prefix");
        d::ElementaryWork expected;
        if(dimension==3)expected.counters[7]=1;
        else if(dimension==2 || (boundary==1 && dimension==1)) {
          expected.counters[5]=2;expected.counters[6]=2;
        } else {
          expected.counters[boundary==1?1:0]=1;expected.counters[5]=1;
        }
        check(q.context.refusal.requested_increment.counters==expected.counters,
              "exact unserved fixed source increment");
      }
    }
  }
  for(unsigned category:{3u,8u}) {
    Request q;q.ledger.served.counters[category]=UINT64_MAX;
    d::prepare_elementary_log(q.context,q.budget,q.scratch);
    check(q.context.status==Status::counter_overflow,"counter overflow");
    check(q.ledger.served.counters[category]==UINT64_MAX,"overflow category unchanged");
  }
  Request bad;bad.ledger.scalar_total=1;
  d::prepare_elementary_log(bad.context,bad.budget,bad.scratch);
  check(bad.context.status==Status::invalid_input,"cached scalar consistency");
  Request missing;
  d::postcheck_log(missing.context,one,zero,missing.budget,missing.scratch,missing.result);
  check(missing.result.status==Status::missing_log_context,"missing context causal refusal");
  Request aliased;
  // Binding an uninitialized object by reference is valid; admission must reject
  // its ownership before any lvalue-to-rvalue read of that object.
  d::postcheck_sqrt(aliased.scratch.scalars[14],one,aliased.budget,aliased.scratch,aliased.result);
  check(aliased.result.status==Status::invalid_ownership && !aliased.result.bound,
        "raw scratch alias refused before read");
}
void environments() {
  const int saved=std::fegetround();
  for(int rounding:{FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}) {
    std::fesetround(rounding);Request q;
    d::postcheck_sqrt(one,one,q.budget,q.scratch,q.result);
    check(q.result.status==Status::unsupported_arithmetic && !q.result.bound,"rounding profile refusal");
  }
  std::fesetround(saved);
#if defined(__x86_64__) && defined(__GNUC__)
  unsigned short x87=0;unsigned mxcsr=0;
  __asm__ volatile("fnstcw %0":"=m"(x87));__asm__ volatile("stmxcsr %0":"=m"(mxcsr));
  const unsigned short narrowed=static_cast<unsigned short>((x87&~0x0300u)|0x0200u);
  __asm__ volatile("fldcw %0"::"m"(narrowed));
  {Request q;d::postcheck_sqrt(one,one,q.budget,q.scratch,q.result);
   check(q.result.status==Status::unsupported_arithmetic,"x87 precision refusal");}
  __asm__ volatile("fldcw %0"::"m"(x87));
  for(unsigned mask:{0x40u,0x8000u}) {
    const unsigned changed=mxcsr|mask;__asm__ volatile("ldmxcsr %0"::"m"(changed));
    Request q;d::postcheck_sqrt(one,one,q.budget,q.scratch,q.result);
    check(q.result.status==Status::unsupported_arithmetic,"FTZ/DAZ refusal");
  }
  __asm__ volatile("ldmxcsr %0"::"m"(mxcsr));
#endif
}
#endif
} // namespace
void *operator new(std::size_t n) {
  if(allocation_active)++new_calls;
  if(void *p=std::malloc(n))return p;throw std::bad_alloc();
}
void *operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void *p) noexcept {std::free(p);}
void operator delete[](void *p) noexcept {std::free(p);}
void operator delete(void *p,std::size_t) noexcept {std::free(p);}
void operator delete[](void *p,std::size_t) noexcept {std::free(p);}
#ifdef IRRED_ELEMENTARY_WRAP_ALLOCATIONS
extern "C" void *__real_malloc(std::size_t);
extern "C" void *__real_calloc(std::size_t,std::size_t);
extern "C" void *__real_realloc(void *,std::size_t);
extern "C" void __real_free(void *);
extern "C" void *__wrap_malloc(std::size_t n){if(allocation_active)++malloc_calls;return __real_malloc(n);}
extern "C" void *__wrap_calloc(std::size_t n,std::size_t s){if(allocation_active)++malloc_calls;return __real_calloc(n,s);}
extern "C" void *__wrap_realloc(void *p,std::size_t n){if(allocation_active)++malloc_calls;return __real_realloc(p,n);}
extern "C" void __wrap_free(void *p){__real_free(p);}
#endif
int main(int argc,char **argv) {
  elementary_test::Wire wire;
#ifdef IRRED_ELEMENTARY_SOURCE_SHA256
  if(argc==3 && std::strcmp(argv[1],"--source-sha")==0)
    return std::strlen(argv[2])==64 && std::strcmp(argv[2],IRRED_ELEMENTARY_SOURCE_SHA256)==0?0:2;
#endif
#ifndef IRRED_ELEMENTARY_FIXTURE_ONLY
  if(argc==3 && std::strcmp(argv[1],"--case")==0) {
    unsigned id=0;
    if(!*argv[2])return 2;
    for(const char *p=argv[2];*p;++p){if(*p<'0'||*p>'9'||id>=4)return 2;id=10*id+static_cast<unsigned>(*p-'0');}
    if(id>=case_count)return 2;
    return run_case(id,&wire)&&!failures?0:1;
  }
  if(argc>2 || (argc==2 && std::strcmp(argv[1],"--facts")!=0))return 2;
  if(!fixture(wire))return 1;
  if(argc==2 && std::strcmp(argv[1],"--facts")==0)return failures?1:0;
  for(unsigned i=0;i<case_count;++i)run_case(i);
  prefixes();environments();
  check(case_count==38 && requests_started==63,"complete mandatory campaign inventory");
#else
  if(argc>2 || (argc==2 && std::strcmp(argv[1],"--facts")!=0))return 2;
  if(!fixture(wire))return 1;
  check(requests_started==1,"original fixture request identity");
#endif
  return failures?1:0;
}
