#pragma once
#include "../src/thermal_elementary_postcheck.hpp"
#include <array>
#include <cfenv>
#include <cstring>
#include <unistd.h>

namespace elementary_test {
namespace d=irred::cosmology::detail;
// Integer-only encoding; no allocating/floating formatter or by-value Wide.
struct Wire {
  std::array<char,1024> bytes;
  std::array<unsigned char,10> significant;
  std::size_t position=0;
  std::uint64_t io_calls=0,io_bytes=0,records=0;
  std::uint32_t byte_reads=0,byte_writes=0,integer_steps=0,branch_upper=0;
  std::uint32_t maximum_reads=0,maximum_writes=0,maximum_steps=0,maximum_branches=0;
  std::array<std::uint32_t,4> whole_control{};
  bool good=true;
  void add_control(const d::ElementaryControlWork &c) noexcept {
    const std::array<std::uint32_t,4> delta={c.transitions,c.predicates,c.iterations,c.integer_operation_upper};
    for(unsigned i=0;i<4;++i)if(delta[i]>UINT32_MAX-whole_control[i]){good=false;return;}
    for(unsigned i=0;i<4;++i)whole_control[i]+=delta[i];
  }
  bool reserve(unsigned reads,unsigned writes,unsigned steps,unsigned branches) noexcept {
    // Nonrecursive fixed encoder ledger; no floating work occurs on this path.
    // Bitwise boolean combination avoids unenumerated short-circuit decisions.
    if((!good)|(reads>256-byte_reads)|(writes>1024-byte_writes)|
       (steps>4096-integer_steps)|(branches>512-branch_upper)) {
      good=false;return false;
    }
    byte_reads+=reads;byte_writes+=writes;integer_steps+=steps;branch_upper+=branches;return true;
  }
  void character(char c) noexcept {
    if(!reserve(0,1,8,2))return;
    bytes[position++]=c;
  }
  template<std::size_t N> void text(const char (&p)[N]) noexcept {
    // Only compiled literal/schema text: length known before any byte access.
    if(!reserve(N-1,N-1,8*(N-1)+16,N+2))return;
    for(std::size_t i=0;i+1<N;++i)bytes[position++]=p[i];
  }
  void integer(std::uint64_t n) noexcept {
    // Finite integer-only length planning has <=20 iterations. Its operations
    // and decisions are included in the admitted block, not recursively billed.
    unsigned used=1;std::uint64_t divisor=1,remaining=n;
    while(remaining>=10){remaining/=10;divisor*=10;++used;}
    if(!reserve(0,used,20*used+8,2*used+3))return;
    for(unsigned i=0;i<used;++i){
      bytes[position++]=static_cast<char>('0'+n/divisor);n%=divisor;divisor/=10;
    }
  }
  void wide(const d::ElementaryWide &x) noexcept {
    if(!reserve(20,32,272,13))return;
    std::memcpy(significant.data(),&x,10);bytes[position++]='"';
    for(auto b:significant){
      const unsigned h=b>>4,l=b&15;
      bytes[position++]=static_cast<char>('0'+h+((h+6)>>4)*39);
      bytes[position++]=static_cast<char>('0'+l+((l+6)>>4)*39);
    }
    bytes[position++]='"';
  }
  void work(const d::ElementaryWork &w) noexcept {
    character('[');integer(w.counters[0]);
#define IRRED_WIRE_FIELD(i) character(',');integer(w.counters[i])
    IRRED_WIRE_FIELD(1);IRRED_WIRE_FIELD(2);IRRED_WIRE_FIELD(3);IRRED_WIRE_FIELD(4);
    IRRED_WIRE_FIELD(5);IRRED_WIRE_FIELD(6);IRRED_WIRE_FIELD(7);IRRED_WIRE_FIELD(8);
#undef IRRED_WIRE_FIELD
    character(']');
  }
  void control(const d::ElementaryControlWork &c) noexcept {
    character('[');integer(c.transitions);character(',');integer(c.predicates);
    character(',');integer(c.iterations);character(',');integer(c.integer_operation_upper);character(']');
  }
  void aggregate_control() noexcept {
    character('[');integer(whole_control[0]);character(',');integer(whole_control[1]);
    character(',');integer(whole_control[2]);character(',');integer(whole_control[3]);character(']');
  }
  bool flush(bool final_dominated=false) noexcept {
    // Fixed successful/failing leaf includes I/O disposition and four extrema.
    if(!reserve(0,1,64,16))return false;
    bytes[position++]='\n';
    if(final_dominated && ((byte_reads>maximum_reads)|(byte_writes>maximum_writes)|
       (integer_steps>maximum_steps)|(branch_upper>maximum_branches))){good=false;return false;}
    ++io_calls;
    const auto n=::write(STDOUT_FILENO,bytes.data(),position);
    if(n<0 || static_cast<std::size_t>(n)!=position){good=false;return false;}
    io_bytes+=position;++records;position=0;
    if(byte_reads>maximum_reads)maximum_reads=byte_reads;
    if(byte_writes>maximum_writes)maximum_writes=byte_writes;
    if(integer_steps>maximum_steps)maximum_steps=integer_steps;
    if(branch_upper>maximum_branches)maximum_branches=branch_upper;
    byte_reads=0;byte_writes=0;integer_steps=0;branch_upper=0;return true;
  }
  bool result(unsigned id,unsigned operation,const d::ElementaryWide &x,
              const d::ElementaryWide *candidate,const d::ElementaryResult &r,
              const d::ElementaryLedger &ledger) noexcept {
    character('[');integer(id);character(',');integer(operation);character(',');
    integer(static_cast<unsigned>(r.status));character(',');integer(static_cast<unsigned>(r.stage));
    character(',');wide(x);character(',');if(candidate)wide(*candidate);else text("null");
    character(',');
    if(r.bound){character('[');wide(r.bound->lower);character(',');wide(r.bound->upper);
      character(',');wide(r.bound->absolute_error);character(']');}else text("null");
    character(',');work(r.call_work);character(',');work(ledger.served);character(',');
    character('[');integer(r.refusal.occurred);character(',');integer(static_cast<unsigned>(r.refusal.stage));
    character(',');work(r.refusal.requested_increment);character(']');character(',');
    control(r.control_work);character(']');return flush();
  }
  bool refused(unsigned id,unsigned operation,const double &original,
               const d::ElementaryWide *x,const d::ElementaryWide *candidate,
               const d::ElementaryRefusal &r,const d::ElementaryControlWork &c,
               const d::ElementaryLedger &ledger) noexcept {
    text("[\"caller-refusal\",");integer(id);character(',');integer(operation);character(',');
    if(!reserve(16,26,224,11))return false;
    std::memcpy(significant.data(),&original,8);bytes[position++]='"';
    for(unsigned i=0;i<8;++i){
      const unsigned h=significant[i]>>4,l=significant[i]&15;
      bytes[position++]=static_cast<char>('0'+h+((h+6)>>4)*39);
      bytes[position++]=static_cast<char>('0'+l+((l+6)>>4)*39);}
    bytes[position++]='"';character(',');if(x)wide(*x);else text("null");
    character(',');if(candidate)wide(*candidate);else text("null");character(',');
    integer(static_cast<unsigned>(r.stage));character(',');work(r.requested_increment);
    character(',');work(ledger.served);character(',');control(c);character(']');return flush();
  }
  bool gate(unsigned id,unsigned operation,bool accepted,const d::ElementaryRefusal &r,
            const d::ElementaryControlWork &c,const d::ElementaryLedger &ledger) noexcept {
    text("[\"gate\",");integer(id);character(',');integer(operation);character(',');integer(accepted);
    character(',');integer(static_cast<unsigned>(r.stage));character(',');work(r.requested_increment);
    character(',');work(ledger.served);character(',');control(c);character(']');return flush();
  }
  bool request(unsigned id,unsigned kind,unsigned descriptor,const d::ElementaryBudget &b) noexcept {
    text("[\"request\",");integer(id);character(',');integer(kind);character(',');integer(descriptor);
    character(',');character('[');integer(b.maximum_scalar);character(',');integer(b.maximum_writes);
    character(',');integer(b.maximum_copies);character(',');integer(b.maximum_guards);character(']');
    character(',');work(b.ledger.served);character(',');integer(b.ledger.scalar_total);character(',');
    if((kind==8)|(kind==9)|(kind==10)|(kind==14)) {
      unsigned short x87=0;unsigned mxcsr=0;
#if defined(__x86_64__) && defined(__GNUC__)
      __asm__ volatile("fnstcw %0":"=m"(x87));__asm__ volatile("stmxcsr %0":"=m"(mxcsr));
#endif
      character('[');integer(static_cast<unsigned>(std::fegetround()));character(',');integer(x87);
      character(',');integer(mxcsr);character(']');
    } else text("null");
    character(']');return flush();
  }
  bool context(unsigned id,unsigned phase,const d::ElementaryLogContext &c,
               const d::ElementaryLedger &ledger) noexcept {
    text("[\"context\",");integer(id);character(',');integer(phase);character(',');
    integer(static_cast<unsigned>(c.status));character(',');integer(static_cast<unsigned>(c.stage));
    character(',');
    if(c.log2()){character('[');wide(c.log2()->lower);character(',');wide(c.log2()->upper);character(']');}
    else text("null");
    character(',');work(c.preparation_work);character(',');work(ledger.served);character(',');
    character('[');integer(c.refusal.occurred);character(',');integer(static_cast<unsigned>(c.refusal.stage));
    character(',');work(c.refusal.requested_increment);character(']');character(',');
    control(c.control_work);character(']');return flush();
  }
  bool caller(unsigned id,unsigned entry,bool accepted,const d::ElementaryWide *out,
              const d::ElementaryRefusal &r,const d::ElementaryControlWork &c,
              const d::ElementaryLedger &ledger) noexcept {
    text("[\"caller\",");integer(id);character(',');integer(entry);character(',');integer(accepted);
    character(',');if(out)wide(*out);else text("null");character(',');
    integer(static_cast<unsigned>(r.stage));character(',');work(r.requested_increment);
    character(',');work(ledger.served);character(',');control(c);character(']');return flush();
  }
  bool coordinate(unsigned id,const double &original) noexcept {
    text("[\"coordinate\",");integer(id);character(',');
    if(!reserve(16,26,224,11))return false;
    std::memcpy(significant.data(),&original,8);bytes[position++]='"';
    for(unsigned i=0;i<8;++i){
      const unsigned h=significant[i]>>4,l=significant[i]&15;
      bytes[position++]=static_cast<char>('0'+h+((h+6)>>4)*39);
      bytes[position++]=static_cast<char>('0'+l+((l+6)>>4)*39);}
    bytes[position++]='"';character(']');return flush();
  }
  bool source() noexcept {
    text("[\"control-source\",\"");
#ifdef IRRED_ELEMENTARY_SOURCE_SHA256
    text(IRRED_ELEMENTARY_SOURCE_SHA256);
#else
    text("unconfigured");
#endif
    text("\"]");return flush();
  }
};
static_assert(sizeof(Wire)<=1152);
} // namespace elementary_test
