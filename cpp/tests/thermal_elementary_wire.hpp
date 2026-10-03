#pragma once
#include "../src/thermal_elementary_postcheck.hpp"
#include <array>
#include <cfenv>
#include <cstring>
#include <unistd.h>

namespace elementary_test {
namespace d=irred::cosmology::detail;
// Integer-only encoding; no allocating/floating formatter or by-value Wide.
struct RecordPlan {
  std::array<unsigned char,37> lengths; // Unassigned until the corresponding plan action.
  unsigned fields=0,digits=0,characters=0,wide_values=0,text_bytes=0,text_blocks=0,coordinates=0;
  RecordPlan() noexcept {} // No brace/value initialization of length bytes.
  void integer(const std::uint64_t &n) noexcept {
    std::uint64_t remaining=n;unsigned used=0;
    do{++used;remaining/=10;}while(remaining);
    lengths[fields++]=static_cast<unsigned char>(used);digits+=used;
  }
  void character(char) noexcept {++characters;}
  template<std::size_t N> void text(const char (&)[N]) noexcept {text_bytes+=N-1;++text_blocks;}
  void wide(const d::ElementaryWide &) noexcept {++wide_values;}
  void coordinate_value(const double &) noexcept {++coordinates;}
  void work(const d::ElementaryWork &w) noexcept {
    character('[');integer(w.counters[0]);
#define IRRED_PLAN_FIELD(i) character(',');integer(w.counters[i])
    IRRED_PLAN_FIELD(1);IRRED_PLAN_FIELD(2);IRRED_PLAN_FIELD(3);IRRED_PLAN_FIELD(4);
    IRRED_PLAN_FIELD(5);IRRED_PLAN_FIELD(6);IRRED_PLAN_FIELD(7);IRRED_PLAN_FIELD(8);
#undef IRRED_PLAN_FIELD
    character(']');
  }
  void control(const d::ElementaryControlWork &c) noexcept {
    character('[');integer(c.transitions);character(',');integer(c.predicates);
    character(',');integer(c.iterations);character(',');integer(c.integer_operation_upper);character(']');
  }
  void aggregate_control(const std::array<std::uint32_t,4> &values) noexcept {
    character('[');integer(values[0]);character(',');integer(values[1]);
    character(',');integer(values[2]);character(',');integer(values[3]);character(']');
  }
};
static_assert(sizeof(RecordPlan)<=72);
struct Wire {
  std::array<char,1024> bytes;
  std::array<unsigned char,10> significant;
  std::size_t position=0;
  std::uint64_t io_calls=0,io_bytes=0,records=0;
  std::uint32_t byte_reads=0,byte_writes=0,integer_steps=0,branch_upper=0;
  std::uint32_t maximum_reads=0,maximum_writes=0,maximum_steps=0,maximum_branches=0;
  std::array<std::uint32_t,4> whole_control{};
  bool good=true;
  const RecordPlan *plan=nullptr;unsigned integer_index=0;
  void add_control(const d::ElementaryControlWork &c) noexcept {
    const std::array<std::uint32_t,4> delta={c.transitions,c.predicates,c.iterations,c.integer_operation_upper};
    for(unsigned i=0;i<4;++i)if(delta[i]>UINT32_MAX-whole_control[i]){good=false;return;}
    for(unsigned i=0;i<4;++i)whole_control[i]+=delta[i];
  }
  // Finite compiled schemas ONLY: <=37 integer fields and <=128 aggregate
  // decimal digits are source-derived from this exact fixture/control corpus.
  // Planner, admission, encoder and flush are one nonrecursive resource owner.
  template<class Renderer> bool record(const Renderer &render,bool final_dominated=false) noexcept {
    // This fixed source-corpus planning prefix includes its admission/disposition
    // and <=37 scratch-byte writes. It is not an arbitrary-renderer promise.
    constexpr unsigned planning_prefix_upper=3072;
    static_assert(planning_prefix_upper<=4096);
    if((!good)|(position!=0)|(planning_prefix_upper>4096)){good=false;return false;}
    RecordPlan p;render(p);
    if((!good)|(position!=0)|(p.fields>37)|(p.digits>128)) {good=false;return false;}
    byte_reads=20*p.wide_values+16*p.coordinates+p.text_bytes+p.fields;
    byte_writes=32*p.wide_values+26*p.coordinates+p.text_bytes+p.characters+p.digits+p.fields+1;
    integer_steps=16*p.digits+17*p.fields+288*p.wide_values+8*p.characters+
      8*p.text_bytes+16*p.text_blocks+240*p.coordinates+192;
    branch_upper=2*p.digits+4*p.fields+12*p.wide_values+2*p.text_bytes+2*p.text_blocks+10*p.coordinates+64;
    if((byte_reads>256)|(byte_writes>1024)|(integer_steps>4096)|(branch_upper>512)) {
      good=false;return false;
    }
    plan=&p;integer_index=0;render(*this);
    if(integer_index!=p.fields){good=false;return false;}
    return flush(final_dominated);
  }
  void character(char c) noexcept {bytes[position++]=c;}
  template<std::size_t N> void text(const char (&p)[N]) noexcept {
    for(std::size_t i=0;i+1<N;++i)bytes[position++]=p[i];
  }
  void integer(std::uint64_t n) noexcept {
    const unsigned used=plan->lengths[integer_index++];
    const auto end=position+used;position=end;
    for(unsigned i=0;i<used;++i){
      bytes[end-1-i]=static_cast<char>('0'+n%10);n/=10;
    }
  }
  void wide(const d::ElementaryWide &x) noexcept {
    std::memcpy(significant.data(),&x,10);bytes[position++]='"';
    for(auto b:significant){
      const unsigned h=b>>4,l=b&15;
      bytes[position++]=static_cast<char>('0'+h+((h+6)>>4)*39);
      bytes[position++]=static_cast<char>('0'+l+((l+6)>>4)*39);
    }
    bytes[position++]='"';
  }
  void coordinate_value(const double &original) noexcept {
    std::memcpy(significant.data(),&original,8);bytes[position++]='"';
    for(unsigned i=0;i<8;++i){
      const unsigned h=significant[i]>>4,l=significant[i]&15;
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
  void aggregate_control(const std::array<std::uint32_t,4> &values) noexcept {
    character('[');integer(values[0]);character(',');integer(values[1]);
    character(',');integer(values[2]);character(',');integer(values[3]);character(']');
  }
  bool flush(bool final_dominated=false) noexcept {
    // Whole-record admission owns this complete disposition and its extrema.
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
    byte_reads=0;byte_writes=0;integer_steps=0;branch_upper=0;plan=nullptr;integer_index=0;return true;
  }
  bool result(unsigned id,unsigned operation,const d::ElementaryWide &x,
              const d::ElementaryWide *candidate,const d::ElementaryResult &r,
              const d::ElementaryLedger &ledger) noexcept {
    return record([&](auto &sink){
    sink.character('[');sink.integer(id);sink.character(',');sink.integer(operation);sink.character(',');
    sink.integer(static_cast<unsigned>(r.status));sink.character(',');sink.integer(static_cast<unsigned>(r.stage));
    sink.character(',');sink.wide(x);sink.character(',');if(candidate)sink.wide(*candidate);else sink.text("null");
    sink.character(',');
    if(r.bound){sink.character('[');sink.wide(r.bound->lower);sink.character(',');sink.wide(r.bound->upper);
      sink.character(',');sink.wide(r.bound->absolute_error);sink.character(']');}else sink.text("null");
    sink.character(',');sink.work(r.call_work);sink.character(',');sink.work(ledger.served);sink.character(',');
    sink.character('[');sink.integer(r.refusal.occurred);sink.character(',');sink.integer(static_cast<unsigned>(r.refusal.stage));
    sink.character(',');sink.work(r.refusal.requested_increment);sink.character(']');sink.character(',');
    sink.control(r.control_work);sink.character(']');

    });
  }
  bool refused(unsigned id,unsigned operation,const double &original,
               const d::ElementaryWide *x,const d::ElementaryWide *candidate,
               const d::ElementaryRefusal &r,const d::ElementaryControlWork &c,
               const d::ElementaryLedger &ledger) noexcept {
    return record([&](auto &sink){
    sink.text("[\"caller-refusal\",");sink.integer(id);sink.character(',');sink.integer(operation);sink.character(',');

    sink.coordinate_value(original);sink.character(',');if(x)sink.wide(*x);else sink.text("null");
    sink.character(',');if(candidate)sink.wide(*candidate);else sink.text("null");sink.character(',');
    sink.integer(static_cast<unsigned>(r.stage));sink.character(',');sink.work(r.requested_increment);
    sink.character(',');sink.work(ledger.served);sink.character(',');sink.control(c);sink.character(']');

    });
  }
  bool gate(unsigned id,unsigned operation,bool accepted,const d::ElementaryRefusal &r,
            const d::ElementaryControlWork &c,const d::ElementaryLedger &ledger) noexcept {
    return record([&](auto &sink){
    sink.text("[\"gate\",");sink.integer(id);sink.character(',');sink.integer(operation);sink.character(',');sink.integer(accepted);
    sink.character(',');sink.integer(static_cast<unsigned>(r.stage));sink.character(',');sink.work(r.requested_increment);
    sink.character(',');sink.work(ledger.served);sink.character(',');sink.control(c);sink.character(']');

    });
  }
  bool request(unsigned id,unsigned kind,unsigned descriptor,const d::ElementaryBudget &b) noexcept {
    unsigned short x87=0;unsigned mxcsr=0;const auto rounding=static_cast<unsigned>(std::fegetround());
#if defined(__x86_64__) && defined(__GNUC__)
    __asm__ volatile("fnstcw %0":"=m"(x87));__asm__ volatile("stmxcsr %0":"=m"(mxcsr));
#endif

    return record([&](auto &sink){
    sink.text("[\"request\",");sink.integer(id);sink.character(',');sink.integer(kind);sink.character(',');sink.integer(descriptor);
    sink.character(',');sink.character('[');sink.integer(b.maximum_scalar);sink.character(',');sink.integer(b.maximum_writes);
    sink.character(',');sink.integer(b.maximum_copies);sink.character(',');sink.integer(b.maximum_guards);sink.character(']');
    sink.character(',');sink.work(b.ledger.served);sink.character(',');sink.integer(b.ledger.scalar_total);sink.character(',');
    if((kind==8)|(kind==9)|(kind==10)|(kind==14)) {

      sink.character('[');sink.integer(rounding);sink.character(',');sink.integer(x87);
      sink.character(',');sink.integer(mxcsr);sink.character(']');
    } else sink.text("null");
    sink.character(']');

    });
  }
  bool context(unsigned id,unsigned phase,const d::ElementaryLogContext &c,
               const d::ElementaryLedger &ledger) noexcept {
    return record([&](auto &sink){
    sink.text("[\"context\",");sink.integer(id);sink.character(',');sink.integer(phase);sink.character(',');
    sink.integer(static_cast<unsigned>(c.status));sink.character(',');sink.integer(static_cast<unsigned>(c.stage));
    sink.character(',');
    if(c.log2()){sink.character('[');sink.wide(c.log2()->lower);sink.character(',');sink.wide(c.log2()->upper);sink.character(']');}
    else sink.text("null");
    sink.character(',');sink.work(c.preparation_work);sink.character(',');sink.work(ledger.served);sink.character(',');
    sink.character('[');sink.integer(c.refusal.occurred);sink.character(',');sink.integer(static_cast<unsigned>(c.refusal.stage));
    sink.character(',');sink.work(c.refusal.requested_increment);sink.character(']');sink.character(',');
    sink.control(c.control_work);sink.character(']');

    });
  }
  bool caller(unsigned id,unsigned entry,bool accepted,const d::ElementaryWide *out,
              const d::ElementaryRefusal &r,const d::ElementaryControlWork &c,
              const d::ElementaryLedger &ledger) noexcept {
    return record([&](auto &sink){
    sink.text("[\"caller\",");sink.integer(id);sink.character(',');sink.integer(entry);sink.character(',');sink.integer(accepted);
    sink.character(',');if(out)sink.wide(*out);else sink.text("null");sink.character(',');
    sink.integer(static_cast<unsigned>(r.stage));sink.character(',');sink.work(r.requested_increment);
    sink.character(',');sink.work(ledger.served);sink.character(',');sink.control(c);sink.character(']');

    });
  }
  bool coordinate(unsigned id,const double &original) noexcept {
    return record([&](auto &sink){
    sink.text("[\"coordinate\",");sink.integer(id);sink.character(',');

    sink.coordinate_value(original);sink.character(']');

    });
  }
  bool source() noexcept {
    return record([&](auto &sink){
    sink.text("[\"control-source\",\"");
#ifdef IRRED_ELEMENTARY_SOURCE_SHA256
    sink.text(IRRED_ELEMENTARY_SOURCE_SHA256);
#else
    sink.text("unconfigured");
#endif
    sink.text("\"]");

    });
  }
};
static_assert(sizeof(Wire)<=1152);
} // namespace elementary_test
