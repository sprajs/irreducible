#pragma once
#include "../src/thermal_elementary_postcheck.hpp"
#include <array>
#include <cstring>
#include <unistd.h>

namespace elementary_test {
namespace d=irred::cosmology::detail;
// Integer-only encoding; no allocating/floating formatter or by-value Wide.
struct Wire {
  std::array<char,1024> bytes;
  std::array<char,20> digits;
  std::array<unsigned char,10> significant;
  std::size_t position=0;
  std::uint64_t io_calls=0,io_bytes=0,records=0;
  std::uint32_t byte_reads=0,byte_writes=0,integer_steps=0;
  std::uint32_t last_reads=0,last_writes=0,last_steps=0;
  std::array<std::uint32_t,4> whole_control{};
  bool good=true;
  void add_control(const d::ElementaryControlWork &c) noexcept {
    const std::array<std::uint32_t,4> delta={c.transitions,c.predicates,c.iterations,c.integer_operation_upper};
    for(unsigned i=0;i<4;++i)if(delta[i]>UINT32_MAX-whole_control[i]){good=false;return;}
    for(unsigned i=0;i<4;++i)whole_control[i]+=delta[i];
  }
  bool reserve(unsigned reads,unsigned writes,unsigned steps) noexcept {
    // Nonrecursive fixed encoder ledger; no floating work occurs on this path.
    if(reads>256-byte_reads || writes>1024-byte_writes || steps>4096-integer_steps) {
      good=false;return false;
    }
    byte_reads+=reads;byte_writes+=writes;integer_steps+=steps;return true;
  }
  void character(char c) noexcept {
    if(!good || !reserve(0,1,8))return;
    if(position==bytes.size()){good=false;return;}
    bytes[position++]=c;
  }
  void text(const char *p) noexcept {
    for(;;){if(!good || !reserve(1,0,4))return;if(!*p)break;character(*p++);}
  }
  void integer(std::uint64_t n) noexcept {
    unsigned used=0;
    do{if(!good || !reserve(0,1,8))return;
      digits[used++]=static_cast<char>('0'+n%10);n/=10;}while(n);
    while(used){if(!good || !reserve(1,0,4))return;character(digits[--used]);}
  }
  void wide(const d::ElementaryWide &x) noexcept {
    if(!good || !reserve(10,10,16))return;
    std::memcpy(significant.data(),&x,10);character('"');
    for(auto b:significant){
      if(!good || !reserve(1,0,8))return;
      const unsigned h=b>>4,l=b&15;
      character(static_cast<char>(h<10?'0'+h:'a'+h-10));
      character(static_cast<char>(l<10?'0'+l:'a'+l-10));
    }
    character('"');
  }
  void work(const d::ElementaryWork &w) noexcept {
    character('[');
    for(unsigned i=0;i<9;++i){if(i)character(',');integer(w.counters[i]);}
    character(']');
  }
  void control(const d::ElementaryControlWork &c) noexcept {
    character('[');integer(c.transitions);character(',');integer(c.predicates);
    character(',');integer(c.iterations);character(',');integer(c.integer_operation_upper);character(']');
  }
  void aggregate_control() noexcept {
    character('[');for(unsigned i=0;i<4;++i){if(i)character(',');integer(whole_control[i]);}character(']');
  }
  bool flush() noexcept {
    character('\n');if(!good)return false;
    ++io_calls;
    const auto n=::write(STDOUT_FILENO,bytes.data(),position);
    if(n<0 || static_cast<std::size_t>(n)!=position){good=false;return false;}
    io_bytes+=position;++records;position=0;
    last_reads=byte_reads;last_writes=byte_writes;last_steps=integer_steps;
    byte_reads=0;byte_writes=0;integer_steps=0;return true;
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
    if(!reserve(8,8,16))return false;
    std::memcpy(significant.data(),&original,8);character('"');
    for(unsigned i=0;i<8;++i){if(!reserve(1,0,8))return false;
      const unsigned h=significant[i]>>4,l=significant[i]&15;
      character(static_cast<char>(h<10?'0'+h:'a'+h-10));
      character(static_cast<char>(l<10?'0'+l:'a'+l-10));}
    character('"');character(',');if(x)wide(*x);else text("null");
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
};
static_assert(sizeof(Wire)<=1152);
} // namespace elementary_test
