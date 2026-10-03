#include "irred/windowed_linear_power.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

// Independent test-only Legendre quadrature/table search/bisection. Shared
// working-law ancestry and libm do not constitute independent physical proof.
namespace {
using W=long double;
namespace p=irred::windowed_linear_power;
constexpr std::size_t payload=16*1024*1024;
constexpr W eps=std::numeric_limits<W>::epsilon();
struct Counters {
  std::size_t angular=0,queries=0,roots=0,newton=0,recurrences=0;
  bool probe() { if (angular>=1000000) return false; ++angular; return true; }
};
struct Rule { std::vector<W> nodes,weights; };
std::pair<W,W> legendre(unsigned n,W x,Counters& c) {
  ++c.recurrences;
  W previous=1,current=x;
  for (unsigned j=2;j<=n;++j) {
    const W next=((2*j-1)*x*current-(j-1)*previous)/j;
    previous=current;current=next;
  }
  return {current,n*(x*current-previous)/(x*x-1)};
}
bool make_rule(unsigned n,Rule& rule,Counters& c) {
  rule.nodes.resize(n);rule.weights.resize(n);
  const W pi=std::acos(-1.L);
  for (unsigned j=0;j<n/2;++j) {
    if (++c.roots>56) return false;
    W x=std::cos(pi*(j+.75L)/(n+.5L));
    bool accepted=false;
    const auto before=c.recurrences;
    for (unsigned iteration=0;iteration<64;++iteration) {
      ++c.newton;
      const auto [value,derivative]=legendre(n,x,c);
      if (!std::isfinite(value) || !std::isfinite(derivative) || derivative==0) return false;
      const W step=value/derivative,next=x-step;
      if (!(next>0 && next<1) || !std::isfinite(next)) return false;
      x=next;
      const auto [residual,d]=legendre(n,x,c);
      if (std::abs(step)<=64*eps && std::abs(residual)<=256*n*eps) {
        const W weight=2/((1-x*x)*d*d);
        if (!(weight>0) || !std::isnormal(weight)) return false;
        rule.nodes[j]=-x;rule.nodes[n-1-j]=x;
        rule.weights[j]=rule.weights[n-1-j]=weight;
        accepted=true;break;
      }
    }
    if (!accepted || c.recurrences-before>130) return false;
  }
  for (unsigned j=0;j<n;++j)
    if (!std::isfinite(rule.nodes[j]) || rule.nodes[j]<=-1 || rule.nodes[j]>=1 ||
        (j && rule.nodes[j]<=rule.nodes[j-1])) return false;
  // Every even polynomial moment required by a degree(2n-1) rule.
  for (unsigned power=0;power<2*n;power+=2) {
    W moment=0;
    for (unsigned j=0;j<n;++j) {
      W v=1;
      for (unsigned exponent=0;exponent<power;++exponent) v*=rule.nodes[j];
      moment+=rule.weights[j]*v;
    }
    if (std::abs(moment-2.L/(power+1))>1e-13L) return false;
  }
  return true;
}
W mapped_k(W mu,W k,const p::ModelPoint& m) {
  const W ap=m.alpha_perp,ar=m.alpha_parallel;
  return k*std::sqrt((1-mu*mu)/(ap*ap)+mu*mu/(ar*ar));
}
bool partition(const p::Spectrum& s,W k,const p::ModelPoint& m,
               Counters& c,std::vector<W>& edges) {
  if (!c.probe()) return false;
  const W left=mapped_k(0,k,m);
  if (!c.probe()) return false;
  const W right=mapped_k(1,k,m);
  if (std::min(left,right)<s.k_per_mpc.front() ||
      std::max(left,right)>s.k_per_mpc.back()) return false;
  edges={0};
  if (m.alpha_perp!=m.alpha_parallel) {
    for (double knot:s.k_per_mpc) {
      if (!(knot>std::min(left,right) && knot<std::max(left,right))) continue;
      W a=0,b=1;
      bool resolved=false;
      for (unsigned proposal=0;proposal<128;++proposal) {
        const W mid=a+(b-a)/2;
        if (mid==a || mid==b) {
          if (!c.probe()) return false;
          const W va=mapped_k(a,k,m);
          if (!c.probe()) return false;
          const W vb=mapped_k(b,k,m);
          if (!(knot>=std::min(va,vb) && knot<=std::max(va,vb))) return false;
          const W error=std::min(std::abs(va-knot),std::abs(vb-knot));
          if (error>256*eps*std::abs(static_cast<W>(knot))) return false;
          const W chosen=std::abs(va-knot)<=std::abs(vb-knot) ? a : b;
          if (!(chosen>0 && chosen<1)) return false;
          edges.push_back(chosen);resolved=true;break;
        }
        if (!c.probe()) return false;
        const W value=mapped_k(mid,k,m);
        if (value==knot) { edges.push_back(mid);resolved=true;break; }
        if ((value<knot)==(right>left)) a=mid; else b=mid;
      }
      if (!resolved) return false;
    }
    std::sort(edges.begin()+1,edges.end());
    for (std::size_t j=1;j<edges.size();++j) if (edges[j]<=edges[j-1]) return false;
  }
  edges.push_back(1);return true;
}
bool independent_power(W mu,W k,const p::Spectrum& s,const p::ModelPoint& m,
                       Counters& c,W& out) {
  if (!c.probe()) return false;
  const W ap=m.alpha_perp,ar=m.alpha_parallel;
  const W q=std::sqrt((1-mu*mu)/(ap*ap)+mu*mu/(ar*ar));
  const W kt=k*q,mut=mu/(ar*q);
  if (kt<s.k_per_mpc.front() || kt>s.k_per_mpc.back()) return false;
  std::size_t upper=1;
  while (upper+1<s.k_per_mpc.size() && kt>s.k_per_mpc[upper]) ++upper;
  const W x0=s.k_per_mpc[upper-1],x1=s.k_per_mpc[upper];
  const W value=s.power_mpc3[upper-1]+(kt-x0)/(x1-x0)*
                 (static_cast<W>(s.power_mpc3[upper])-s.power_mpc3[upper-1]);
  const W response=m.b1+static_cast<W>(m.f)*mut*mut;
  out=response*response*value/(ap*ap*ar);
  return std::isfinite(out) && out>=0;
}
bool integrate(const Rule& rule,const std::vector<W>& edges,W k,
               const p::Spectrum& s,const p::ModelPoint& m,Counters& c,
               std::pair<W,W>& result) {
  W p0=0,p2=0;
  for (std::size_t cell=1;cell<edges.size();++cell) {
    const W half=(edges[cell]-edges[cell-1])/2;
    const W center=(edges[cell]+edges[cell-1])/2;
    for (std::size_t j=0;j<rule.nodes.size();++j) {
      const W mu=center+half*rule.nodes[j];W value=0;
      if (!independent_power(mu,k,s,m,c,value)) return false;
      const W weighted=half*rule.weights[j]*value;
      p0+=weighted;p2+=5*weighted*(3*mu*mu-1)/2;
    }
  }
  result={p0,p2};return std::isfinite(p0) && std::isfinite(p2);
}
p::Source source(unsigned type) {
  p::Source s;
  s.identity={"synthetic/peer-dyadic-table","synthetic/peer-identity-window",
      "supplied-matter","supplied-f","typed0then2","unqualified-synthetic"};
  s.spectrum.redshift=.32;
  s.spectrum.k_per_mpc={.015625,.03125,.0625,.125,.25};
  if (type==0) s.spectrum.power_mpc3={3000,3000,3000,3000,3000};
  else if (type==1) s.spectrum.power_mpc3={64,128,256,512,1024};
  else s.spectrum.power_mpc3={1000,1200,900,1600,2000};
  s.window.coordinates=p::Coordinates::physical_mpc;
  s.window.theory_k={.03125,.0625,.125};
  s.window.output0_k=s.window.output2_k=s.window.theory_k;
  s.window.output0_ids={0,1,2};s.window.output2_ids={3,4,5};
  s.window.w00=s.window.w22={1,0,0,0,1,0,0,0,1};
  s.window.w02=s.window.w20=std::vector<double>(9,0);
  return s;
}
#define REQUIRE(x) do { if (!(x)) { std::cerr << "peer line " << __LINE__ << ": " << #x << '\n'; return 1; } } while (false)
}
int main() {
  REQUIRE(std::fegetround()==FE_TONEAREST && std::numeric_limits<W>::digits>=64);
  Counters counts;Rule r16,r32,r64;
  W change16_32=0,change32_64=0;
  REQUIRE(make_rule(16,r16,counts));REQUIRE(make_rule(32,r32,counts));REQUIRE(make_rule(64,r64,counts));
  const p::ModelPoint models[]{ {.32,2,1,1.1,.9}, {.32,1,-1,.9,1.1}, {.32,2,.8,1.03,.97} };
  for (unsigned type=0;type<3;++type) {
    for (const auto& model:models) {
      auto s=source(type);
      std::vector<std::pair<W,W>> reference;
      reference.reserve(s.window.theory_k.size());
      for (double k:s.window.theory_k) {
        REQUIRE(++counts.queries<=64 && s.spectrum.k_per_mpc.size()<=4096);
        std::vector<W> edges;edges.reserve(s.spectrum.k_per_mpc.size()+2);
        REQUIRE(partition(s.spectrum,k,model,counts,edges));
        std::pair<W,W> coarse,medium,fine;
        REQUIRE(integrate(r16,edges,k,s.spectrum,model,counts,coarse));
        REQUIRE(integrate(r32,edges,k,s.spectrum,model,counts,medium));
        REQUIRE(integrate(r64,edges,k,s.spectrum,model,counts,fine));
        change16_32=std::max({change16_32,std::abs(coarse.first-medium.first),std::abs(coarse.second-medium.second)});
        change32_64=std::max({change32_64,std::abs(medium.first-fine.first),std::abs(medium.second-fine.second)});
        REQUIRE(std::abs(medium.first-fine.first)<=1e-7L);
        REQUIRE(std::abs(medium.second-fine.second)<=1e-7L);
        reference.push_back(fine);
      }
      auto owner=p::prepare(std::move(s),{2000000,payload});
      REQUIRE(owner.status()==p::Status::ok);
      const auto actual=owner.evaluate(model,{1e-7,1e-6,8000000,2000000,4096,20,payload,true});
      REQUIRE(actual.status==p::Status::ok && actual.input_multipoles.size()==reference.size());
      for (std::size_t j=0;j<reference.size();++j) {
        REQUIRE(std::abs(static_cast<W>(actual.input_multipoles[j].p0.value)-reference[j].first)<=1e-6L);
        REQUIRE(std::abs(static_cast<W>(actual.input_multipoles[j].p2.value)-reference[j].second)<=1e-6L);
      }
      // Explicit finite fixed-source inventory; not a portable whole RSS bound.
      REQUIRE(owner.retained_payload_bytes() && *owner.retained_payload_bytes()+
          actual.known_live_payload_estimate_bytes+
          reference.capacity()*sizeof(std::pair<W,W>)+
          (r16.nodes.capacity()+r16.weights.capacity()+r32.nodes.capacity()+
           r32.weights.capacity()+r64.nodes.capacity()+r64.weights.capacity())*sizeof(W)<payload);
    }
  }
  REQUIRE(counts.roots==56 && counts.newton<=56*64 && counts.recurrences<=56*130);
  REQUIRE(counts.angular<=1000000 && counts.queries==27);
  std::cout << "windowed linear GL peer PASS; probes=" << counts.angular
            << ", queries=" << counts.queries << ", refinement16-32=" << change16_32
            << ", refinement32-64=" << change32_64 << "; empirical refinement only\n";
}
