#pragma once
#include "finite_opacity_equations.hpp"
#include <algorithm>
#include <array>
#include <limits>

namespace irred::cosmology::detail::finite_opacity {
using S = numerics::Status;
inline constexpr std::array<std::array<W,2>,2> radau_a{{
    {{5.L/12,-1.L/12}},{{3.L/4,1.L/4}}}};
struct Ledger {
  FiniteOpacityWork &work;
  const FiniteOpacityPolicy &policy;
  S status=S::ok;
  bool count(std::size_t &counter, std::size_t count, std::size_t limit) {
    if (count>limit || counter>limit-count) { status=S::work_limit; return false; }
    counter+=count; return true;
  }
  bool writes(std::size_t n) {
    if(work.reserved_refusal_writes>policy.maximum_destination_writes ||
       !count(work.destination_writes,n,policy.maximum_destination_writes-work.reserved_refusal_writes)) {
      work.refused_write_request=n; return false;
    }
    return true;
  }
  bool reserve_refusal_writes(std::size_t n) {
    if(work.destination_writes>policy.maximum_destination_writes ||
       !count(work.reserved_refusal_writes,n,policy.maximum_destination_writes-work.destination_writes)) {
      work.refused_write_request=n;return false;
    }
    return true;
  }
  void commit_refusal_writes(std::size_t n) {
    work.reserved_refusal_writes-=n;work.destination_writes+=n;
  }
  bool background() {
    return count(work.background_clock_calls,1,policy.maximum_background_clock_calls);
  }
  bool step() { return count(work.attempted_steps,1,policy.maximum_attempted_steps); }
  bool solve() { return count(work.coupled_stage_solves,1,policy.maximum_coupled_stage_solves); }
};
using Pair = std::array<W,2>;
using Block = std::array<std::array<W,2>,2>;
inline Pair multiply(const Block &a,const Pair &b) {
  return {a[0][0]*b[0]+a[0][1]*b[1],a[1][0]*b[0]+a[1][1]*b[1]};
}
inline Block multiply(const Block &a,const Block &b) {
  Block out;
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)
    out[i][j]=a[i][0]*b[0][j]+a[i][1]*b[1][j];
  return out;
}
struct TailReduction { Block inverse{}, lower{}; Pair rhs{}; };
struct TransportState {
  Core core{};
  std::array<W,193> temperature{}, polarization{}; // only indices 3..L
};
struct TransportScratch {
  std::array<TailReduction,193> f{}, p{};
  std::array<std::array<W,22>,22> matrix{};
  std::array<W,22> rhs{}, core_solution{};
  std::array<TransportState,2> stages{};
};
// Each elimination is a charged 2x2 inversion. Determinant scaling avoids an
// overflow-prone unscaled test; no regularization or dense hierarchy fallback.
inline bool invert(const Block &m, Block &out, Ledger &ledger, W &minimum_pivot) {
  if (!ledger.count(ledger.work.tail_block_inversions,1,
                    std::numeric_limits<std::size_t>::max())) return false;
  W scale=0; for(const auto &r:m) for(W x:r) scale=std::max(scale,std::abs(x));
  if (!(scale>0) || !std::isfinite(scale)) { ledger.status=S::singular; return false; }
  const W a=m[0][0]/scale,b=m[0][1]/scale,c=m[1][0]/scale,d=m[1][1]/scale;
  const W determinant=a*d-b*c;
  minimum_pivot=std::min(minimum_pivot,std::abs(determinant));
  if (!std::isfinite(determinant) || std::abs(determinant)<1024*std::numeric_limits<W>::epsilon()) {
    ledger.status=S::conditioning_budget_exceeded; return false;
  }
  if (!ledger.writes(4)) return false;
  out[0][0]=d/(determinant*scale);out[0][1]=-b/(determinant*scale);
  out[1][0]=-c/(determinant*scale);out[1][1]=a/(determinant*scale);
  for(const auto &r:out) for(W x:r) if(!std::isfinite(x)) { ledger.status=S::overflow; return false; }
  return true;
}
inline bool reduce_tail(unsigned lmax,W eta,W h,W k,const Pair &opacity,
                        const std::array<W,193> &old,
                        std::array<TailReduction,193> &tail,
                        Ledger &ledger,W &minimum_pivot) {
  for(unsigned l=lmax;;--l) {
    Block diagonal,lower,upper; // all twelve entries assigned below
    if(!ledger.writes(14)) return false; // twelve block entries plus rhs pair
    for(unsigned i=0;i<2;++i) for(unsigned j=0;j<2;++j) {
      const W stage_eta=eta+(j==0?h/3:h);
      const W d=-opacity[j]-(l==lmax?(l+1)/stage_eta:0);
      const W low=l==lmax?k:k*l/(2*l+1);
      const W up=l==lmax?0:-k*(l+1)/(2*l+1);
      diagonal[i][j]=(i==j?1:0)-h*radau_a[i][j]*d;
      lower[i][j]=-h*radau_a[i][j]*low;
      upper[i][j]=-h*radau_a[i][j]*up;
    }
    Pair r{old[l],old[l]};
    if(l<lmax) {
      if(!ledger.writes(10)) return false; // two block products and one pair
      const auto v=multiply(upper,tail[l+1].inverse);
      const auto q=multiply(v,tail[l+1].lower); const auto z=multiply(v,tail[l+1].rhs);
      if(!ledger.writes(6)) return false;
      for(unsigned i=0;i<2;++i) {
        r[i]-=z[i]; for(unsigned j=0;j<2;++j) diagonal[i][j]-=q[i][j];
      }
    }
    if(!invert(diagonal,tail[l].inverse,ledger,minimum_pivot)) return false;
    if(!ledger.writes(6)) return false;
    tail[l].lower=lower; tail[l].rhs=r;
    if(l==3) break;
  }
  return true;
}
inline bool factor_core(TransportScratch &scratch,Ledger &ledger,W &minimum_pivot) {
  if(!ledger.count(ledger.work.core_factorizations,1,
                   std::numeric_limits<std::size_t>::max())) return false;
  auto &a=scratch.matrix; auto &b=scratch.rhs;
  for(unsigned col=0;col<22;++col) {
    unsigned pivot=col; W norm=0;
    for(unsigned row=col;row<22;++row) {
      if(std::abs(a[row][col])>std::abs(a[pivot][col])) pivot=row;
      for(unsigned j=col;j<22;++j) norm=std::max(norm,std::abs(a[row][j]));
    }
    const W relative=norm>0?std::abs(a[pivot][col])/norm:0;
    minimum_pivot=std::min(minimum_pivot,relative);
    if(!std::isfinite(relative) || relative<1024*std::numeric_limits<W>::epsilon()) {
      ledger.status=S::conditioning_budget_exceeded; return false;
    }
    if(pivot!=col) {
      if(!ledger.writes(46)) return false;
      std::swap(a[pivot],a[col]); std::swap(b[pivot],b[col]);
    }
    for(unsigned row=col+1;row<22;++row) {
      const W multiplier=a[row][col]/a[col][col];
      if(!ledger.writes(23-col)) return false;
      a[row][col]=0;
      for(unsigned j=col+1;j<22;++j) a[row][j]-=multiplier*a[col][j];
      b[row]-=multiplier*b[col];
    }
  }
  for(unsigned row=22;row-->0;) {
    W v=b[row]; for(unsigned j=row+1;j<22;++j) v-=a[row][j]*scratch.core_solution[j];
    if(!ledger.writes(1)) return false;
    scratch.core_solution[row]=v/a[row][row];
    if(!std::isfinite(scratch.core_solution[row])) { ledger.status=S::overflow; return false; }
  }
  return true;
}
inline bool back_substitute(unsigned lmax,const std::array<TailReduction,193> &tail,
                            Pair previous,std::array<TransportState,2> &stages,
                            bool temperature,Ledger &ledger) {
  for(unsigned l=3;l<=lmax;++l) {
    if(!ledger.writes(8)) return false; // temporary pairs and previous update
    const auto q=multiply(tail[l].lower,previous);
    const Pair r{tail[l].rhs[0]-q[0],tail[l].rhs[1]-q[1]};
    previous=multiply(tail[l].inverse,r);
    if(!ledger.writes(2)) return false;
    for(unsigned s=0;s<2;++s) {
      if(!std::isfinite(previous[s])) { ledger.status=S::overflow; return false; }
      (temperature?stages[s].temperature:stages[s].polarization)[l]=previous[s];
    }
  }
  return true;
}
inline W tail_derivative(unsigned l,unsigned lmax,W eta,W k,W opacity,
                         const TransportState &y,bool temperature) {
  const auto &tail=temperature?y.temperature:y.polarization;
  const W previous=l==3?(temperature?2*y.core[sg]:y.core[g2]):tail[l-1];
  if(l==lmax) return k*previous-(l+1)*tail[l]/eta-opacity*tail[l];
  return k*(l*previous-(l+1)*tail[l+1])/(2*l+1)-opacity*tail[l];
}
// One actual full linear Radau stage solve. Tails are reduced from L down to 3,
// the 22 core is pivoted, then both tails are recovered. Stiff accuracy returns
// stage 2; residual checks use the original full operator, not the Schur rows.
inline bool radau_step(const std::array<Epoch,2> &epoch,const Pair &opacity,
                       unsigned lmax,W eta,W h,W k,TransportState &state,
                       TransportScratch &scratch,Ledger &ledger,
                       FiniteOpacityAttemptReceipt &receipt) {
  if(!ledger.solve()) return false;
  if(!reduce_tail(lmax,eta,h,k,opacity,state.temperature,scratch.f,ledger,receipt.minimum_scaled_pivot) ||
     !reduce_tail(lmax,eta,h,k,opacity,state.polarization,scratch.p,ledger,receipt.minimum_scaled_pivot)) return false;
  if(!ledger.writes(12)) return false; // two reduced pairs and two blocks
  const Pair fr=multiply(scratch.f[3].inverse,scratch.f[3].rhs);
  const Pair pr=multiply(scratch.p[3].inverse,scratch.p[3].rhs);
  const Block fq=multiply(scratch.f[3].inverse,scratch.f[3].lower);
  const Block pq=multiply(scratch.p[3].inverse,scratch.p[3].lower);
  if(!ledger.writes(506)) return false;
  for(unsigned i=0;i<22;++i) {
    scratch.rhs[i]=state.core[i%11];
    for(unsigned j=0;j<22;++j) scratch.matrix[i][j]=i==j?1:0;
  }
  for(unsigned stage=0;stage<2;++stage) for(unsigned col=0;col<11;++col) {
    if(!ledger.writes(23)) return false; // unit initialization/update and force
    Core unit{}; unit[col]=1;
    const auto d=core_derivative(epoch[stage],k,opacity[stage],unit,0,0);
    if(!ledger.writes(22)) return false;
    for(unsigned rowstage=0;rowstage<2;++rowstage) for(unsigned row=0;row<11;++row)
      scratch.matrix[rowstage*11+row][stage*11+col]-=h*radau_a[rowstage][stage]*d[row];
  }
  for(unsigned i=0;i<2;++i) for(unsigned j=0;j<2;++j) {
    const W fforce=-3*k/10,pforce=-3*k/5,ha=h*radau_a[i][j];
    if(!ledger.writes(6)) return false;
    scratch.rhs[i*11+sg]+=ha*fforce*fr[j];
    scratch.rhs[i*11+g2]+=ha*pforce*pr[j];
    for(unsigned s=0;s<2;++s) {
      // F2=2 sigma supplies the lower boundary of the temperature tail.
      scratch.matrix[i*11+sg][s*11+sg]+=ha*fforce*2*fq[j][s];
      scratch.matrix[i*11+g2][s*11+g2]+=ha*pforce*pq[j][s];
    }
  }
  if(!factor_core(scratch,ledger,receipt.minimum_scaled_pivot)) return false;
  if(!ledger.writes(22)) return false;
  for(unsigned s=0;s<2;++s) for(unsigned j=0;j<11;++j)
    scratch.stages[s].core[j]=scratch.core_solution[s*11+j];
  if(!back_substitute(lmax,scratch.f,{2*scratch.stages[0].core[sg],2*scratch.stages[1].core[sg]},scratch.stages,true,ledger) ||
     !back_substitute(lmax,scratch.p,{scratch.stages[0].core[g2],scratch.stages[1].core[g2]},scratch.stages,false,ledger)) return false;
  if(!ledger.writes(22)) return false;
  const std::array<Core,2> derivative{
      core_derivative(epoch[0],k,opacity[0],scratch.stages[0].core,scratch.stages[0].temperature[3],scratch.stages[0].polarization[3]),
      core_derivative(epoch[1],k,opacity[1],scratch.stages[1].core,scratch.stages[1].temperature[3],scratch.stages[1].polarization[3])};
  W maximum=0;
  for(unsigned s=0;s<2;++s) {
    for(unsigned j=0;j<11;++j) {
      const W term=h*(radau_a[s][0]*derivative[0][j]+radau_a[s][1]*derivative[1][j]);
      const W residual=scratch.stages[s].core[j]-state.core[j]-term;
      const W scale=1+std::abs(state.core[j])+std::abs(scratch.stages[s].core[j])+std::abs(term);
      maximum=std::max(maximum,std::abs(residual)/scale);
    }
    for(bool temperature:{true,false}) for(unsigned l=3;l<=lmax;++l) {
      const auto &old=temperature?state.temperature:state.polarization;
      const auto &y=temperature?scratch.stages[s].temperature:scratch.stages[s].polarization;
      const W d0=tail_derivative(l,lmax,eta+h/3,k,opacity[0],scratch.stages[0],temperature);
      const W d1=tail_derivative(l,lmax,eta+h,k,opacity[1],scratch.stages[1],temperature);
      const W term=h*(radau_a[s][0]*d0+radau_a[s][1]*d1);
      maximum=std::max(maximum,std::abs(y[l]-old[l]-term)/(1+std::abs(y[l])+std::abs(old[l])+std::abs(term)));
    }
  }
  receipt.maximum_stage_residual=std::max(receipt.maximum_stage_residual,maximum);
  if(!std::isfinite(maximum) || maximum>8192*std::numeric_limits<W>::epsilon()) {
    ledger.status=S::conditioning_budget_exceeded; return false;
  }
  if(!ledger.writes(11+2*(lmax-2))) return false;
  state.core=scratch.stages[1].core;
  for(unsigned l=3;l<=lmax;++l) {
    state.temperature[l]=scratch.stages[1].temperature[l];
    state.polarization[l]=scratch.stages[1].polarization[l];
  }
  return true;
}
} // namespace irred::cosmology::detail::finite_opacity
