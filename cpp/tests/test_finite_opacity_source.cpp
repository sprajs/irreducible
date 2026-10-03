#include "irred/finite_opacity_source.hpp"
#include "../src/finite_opacity_transport.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <type_traits>
#include <vector>

namespace {
namespace f=irred::cosmology::detail::finite_opacity;
using W=long double;
void need(bool ok,const char *why) {
  if(!ok){std::cerr<<"FAIL "<<why<<'\n';std::exit(1);}
}
f::Epoch coefficients() {
  f::Epoch e;e.status=irred::numerics::Status::ok;e.hcal=.05L;
  e.fb=.1L;e.fc=.2L;e.fg=.7L;e.fr=e.fg;e.x2=4;e.enthalpy_fraction=e.fb+e.fc+4*e.fg/3;
  return e;
}
void limits() {
  const auto e=coefficients();const W k=.1L,opacity=17;
  const auto seed=f::finite_seed(e,k);W scale=0;
  need(std::abs(f::hamiltonian(e,k,seed,scale))<1e-17L*scale,"finite seed Hamiltonian");
  need(std::abs(-seed[f::phi]+seed[f::db]/3-1)<1e-17L,"local zeta one");
  need(seed[f::db]==seed[f::dc] && seed[f::dg]==4*seed[f::db]/3 &&
       seed[f::ug]==seed[f::ub] && seed[f::ug]==seed[f::uc],"finite entropy/slip seed");
  for(unsigned i:{f::sg,f::g0,f::g1,f::g2})need(seed[i]==0,"explicit initial higher moments zero");
  auto y=seed;y[f::ug]=.2L;y[f::ub]=-.3L;
  const auto scattered=f::core_derivative(e,k,opacity,y,0,0);
  const auto free=f::core_derivative(e,k,0,y,0,0);
  need(scattered[f::dg]==free[f::dg],"zero monopole Thomson energy exchange");
  need(std::abs(4*e.fg*(scattered[f::ug]-free[f::ug])/3+
                e.fb*(scattered[f::ub]-free[f::ub]))<1e-16L,"enthalpy weighted momentum exchange");
  f::Core fast{};fast[f::sg]=.05L;fast[f::g0]=.5L;fast[f::g2]=.1L;
  const auto decay=f::core_derivative(e,k,opacity,fast,0,0);
  need(std::abs(decay[f::sg]+.3L*opacity*fast[f::sg])<1e-17L &&
       std::abs(decay[f::g0]+.3L*opacity*fast[f::g0])<1e-17L &&
       std::abs(decay[f::g2]+.3L*opacity*fast[f::g2])<1e-17L,"polarized fast eigenvalue -0.3K");
  f::Core slow{};slow[f::ug]=.2L;slow[f::sg]=16*k*slow[f::ug]/(45*opacity);
  slow[f::g0]=2.5L*slow[f::sg];slow[f::g2]=.5L*slow[f::sg];
  const auto stationary=f::core_derivative(e,k,opacity,slow,0,0);
  need(std::abs(stationary[f::sg])<1e-18L && std::abs(stationary[f::g0])<1e-18L &&
       std::abs(stationary[f::g2])<1e-18L,"stationary leading polarized 16/45 shear");
  need(f::potential_psi(e,k,slow)!=slow[f::phi],"finite photon shear keeps phi != psi");
  const W w=.4L;const auto source=f::raw_source(e,k,opacity,w,slow);
  const W pi=2*slow[f::sg]+slow[f::g0]+slow[f::g2];
  need(source[2]==opacity*w*pi/8 && source[3]==std::sqrt(6.L)*source[2],"Pi=2sigma+G0+G2 split normalization");
}
// Independent dense stage algorithm for a small SYNTHETIC finite hierarchy.
// Independently assembled equations use theta and F2, converting only at the
// boundary. This controls the native Schur implementation; it is not the
// angular/precision reference or full physical/error qualification.
std::vector<W> dense_derivative(const f::Epoch &e,W k,W opacity,W eta,
                               unsigned lmax,const std::vector<W> &x) {
  const unsigned tail=lmax-2,n=11+2*tail;std::vector<W> d(n);
  const W theta_g=k*x[1],theta_b=k*x[7],theta_c=k*x[9],F2=2*x[2];
  const W psi=x[10]-3*e.hcal*e.hcal*e.fg*F2/(k*k);
  const W phip=-e.hcal*psi+1.5L*e.hcal*e.hcal*(e.fc*theta_c+e.fb*theta_b+4*e.fg*theta_g/3)/(k*k);
  const W pi=F2+x[3]+x[5];
  d[0]=-4*theta_g/3+4*phip;
  d[1]=(k*k*(x[0]/4-F2/2+psi)+opacity*(theta_b-theta_g))/k;
  d[2]=(8*theta_g/15-3*k*x[11]/5-opacity*F2+opacity*pi/10)/2;
  d[3]=-k*x[4]+opacity*(-x[3]+pi/2);
  d[4]=k*(x[3]-2*x[5])/3-opacity*x[4];
  d[5]=k*(2*x[4]-3*x[11+tail])/5+opacity*(-x[5]+pi/10);
  d[6]=-theta_b+3*phip;
  d[7]=(-e.hcal*theta_b+k*k*psi+4*e.fg*opacity*(theta_g-theta_b)/(3*e.fb))/k;
  d[8]=-theta_c+3*phip;d[9]=(-e.hcal*theta_c+k*k*psi)/k;d[10]=phip;
  for(unsigned role=0;role<2;++role) for(unsigned l=3;l<=lmax;++l) {
    const unsigned index=11+role*tail+l-3;
    const W previous=l==3?(role?x[5]:F2):x[index-1];
    d[index]=l==lmax?k*previous-(l+1)*x[index]/eta-opacity*x[index]:
        k*(l*previous-(l+1)*x[index+1])/(2*l+1)-opacity*x[index];
  }
  return d;
}
std::vector<W> dense_step(const std::array<f::Epoch,2> &e,const f::Pair &opacity,
                          W eta,W h,W k,unsigned lmax,const f::TransportState &old) {
  const unsigned tail=lmax-2,n=11+2*tail,total=2*n;
  std::vector<W> initial(n),rhs(total),matrix(total*total);
  for(unsigned j=0;j<11;++j)initial[j]=old.core[j];
  for(unsigned l=3;l<=lmax;++l){initial[11+l-3]=old.temperature[l];initial[11+tail+l-3]=old.polarization[l];}
  for(unsigned i=0;i<total;++i){rhs[i]=initial[i%n];matrix[i*total+i]=1;}
  const W a[2][2]={{5.L/12,-1.L/12},{3.L/4,1.L/4}};
  for(unsigned s=0;s<2;++s) for(unsigned j=0;j<n;++j) {
    std::vector<W> unit(n);unit[j]=1;
    const auto d=dense_derivative(e[s],k,opacity[s],eta+(s?h:h/3),lmax,unit);
    for(unsigned r=0;r<2;++r) for(unsigned i=0;i<n;++i)
      matrix[(r*n+i)*total+s*n+j]-=h*a[r][s]*d[i];
  }
  // Gauss-Jordan, distinct from the native LU/block elimination.
  for(unsigned j=0;j<total;++j) {
    unsigned p=j;for(unsigned i=j+1;i<total;++i)if(std::abs(matrix[i*total+j])>std::abs(matrix[p*total+j]))p=i;
    for(unsigned c=0;c<total;++c)std::swap(matrix[j*total+c],matrix[p*total+c]);
    std::swap(rhs[j],rhs[p]);
    const W pivot=matrix[j*total+j];need(std::abs(pivot)>1e-24L,"dense synthetic pivot");
    for(unsigned c=0;c<total;++c)matrix[j*total+c]/=pivot;
    rhs[j]/=pivot;
    for(unsigned i=0;i<total;++i)if(i!=j) {
      const W q=matrix[i*total+j];
      for(unsigned c=0;c<total;++c)matrix[i*total+c]-=q*matrix[j*total+c];
      rhs[i]-=q*rhs[j];
    }
  }
  return {rhs.begin()+n,rhs.end()};
}
void schur_and_refusals() {
  const auto e=coefficients();const W eta=20,h=.08L,k=.1L;const unsigned lmax=8;
  const std::array<f::Epoch,2> stages{e,e};
  for(W stiffness:{0.L,3.L,100000.L}) {
    const f::Pair opacity{stiffness,1.03L*stiffness};f::TransportState state;
    state.core=f::finite_seed(e,k);
    for(unsigned l=3;l<=lmax;++l){state.temperature[l]=.001L/l;state.polarization[l]=-.002L/l;}
    const auto expected=dense_step(stages,opacity,eta,h,k,lmax,state);
    irred::cosmology::FiniteOpacityWork work;irred::cosmology::FiniteOpacityPolicy policy;
    f::Ledger ledger{work,policy};f::TransportScratch scratch;irred::cosmology::FiniteOpacityAttemptReceipt receipt;
    need(f::radau_step(stages,opacity,lmax,eta,h,k,state,scratch,ledger,receipt),"synthetic Schur step");
    for(unsigned j=0;j<11;++j)need(std::abs(state.core[j]-expected[j])<2e-14L,"core vs independent dense step");
    for(unsigned l=3;l<=lmax;++l) {
      need(std::abs(state.temperature[l]-expected[11+l-3])<2e-14L &&
           std::abs(state.polarization[l]-expected[11+lmax-2+l-3])<2e-14L,"both tails vs dense step");
    }
    need(work.coupled_stage_solves==1 && work.tail_block_inversions==2*(lmax-2) &&
         work.core_factorizations==1 && work.destination_writes>1000,"causal solve/block/write counters");
  }
  f::TransportState state;state.core=f::finite_seed(e,k);const auto saved=state.core;
  irred::cosmology::FiniteOpacityWork work;irred::cosmology::FiniteOpacityPolicy policy;policy.maximum_destination_writes=1;
  f::Ledger ledger{work,policy};f::TransportScratch scratch;irred::cosmology::FiniteOpacityAttemptReceipt receipt;
  need(!f::radau_step(stages,{1,1},lmax,eta,h,k,state,scratch,ledger,receipt) &&
       ledger.status==irred::numerics::Status::work_limit && state.core==saved &&
       work.coupled_stage_solves==1 && work.refused_write_request>0,"write refusal preserves original state");
}
void record_prefix_and_nonfinite_residuals() {
  using namespace irred::cosmology;
  FiniteOpacityPolicy policy;
  const std::size_t node_total=f::node_local_write_allowance+f::node_owned_copy_allowance;
  for(bool deny_local:{true,false}) {
    FiniteOpacityWork work;policy.maximum_destination_writes=deny_local?node_total:2*node_total-1;
    f::Ledger ledger{work,policy};std::vector<FiniteOpacityNode> nodes;nodes.reserve(2);
    need(ledger.writes(f::node_local_write_allowance),"first node local allowance");
    FiniteOpacityNode first;first.eta_index=7;
    need(f::publish_node(std::move(first),nodes,ledger),"first owned node publication");
    const bool local=ledger.writes(f::node_local_write_allowance);
    if(local) {FiniteOpacityNode second;second.eta_index=8;need(!f::publish_node(std::move(second),nodes,ledger),"node copy denied");}
    need(local!=deny_local && nodes.size()==1 && nodes[0].eta_index==7 &&
         ledger.status==irred::numerics::Status::work_limit && work.refused_write_request==
             (deny_local?f::node_local_write_allowance:f::node_owned_copy_allowance),"node prefix and exact denied allowance");
  }
  const std::size_t diagnostic_total=f::diagnostic_local_write_allowance+f::diagnostic_owned_copy_allowance;
  for(bool deny_local:{true,false}) {
    FiniteOpacityWork work;policy.maximum_destination_writes=deny_local?diagnostic_total:2*diagnostic_total-1;
    f::Ledger ledger{work,policy};std::vector<FiniteOpacityChannelDiagnostic> records;records.reserve(2);
    need(ledger.writes(f::diagnostic_local_write_allowance),"first diagnostic local allowance");
    FiniteOpacityChannelDiagnostic first;first.k_index=3;
    need(f::publish_diagnostic(std::move(first),records,ledger),"first owned diagnostic publication");
    const bool local=ledger.writes(f::diagnostic_local_write_allowance);
    if(local) {FiniteOpacityChannelDiagnostic second;second.k_index=4;need(!f::publish_diagnostic(std::move(second),records,ledger),"diagnostic copy denied");}
    need(local!=deny_local && records.size()==1 && records[0].k_index==3 &&
         ledger.status==irred::numerics::Status::work_limit && work.refused_write_request==
             (deny_local?f::diagnostic_local_write_allowance:f::diagnostic_owned_copy_allowance),"diagnostic prefix and denied allowance");
  }
  for(W hostile:{std::numeric_limits<W>::quiet_NaN(),std::numeric_limits<W>::infinity()}) {
    FiniteOpacityWork work;policy=FiniteOpacityPolicy{};f::Ledger ledger{work,policy};W maximum=.25L;
    need(!f::accumulate_stage_residual(1,1,hostile,0,1,.5L,.5L,maximum,ledger) &&
         ledger.status==irred::numerics::Status::overflow && maximum==.25L,
         "nonfinite derivative refused before maximum aggregation");
  }
  FiniteOpacityWork work;f::Ledger ledger{work,policy};W maximum=.25L;
  need(!f::accumulate_stage_residual(1,1,std::numeric_limits<W>::max(),0,2,1,0,maximum,ledger) &&
       ledger.status==irred::numerics::Status::overflow && maximum==.25L,"finite input overflowing term refused");
}
} // namespace
int main() {
  static_assert(!std::is_copy_constructible_v<irred::cosmology::FiniteOpacitySourceProducer>);
  static_assert(!std::is_copy_constructible_v<irred::cosmology::FiniteOpacitySourceResult>);
  limits();schur_and_refusals();record_prefix_and_nonfinite_residuals();
  std::cout<<"finite opacity analytic/source/Schur contract passed\n";
}
