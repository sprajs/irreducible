// Frozen original Decimal110/150 exact binary inputs and Machin pi;
// charge neutrality solved by direct bisection, independent of native Newton.
// Shared SI/CODATA/ASD central facts; supplied synthetic effective masses.
// Density2e-15rel, LTE5e-13rel, reference refinement4.597e-108 (<5%).
#include "irred/baryon_abundance.hpp"
#include "baryon_abundance_peer_facts.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace {
using namespace irred::cosmology;
namespace f=abundance_peer_facts;
using S=irred::numerics::Status;using W=long double;
unsigned checks=0;
void need(bool b,const char* why){++checks;if(!b)throw std::runtime_error(why);}
void near(double a,W b,W budget,const char* why){need(std::abs(W(a)/b-1)<=budget,why);}
BaryonAbundanceSource source(){return {f::source[0],f::source[1],f::source[2],f::source[3],"original synthetic omega/Y","supplied neutral effective synthetic masses"};}
void reference(){
 auto op=prepare_baryon_abundance(source());need(op.status()==S::ok,"mapped owner");
 std::array<double,5> scale{};std::array<BaryonAbundanceLteQuery,5> queries{};
 for(size_t i=0;i<5;++i){scale[i]=f::queries[i][0];queries[i]={scale[i],f::queries[i][1]};}
 auto mapped=op.evaluate(scale);need(mapped.status==S::ok&&mapped.rows.size()==5,"mapped batch");
 auto lte=op.evaluate_equilibrium(queries,"independent supplied synthetic matter temperatures");
 need(lte.status==S::ok&&lte.rows.size()==5&&lte.solves==5,"one coarse LTE batch");
 for(size_t i=0;i<5;++i){auto &r=mapped.rows[i];need(r.hydrogen_nuclei.status==S::ok&&r.helium_nuclei.status==S::ok,"independent mapping groups");
  near(*r.hydrogen_nuclei.value,f::values[i][0],2e-15L,"Decimal hydrogen mapping");near(*r.helium_nuclei.value,f::values[i][1],2e-15L,"Decimal helium mapping");
  need(r.hydrogen_nuclei.relative_arithmetic_estimate<=1e-15&&r.helium_nuclei.relative_arithmetic_estimate<=1e-15,"mapping diagnostic budget");
  const auto &e=lte.rows[i].equilibrium;std::array<const irred::atomic::HydrogenHeliumValue*,6> v{&e.hydrogen_neutral,&e.hydrogen_ionized,&e.helium_neutral,&e.helium_singly_ionized,&e.helium_doubly_ionized,&e.electron_density};
  for(size_t j=0;j<6;++j){need(v[j]->status==S::ok&&v[j]->value,"native LTE output admitted");near(*v[j]->value,f::values[i][j+2],5e-13L,"independent charge root fractions/electrons");need(v[j]->relative_arithmetic_estimate<=1e-13,"combined LTE diagnostic");}
  W rho=W(*r.hydrogen_nuclei.value)*f::source[2]+W(*r.helium_nuclei.value)*f::source[3];
  near(double(W(*r.helium_nuclei.value)*f::source[3]/rho),f::source[1],2e-15L,"inverse mass fraction");
  W charge=W(*r.hydrogen_nuclei.value)*(*e.hydrogen_ionized.value)+W(*r.helium_nuclei.value)*(*e.helium_singly_ionized.value+2*W(*e.helium_doubly_ionized.value));
  near(*e.electron_density.value,charge,5e-13L,"independent charge conservation");
 }
 auto copy=op;auto moved=std::move(op);need(op.status()!=S::ok,"moved owner invalid");auto repeated=moved.evaluate(scale);
 need(*repeated.rows[2].hydrogen_nuclei.value==*mapped.rows[2].hydrogen_nuclei.value&&copy.source()->mass_origin==source().mass_origin,"retained source and bits");
}
void endpoints(){
 std::array<double,2> a{1,.5};
 for(double Y:{0.,1.,std::nextafter(1.,0.)}){auto s=source();s.helium4_mass_fraction=Y;auto op=prepare_baryon_abundance(s);auto r=op.evaluate(a);need(r.status==S::ok,"endpoint mapping admitted");
  for(size_t i=0;i<2;++i){need(r.rows[i].hydrogen_nuclei.value&&r.rows[i].helium_nuclei.value,"endpoint species values");if(Y==0)need(*r.rows[i].helium_nuclei.value==0,"only exact absent helium zero");if(Y==1)need(*r.rows[i].hydrogen_nuclei.value==0,"only exact absent hydrogen zero");}
  if(Y==0||Y==1){std::array<BaryonAbundanceLteQuery,1> q{{{1,1e4}}};auto e=op.evaluate_equilibrium(q,"supplied");need(e.solves==0&&e.rows[0].equilibrium.admission_status!=S::ok,"absent species LTE explicit refusal");}
 }
 auto s=source();s.helium4_mass_fraction=.5;s.helium4_effective_mass_kg=3*s.hydrogen1_effective_mass_kg;auto r=prepare_baryon_abundance(s).evaluate(a);
 near(*r.rows[0].hydrogen_nuclei.value,W(*r.rows[0].helium_nuclei.value)*3,2e-15L,"unequal masses are not automatic factor four");near(*r.rows[1].hydrogen_nuclei.value,W(*r.rows[0].hydrogen_nuclei.value)*8,2e-15L,"a^-3 scaling");
}
}
int main(){try{reference();endpoints();std::cout<<"PASS "<<checks<<" independent abundance controls\n";}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
