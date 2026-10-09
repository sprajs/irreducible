#include "irred/plummer_population.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <array>
#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <utility>
using namespace irred::gravity;
using S=irred::numerics::Status;
using W=long double;
void need(bool condition,const char* message) {
  if(!condition) {std::cerr<<"FAIL "<<message<<'\n';std::abort();}
}
void close(W value,W reference,W relative=2e-9L) {
  if(!std::isfinite(value)||!std::isfinite(reference)||std::abs(value-reference)>relative*std::abs(reference)) {
    std::cerr<<std::setprecision(20)<<"comparison failure value="<<value<<" reference="<<reference<<" relative="<<relative<<'\n';
    std::abort();
  }
}
struct Integrals { W normalization=0,density=0,second=0,fourth=0,error=0; };
Integrals velocity_integrals(const PlummerPopulation& population,std::size_t index,std::size_t count) {
  const auto& state=population.radii()[index];
  W escape=std::sqrt(2*W(state.relative_potential_m2_s2.value));
  W step=std::numbers::pi_v<W>/(2*count);
  std::vector<PlummerVelocityRequest> requests;
  for(std::size_t i=0;i<count;++i) requests.push_back({index,double(escape*std::sin((i+0.5L)*step))});
  // Interior-only transformed midpoint rule; no unresolved endpoint is dropped
  // or replaced by zero. Near-tail point allocation is declared before running.
  PlummerPopulationEvaluationPolicy policy;policy.relative_tolerance=1e-7;
  auto batch=population.evaluate(requests,policy);need(batch.status==S::ok,"velocity integration all nodes");
  Integrals out;
  for(std::size_t i=0;i<count;++i) {
    const auto& row=batch.rows[i];need(row.support==PlummerEnergySupport::bound,"integration node bound");
    W theta=(i+0.5L)*step,dv=escape*std::cos(theta)*step,v=row.request.speed_m_s;
    out.normalization+=dv*row.speed_density_s_m.value;
    out.second+=dv*row.speed_density_s_m.value*v*v;
    out.fourth+=dv*row.speed_density_s_m.value*v*v*v*v;
    out.density+=dv*4*std::numbers::pi_v<W>*v*v*row.distribution_function_kg_s3_m6.value;
    out.error+=dv*row.speed_density_s_m.error_estimate;
  }
  need(batch.velocity_evaluations==count,"velocity work count");
  return out;
}
W projected_variance(PlummerSphere model,double radius,std::size_t count) {
  W scale=std::hypot(W(radius),W(model.scale_radius_metres));
  W step=std::numbers::pi_v<W>/(2*count);
  std::vector<double> radii;
  for(std::size_t i=0;i<count;++i) radii.push_back(double(std::hypot(W(radius),scale*std::tan((i+0.5L)*step))));
  auto prepared=prepare_plummer_population(model,radii);need(prepared.status()==S::ok,"LOS all native radii");
  W numerator=0;
  for(std::size_t i=0;i<count;++i) {
    const auto& state=prepared.radii()[i];W c=std::cos((i+0.5L)*step);
    numerator+=2*state.density_kg_m3.value*state.one_axis_variance_m2_s2.value*scale/(c*c)*step;
  }
  auto projection=project_plummer_sphere(model,radius);need(projection.surface_density_kg_m2.status==S::ok,"native LOS surface density");
  return numerator/projection.surface_density_kg_m2.value;
}
int main() {
  PlummerSphere source{1,1,1};
  std::array<double,7> radii{0,0.1,0.5,std::sqrt(0.5),1,2,10};
  auto population=prepare_plummer_population(source,radii);need(population.status()==S::ok,"prepare required native fields");
  need(population.native_evaluations()==2*radii.size(),"prepare native count");
  need(population.radii()[3].status==S::ok,"irrelevant native tidal zero ignored");
  for(std::size_t i=0;i<radii.size();++i) {
    const auto& state=population.radii()[i];
    auto a=velocity_integrals(population,i,256),b=velocity_integrals(population,i,512);
    close(a.normalization,b.normalization,1e-10L);close(a.density,b.density,1e-10L);
    close(a.second,b.second,1e-10L);close(a.fourth,b.fourth,1e-10L);
    close(b.normalization,1);close(b.density,state.density_kg_m3.value);
    close(b.second,3*W(state.one_axis_variance_m2_s2.value));
    close(b.fourth,5*W(state.relative_potential_m2_s2.value)*state.relative_potential_m2_s2.value/14);
    need(b.error<1e-10L,"integrated point arithmetic allocation");
    W los1=projected_variance(source,radii[i],128),los2=projected_variance(source,radii[i],256);
    close(los1,los2,1e-10L);close(los2,state.projected_los_variance_m2_s2.value);
    if(radii[i]>0) {
      double dr=radii[i]*1e-5;
      std::array<double,2> neighbors{radii[i]-dr,radii[i]+dr};
      auto p=prepare_plummer_population(source,neighbors);need(p.status()==S::ok,"Jeans neighbors");
      W derivative=(W(p.radii()[1].density_kg_m3.value)*p.radii()[1].one_axis_variance_m2_s2.value-
                    W(p.radii()[0].density_kg_m3.value)*p.radii()[0].one_axis_variance_m2_s2.value)/(2*dr);
      // Shell-force RHS uses retained native enclosed mass, not copied gravity.
      W rhs=-W(state.density_kg_m3.value)*source.gravitational_coupling_m3_kg_s2*
            state.enclosed_mass_kg.value/(W(radii[i])*radii[i]);
      close(derivative,rhs,2e-7L);
    }
  }
  // Literal oracle: 70-/90-digit Eddington Beta inversion, independent density
  // shell potential and LOS integration, agreeing to 1e-50 relative.
  struct Fixture {double radius,speed;std::array<double,8> value;};
  const Fixture fixtures[] = {
{0,0, {1.0,0.23873241463784300365332564505877,0.15637905395236658903782592142322,0.65503904942943572898959330232314,0.0,1.4142135623730950488016887242097,0.16666666666666666666666666666667,0.14726215563702155805293640859123}},
{0,.5, {1.0,0.23873241463784300365332564505877,0.097995643031569557850877044596647,0.41048318964238231688212316599607,1.2895709730026141938876520656777,1.4142135623730950488016887242097,0.16666666666666666666666666666667,0.14726215563702155805293640859123}},
{0,1.3, {1.0,0.23873241463784300365332564505877,0.00022926593196871197896384465614608,0.00096034689012176383396118352477473,0.020395046647953322436571024144883,1.4142135623730950488016887242097,0.16666666666666666666666666666667,0.14726215563702155805293640859123}},
{1,0, {0.70710678118654752440084436210485,0.042202327319864347010399962078441,0.046491770899387184152771876709998,1.1016399770328264661475853586497,0.0,1.1892071150027210667174999705605,0.11785113019775792073347406035081,0.10413006886308670891443471070455}},
{1,.5, {0.70710678118654752440084436210485,0.042202327319864347010399962078441,0.023533573934612863874039731432654,0.55763687524255970306230988277142,1.7518679106327936145899126622051,1.1892071150027210667174999705605,0.11785113019775792073347406035081,0.10413006886308670891443471070455}},
{1,1, {0.70710678118654752440084436210485,0.042202327319864347010399962078441,0.00063220612331588765931031372179597,0.014980361592956807343810912347914,0.18824877571420719389961682328399,1.1892071150027210667174999705605,0.11785113019775792073347406035081,0.10413006886308670891443471070455}},
{10,0, {0.099503719020998913566527375373857,0.0000023286700428711438429707801611839,0.000048597751651255802860604985707646,20.869316286362748883910041820756,0.0,0.44610249723802020780484572295978,0.016583953170166485594421229228976,0.014653132156932804387283889326947}},
{10,.1, {0.099503719020998913566527375373857,0.0000023286700428711438429707801611839,0.000040574178028113543669618491650785,17.423755740889522782170012504592,2.1895337213408603912697094651487,0.44610249723802020780484572295978,0.016583953170166485594421229228976,0.014653132156932804387283889326947}},
{10,.3, {0.099503719020998913566527375373857,0.0000023286700428711438429707801611839,0.000005911123881751578122035322702301,2.5384119574379169062157219622748,2.8708762886177207994437071831241,0.44610249723802020780484572295978,0.016583953170166485594421229228976,0.014653132156932804387283889326947}},
  };
  for(auto fixture:fixtures) {
    std::array<double,1> radius{fixture.radius};
    auto p=prepare_plummer_population(source,radius);need(p.status()==S::ok,"fixture preparation");
    std::array<PlummerVelocityRequest,1> request{{{0,fixture.speed}}};
    auto b=p.evaluate(request);need(b.status==S::ok,"fixture velocity");
    const auto& state=p.radii()[0];const auto& row=b.rows[0];
    std::array values{state.relative_potential_m2_s2,state.density_kg_m3,
      row.distribution_function_kg_s3_m6,row.vector_velocity_density_s3_m3,
      row.speed_density_s_m,state.escape_speed_m_s,state.one_axis_variance_m2_s2,
      state.projected_los_variance_m2_s2};
    for(unsigned i=0;i<values.size();++i) {
      if(fixture.value[i]==0) need(values[i].value==0,"fixture exact speed zero");
      else close(values[i].value,fixture.value[i],1e-13L);
      need(std::abs(values[i].value-fixture.value[i])<=values[i].error_estimate+
           2*std::numeric_limits<double>::epsilon()*std::abs(fixture.value[i]),"fixture within inherited arithmetic diagnostic");
    }
  }
  std::array<double,1> base_radius{1},scaled_radius{5};
  auto original=prepare_plummer_population(source,base_radius);
  auto scaled=prepare_plummer_population({3,5,7},scaled_radius);
  need(original.status()==S::ok && scaled.status()==S::ok,"SI rescaled preparation");
  const auto& a=original.radii()[0];const auto& b=scaled.radii()[0];
  W speed_scale=std::sqrt(21.L/5),phase_scale=(3.L/125)/std::pow(speed_scale,3);
  close(b.density_kg_m3.value,a.density_kg_m3.value*3.L/125,1e-13L);
  close(b.projected_surface_density_kg_m2.value,a.projected_surface_density_kg_m2.value*3.L/25,1e-13L);
  close(b.enclosed_mass_kg.value,a.enclosed_mass_kg.value*3,1e-13L);
  close(b.relative_potential_m2_s2.value,a.relative_potential_m2_s2.value*21.L/5,1e-13L);
  close(b.escape_speed_m_s.value,a.escape_speed_m_s.value*speed_scale,1e-13L);
  close(b.one_axis_variance_m2_s2.value,a.one_axis_variance_m2_s2.value*21.L/5,1e-13L);
  close(b.projected_los_variance_m2_s2.value,a.projected_los_variance_m2_s2.value*21.L/5,1e-13L);
  std::array<PlummerVelocityRequest,1> base_request{{{0,.5}}},scaled_request{{{0,double(.5L*speed_scale)}}};
  auto old=original.evaluate(base_request),rescaled=scaled.evaluate(scaled_request);
  need(old.status==S::ok && rescaled.status==S::ok,"SI rescaled velocity law");
  close(rescaled.rows[0].distribution_function_kg_s3_m6.value,old.rows[0].distribution_function_kg_s3_m6.value*phase_scale,1e-13L);
  close(rescaled.rows[0].vector_velocity_density_s3_m3.value,old.rows[0].vector_velocity_density_s3_m3.value/std::pow(speed_scale,3),1e-13L);
  close(rescaled.rows[0].speed_density_s_m.value,old.rows[0].speed_density_s_m.value/speed_scale,1e-13L);
  // Tiny but normal native spatial density: bound DF underflows while the
  // normalized PDFs remain representable. No positive DF becomes support zero.
  std::array<double,1> at_center{0};
  auto low_mass=prepare_plummer_population({1e-300,1,1e300},at_center);
  need(low_mass.status()==S::ok,"small mass native state");
  std::array<PlummerVelocityRequest,1> near_escape{{{0,low_mass.radii()[0].escape_speed_m_s.value*.999}}};
  auto underflow=low_mass.evaluate(near_escape);
  need(underflow.rows[0].support==PlummerEnergySupport::bound && underflow.rows[0].distribution_function_kg_s3_m6.status==S::outside_domain,"positive bound DF underflow refused");
  need(underflow.rows[0].vector_velocity_density_s3_m3.status==S::ok && underflow.rows[0].vector_velocity_density_s3_m3.value>0,"normalized velocity law retains positive value");
  const auto& center=population.radii()[0];
  std::vector<PlummerVelocityRequest> requests{{0,0},{0,0.5},{0,center.escape_speed_m_s.value},
      {0,std::nextafter(center.escape_speed_m_s.value,0.0)},
      {0,std::nextafter(center.escape_speed_m_s.value,std::numeric_limits<double>::infinity())},
      {0,center.escape_speed_m_s.value*1.001},{0,1e308},{100,1},{0,-1},{0,std::numeric_limits<double>::quiet_NaN()}};
  auto batch=population.evaluate(requests);
  need(batch.rows.size()==requests.size(),"preserve every input row");
  need(batch.rows[0].status==S::ok && batch.rows[0].speed_density_s_m.value==0,"zero speed has zero speed density");
  need(batch.rows[0].distribution_function_kg_s3_m6.value>0 && batch.rows[0].vector_velocity_density_s3_m3.value>0,"zero-speed positive DF/vector density");
  for(unsigned i:{2u,3u,4u}) {
    need(batch.rows[i].status==S::conditioning_budget_exceeded && batch.rows[i].support==PlummerEnergySupport::ambiguous,"rounded escape unresolved support");
    need(std::abs(batch.rows[i].binding_energy_m2_s2)<=batch.rows[i].binding_energy_absolute_error_m2_s2,"ambiguous energy diagnostics");
  }
  for(unsigned i:{5u,6u}) {
    need(batch.rows[i].status==S::ok && batch.rows[i].support==PlummerEnergySupport::outside,"outside empirical energy interval");
    need(batch.rows[i].distribution_function_kg_s3_m6.value==0 && batch.rows[i].vector_velocity_density_s3_m3.value==0 && batch.rows[i].speed_density_s_m.value==0,"outside support exact zeros");
  }
  for(unsigned i:{7u,8u,9u}) need(batch.rows[i].status==S::invalid_input,"invalid velocity row retained");
  std::array<PlummerVelocityRequest,1> single{{{0,0.5}}};
  auto strict=PlummerPopulationEvaluationPolicy{};strict.relative_tolerance=1e-30;
  need(population.evaluate(single,strict).status==S::conditioning_budget_exceeded,"tight velocity budget refusal");
  auto tiny=PlummerPopulationPreparationPolicy{};tiny.relative_tolerance=std::numeric_limits<double>::denorm_min();
  need(prepare_plummer_population(source,radii,tiny).status()==S::conditioning_budget_exceeded,"native allocation underflow refusal");
  auto prep=PlummerPopulationPreparationPolicy{};
  prep.maximum_radii=radii.size()-1;need(prepare_plummer_population(source,radii,prep).status()==S::work_limit,"radius cap");
  prep={};prep.maximum_native_evaluations=2*radii.size()-1;need(prepare_plummer_population(source,radii,prep).status()==S::work_limit,"native cap");
  prep={};prep.maximum_payload_bytes=*plummer_population_payload_bound(radii.size())-1;
  need(prepare_plummer_population(source,radii,prep).status()==S::work_limit,"prepare byte cap");
  prep.maximum_payload_bytes+=1;need(prepare_plummer_population(source,radii,prep).status()==S::ok,"exact prepare byte cap");
  auto eval=PlummerPopulationEvaluationPolicy{};eval.maximum_rows=single.size()-1;
  need(population.evaluate(single,eval).status==S::invalid_input,"zero policy row cap invalid");
  std::array<PlummerVelocityRequest,2> pair{{{0,.5},{0,.25}}};
  eval={};eval.maximum_rows=1;need(population.evaluate(pair,eval).status==S::work_limit,"nonzero row cap");
  eval={};eval.maximum_velocity_evaluations=1;need(population.evaluate(pair,eval).status==S::work_limit,"velocity work cap");
  eval={};eval.maximum_payload_bytes=*plummer_population_batch_payload_bound(single.size())-1;
  need(population.evaluate(single,eval).status==S::work_limit,"batch byte cap");
  eval.maximum_payload_bytes+=1;need(population.evaluate(single,eval).status==S::ok,"exact batch byte cap");
  need(!plummer_population_payload_bound(std::numeric_limits<std::size_t>::max()),"prepare byte overflow");
  need(!plummer_population_batch_payload_bound(std::numeric_limits<std::size_t>::max()),"batch byte overflow");
  std::array<double,2> invalid{0,-1};auto failed=prepare_plummer_population(source,invalid);
  need(failed.status()==S::invalid_input && failed.radii().size()==2 && failed.radii()[0].status==S::ok && failed.radii()[1].status==S::invalid_input,"failed preparation preserves radius diagnostics");
  need(failed.evaluate(single).status==S::invalid_input && failed.evaluate(single).rows.empty(),"failed preparation withholds batch");
  need(prepare_plummer_population({-1,1,1},radii).status()==S::invalid_input,"invalid source");
  need(prepare_plummer_population(source,{}).status()==S::invalid_input,"empty radius grid");
  auto copy=population;auto moved=std::move(copy);
  need(moved.status()==S::ok && copy.status()==S::invalid_input && copy.source()==nullptr,"move invalidates source owner");
  auto* alias=&moved;moved=std::move(*alias);need(moved.evaluate(single).status==S::ok,"self move keeps owner");
  radii[0]=999;source.total_mass_kg=999;
  need(population.radii()[0].radius_metres==0 && population.source()->total_mass_kg==1,"owns radius and source inputs");
  need(population.evaluate({}).status==S::ok,"empty velocity batch");
  std::fesetround(FE_UPWARD);need(population.evaluate(single).status==S::invalid_input,"runtime arithmetic mode refusal");std::fesetround(FE_TONEAREST);
}
