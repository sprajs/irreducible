// Original finite synthetic controls. Independent positive GL8 integration
// and 8/16-panel refinement are portable empirical controls, not a certified
// interval oracle. The separate high-precision enclosure gate remains open.
#include "irred/two_deflector_forward.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
using namespace irred::lensing;
using S = irred::numerics::Status;
using W = long double;
namespace {
unsigned checks=0;
std::size_t reference_field_starts=0;
constexpr std::size_t reference_cap=8000000;
void need(bool ok,const char *why) {
  ++checks; if(!ok) throw std::runtime_error(why);
}
W value(const ForwardPixelMean &r) {
  need(r.status==S::ok&&r.electrons&&r.errors,"available native mean/diagnostics");
  const auto &e=*r.errors;
  need(std::isfinite(e.total_electrons)&&e.total_electrons>=0&&
       W(e.total_electrons)>=W(e.inner_quadrature_electrons)+e.outer_quadrature_electrons+
            e.field_psf_arithmetic_electrons+e.projection_area_electrons,
       "stored aggregate contains stored contributors");
  return *r.electrons;
}
W budget(W ref) { return 1e-10L+1e-10L*std::abs(ref); }
W portable_reference_error(W coarse,W fine,std::size_t states) {
  // Counted empirical roundoff proposal for all8/16-panel positive GL8 terms
  // and original principal-coordinate primitive paths. This is not a libm
  // or quadrature certificate; the directed interval reference is separate.
  const W operations=512.L*states*64*(8*8+16*16);
  const W gamma=operations*std::numeric_limits<W>::epsilon()/
      (1-operations*std::numeric_limits<W>::epsilon());
  return std::abs(fine-coarse)+gamma*(std::abs(coarse)+std::abs(fine));
}
void compare(W actual,W ref,W reference_error,const char *why) {
  const W total=budget(ref);
  need(reference_error>=0&&reference_error<=.05L*total,"reference within original reservation");
  if(std::abs(actual-ref)+reference_error>total)
    std::cerr<<why<<" native="<<double(actual)<<" ref="<<double(ref)
             <<" ref_error="<<double(reference_error)<<'\n';
  need(std::abs(actual-ref)+reference_error<=total,why);
}
TwoDeflectorScene scene() {
  TwoDeflectorScene s;
  s.deflectors={SoftenedPotentialComponent{{-.3,.1},.6,.3,.8,.125},
                SoftenedPotentialComponent{{.4,-.15},.4,.2,.9,-.25}};
  s.source={{.0625,.125},1.5,1,.25,1e13};
  s.shear_1=.04; s.shear_2=.03; s.theta_scale_radians=1e-5;
  s.exposure_seconds=600; s.uniform_sky_electrons_per_second_per_radian_squared=1e9;
  s.origin="synthetic-fixed-two-potential/Gaussian-source/affine-electron-response";
  return s;
}
AffineDetectorCutout cutout() { return {{0,0},{.02,.004,.002,.015}}; }
std::array<ShiftPSFComponent,3> psf() { return {{{{-.02,0},.25},{{0,0},.5},{{.02,0},.25}}}; }
std::array<ForwardPixelRectangle,16> pixels() {
  std::array<ForwardPixelRectangle,16> p{};
  constexpr double centers[]{-48,-16,16,48};
  for(unsigned iy=0;iy<4;++iy) for(unsigned ix=0;ix<4;++ix) {
    const auto i=4*iy+ix; const auto x=centers[ix],y=centers[iy];
    p[i]={x-.05,x+.05,y-.05,y+.05,i};
  }
  return p;
}
// Independent principal-coordinate equations: no prepared Q/P/source array
// or production map/getter is used as a reference.
W brightness(const TwoDeflectorScene &s,W u,W v) {
  need(reference_field_starts<reference_cap,"original reference cap before next start");
  ++reference_field_starts;
  W alpha_x=0,alpha_y=0;
  for(const auto &d:s.deflectors) {
    const W c=std::cos(W(d.angle_radians)),t=std::sin(W(d.angle_radians));
    const W x=u-d.center.x,y=v-d.center.y;
    const W principal_x=c*x+t*y,principal_y=-t*x+c*y;
    const W q2=W(d.axis_ratio)*d.axis_ratio;
    const W radius=std::sqrt(W(d.core)*d.core+principal_x*principal_x+principal_y*principal_y/q2);
    const W px=d.strength*principal_x/radius,py=d.strength*principal_y/(q2*radius);
    alpha_x+=c*px-t*py; alpha_y+=t*px+c*py;
  }
  const W bx=s.mass_scale_lambda*((1-W(s.shear_1))*u-s.shear_2*v-alpha_x)-s.source.center.x;
  const W by=s.mass_scale_lambda*((1+W(s.shear_1))*v-s.shear_2*u-alpha_y)-s.source.center.y;
  const W c=std::cos(W(s.source.angle_radians)),t=std::sin(W(s.source.angle_radians));
  const W x=(c*bx+t*by)/s.source.major_width,y=(-t*bx+c*by)/s.source.minor_width;
  return s.source.peak_electrons_per_second_per_radian_squared*std::exp(-(x*x+y*y)/2);
}
W reference(const TwoDeflectorScene &s,const AffineDetectorCutout &a,
            std::span<const ShiftPSFComponent> h,const ForwardPixelRectangle &p,unsigned panels) {
  constexpr W nodes[]{.183434642495649804939476142360184L,.525532409916328985817739049189246L,
                     .796666477413626739591553936475831L,.960289856497536231683560868569473L};
  constexpr W weights[]{.362683783378361982965150449277195L,.313706645877887287337962201986601L,
                       .222381034453374470544355994426240L,.101228536290376259152531354309962L};
  const W hx=(W(p.x_high)-p.x_low)/panels,hy=(W(p.y_high)-p.y_low)/panels;
  W sum=0;
  for(unsigned iy=0;iy<panels;++iy) for(unsigned ix=0;ix<panels;++ix)
    for(unsigned j=0;j<4;++j) for(unsigned k=0;k<4;++k)
      for(W sy:{-1.L,1.L}) for(W sx:{-1.L,1.L}) {
        const W x=p.x_low+(ix+.5L)*hx+sx*nodes[k]*hx/2;
        const W y=p.y_low+(iy+.5L)*hy+sy*nodes[j]*hy/2;
        W field=0;
        for(const auto &state:h)
          field+=state.probability*brightness(s,a.origin.x+W(a.matrix[0])*x+W(a.matrix[1])*y-state.shift.x,
                                              a.origin.y+W(a.matrix[2])*x+W(a.matrix[3])*y-state.shift.y);
        sum+=weights[k]*weights[j]*field*hx*hy/4;
      }
  const W area=(W(p.x_high)-p.x_low)*(W(p.y_high)-p.y_low);
  const W det=W(a.matrix[0])*a.matrix[3]-W(a.matrix[1])*a.matrix[2];
  return s.exposure_seconds*W(s.theta_scale_radians)*s.theta_scale_radians*std::abs(det)*
      (sum+s.uniform_sky_electrons_per_second_per_radian_squared*area);
}
void complete_work(const ForwardPixelMeans &r,std::size_t states) {
  ForwardWork sum{};
  for(const auto &row:r.rows()) {
    sum.pixels_started+=row.work.pixels_started;
    sum.psf_preflight_started+=row.work.psf_preflight_started;
    sum.field_samples_started+=row.work.field_samples_started;
    sum.deflector_evaluations_started+=row.work.deflector_evaluations_started;
    sum.inner_callbacks_started+=row.work.inner_callbacks_started;
    sum.outer_callbacks_started+=row.work.outer_callbacks_started;
    sum.inner_integrations_started+=row.work.inner_integrations_started;
    sum.outer_integrations_started+=row.work.outer_integrations_started;
    sum.failed_starts+=row.work.failed_starts;
    if(row.status==S::ok&&row.work.field_samples_started)
      need(row.work.field_samples_started==row.work.psf_preflight_started+
           states*row.work.inner_callbacks_started,"preflight+all PSF field starts");
  }
  need(sum.pixels_started==r.work.pixels_started&&sum.psf_preflight_started==r.work.psf_preflight_started&&
       sum.field_samples_started==r.work.field_samples_started&&
       sum.deflector_evaluations_started==r.work.deflector_evaluations_started&&
       sum.inner_callbacks_started==r.work.inner_callbacks_started&&
       sum.outer_callbacks_started==r.work.outer_callbacks_started&&
       sum.inner_integrations_started==r.work.inner_integrations_started&&
       sum.outer_integrations_started==r.work.outer_integrations_started&&
       sum.failed_starts==r.work.failed_starts,"whole actual work contains every row");
  need(r.work.deflector_evaluations_started==2*r.work.field_samples_started,"both required deflectors count");
}
void nonlinear_and_invariance() {
  const auto s=scene(); const auto a=cutout(); const auto h=psf(); const auto p=pixels();
  auto lens=prepare_two_deflector_forward(s,a,h);
  need(lens.status()==S::ok&&lens.preparation_work().preparation_started==1,"retained synthetic scene");
  auto base=lens.means(p); need(base.status==S::ok&&base.rows().size()==16,"coarse ordered means");
  complete_work(base,3);
  for(std::size_t i=0;i<16;++i) {
    const W coarse=reference(s,a,h,p[i],8),fine=reference(s,a,h,p[i],16);
    const W error=portable_reference_error(coarse,fine,h.size());
    compare(value(base.rows()[i]),fine,error,"independent true two-offset GL8 pixel");
    need(base.rows()[i].original_index==i,"supplied original order");
  }
  for(double lambda:{.75,.5,.25}) {
    auto transformed=s; transformed.mass_scale_lambda=lambda;
    transformed.source.center.x*=lambda; transformed.source.center.y*=lambda;
    transformed.source.major_width*=lambda; transformed.source.minor_width*=lambda;
    auto retained=prepare_two_deflector_forward(std::move(transformed),a,h);
    auto result=retained.means(p); need(result.status==S::ok,"MST retains admitted domain");
    for(std::size_t i=0;i<16;++i) compare(value(result.rows()[i]),value(base.rows()[i]),0,"full PSF MST invariance");
  }
  auto swapped=s; std::swap(swapped.deflectors[0],swapped.deflectors[1]);
  auto exchange=prepare_two_deflector_forward(std::move(swapped),a,h).means(p);
  need(exchange.status==S::ok,"component exchange");
  auto doubled=s; doubled.exposure_seconds=1200;
  auto twice=prepare_two_deflector_forward(std::move(doubled),a,h).means(p);
  need(twice.status==S::ok,"full exposure doubled");
  auto source=s; source.uniform_sky_electrons_per_second_per_radian_squared=0;
  auto source_means=prepare_two_deflector_forward(source,a,h).means(p);
  source.source.peak_electrons_per_second_per_radian_squared*=2;
  auto doubled_source=prepare_two_deflector_forward(std::move(source),a,h).means(p);
  need(source_means.status==S::ok&&doubled_source.status==S::ok,"source-amplitude conditional means");
  auto reflected=s; reflected.source.center.y=-reflected.source.center.y;
  reflected.source.angle_radians=-reflected.source.angle_radians; reflected.shear_2=-reflected.shear_2;
  for(auto &d:reflected.deflectors) { d.center.y=-d.center.y; d.angle_radians=-d.angle_radians; }
  auto reflected_a=a; reflected_a.origin.y=-reflected_a.origin.y;
  reflected_a.matrix[2]=-reflected_a.matrix[2]; reflected_a.matrix[3]=-reflected_a.matrix[3];
  auto reflected_h=h; for(auto &state:reflected_h) state.shift.y=-state.shift.y;
  auto mirror=prepare_two_deflector_forward(std::move(reflected),reflected_a,reflected_h).means(p);
  need(mirror.status==S::ok,"explicit complete reflection request");
  std::array<ForwardPixelMeans,3> singleton;
  for(std::size_t k=0;k<3;++k) {
    const std::array<ShiftPSFComponent,1> single{{{h[k].shift,1}}};
    singleton[k]=prepare_two_deflector_forward(s,a,single).means(p);
    need(singleton[k].status==S::ok,"normalized singleton PSF mass1");
  }
  for(std::size_t i=0;i<16;++i) {
    compare(value(exchange.rows()[i]),value(base.rows()[i]),0,"deflector-label invariance");
    compare(value(twice.rows()[i]),2*value(base.rows()[i]),0,"full exposure linearity");
    compare(value(doubled_source.rows()[i]),2*value(source_means.rows()[i]),0,"source peak linearity");
    compare(value(mirror.rows()[i]),value(base.rows()[i]),0,"complete reflection invariance");
    W mixed=0; for(std::size_t k=0;k<3;++k) mixed+=h[k].probability*value(singleton[k].rows()[i]);
    compare(value(base.rows()[i]),mixed,0,"original-weight PSF linearity");
  }
}
void analytic_limits() {
  auto s=scene(); for(auto &d:s.deflectors) d.strength=0;
  s.shear_1=s.shear_2=0; s.source.angle_radians=0;
  const AffineDetectorCutout a{{0,0},{.02,0,0,.015}}; const auto h=psf();
  const std::array<ForwardPixelRectangle,4> p{{{-16.05,-15.95,-16.05,-15.95,9},
      {15.95,16.05,-16.05,-15.95,2},{-16.05,-15.95,15.95,16.05,7},{15.95,16.05,15.95,16.05,1}}};
  auto means=prepare_two_deflector_forward(s,a,h).means(p);
  need(means.status==S::ok,"no-deflection diagonal Gaussian request");
  const W total_source=s.exposure_seconds*W(s.theta_scale_radians)*s.theta_scale_radians*
      s.source.peak_electrons_per_second_per_radian_squared*2*std::numbers::pi_v<W>*
      s.source.major_width*s.source.minor_width;
  for(std::size_t i=0;i<p.size();++i) {
    W probability=0,error=0;
    for(const auto &state:h) {
      W px[2]{},py[2]{};
      for(unsigned k=0;k<2;++k) {
        px[k]=std::erf((W(a.matrix[0])*(k?p[i].x_high:p[i].x_low)-state.shift.x-s.source.center.x)/
                       (std::sqrt(2.L)*s.source.major_width));
        py[k]=std::erf((W(a.matrix[3])*(k?p[i].y_high:p[i].y_low)-state.shift.y-s.source.center.y)/
                       (std::sqrt(2.L)*s.source.minor_width));
      }
      const W dx=(px[1]-px[0])/2,dy=(py[1]-py[0])/2;
      probability+=state.probability*dx*dy;
      error+=state.probability*128*std::numeric_limits<W>::epsilon()*
          ((std::abs(px[0])+std::abs(px[1]))*std::abs(dy)+
           (std::abs(py[0])+std::abs(py[1]))*std::abs(dx));
    }
    const W area=(W(p[i].x_high)-p[i].x_low)*(W(p[i].y_high)-p[i].y_low);
    const W sky=s.exposure_seconds*W(s.theta_scale_radians)*s.theta_scale_radians*
        a.matrix[0]*a.matrix[3]*s.uniform_sky_electrons_per_second_per_radian_squared*area;
    compare(value(means.rows()[i]),total_source*probability+sky,total_source*error,"independent no-deflection CDF pixel");
    need(means.rows()[i].original_index==p[i].original_index,"nonmonotone IDs preserve supplied order");
  }
  // The identical-Q, same-core, zero-shear, Lambda1 centered-source limit
  // has FIVE stationary points, including its central one; not exactly4.
  auto five=scene(); five.deflectors={SoftenedPotentialComponent{{0,0},.7,.2,.8,0},
                                    SoftenedPotentialComponent{{0,0},.3,.2,.8,0}};
  five.shear_1=five.shear_2=0; five.source={{0,0},.1,.1,0,1e13};
  five.uniform_sky_electrons_per_second_per_radian_squared=0;
  const W strength=W(five.deflectors[0].strength)+five.deflectors[1].strength;
  const W core=five.deflectors[0].core,q=five.deflectors[0].axis_ratio;
  const double x=static_cast<double>(std::sqrt(strength*strength-core*core));
  const double y=static_cast<double>(std::sqrt(strength*strength/(q*q)-core*core*q*q));
  const std::array<ForwardPoint2,5> centers{{{0,0},{x,0},{-x,0},{0,y},{0,-y}}};
  std::array<ForwardPixelRectangle,5> controls{};
  for(std::size_t i=0;i<5;++i) {
    const auto c=centers[i]; const W radius=std::sqrt(core*core+W(c.x)*c.x+W(c.y)*c.y/(q*q));
    need(std::hypot(W(c.x)*(1-strength/radius),W(c.y)*(1-strength/(q*q*radius)))<=1e-14L,
         "rounded analytic coordinates retain small stationarity residual");
    controls[i]={c.x-.00005,c.x+.00005,c.y-.00005,c.y+.00005,i};
  }
  const std::array<ShiftPSFComponent,1> single{{{{0,0},1}}};
  const AffineDetectorCutout identity;
  auto all=prepare_two_deflector_forward(five,identity,single).means(controls);
  need(all.status==S::ok&&all.rows().size()==5,"central and four outer brightness pixels retained");
  for(std::size_t i=0;i<5;++i) {
    const W coarse=reference(five,identity,single,controls[i],8),fine=reference(five,identity,single,controls[i],16);
    compare(value(all.rows()[i]),fine,portable_reference_error(coarse,fine,1),
            "independent five-stationary-point extended pixel");
  }
}
void refusal_and_ownership() {
  const auto a=cutout(); const auto h=psf(); const auto p=pixels();
  for(unsigned field=0;field<24;++field) {
    auto nonfinite=scene();
    const std::array<double *,24> fields{{
        &nonfinite.source.center.x,&nonfinite.source.center.y,&nonfinite.source.major_width,
        &nonfinite.source.minor_width,&nonfinite.source.angle_radians,
        &nonfinite.source.peak_electrons_per_second_per_radian_squared,
        &nonfinite.shear_1,&nonfinite.shear_2,&nonfinite.mass_scale_lambda,
        &nonfinite.theta_scale_radians,&nonfinite.exposure_seconds,
        &nonfinite.uniform_sky_electrons_per_second_per_radian_squared,
        &nonfinite.deflectors[0].center.x,&nonfinite.deflectors[0].center.y,
        &nonfinite.deflectors[0].strength,&nonfinite.deflectors[0].core,
        &nonfinite.deflectors[0].axis_ratio,&nonfinite.deflectors[0].angle_radians,
        &nonfinite.deflectors[1].center.x,&nonfinite.deflectors[1].center.y,
        &nonfinite.deflectors[1].strength,&nonfinite.deflectors[1].core,
        &nonfinite.deflectors[1].axis_ratio,&nonfinite.deflectors[1].angle_radians}};
    *fields[field]=std::numeric_limits<double>::quiet_NaN();
    const auto refused=prepare_two_deflector_forward(std::move(nonfinite),a,h);
    need(refused.status()==S::nonfinite_input&&!refused.scene()&&
         refused.preparation_work().failed_starts==1,"every scene scalar finite gate before availability");
  }
  for(unsigned field=0;field<6;++field) {
    auto affine=a;
    const std::array<double *,6> fields{{&affine.origin.x,&affine.origin.y,
        &affine.matrix[0],&affine.matrix[1],&affine.matrix[2],&affine.matrix[3]}};
    *fields[field]=std::numeric_limits<double>::infinity();
    need(prepare_two_deflector_forward(scene(),affine,h).status()==S::nonfinite_input,
         "every affine scalar finite gate");
  }
  for(unsigned field=0;field<9;++field) {
    auto states=h;
    const std::array<double *,9> fields{{&states[0].shift.x,&states[0].shift.y,&states[0].probability,
        &states[1].shift.x,&states[1].shift.y,&states[1].probability,
        &states[2].shift.x,&states[2].shift.y,&states[2].probability}};
    *fields[field]=std::numeric_limits<double>::infinity();
    need(prepare_two_deflector_forward(scene(),a,states).status()==S::nonfinite_input,
         "every required PSF scalar finite gate");
  }
  auto bad=scene(); bad.deflectors[0].core=0;
  need(prepare_two_deflector_forward(std::move(bad),a,h).status()==S::outside_domain,"zero core refused");
  bad=scene(); bad.source.major_width=std::numeric_limits<double>::quiet_NaN();
  need(prepare_two_deflector_forward(std::move(bad),a,h).status()==S::nonfinite_input,"nonfinite source refused");
  bad=scene(); bad.origin.clear();
  need(prepare_two_deflector_forward(std::move(bad),a,h).status()==S::invalid_input,"missing physical provenance refused");
  bad=scene(); bad.origin.assign(65537,'x');
  need(prepare_two_deflector_forward(std::move(bad),a,h).status()==S::work_limit,"original origin byte cap");
  bad=scene(); bad.mass_scale_lambda=1.01;
  need(prepare_two_deflector_forward(std::move(bad),a,h).status()==S::outside_domain,"unsupported negative-sheet scale");
  bad=scene(); bad.source.minor_width=2;
  need(prepare_two_deflector_forward(std::move(bad),a,h).status()==S::outside_domain,"minor exceeds major refused");
  bad=scene(); bad.exposure_seconds=0;
  need(prepare_two_deflector_forward(std::move(bad),a,h).status()==S::outside_domain,"zero full exposure refused");
  const std::array<ShiftPSFComponent,10> too_many_states{};
  need(prepare_two_deflector_forward(scene(),a,{}).status()==S::invalid_input&&
       prepare_two_deflector_forward(scene(),a,too_many_states).status()==S::invalid_input,
       "empty and excessive PSF shapes before copies");
  auto invalid_h=h; invalid_h[0].probability=.125;
  need(prepare_two_deflector_forward(scene(),a,invalid_h).status()==S::invalid_input,"nonunit PSF never normalized");
  const std::array<ShiftPSFComponent,1> third{{{{0,0},1./3}}};
  need(prepare_two_deflector_forward(scene(),a,third).status()==S::outside_domain,"non-dyadic PSF domain refused");
  invalid_h=h; invalid_h[2].probability=0;
  need(prepare_two_deflector_forward(scene(),a,invalid_h).status()==S::outside_domain,"zero required PSF state refused");
  auto singular=a; singular.matrix={1,1,1,1};
  need(prepare_two_deflector_forward(scene(),singular,h).status()==S::singular,"singular affine refused");
  singular.matrix={1,0,0,.0001};
  need(prepare_two_deflector_forward(scene(),singular,h).status()==S::outside_domain,"ill-conditioned affine refused");
  auto input=scene(); const auto prefix=two_deflector_payload_bound(0,input.origin.capacity());
  need(prefix.has_value(),"preparation charge exists");
  TwoDeflectorPreparationPolicy prep_policy; prep_policy.maximum_native_bytes=*prefix-1;
  auto prep_refusal=prepare_two_deflector_forward(std::move(input),a,h,prep_policy);
  need(prep_refusal.status()==S::work_limit&&!prep_refusal.scene()&&
       prep_refusal.preparation_work().preparation_started==1&&prep_refusal.preparation_work().failed_starts==1,
       "preallocation preparation failure receipt");
  auto lens=prepare_two_deflector_forward(scene(),a,h);
  auto empty=lens.means({}); need(empty.status==S::invalid_input&&empty.rows().empty(),"causal empty batch");
  const std::array<ForwardPixelRectangle,4097> too_many_pixels{};
  need(lens.means(too_many_pixels).status==S::work_limit,"original4096 output shape cap before allocation");
  auto repeated=p; repeated[1].original_index=repeated[0].original_index;
  auto duplicate=lens.means(repeated);
  need(duplicate.status==S::invalid_input&&duplicate.rows().empty()&&duplicate.work.index_comparisons_started==1,
       "duplicate index refused before output allocation");
  ForwardPixelPolicy policy; const auto charge=two_deflector_payload_bound(p.size(),lens.scene()->origin.capacity());
  need(charge.has_value(),"evaluation source charge"); policy.maximum_native_bytes=*charge-1;
  auto payload=lens.means(p,policy);
  need(payload.status==S::work_limit&&payload.rows().empty()&&!payload.admitted_payload_bytes,"prefix bytes rejected before allocation");
  policy={}; policy.maximum_field_samples=3;
  auto capped=lens.means(p,policy);
  need(capped.status==S::work_limit&&capped.rows().size()==16&&capped.work.field_samples_started==3,
       "one original whole cap across all attempted rows"); complete_work(capped,3);
  for(const auto &r:capped.rows()) need(!r.electrons&&r.status==S::work_limit,"cap keeps every required failed row");
  policy={}; policy.maximum_field_samples_per_pixel=2;
  auto row_cap=lens.means(std::span(p).first(1),policy);
  need(row_cap.status==S::work_limit&&row_cap.work.field_samples_started==2,"row cap includes preflight states");
  auto malformed=p; malformed[0].x_high=malformed[0].x_low;
  auto structural=lens.means(malformed);
  need(structural.status==S::invalid_input&&!structural.rows()[0].electrons&&structural.rows()[1].electrons,
       "bad row retains later independent required attempt");
  malformed=p; malformed[0].x_low=std::numeric_limits<double>::infinity();
  need(lens.means(malformed).rows()[0].status==S::nonfinite_input,"nonfinite row cause retained");
  malformed=p; malformed[0]={1e8,1e8+1,0,1,0};
  need(lens.means(malformed).rows()[0].status==S::outside_domain,"footprint domain refusal");
  policy={}; policy.absolute_tolerance_electrons=1e-30; policy.relative_tolerance=0; policy.maximum_depth=0;
  auto tight=lens.means(std::span(p).first(1),policy);
  need(tight.status!=S::ok&&!tight.rows()[0].electrons,"original tight budget/refinement refusal not weakened");
  // F01 preserved source counterexample: positive inner values may remain
  // binary64, but the returned outer integral can underflow before sky is
  // added. The positive contribution is refused even when sky hides it.
  auto tiny=scene(); for(auto &d:tiny.deflectors) d.strength=0;
  tiny.shear_1=tiny.shear_2=0; tiny.source={{0,0},2,2,0,1e-310};
  const std::array<ShiftPSFComponent,1> single{{{{0,0},1}}};
  const std::array<ForwardPixelRectangle,1> needle{{{0,1e-20,-.005,.005,17}}};
  const auto positive_underflow=prepare_two_deflector_forward(std::move(tiny),AffineDetectorCutout{},single).means(needle);
  need(positive_underflow.status==S::overflow&&!positive_underflow.rows()[0].electrons&&
       positive_underflow.work.outer_integrations_started==1&&positive_underflow.work.inner_integrations_started>0,
       "positive outer-integral zero never hidden under sky");
  auto broad=scene(); broad.source.minor_width=broad.source.major_width=.02;
  auto unresolved=prepare_two_deflector_forward(std::move(broad),a,h).means(std::span(p).first(1));
  need(unresolved.status==S::outside_domain&&!unresolved.rows()[0].electrons,"unresolved smooth-cell domain refused");
  auto zero=scene(); zero.source.peak_electrons_per_second_per_radian_squared=0;
  zero.uniform_sky_electrons_per_second_per_radian_squared=0;
  auto zero_lens=prepare_two_deflector_forward(zero,a,h); policy={}; policy.maximum_field_samples=0;
  auto exact=zero_lens.means(p,policy);
  need(exact.status==S::ok&&exact.work.field_samples_started==0&&exact.work.outer_integrations_started==0,
       "exact source-zero omits unused physical/quadrature work");
  for(const auto &r:exact.rows()) need(value(r)==0&&r.errors->total_electrons==0,"complete exact-zero law");
  auto invalid_zero=scene(); invalid_zero.source.peak_electrons_per_second_per_radian_squared=0;
  invalid_zero.deflectors[1].core=0;
  need(prepare_two_deflector_forward(std::move(invalid_zero),a,h).status()==S::outside_domain,
       "zero source does not drop invalid required model state");
  zero.uniform_sky_electrons_per_second_per_radian_squared=1e9;
  auto sky=prepare_two_deflector_forward(std::move(zero),a,h).means(p,policy);
  need(sky.status==S::ok&&sky.work.field_samples_started==0,"direct sky with zero work cap");
  const W det=W(a.matrix[0])*a.matrix[3]-W(a.matrix[1])*a.matrix[2];
  for(std::size_t i=0;i<16;++i)
    compare(value(sky.rows()[i]),600*W(1e-5)*1e-5*std::abs(det)*1e9*
      (W(p[i].x_high)-p[i].x_low)*(W(p[i].y_high)-p[i].y_low),0,"exact sky-area convention");
  auto moved=std::move(lens); need(lens.status()!=S::ok&&!lens.scene()&&moved.status()==S::ok,"move invalidates old owner");
  moved=std::move(moved); need(moved.status()==S::ok,"owner self move retains state");
  auto result=moved.means(std::span(p).first(1)); const auto *buffer=result.rows().data();
  auto result_moved=std::move(result);
  need(result.status!=S::ok&&result.rows().empty()&&result_moved.rows().data()==buffer,"return buffer transfer without copy");
  result_moved=std::move(result_moved); need(result_moved.status==S::ok,"result self move retains buffer");
  moved={}; need(value(result_moved.rows()[0])>0,"mean result outlives scene owner");
  need(!two_deflector_payload_bound(1,std::numeric_limits<std::size_t>::max()),"capacity overflow refused");
  const int saved=std::fegetround(); need(std::fesetround(FE_DOWNWARD)==0,"rounding control available");
  const auto arithmetic_refusal=prepare_two_deflector_forward(scene(),a,h).status();
  std::fesetround(saved); need(arithmetic_refusal==S::outside_domain,"unsupported rounding profile refused");
}
} // namespace
int main() {
  static_assert(!std::is_copy_constructible_v<PreparedTwoDeflectorForward>);
  static_assert(!std::is_copy_constructible_v<ForwardPixelMeans>);
  try { nonlinear_and_invariance(); analytic_limits(); refusal_and_ownership();
    std::cout<<"PASS two-deflector synthetic controls checks="<<checks
             <<" independent_GL8_field_starts="<<reference_field_starts
             <<" (empirical refinement; measured/interval gates open)\n";
    return 0;
  } catch(const std::exception &e) { std::cerr<<"FAIL "<<e.what()<<'\n'; return 1; }
}
