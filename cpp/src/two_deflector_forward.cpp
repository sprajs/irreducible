#include "irred/two_deflector_forward.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
#include <new>
#include <utility>

namespace irred::lensing {
namespace {
using S = numerics::Status;
using W = long double;
constexpr std::size_t hard_pixels = 4096, hard_samples = 8000000;
constexpr std::size_t hard_row_samples = 500000, hard_bytes = 16 * 1024 * 1024;
constexpr std::size_t hard_origin = 65536;
constexpr W eps = std::numeric_limits<W>::epsilon();
struct Point { W x = 0, y = 0; };
using Matrix = std::array<W, 4>;
W norm(Point p) noexcept { return std::hypot(p.x, p.y); }
Point subtract(Point a, Point b) noexcept { return {a.x-b.x, a.y-b.y}; }
Point multiply(const Matrix &q, Point v) noexcept {
  return {q[0]*v.x+q[1]*v.y, q[2]*v.x+q[3]*v.y};
}
Point point(ForwardPoint2 p) noexcept { return {p.x,p.y}; }
bool finite(ForwardPoint2 p) noexcept {
  return std::isfinite(p.x) && std::isfinite(p.y);
}
bool arithmetic() noexcept {
  return std::numeric_limits<double>::is_iec559 &&
         std::numeric_limits<double>::digits == 53 &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >
             std::numeric_limits<double>::max_exponent &&
         std::fegetround() == FE_TONEAREST;
}
bool within(double v, W low, W high) noexcept { return v >= low && v <= high; }
bool center_domain(ForwardPoint2 p, W bound) noexcept {
  return std::abs(p.x) <= bound && std::abs(p.y) <= bound;
}
Matrix rotate_diagonal(W angle, W first, W second) noexcept {
  const W c=std::cos(angle), s=std::sin(angle), off=(first-second)*c*s;
  return {first*c*c+second*s*s, off, off, first*s*s+second*c*c};
}
// Outward conversion of a nonnegative empirical diagnostic. The primitives
// themselves remain empirical; nextafter does not turn them into certificates.
double diagnostic(W x) noexcept {
  if (!(x >= 0) || !std::isfinite(x) ||
      x > std::numeric_limits<double>::max())
    return std::numeric_limits<double>::infinity();
  if (x == 0) return 0;
  const double d=static_cast<double>(x);
  return std::nextafter(d,std::numeric_limits<double>::infinity());
}
bool increment(std::size_t &a, std::size_t n=1) noexcept {
  if (n > std::numeric_limits<std::size_t>::max()-a) return false;
  a+=n; return true;
}
bool add_work(ForwardWork &a,const ForwardWork &b) noexcept {
#define IRRED_FORWARD_ADD(field) if (!increment(a.field,b.field)) return false
  IRRED_FORWARD_ADD(preparation_started); IRRED_FORWARD_ADD(pixels_started);
  IRRED_FORWARD_ADD(index_comparisons_started); IRRED_FORWARD_ADD(psf_preflight_started);
  IRRED_FORWARD_ADD(field_samples_started); IRRED_FORWARD_ADD(deflector_evaluations_started);
  IRRED_FORWARD_ADD(outer_integrations_started); IRRED_FORWARD_ADD(inner_integrations_started);
  IRRED_FORWARD_ADD(outer_callbacks_started); IRRED_FORWARD_ADD(inner_callbacks_started);
  IRRED_FORWARD_ADD(failed_starts);
#undef IRRED_FORWARD_ADD
  return true;
}
struct Context {
  const TwoDeflectorScene &scene;
  const AffineDetectorCutout &cutout;
  const std::array<Matrix,3> &derived;
  std::span<const ShiftPSFComponent> psf;
  const ForwardPixelRectangle &pixel;
  const ForwardPixelPolicy &policy;
  ForwardWork &batch_work;
  ForwardWork &row_work;
  S cause=S::ok;
  W x=0, maximum_inner_error=0, maximum_field_error=0, latest_field_error=0;
  double inner_tolerance=0;

  void fail(S s) noexcept { if(cause==S::ok) { cause=s; increment(row_work.failed_starts); } }
  bool count(std::size_t ForwardWork::*member) noexcept {
    if(!increment(row_work.*member)) { fail(S::overflow); return false; }
    return true;
  }
  bool begin_field() noexcept {
    if (cause!=S::ok) return false;
    if (batch_work.field_samples_started >= policy.maximum_field_samples ||
        row_work.field_samples_started >= policy.maximum_field_samples_per_pixel) {
      fail(S::work_limit); return false;
    }
    if(!increment(batch_work.field_samples_started) ||
       !increment(row_work.field_samples_started)) { fail(S::overflow); return false; }
    return true;
  }
  Point detector(W xx,W yy,ForwardPoint2 shift={}) const noexcept {
    return {cutout.origin.x+W(cutout.matrix[0])*xx+W(cutout.matrix[1])*yy-shift.x,
            cutout.origin.y+W(cutout.matrix[2])*xx+W(cutout.matrix[3])*yy-shift.y};
  }
  Point map(Point u,W &map_error) noexcept {
    Point deflection{}; W absolute=1+norm(u);
    for(std::size_t j=0;j<2;++j) {
      if(!count(&ForwardWork::deflector_evaluations_started)) return {};
      const auto &lens=scene.deflectors[j];
      const auto v=subtract(u,point(lens.center));
      const auto p=multiply(derived[j],v);
      const W d2=W(lens.core)*lens.core+v.x*p.x+v.y*p.y;
      if (!(d2>0) || !std::isfinite(d2)) { fail(S::overflow); return {}; }
      const W f=W(lens.strength)/std::sqrt(d2);
      const Point alpha{f*p.x,f*p.y};
      deflection.x+=alpha.x; deflection.y+=alpha.y; absolute+=norm(alpha);
    }
    const W l=scene.mass_scale_lambda;
    Point b{l*((1-W(scene.shear_1))*u.x-W(scene.shear_2)*u.y-deflection.x),
            l*((1+W(scene.shear_1))*u.y-W(scene.shear_2)*u.x-deflection.y)};
    if(!std::isfinite(b.x)||!std::isfinite(b.y)) fail(S::overflow);
    // Empirical operation/primitive allowance, including prepared rotations,
    // radius, affine/map cancellation and cast-independent wide accumulation.
    // 512 exceeds the per-map arithmetic count; it is not a libm proof.
    map_error=512*eps*l*absolute;
    return b;
  }
  W field(Point u) noexcept {
    if(!begin_field()) return 0;
    W map_error=0; const auto b=map(u,map_error);
    if(cause!=S::ok) return 0;
    const auto r=subtract(b,point(scene.source.center));
    const auto p=multiply(derived[2],r);
    const W quadratic=r.x*p.x+r.y*p.y;
    if(!(quadratic>=0)||!std::isfinite(quadratic)) { fail(S::overflow); return 0; }
    const W intensity=W(scene.source.peak_electrons_per_second_per_radian_squared)*
                       std::exp(-quadratic/2);
    const double cast=static_cast<double>(intensity);
    if(!(intensity>0)||!std::isfinite(intensity)||!std::isfinite(cast)||cast==0) {
      fail(S::overflow); return 0;
    }
    const W pnorm=1/(W(scene.source.minor_width)*scene.source.minor_width);
    const W exponent_error=128*eps*(1+quadratic)+
      pnorm*(norm(r)*map_error+map_error*map_error/2);
    const W error=intensity*(std::expm1(exponent_error)+64*eps)+
                  std::abs(intensity-W(cast));
    if(!(error>=0)||!std::isfinite(error)) { fail(S::overflow); return 0; }
    latest_field_error=error;
    return intensity;
  }
  bool preflight() noexcept {
    const W cx=(W(pixel.x_low)+pixel.x_high)/2,
            cy=(W(pixel.y_low)+pixel.y_high)/2;
    for(const auto &h:psf) {
      if(!count(&ForwardWork::psf_preflight_started)) return false;
      const auto u=detector(cx,cy,h.shift); W radius=0;
      for(W xx:{W(pixel.x_low),W(pixel.x_high)})
        for(W yy:{W(pixel.y_low),W(pixel.y_high)}) {
          const auto corner=detector(xx,yy,h.shift);
          const W affine_error=64*eps*(1+std::abs(cutout.origin.x)+
              std::abs(cutout.origin.y)+std::abs(W(cutout.matrix[0])*xx)+
              std::abs(W(cutout.matrix[1])*yy)+std::abs(W(cutout.matrix[2])*xx)+
              std::abs(W(cutout.matrix[3])*yy)+norm(point(h.shift)));
          if(!std::isfinite(corner.x)||!std::isfinite(corner.y) ||
             std::abs(corner.x)+affine_error>8 ||
             std::abs(corner.y)+affine_error>8) { fail(S::outside_domain); return false; }
          radius=std::max(radius,norm(subtract(corner,u))+affine_error);
        }
      // Structural zero omits unused field work, while all footprints/states
      // and all source/model fields remain mandatory and validated.
      if(scene.source.peak_electrons_per_second_per_radian_squared==0) continue;
      if(!begin_field()) return false;
      W map_error=0; const auto b=map(u,map_error);
      if(cause!=S::ok) return false;
      W lipschitz=1+std::hypot(W(scene.shear_1),W(scene.shear_2));
      for(const auto &lens:scene.deflectors) {
        const W rmin=std::max(W(0),norm(subtract(u,point(lens.center)))-radius);
        const W dmin=std::sqrt(W(lens.core)*lens.core+rmin*rmin);
        lipschitz+=W(lens.strength)/(W(lens.axis_ratio)*lens.axis_ratio*dmin);
      }
      lipschitz*=scene.mass_scale_lambda;
      const W delta=lipschitz*radius+map_error;
      const W variation=(norm(subtract(b,point(scene.source.center)))*delta+
                         delta*delta/2)/(W(scene.source.minor_width)*scene.source.minor_width);
      if(!std::isfinite(variation)||variation*(1+128*eps)>W(.25)) {
        fail(S::outside_domain); return false;
      }
    }
    return true;
  }
};
double inner_callback(double y,const void *opaque) noexcept {
  auto &c=*static_cast<Context *>(const_cast<void *>(opaque));
  if(c.cause!=S::ok) return std::numeric_limits<double>::quiet_NaN();
  if(!c.count(&ForwardWork::inner_callbacks_started))
    return std::numeric_limits<double>::quiet_NaN();
  W value=0, absolute=0, field_error=0;
  for(const auto &h:c.psf) {
    const W field=c.field(c.detector(c.x,y,h.shift));
    if(c.cause!=S::ok) return std::numeric_limits<double>::quiet_NaN();
    const W term=W(h.probability)*field;
    if(!(term>0)||!std::isfinite(term)) { c.fail(S::overflow); return std::numeric_limits<double>::quiet_NaN(); }
    value+=term; absolute+=term;
    field_error+=W(h.probability)*c.latest_field_error;
  }
  const double cast=static_cast<double>(value);
  if(!(cast>0)||!std::isfinite(cast)) { c.fail(S::overflow); return std::numeric_limits<double>::quiet_NaN(); }
  // The positive weights sum to exactly one in admitted wide arithmetic.
  c.maximum_field_error=std::max(c.maximum_field_error,
      field_error+32*eps*absolute+std::abs(value-W(cast)));
  return cast;
}
double outer_callback(double x,const void *opaque) noexcept {
  auto &c=*static_cast<Context *>(const_cast<void *>(opaque));
  if(c.cause!=S::ok) return std::numeric_limits<double>::quiet_NaN();
  if(!c.count(&ForwardWork::outer_callbacks_started)||
     !c.count(&ForwardWork::inner_integrations_started))
    return std::numeric_limits<double>::quiet_NaN();
  c.x=x;
  const auto inner=numerics::integrate(inner_callback,&c,c.pixel.y_low,c.pixel.y_high,
      {c.inner_tolerance,0,std::max(std::size_t(3),c.policy.maximum_field_samples_per_pixel),c.policy.maximum_depth});
  if(c.cause==S::ok&&inner.status!=S::ok) c.fail(inner.status);
  if(c.cause!=S::ok) return std::numeric_limits<double>::quiet_NaN();
  if(!(inner.value>0)||!std::isfinite(inner.error_estimate)||inner.error_estimate<0) {
    c.fail(S::overflow); return std::numeric_limits<double>::quiet_NaN();
  }
  c.maximum_inner_error=std::max(c.maximum_inner_error,W(inner.error_estimate));
  return inner.value;
}
} // namespace

ForwardPixelMeans::ForwardPixelMeans(ForwardPixelMeans &&other) noexcept
    : status(std::exchange(other.status,S::invalid_input)), work(other.work),
      admitted_payload_bytes(std::exchange(other.admitted_payload_bytes,{})),
      rows_(std::move(other.rows_)), row_count_(std::exchange(other.row_count_,0)) {
  other.work={};
}
ForwardPixelMeans &ForwardPixelMeans::operator=(ForwardPixelMeans &&other) noexcept {
  if(this!=&other) {
    status=std::exchange(other.status,S::invalid_input); work=other.work; other.work={};
    admitted_payload_bytes=std::exchange(other.admitted_payload_bytes,{});
    rows_=std::move(other.rows_); row_count_=std::exchange(other.row_count_,0);
  }
  return *this;
}

void PreparedTwoDeflectorForward::swap(PreparedTwoDeflectorForward &other) noexcept {
  using std::swap;
  swap(status_,other.status_); swap(scene_,other.scene_); swap(cutout_,other.cutout_);
  swap(psf_,other.psf_); swap(psf_count_,other.psf_count_);
  swap(derived_,other.derived_); swap(determinant_,other.determinant_); swap(work_,other.work_);
}
PreparedTwoDeflectorForward::PreparedTwoDeflectorForward(PreparedTwoDeflectorForward &&other) noexcept {
  swap(other);
}
PreparedTwoDeflectorForward &PreparedTwoDeflectorForward::operator=(PreparedTwoDeflectorForward &&other) noexcept {
  if(this!=&other) { PreparedTwoDeflectorForward moved(std::move(other)); swap(moved); }
  return *this;
}
std::optional<std::size_t> PreparedTwoDeflectorForward::retained_payload_bound() const noexcept {
  detail::PayloadAccounting bytes(sizeof(*this));
  if(scene_) bytes.string(scene_->origin);
  return bytes.result();
}
std::optional<std::size_t> two_deflector_payload_bound(std::size_t pixels,
                                                     std::size_t capacity) noexcept {
  if(capacity==std::numeric_limits<std::size_t>::max()) return {};
  // Both no-elision owner/result headers plus incoming scene and two string
  // capacities (including conservative moved-from SSO) simultaneously.
  detail::PayloadAccounting bytes(2*sizeof(PreparedTwoDeflectorForward)+
      sizeof(TwoDeflectorScene)+2*sizeof(ForwardPixelMeans));
  bytes.add(2,capacity+1); bytes.add(pixels,sizeof(ForwardPixelMean));
  return bytes.result();
}

PreparedTwoDeflectorForward prepare_two_deflector_forward(
    TwoDeflectorScene scene,const AffineDetectorCutout &cutout,
    std::span<const ShiftPSFComponent> psf,TwoDeflectorPreparationPolicy policy) {
  PreparedTwoDeflectorForward out; out.work_.preparation_started=1;
  auto refuse=[&](S cause) { out.status_=cause; ++out.work_.failed_starts; return std::move(out); };
  if(!arithmetic()) return refuse(S::outside_domain);
  if(policy.maximum_origin_bytes>hard_origin||policy.maximum_native_bytes>hard_bytes)
    return refuse(S::invalid_input);
  if(scene.origin.empty()||psf.empty()||psf.size()>9) return refuse(S::invalid_input);
  if(scene.origin.size()>policy.maximum_origin_bytes) return refuse(S::work_limit);
  const auto payload=two_deflector_payload_bound(0,scene.origin.capacity());
  if(!payload||*payload>policy.maximum_native_bytes) return refuse(S::work_limit);
  const auto &source=scene.source;
  if(!finite(source.center)||!std::isfinite(source.major_width)||!std::isfinite(source.minor_width)||
     !std::isfinite(source.angle_radians)||!std::isfinite(source.peak_electrons_per_second_per_radian_squared)||
     !std::isfinite(scene.shear_1)||!std::isfinite(scene.shear_2)||!std::isfinite(scene.mass_scale_lambda)||
     !std::isfinite(scene.theta_scale_radians)||!std::isfinite(scene.exposure_seconds)||
     !std::isfinite(scene.uniform_sky_electrons_per_second_per_radian_squared)||!finite(cutout.origin))
    return refuse(S::nonfinite_input);
  const W pi=std::numbers::pi_v<W>;
  if(!center_domain(source.center,4)||!within(source.major_width,.02,4)||
     !within(source.minor_width,.02,source.major_width)||!within(source.angle_radians,-pi,pi)||
     !within(source.peak_electrons_per_second_per_radian_squared,0,1e15L)||
     std::hypot(W(scene.shear_1),W(scene.shear_2))>W(.3)||
     !within(scene.mass_scale_lambda,.25,1)||!within(scene.theta_scale_radians,1e-6,1e-4)||
     !(scene.exposure_seconds>0)||scene.exposure_seconds>3600||
     !within(scene.uniform_sky_electrons_per_second_per_radian_squared,0,1e15L)||
     !center_domain(cutout.origin,8)) return refuse(S::outside_domain);
  for(std::size_t j=0;j<2;++j) {
    const auto &lens=scene.deflectors[j];
    if(!finite(lens.center)||!std::isfinite(lens.strength)||!std::isfinite(lens.core)||
       !std::isfinite(lens.axis_ratio)||!std::isfinite(lens.angle_radians))
      return refuse(S::nonfinite_input);
    if(!center_domain(lens.center,4)||!within(lens.strength,0,4)||!within(lens.core,.02,1)||
       !within(lens.axis_ratio,.5,1)||!within(lens.angle_radians,-pi,pi))
      return refuse(S::outside_domain);
    out.derived_[j]=rotate_diagonal(lens.angle_radians,1,1/(W(lens.axis_ratio)*lens.axis_ratio));
  }
  for(double v:cutout.matrix) {
    if(!std::isfinite(v)) return refuse(S::nonfinite_input);
    if(std::abs(v)>16) return refuse(S::outside_domain);
  }
  const auto &a=cutout.matrix;
  out.determinant_=W(a[0])*a[3]-W(a[1])*a[2];
  if(out.determinant_==0) return refuse(S::singular);
  const W trace=W(a[0])*a[0]+W(a[1])*a[1]+W(a[2])*a[2]+W(a[3])*a[3];
  const W largest=(trace+std::sqrt(std::max(W(0),trace*trace-4*out.determinant_*out.determinant_)))/2;
  if(!std::isfinite(largest)||largest*(1+64*eps)>100*std::abs(out.determinant_))
    return refuse(S::outside_domain);
  // Exact finite dyadic law: each supplied mass is an integer multiple of
  // 2^-32. This keeps nine-state normalization exact without a tolerance or
  // a rounded wide sum that could accept a nonunit real mass.
  std::uint64_t total=0;
  for(const auto &h:psf) {
    if(!finite(h.shift)||!std::isfinite(h.probability)) return refuse(S::nonfinite_input);
    if(!center_domain(h.shift,1)||!(h.probability>0)||h.probability>1)
      return refuse(S::outside_domain);
    const W mass=W(h.probability)*4294967296.L;
    if(std::floor(mass)!=mass) return refuse(S::outside_domain);
    total+=static_cast<std::uint64_t>(mass);
  }
  if(total!=4294967296ULL) return refuse(S::invalid_input);
  out.derived_[2]=rotate_diagonal(source.angle_radians,
      1/(W(source.major_width)*source.major_width),1/(W(source.minor_width)*source.minor_width));
  out.cutout_=cutout; out.psf_count_=psf.size();
  std::copy(psf.begin(),psf.end(),out.psf_.begin());
  out.scene_.emplace(std::move(scene));
  const auto actual=out.retained_payload_bound();
  if(!actual||*actual>policy.maximum_native_bytes) {
    out.scene_.reset(); return refuse(S::work_limit);
  }
  out.status_=S::ok; return out;
}

ForwardPixelMeans PreparedTwoDeflectorForward::means(
    std::span<const ForwardPixelRectangle> pixels,ForwardPixelPolicy policy) const {
  ForwardPixelMeans out;
  auto refuse=[&](S cause) { out.status=cause; increment(out.work.failed_starts); return std::move(out); };
  if(status_!=S::ok||!scene_) return refuse(status_);
  if(!arithmetic()) return refuse(S::outside_domain);
  if(!std::isfinite(policy.absolute_tolerance_electrons)||!std::isfinite(policy.relative_tolerance)) {
    return refuse(S::nonfinite_input);
  }
  if(!(policy.absolute_tolerance_electrons>0)||policy.relative_tolerance<0||
     policy.maximum_pixels>hard_pixels||policy.maximum_field_samples>hard_samples||
     policy.maximum_field_samples_per_pixel>hard_row_samples||policy.maximum_depth>20||
     policy.maximum_native_bytes>hard_bytes) return refuse(S::invalid_input);
  if(pixels.empty()) return refuse(S::invalid_input);
  if(pixels.size()>hard_pixels||pixels.size()>policy.maximum_pixels) {
    return refuse(S::work_limit);
  }
  for(std::size_t i=0;i<pixels.size();++i)
    for(std::size_t j=0;j<i;++j) {
      if(!increment(out.work.index_comparisons_started)) return refuse(S::overflow);
      if(pixels[i].original_index==pixels[j].original_index) return refuse(S::invalid_input);
    }
  const auto payload=two_deflector_payload_bound(pixels.size(),scene_->origin.capacity());
  if(!payload||*payload>policy.maximum_native_bytes) return refuse(S::work_limit);
  out.admitted_payload_bytes=payload;
  // A single exact-count array allocation follows admission; no vector growth
  // capacity or copied return buffer. The allocator/cookie are excluded.
  try { out.rows_=std::make_unique<ForwardPixelMean[]>(pixels.size()); }
  catch(const std::bad_alloc &) { return refuse(S::work_limit); }
  out.row_count_=pixels.size(); out.status=S::ok;
  const W prefactor=W(scene_->exposure_seconds)*scene_->theta_scale_radians*
      scene_->theta_scale_radians*std::abs(determinant_);
  for(std::size_t i=0;i<pixels.size();++i) {
    const auto &p=pixels[i]; auto &r=out.rows_[i]; r.original_index=p.original_index;
    r.work.pixels_started=1;
    Context context{*scene_,cutout_,derived_,{psf_.data(),psf_count_},p,policy,out.work,r.work};
    if(!std::isfinite(p.x_low)||!std::isfinite(p.x_high)||!std::isfinite(p.y_low)||!std::isfinite(p.y_high))
      context.fail(S::nonfinite_input);
    else if(!(p.x_low<p.x_high)||!(p.y_low<p.y_high)) context.fail(S::invalid_input);
    else if(std::abs(p.x_low)>1e9||std::abs(p.x_high)>1e9||
            std::abs(p.y_low)>1e9||std::abs(p.y_high)>1e9) context.fail(S::outside_domain);
    const W width=W(p.x_high)-p.x_low, height=W(p.y_high)-p.y_low, area=width*height;
    if(context.cause==S::ok&&(!(area>0)||!std::isfinite(area)||!(prefactor>0)||!std::isfinite(prefactor)))
      context.fail(S::overflow);
    if(context.cause==S::ok) context.preflight();
    W source_integral=0,inner_error=0,outer_error=0;
    if(context.cause==S::ok&&scene_->source.peak_electrons_per_second_per_radian_squared>0) {
      const W inner_tolerance=.35L*policy.absolute_tolerance_electrons/(prefactor*width);
      const W outer_tolerance=.35L*policy.absolute_tolerance_electrons/prefactor;
      context.inner_tolerance=static_cast<double>(inner_tolerance);
      const double ot=static_cast<double>(outer_tolerance);
      if(!(context.inner_tolerance>0)||!std::isfinite(context.inner_tolerance)||!(ot>0)||!std::isfinite(ot))
        context.fail(S::conditioning_budget_exceeded);
      else {
        context.count(&ForwardWork::outer_integrations_started);
        const auto integral=numerics::integrate(outer_callback,&context,p.x_low,p.x_high,
            {ot,0,std::max(std::size_t(3),policy.maximum_field_samples_per_pixel),policy.maximum_depth});
        if(context.cause==S::ok&&integral.status!=S::ok) context.fail(integral.status);
        if(context.cause==S::ok) {
          source_integral=integral.value;
          inner_error=prefactor*width*context.maximum_inner_error;
          outer_error=prefactor*integral.error_estimate;
        }
      }
    }
    if(context.cause==S::ok) {
      const W sky=scene_->uniform_sky_electrons_per_second_per_radian_squared*area;
      const W mean=prefactor*(source_integral+sky);
      const double cast=static_cast<double>(mean);
      if(!(mean>=0)||!std::isfinite(mean)||!std::isfinite(cast)||(mean>0&&cast==0)) context.fail(S::overflow);
      else {
        r.errors.emplace(); auto &e=*r.errors;
        e.inner_quadrature_electrons=diagnostic(inner_error);
        e.outer_quadrature_electrons=diagnostic(outer_error);
        e.field_psf_arithmetic_electrons=diagnostic(prefactor*area*context.maximum_field_error);
        e.projection_area_electrons=diagnostic(64*eps*std::abs(mean)+std::abs(mean-W(cast)));
        const W total=W(e.inner_quadrature_electrons)+e.outer_quadrature_electrons+
            e.field_psf_arithmetic_electrons+e.projection_area_electrons;
        e.total_electrons=diagnostic(total);
        const W budget=W(policy.absolute_tolerance_electrons)+policy.relative_tolerance*std::abs(W(cast));
        if(!std::isfinite(budget)||!std::isfinite(e.total_electrons)||
           e.inner_quadrature_electrons>.35L*budget||e.outer_quadrature_electrons>.35L*budget||
           e.field_psf_arithmetic_electrons>.20L*budget||e.projection_area_electrons>.05L*budget||
           e.total_electrons>.95L*budget) context.fail(S::conditioning_budget_exceeded);
        else r.electrons=cast;
      }
    }
    r.status=context.cause;
    if(out.status==S::ok&&r.status!=S::ok) out.status=r.status;
    // Field samples were charged to the original whole scope before starts.
    auto row_other=r.work; row_other.field_samples_started=0;
    if(!add_work(out.work,row_other)) { out.status=S::overflow; r.status=S::overflow; r.electrons.reset(); }
  }
  return out;
}
} // namespace irred::lensing
