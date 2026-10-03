#include "irred/windowed_linear_power.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>

// Requested C++ heap only: excludes malloc users, allocator overhead and RSS.
namespace {
struct alignas(std::max_align_t) Allocation { std::size_t bytes; };
std::size_t heap_live=0, heap_peak=0;
}
void* operator new(std::size_t bytes) {
  if (bytes>std::numeric_limits<std::size_t>::max()-sizeof(Allocation) ||
      bytes>std::numeric_limits<std::size_t>::max()-heap_live) throw std::bad_alloc();
  auto* p=static_cast<Allocation*>(std::malloc(sizeof(Allocation)+(bytes ? bytes : 1)));
  if (!p) throw std::bad_alloc();
  p->bytes=bytes; heap_live+=bytes;
  heap_peak=std::max(heap_peak,heap_live);
  return p+1;
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* p) noexcept {
  if (!p) return;
  auto* h=static_cast<Allocation*>(p)-1; heap_live-=h->bytes; std::free(h);
}
void operator delete[](void* p) noexcept { ::operator delete(p); }
void operator delete(void* p,std::size_t) noexcept { ::operator delete(p); }
void operator delete[](void* p,std::size_t) noexcept { ::operator delete(p); }

namespace {
namespace p=irred::windowed_linear_power;
constexpr std::size_t bytes=16*1024*1024;
const p::PreparationPolicy preparation{2000000,bytes};
p::EvaluationPolicy policy(bool multipoles=true) {
  return {1e-7,1e-6,8000000,2000000,4096,20,bytes,multipoles};
}
p::Source fixture(double amplitude=3000) {
  p::Source s;
  s.identity={"synthetic/dyadic-spectrum","synthetic/identity-window",
              "supplied-matter","supplied-scalar-f","ordered-0-then-2",
              "synthetic-fixed-not-observed"};
  s.spectrum={.32,{.015625,.25},{amplitude,amplitude}};
  s.window.coordinates=p::Coordinates::physical_mpc;
  s.window.theory_k={.03125,.0625,.125};
  s.window.output0_k=s.window.output2_k=s.window.theory_k;
  s.window.output0_ids={10,11,12};s.window.output2_ids={20,21,22};
  s.window.w00=s.window.w22={1,0,0,0,1,0,0,0,1};
  s.window.w02=s.window.w20=std::vector<double>(9,0);
  return s;
}
bool close(double a,long double b,double tolerance=1e-6) {
  return std::isfinite(a) && std::abs(static_cast<long double>(a)-b)<=tolerance;
}
bool empty_failure(const p::Result& r,p::Status status) {
  return r.status==status && r.mean.empty() && r.input_multipoles.empty() &&
      r.mean.capacity()==0 && r.input_multipoles.capacity()==0;
}
#define REQUIRE(x) do { if (!(x)) { std::cerr << "line " << __LINE__ << ": " << #x << '\n'; return 1; } } while (false)
}

int main() {
  static_assert(!std::is_copy_constructible_v<p::Prepared>);
  static_assert(!std::is_copy_assignable_v<p::Prepared>);
  static_assert(std::is_nothrow_move_constructible_v<p::Prepared>);
  const p::ModelPoint standard{.32,2,1,1,1};
  auto source=fixture();
  auto owner=p::prepare(std::move(source),preparation);
  REQUIRE(owner.status()==p::Status::ok && owner.owns_source());
  REQUIRE(owner.retained_payload_bytes() && *owner.retained_payload_bytes()<bytes);
  auto r=owner.evaluate(standard,policy());
  REQUIRE(r.status==p::Status::ok && r.mean.size()==6 && r.input_multipoles.size()==3);
  REQUIRE(r.identity.spectrum=="synthetic/dyadic-spectrum" && !r.observational_qualified);
  for (std::size_t j=0;j<3;++j) {
    REQUIRE(close(r.input_multipoles[j].p0.value,3000.L*83/15));
    REQUIRE(close(r.input_multipoles[j].p2.value,3000.L*68/21));
    REQUIRE(r.mean[j].multipole==0 && r.mean[j].source_row_id==10+j);
    REQUIRE(r.mean[j+3].multipole==2 && r.mean[j+3].source_row_id==20+j);
    REQUIRE(close(r.mean[j].power.value,r.input_multipoles[j].p0.value));
    REQUIRE(close(r.mean[j+3].power.value,r.input_multipoles[j].p2.value));
  }
  REQUIRE(r.work.angular_callbacks>0 && r.work.bracket_comparisons>0);
  REQUIRE(r.work.interpolation_calls==r.work.angular_callbacks);
  REQUIRE(r.work.window_products==36 && r.work.output_writes==9);
  // A negative supplied response coefficient is valid; P2 must stay signed.
  const auto signed_r=owner.evaluate({.32,1,-1,1,1},policy());
  REQUIRE(signed_r.status==p::Status::ok);
  REQUIRE(close(signed_r.mean[0].power.value,3000.L*8/15));
  REQUIRE(close(signed_r.mean[3].power.value,-3000.L*16/21));
  const auto isotropic=owner.evaluate({.32,2,1,2,2},policy());
  REQUIRE(isotropic.status==p::Status::ok);
  REQUIRE(close(isotropic.mean[0].power.value,3000.L*83/120));
  REQUIRE(close(isotropic.mean[3].power.value,3000.L*68/168));
  const auto anisotropic=owner.evaluate({.32,2,0,1.1,.9},policy());
  REQUIRE(anisotropic.status==p::Status::ok);
  REQUIRE(close(anisotropic.mean[0].power.value,12000.L/(1.1L*1.1L*.9L)));
  REQUIRE(close(anisotropic.mean[3].power.value,0));
  auto affine=fixture();
  affine.spectrum.power_mpc3={64,1024}; // exactly P=4096*k on dyadic nodes
  auto affine_owner=p::prepare(std::move(affine),preparation);
  auto affine_r=affine_owner.evaluate({.32,2,1,2,2},policy());
  REQUIRE(affine_r.status==p::Status::ok);
  for (std::size_t j=0;j<3;++j)
    REQUIRE(close(affine_r.mean[j].power.value,4096.L*affine_r.mean[j].source_k*83/(15*16)));

  auto mixed=fixture();
  mixed.window.w00={.5,.5,0,0,.25,.75,-1,0,0};
  mixed.window.w02={.02,-.01,0,0,.01,-.02,0,0,0};
  mixed.window.w20={.03,-.02,0,0,-.01,.02,1,0,0};
  mixed.window.w22={.5,.5,0,0,.25,.75,0,0,0};
  auto mixed_owner=p::prepare(std::move(mixed),preparation);
  const auto mix=mixed_owner.evaluate(standard,policy());
  REQUIRE(mix.status==p::Status::ok);
  const long double mono=3000.L*83/15,quad=3000.L*68/21;
  REQUIRE(close(mix.mean[0].power.value,mono+.01L*quad));
  REQUIRE(close(mix.mean[1].power.value,mono-.01L*quad));
  REQUIRE(close(mix.mean[2].power.value,-mono));
  REQUIRE(close(mix.mean[3].power.value,.01L*mono+quad));
  REQUIRE(close(mix.mean[4].power.value,.01L*mono+quad));
  REQUIRE(close(mix.mean[5].power.value,mono));
  // f=0 cannot remove a quadrupole produced by W20 mixing the monopole.
  const auto mixed_zero_f=mixed_owner.evaluate({.32,2,0,1,1},policy());
  REQUIRE(mixed_zero_f.status==p::Status::ok && mixed_zero_f.mean[5].power.value>0);

  // Unequal column powers expose a unilateral W02/W20 column permutation.
  // The emitted dyadic table is exactly P(k)=4096*k, giving128,256,512.
  auto affine_cross=fixture();
  affine_cross.spectrum.power_mpc3={64,1024};
  affine_cross.window.w02={.125,.5,.25,0,0,0,0,0,0};
  affine_cross.window.w20={.5,-.25,.125,0,0,0,0,0,0};
  auto affine_cross_owner=p::prepare(std::move(affine_cross),preparation);
  const auto affine_cross_result=affine_cross_owner.evaluate(standard,policy());
  REQUIRE(affine_cross_result.status==p::Status::ok);
  const long double c0=83.L/15,c2=68.L/21;
  REQUIRE(close(affine_cross_result.mean[0].power.value,c0*128+c2*272));
  REQUIRE(close(affine_cross_result.mean[3].power.value,c2*128+c0*64));

  auto h_source=fixture();
  constexpr double h=.6711;
  h_source.window.coordinates=p::Coordinates::fixed_h_reference;
  h_source.window.h_reference=h;
  for (auto* axis:{&h_source.window.theory_k,&h_source.window.output0_k,&h_source.window.output2_k})
    for (auto& k:*axis) k/=h;
  auto h_owner=p::prepare(std::move(h_source),preparation);
  const auto h_r=h_owner.evaluate(standard,policy());
  REQUIRE(h_r.status==p::Status::ok && h_r.h_reference==h);
  for (std::size_t j=0;j<6;++j)
    REQUIRE(close(h_r.mean[j].power.value,static_cast<long double>(h)*h*h*r.mean[j].power.value));
  REQUIRE(h_r.mean[0].power.volume_projection_estimate>0);
  // Atomic outputs and acquired metadata survive required late failure.
  auto narrow=policy();narrow.maximum_actions=30;
  const auto limited=owner.evaluate(standard,narrow);
  REQUIRE(empty_failure(limited,p::Status::work_limit));
  REQUIRE(limited.identity.spectrum==r.identity.spectrum);
  REQUIRE(limited.work.angular_callbacks>0 && limited.work.bracket_comparisons>0);
  auto callback_limit=policy();callback_limit.maximum_angular_callbacks=3;
  const auto cb=owner.evaluate(standard,callback_limit);
  REQUIRE(empty_failure(cb,p::Status::work_limit) && cb.work.angular_callbacks==3);
  auto bad_budget=policy();bad_budget.maximum_actions=std::numeric_limits<std::size_t>::max();
  REQUIRE(empty_failure(owner.evaluate(standard,bad_budget),p::Status::invalid_input));
  auto small_payload=policy();small_payload.maximum_live_payload_bytes=1;
  REQUIRE(empty_failure(owner.evaluate(standard,small_payload),p::Status::payload_limit));
  REQUIRE(empty_failure(owner.evaluate({.31,2,1,1,1},policy()),p::Status::incompatible_epoch));
  REQUIRE(empty_failure(owner.evaluate({.32,2,1,0,1},policy()),p::Status::invalid_input));
  REQUIRE(empty_failure(owner.evaluate({.32,2,1,.1,1},policy()),p::Status::outside_support));
  REQUIRE(empty_failure(owner.evaluate({.32,std::numeric_limits<double>::infinity(),1,1,1},policy()),p::Status::nonfinite_input));
  auto exact_zero=owner.evaluate({.32,0,0,1,1},policy());
  REQUIRE(exact_zero.status==p::Status::ok);
  for (const auto& row:exact_zero.mean) REQUIRE(row.power.value==0 && row.power.combined_estimate==0);
  auto zero_source=fixture(0);auto zero_owner=p::prepare(std::move(zero_source),preparation);
  auto zero=zero_owner.evaluate(standard,policy());
  REQUIRE(zero.status==p::Status::ok);
  for (const auto& row:zero.mean) REQUIRE(row.power.value==0 && row.power.combined_estimate==0);
  auto tiny_source=fixture(1e-300);
  auto tiny_owner=p::prepare(std::move(tiny_source),preparation);
  REQUIRE(empty_failure(tiny_owner.evaluate({.32,1e-300,0,1,1},policy()),p::Status::unresolved_projection));
  const int original_rounding=std::fegetround();
  REQUIRE(std::fesetround(FE_UPWARD)==0);
  const auto wrong_rounding=owner.evaluate(standard,policy());
  REQUIRE(std::fesetround(original_rounding)==0);
  REQUIRE(empty_failure(wrong_rounding,p::Status::unsupported_arithmetic));

  auto unsorted=fixture();std::swap(unsorted.window.theory_k[0],unsorted.window.theory_k[1]);
  auto rejected=p::prepare(std::move(unsorted),preparation);
  REQUIRE(rejected.status()==p::Status::invalid_input && !rejected.owns_source());
  REQUIRE(unsorted.spectrum.power_mpc3.size()==2); // not consumed
  auto negative=fixture();negative.spectrum.power_mpc3[0]=-1;
  auto negative_owner=p::prepare(std::move(negative),preparation);
  REQUIRE(negative_owner.status()==p::Status::invalid_input && !negative_owner.owns_source());
  auto bad_shape=fixture();bad_shape.window.w02.pop_back();
  REQUIRE(p::prepare(std::move(bad_shape),preparation).status()==p::Status::invalid_input);
  auto invalid_h=fixture();invalid_h.window.h_reference=.7;
  REQUIRE(p::prepare(std::move(invalid_h),preparation).status()==p::Status::invalid_input);
  auto duplicate=fixture();duplicate.window.output2_ids[1]=20;
  REQUIRE(p::prepare(std::move(duplicate),preparation).status()==p::Status::invalid_input);
  auto prep_limit_source=fixture();
  const auto prep_limited=p::prepare(std::move(prep_limit_source),{10,bytes});
  REQUIRE(prep_limited.status()==p::Status::work_limit && !prep_limited.owns_source());
  REQUIRE(prep_limit_source.window.w00.size()==9);
  auto byte_source=fixture();
  REQUIRE(p::prepare(std::move(byte_source),{2000000,1}).status()==p::Status::payload_limit);
  // Physical-axis endpoint touch is valid; mapping diagnostic touch can refuse.
  auto boundary=fixture();boundary.spectrum.k_per_mpc={.03125,.125};
  auto boundary_owner=p::prepare(std::move(boundary),preparation);
  REQUIRE(boundary_owner.evaluate(standard,policy()).status==p::Status::ok);
  auto displaced=fixture();displaced.window.coordinates=p::Coordinates::fixed_h_reference;
  displaced.window.h_reference=1;
  displaced.spectrum.k_per_mpc={.03125,.125};
  auto displaced_owner=p::prepare(std::move(displaced),preparation);
  REQUIRE(empty_failure(displaced_owner.evaluate(standard,policy()),p::Status::outside_support));

  // Optional output gate must still run when the window annihilates everything.
  auto annihilate=fixture();
  annihilate.window.coordinates=p::Coordinates::fixed_h_reference;
  annihilate.window.h_reference=h;
  for (auto& k:annihilate.window.theory_k) k/=h;
  for (auto* w:{&annihilate.window.w00,&annihilate.window.w02,&annihilate.window.w20,&annihilate.window.w22})
    for (auto& coefficient:*w) coefficient=0;
  auto annihilate_owner=p::prepare(std::move(annihilate),preparation);
  const p::ModelPoint constant_response{.32,2,0,1,1};
  auto projected=annihilate_owner.evaluate(constant_response,policy());
  auto unprojected=owner.evaluate(constant_response,policy());
  REQUIRE(projected.status==p::Status::ok && unprojected.status==p::Status::ok);
  const long double factor=static_cast<long double>(h)*h*h;
  const long double physical_max=std::max(unprojected.input_multipoles[0].p0.combined_estimate,
                                         unprojected.input_multipoles[0].p2.combined_estimate);
  const long double serialized_max=std::max(projected.input_multipoles[0].p0.combined_estimate,
                                            projected.input_multipoles[0].p2.combined_estimate)/factor;
  REQUIRE(serialized_max>physical_max);
  auto output_budget=policy();
  output_budget.maximum_multipole_absolute_error_estimate_mpc3=
      static_cast<double>((physical_max+serialized_max)/2);
  REQUIRE(empty_failure(annihilate_owner.evaluate(constant_response,output_budget),p::Status::numerical_budget_exceeded));
  output_budget.retain_input_multipoles=false;
  REQUIRE(annihilate_owner.evaluate(constant_response,output_budget).status==p::Status::ok);

  auto moved=std::move(owner);
  REQUIRE(owner.status()==p::Status::invalid_owner && !owner.owns_source());
  REQUIRE(empty_failure(owner.evaluate(standard,policy()),p::Status::invalid_owner));
  moved=std::move(moved);
  REQUIRE(moved.status()==p::Status::ok);
  p::Prepared assigned;assigned=std::move(moved);
  REQUIRE(moved.status()==p::Status::invalid_owner && assigned.evaluate(standard,policy()).status==p::Status::ok);
  p::Result survives;
  { auto temporary_source=fixture();auto temporary=p::prepare(std::move(temporary_source),preparation);
    survives=temporary.evaluate(standard,policy()); }
  REQUIRE(survives.status==p::Status::ok && survives.identity.spectrum=="synthetic/dyadic-spectrum");
  REQUIRE(close(survives.mean[0].power.value,mono));
  REQUIRE(heap_peak<bytes);
  const auto& prep=assigned.preparation_work();
  const auto& work=r.work;
  std::cout << "windowed linear controls PASS; requested-new peak=" << heap_peak
            << ", owner-known-payload=" << r.known_live_payload_estimate_bytes
            << ", preparation-actions=" << prep.fields+prep.spectrum_nodes+
                prep.axis_nodes+prep.window_coefficients+prep.mapped_axis_writes
            << ", evaluation-actions=" << work.model_fields+work.support_columns+
                work.split_candidates+work.angular_callbacks+work.bracket_comparisons+
                work.interpolation_calls+work.window_products+work.output_writes
            << ", callbacks=" << work.angular_callbacks
            << ", refused-actions=" << limited.work.model_fields+limited.work.support_columns+
                limited.work.split_candidates+limited.work.angular_callbacks+
                limited.work.bracket_comparisons+limited.work.interpolation_calls+
                limited.work.window_products+limited.work.output_writes
            << ", refused-callbacks=" << cb.work.angular_callbacks
            << " (not RSS/opaque allocation); native fit is case-specific\n";
}
