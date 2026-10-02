#include "irred/baryon_abundance.hpp"
#include "payload_accounting.hpp"
#include "thermal_constants.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
constexpr std::size_t hard_rows = 65536, hard_bytes = 1ull << 30;
constexpr W eps = std::numeric_limits<W>::epsilon();
bool profile() {
  return std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384 &&
         std::fegetround() == FE_TONEAREST;
}
bool policy(BaryonAbundancePolicy p) {
  return p.maximum_rows <= hard_rows && p.maximum_native_bytes <= hard_bytes;
}
// Relative operation roundoff diagnostic, including explicit cast loss. A
// wide subnormal result cannot use this normal-arithmetic estimate.
BaryonAbundanceValue project(W v, W error) {
  BaryonAbundanceValue out;
  out.status = S::conditioning_budget_exceeded;
  if (v == 0) {
    out.status = S::ok; out.value = 0; return out;
  }
  if (!std::isnormal(v) || !(v > 0) || !std::isfinite(error)) return out;
  const double d = static_cast<double>(v);
  if (!std::isfinite(d) || !(d > 0)) return out;
  W e = error + std::abs(W(d)-v)/v + 4*eps;
  if (!(e <= 1e-15L)) return out;
  double diagnostic = static_cast<double>(e);
  if (W(diagnostic) < e) diagnostic = std::nextafter(diagnostic, INFINITY);
  if (!(diagnostic > 0) || diagnostic > 1e-15) return out;
  out.status = S::ok; out.value = d; out.relative_arithmetic_estimate = diagnostic;
  return out;
}
void fail_map(BaryonAbundanceRow &r, S s) {
  r.admission_status = s;
  r.hydrogen_nuclei.status = r.helium_nuclei.status = s;
}
void fail_lte(atomic::HydrogenHeliumRow &r, S s, unsigned mask) {
  r.admission_status = s;
  std::array<atomic::HydrogenHeliumValue *,6> values{
    &r.hydrogen_neutral,&r.hydrogen_ionized,&r.helium_neutral,
    &r.helium_singly_ionized,&r.helium_doubly_ionized,&r.electron_density};
  for (std::size_t i=0;i<values.size();++i) if (mask & (1u<<i)) {
    *values[i] = {}; values[i]->status = s;
  }
}
}
BaryonAbundance &BaryonAbundance::operator=(const BaryonAbundance &o) {
  if (this != &o) {
    if (o.source_) { BaryonAbundanceSource copy=*o.source_; source_=std::move(copy); }
    else source_.reset();
    status_=o.status_; rho0_=o.rho0_; rho_error_=o.rho_error_;
  }
  return *this;
}
BaryonAbundance::BaryonAbundance(BaryonAbundance &&o) noexcept
    : status_(o.status_), source_(std::move(o.source_)), rho0_(o.rho0_), rho_error_(o.rho_error_) {
  o.status_=S::invalid_input; o.source_.reset(); o.rho0_=o.rho_error_=0;
}
BaryonAbundance &BaryonAbundance::operator=(BaryonAbundance &&o) noexcept {
  if (this != &o) {
    status_=o.status_; source_=std::move(o.source_); rho0_=o.rho0_; rho_error_=o.rho_error_;
    o.status_=S::invalid_input; o.source_.reset(); o.rho0_=o.rho_error_=0;
  }
  return *this;
}
std::optional<std::size_t> BaryonAbundance::payload() const noexcept {
  irred::detail::PayloadAccounting b(sizeof(BaryonAbundance));
  if (source_) { b.string(source_->source_origin); b.string(source_->mass_origin); }
  return b.result();
}
BaryonAbundance prepare_baryon_abundance(const BaryonAbundanceSource &s, BaryonAbundancePolicy p) {
  BaryonAbundance out;
  if (!policy(p)) return out;
  irred::detail::PayloadAccounting bytes(sizeof(BaryonAbundance));
  bytes.string(s.source_origin); bytes.string(s.mass_origin);
  auto size=bytes.result();
  if (!size || *size>p.maximum_native_bytes) { out.status_=S::work_limit; return out; }
  // Structural admission precedes scans and immutable source copy.
  if (s.source_origin.empty() || s.mass_origin.empty()) return out;
  out.source_=s;
  for (double x : {s.physical_baryon_density,s.helium4_mass_fraction,
                   s.hydrogen1_effective_mass_kg,s.helium4_effective_mass_kg})
    if (!std::isfinite(x)) { out.status_=S::nonfinite_input; return out; }
  if (s.physical_baryon_density<0 || s.helium4_mass_fraction<0 || s.helium4_mass_fraction>1 ||
      !std::isnormal(s.hydrogen1_effective_mass_kg) || !(s.hydrogen1_effective_mass_kg>0) ||
      !std::isnormal(s.helium4_effective_mass_kg) || !(s.helium4_effective_mass_kg>0)) {
    out.status_=S::outside_domain; return out;
  }
  if (!profile()) return out;
  out.rho0_=detail::critical_mass_density_si(100)*W(s.physical_baryon_density);
  // 32 epsilon covers the original critical-density helper's conversions,
  // products/divisions, pi representation and multiplication by supplied omega.
  out.rho_error_=32*eps;
  if (s.physical_baryon_density>0 && !std::isnormal(out.rho0_)) {
    out.status_=S::conditioning_budget_exceeded; return out;
  }
  out.status_=S::ok; return out;
}
BaryonAbundanceRow BaryonAbundance::map(double a) const {
  BaryonAbundanceRow r; r.scale_factor=a;
  if (!std::isfinite(a)) { fail_map(r,S::nonfinite_input); return r; }
  if (!(a>0) || a>1) { fail_map(r,S::outside_domain); return r; }
  if (!profile()) { fail_map(r,S::invalid_input); return r; }
  const auto &s=*source_;
  // Wide intermediates preserve each species even if its partner is not
  // representable in binary64. No total-density binary64 bottleneck.
  W aa=a, volume=aa*aa*aa, rho=rho0_/volume;
  if (!(volume>0) || !std::isfinite(rho)) {
    fail_map(r,S::conditioning_budget_exceeded); return r;
  }
  W y=s.helium4_mass_fraction, h=1-y;
  W e=rho_error_+8*eps;
  r.hydrogen_nuclei=project(rho*h/W(s.hydrogen1_effective_mass_kg),e+2*eps);
  r.helium_nuclei=project(rho*y/W(s.helium4_effective_mass_kg),e+2*eps);
  r.admission_status=S::ok;
  return r;
}
BaryonAbundanceBatch BaryonAbundance::evaluate(std::span<const double> a, BaryonAbundancePolicy p) const {
  BaryonAbundanceBatch out;
  if (!policy(p)) return out;
  if (status_!=S::ok) { out.status=status_; return out; }
  auto n=payload();
  if (a.size()>p.maximum_rows || !n ||
      !irred::detail::checked_payload_add(*n,1,sizeof(out)) ||
      !irred::detail::checked_payload_add(*n,a.size(),sizeof(BaryonAbundanceRow)) ||
      *n>p.maximum_native_bytes) { out.status=S::work_limit; return out; }
  out.rows.reserve(a.size());
  for (double v:a) out.rows.push_back(map(v));
  out.status=S::ok; return out;
}
BaryonAbundanceLteBatch BaryonAbundance::evaluate_equilibrium(
    std::span<const BaryonAbundanceLteQuery> q, std::string_view origin,
    atomic::HydrogenHeliumPolicy h, BaryonAbundancePolicy p) const {
  BaryonAbundanceLteBatch out;
  if (!policy(p) || origin.empty() || !h.requested_outputs || (h.requested_outputs&~63u) ||
      h.maximum_rows>hard_rows || h.maximum_solves>hard_rows || h.maximum_root_iterations>256 ||
      h.maximum_charge_evaluations>4000000 || h.maximum_native_bytes>hard_bytes) return out;
  if (status_!=S::ok) { out.status=status_; return out; }
  auto n=payload();
  // Charge peak owner + final rows + temporary native inputs + temporary native
  // LTE result. The called owner's own bound is additionally checked below.
  if (q.size()>p.maximum_rows || q.size()>h.maximum_rows || !n ||
      !irred::detail::checked_payload_add(*n,1,sizeof(out)+sizeof(atomic::HydrogenHeliumBatch)) ||
      !irred::detail::checked_payload_add(*n,std::max(origin.size(),std::string{}.capacity()),1) ||
      !irred::detail::checked_payload_add(*n,1,1) ||
      !irred::detail::checked_payload_add(*n,q.size(),sizeof(BaryonAbundanceLteRow)+
           sizeof(atomic::HydrogenHeliumState)+sizeof(atomic::HydrogenHeliumRow)) ||
      *n>p.maximum_native_bytes) { out.status=S::work_limit; return out; }
  std::size_t lte_bytes=sizeof(atomic::HydrogenHeliumBatch);
  if (!irred::detail::checked_payload_add(lte_bytes,q.size(),sizeof(atomic::HydrogenHeliumRow)) ||
      lte_bytes>h.maximum_native_bytes) { out.status=S::work_limit; return out; }
  out.temperature_origin=origin;
  out.rows.resize(q.size());
  std::vector<atomic::HydrogenHeliumState> states(q.size());
  const bool interior=source_->physical_baryon_density>0 && source_->helium4_mass_fraction>0 &&
      source_->helium4_mass_fraction<1;
  for (std::size_t i=0;i<q.size();++i) {
    auto &r=out.rows[i]; r.source=q[i]; r.nuclei=map(q[i].scale_factor);
    states[i].temperature_kelvin=q[i].temperature_kelvin;
    if (interior && r.nuclei.admission_status==S::ok &&
        r.nuclei.hydrogen_nuclei.status==S::ok && r.nuclei.helium_nuclei.status==S::ok) {
      states[i].hydrogen_nuclei_per_cubic_metre=*r.nuclei.hydrogen_nuclei.value;
      states[i].helium_nuclei_per_cubic_metre=*r.nuclei.helium_nuclei.value;
    }
  }
  auto batch=atomic::evaluate_hydrogen_helium_equilibrium(states,h);
  out.status=batch.status;
  out.solves=batch.solves; out.root_iterations=batch.root_iterations;
  out.charge_evaluations=batch.charge_evaluations;
  if (batch.status!=S::ok || batch.rows.size()!=q.size()) { out.rows.clear(); return out; }
  for (std::size_t i=0;i<q.size();++i) {
    auto &r=out.rows[i]; r.equilibrium=std::move(batch.rows[i]);
    S cause = !interior ? S::outside_domain : r.nuclei.admission_status;
    if (cause==S::ok && r.nuclei.hydrogen_nuclei.status!=S::ok) cause=r.nuclei.hydrogen_nuclei.status;
    if (cause==S::ok && r.nuclei.helium_nuclei.status!=S::ok) cause=r.nuclei.helium_nuclei.status;
    if (cause!=S::ok) { fail_lte(r.equilibrium,cause,h.requested_outputs); continue; }
    W d=std::max(r.nuclei.hydrogen_nuclei.relative_arithmetic_estimate,
                 r.nuclei.helium_nuclei.relative_arithmetic_estimate);
    // Relative density perturbation -> conservative absolute log perturbation.
    W log_error=d/(1-d);
    std::array<atomic::HydrogenHeliumValue *,6> values{
      &r.equilibrium.hydrogen_neutral,&r.equilibrium.hydrogen_ionized,
      &r.equilibrium.helium_neutral,&r.equilibrium.helium_singly_ionized,
      &r.equilibrium.helium_doubly_ionized,&r.equilibrium.electron_density};
    for (std::size_t j=0;j<values.size();++j) {
      auto &v=*values[j]; if (v.status!=S::ok || !v.value) continue;
      W combined=W(v.relative_arithmetic_estimate)+(j==5?1:2)*log_error+8*eps;
      double e=static_cast<double>(combined);
      if (W(e)<combined) e=std::nextafter(e,INFINITY);
      if (!(e>0) || e>1e-13) { v.value.reset(); v.status=S::conditioning_budget_exceeded;
        v.relative_arithmetic_estimate=0; }
      else v.relative_arithmetic_estimate=e;
    }
  }
  return out;
}
} // namespace irred::cosmology
