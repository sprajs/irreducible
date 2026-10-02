#include "irred/thermal_neutrino.hpp"
#include "irred/quantities.hpp"
#include "payload_accounting.hpp"
#include "thermal_constants.hpp"
#include "thermal_cc_sample.hpp"
#include "thermal_cc_constants.hpp"
#include "thermal_fd_sample.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>
#include <utility>

namespace irred::cosmology {
namespace {
using S = numerics::Status;
constexpr long double pi = std::numbers::pi_v<long double>;
constexpr long double ev_joule = electron_volt_joule;
constexpr long double boltzmann_ev_kelvin = boltzmann_constant_joule_per_kelvin / ev_joule;
constexpr long double arithmetic_relative =
    64.L * std::numeric_limits<double>::epsilon();
bool arithmetic_supported() {
  return std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384 &&
         std::fegetround() == FE_TONEAREST;
}
bool valid_policy(const ThermalPolicy &p) {
  return std::isfinite(p.absolute_tolerance) &&
         std::isfinite(p.relative_tolerance) && p.absolute_tolerance >= 0 &&
         p.relative_tolerance >= 0 &&
         (p.absolute_tolerance > 0 || p.relative_tolerance > 0) &&
         p.maximum_depth <= 60 &&
         (p.momentum_method == ThermalMomentumMethod::direct_adaptive ||
          p.momentum_method == ThermalMomentumMethod::nested_clenshaw_curtis);
}
bool normal_positive(long double x) {
  return std::isfinite(x) && x >= std::numeric_limits<double>::min() &&
         x <= std::numeric_limits<double>::max();
}
// Every positive estimate remains positive and normal after storage. Refuse an
// unrepresentable diagnostic rather than silently reporting exact zero.
bool store(long double x, long double error, std::optional<double> &value,
           double &diagnostic) {
  const double rounded = static_cast<double>(x);
  error += std::abs(x - static_cast<long double>(rounded));
  if (!normal_positive(x) || !normal_positive(error)) return false;
  value = rounded;
  diagnostic = static_cast<double>(error);
  return true;
}
S species_status(const ThermalSpecies &s, double a) {
  for (double x : {s.mass_ev, s.temperature_today_ev, s.statistical_weight, a})
    if (!std::isfinite(x)) return S::nonfinite_input;
  if (s.mass_ev < 0 || !(s.temperature_today_ev > 0) ||
      !(s.statistical_weight > 0) || !(a > 0) || a > 1)
    return S::outside_domain;
  return S::ok;
}
struct MomentState { long double y, scale; bool pressure; };
double momentum(double q, const void *v) {
  const auto &s = *static_cast<const MomentState *>(v);
  return detail::thermal_fd_value(detail::thermal_fd_state(q,s.y,s.scale),s.pressure);
}
long double exponential_moment(unsigned n, long double q) {
  long double power = 1, polynomial = 1;
  for (unsigned k = 1; k <= n; ++k) {
    power *= q;
    polynomial = power + k * polynomial;
  }
  return std::exp(-q) * polynomial;
}
struct WideMoments {
  S status = S::invalid_input;
  long double rho = 0, pressure = 0, rho_error = 0, pressure_error = 0;
  std::size_t callbacks = 0;
};
struct CCPanel { long double rho=0, pressure=0, rho_error=0, pressure_error=0; S status=S::ok; };
CCPanel cc_panel(double lower, double upper, long double y, long double scale,
                 bool need_pressure, long double absolute, long double relative,
                 unsigned depth, const ThermalPolicy &p, std::size_t limit,
                 std::size_t &callbacks) {
  CCPanel out;
  std::array<detail::ThermalCCSample,65> samples{};
  std::array<long double,3> rho{}, pressure{};
  const double width=upper-lower;
  const unsigned charge=need_pressure?2:1;
  for (unsigned level=0; level<3; ++level) {
    const unsigned n=16u<<level, stride=64/n;
    for (unsigned j=0; j<=n; ++j) {
      const unsigned index=j*stride;
      if (level==0 || j%2) {
        if (limit-callbacks < charge) { out.status=S::work_limit; return out; }
        callbacks+=charge;
        const double q=j==0?lower:(j==n?upper:lower+width*detail::cc_nodes[index]);
        samples[index]=detail::thermal_cc_sample(q,y,scale,need_pressure);
        if (!std::isfinite(samples[index].rho) ||
            (need_pressure && !std::isfinite(samples[index].pressure))) {
          out.status=S::conditioning_budget_exceeded; return out;
        }
      }
      const double weight=level==0?detail::cc_weights_16[j]:
                          level==1?detail::cc_weights_32[j]:detail::cc_weights_64[j];
      rho[level]+=static_cast<long double>(weight)*samples[index].rho*width;
      if (need_pressure) pressure[level]+=static_cast<long double>(weight)*samples[index].pressure*width;
    }
  }
  out.rho=rho[2]; out.pressure=pressure[2];
  out.rho_error=8*std::max(std::abs(rho[2]-rho[1]),std::abs(rho[1]-rho[0]));
  if (need_pressure) out.pressure_error=8*std::max(std::abs(pressure[2]-pressure[1]),
                                                std::abs(pressure[1]-pressure[0]));
  if (out.rho_error<=absolute+relative*out.rho &&
      (!need_pressure || out.pressure_error<=absolute+relative*out.pressure)) return out;
  if (depth>=p.maximum_depth) { out.status=S::conditioning_budget_exceeded; return out; }
  const double middle=lower+width/2;
  if (!(middle>lower && middle<upper)) { out.status=S::conditioning_budget_exceeded; return out; }
  auto left=cc_panel(lower,middle,y,scale,need_pressure,absolute/2,relative,depth+1,p,limit,callbacks);
  if (left.status!=S::ok) return left;
  auto right=cc_panel(middle,upper,y,scale,need_pressure,absolute/2,relative,depth+1,p,limit,callbacks);
  if (right.status!=S::ok) return right;
  out.rho=left.rho+right.rho; out.pressure=left.pressure+right.pressure;
  out.rho_error=left.rho_error+right.rho_error;
  out.pressure_error=left.pressure_error+right.pressure_error;
  return out;
}
WideMoments moments(long double y, const ThermalPolicy &p,
                    std::size_t remaining, bool need_pressure = true) {
  WideMoments out;
  if (!std::isfinite(y) || y < 0) { out.status = S::outside_domain; return out; }
  const long double scale = std::hypot(1.L, y);
  // On q in [1,2], the normalized rho integral is at least
  // 1/(exp(2)+1), and normalized pressure at least 1/[6*(exp(2)+1)].
  // Thus .01 is a conservative positive lower control for either integral.
  // Their exponential upper controls are below 16; these constants select
  // tail/arithmetic work, while final admission uses the evaluated moment.
  const long double requested_lower = p.absolute_tolerance +
                                     p.relative_tolerance * .01L;
  if (p.absolute_tolerance + 16.L*p.relative_tolerance <
      arithmetic_relative * .01L) {
    out.status = S::conditioning_budget_exceeded; return out;
  }
  if (y == 0) {
    out.rho = 7*pi*pi*pi*pi/120;
    out.pressure = out.rho/3;
    out.rho_error = arithmetic_relative*out.rho;
    out.pressure_error = arithmetic_relative*out.pressure;
    out.status = S::ok;
  } else {
    // f(q)<=exp(-q), sqrt(q^2+y^2)<=q+y, and 1/sqrt(q^2+y^2)
    // <=min(1/q,1/y). These bound only the omitted mathematical tail.
    unsigned cutoff = 0;
    long double tails[2]{};
    for (unsigned q = 16; q <= 256; q += 8) {
      tails[0] = (exponential_moment(3,q) + y*exponential_moment(2,q))/scale;
      tails[1] = y >= 1 ? scale/y*exponential_moment(4,q)/3
                          : scale*exponential_moment(3,q)/3;
      if (std::max(tails[0],tails[1]) <= requested_lower/8) {
        cutoff = q; break;
      }
    }
    if (!cutoff) { out.status = S::conditioning_budget_exceeded; return out; }
    const unsigned panels = cutoff/4;
    long double sums[2]{}, errors[2]{};
    const std::size_t limit = std::min(p.maximum_callbacks_per_evaluation,
                                      remaining);
    bool cc_accepted=false;
    if (p.momentum_method==ThermalMomentumMethod::nested_clenshaw_curtis) {
      S cc_status=S::ok;
      for (unsigned panel=0; panel<panels; ++panel) {
        auto part=cc_panel(4*panel,4*(panel+1),y,scale,need_pressure,
                           p.absolute_tolerance/(4.L*panels),p.relative_tolerance/4.L,
                           0,p,limit,out.callbacks);
        if (part.status!=S::ok) { cc_status=part.status; break; }
        sums[0]+=part.rho; sums[1]+=part.pressure;
        errors[0]+=part.rho_error; errors[1]+=part.pressure_error;
      }
      if (cc_status==S::work_limit) { out.status=cc_status; return out; }
      cc_accepted=cc_status==S::ok;
      if (!cc_accepted) { sums[0]=sums[1]=errors[0]=errors[1]=0; }
    }
    for (unsigned kind = 0; kind < (need_pressure ? 2u : 1u); ++kind) {
      MomentState state{y, scale, kind == 1};
      for (unsigned panel = 0; !cc_accepted && panel < panels; ++panel) {
        const auto available = limit - out.callbacks;
        if (available < 5) { out.status = S::work_limit; return out; }
        const double absolute = p.absolute_tolerance/(4*panels);
        const double relative = p.relative_tolerance/4;
        if (absolute == 0 && relative == 0) {
          out.status = S::conditioning_budget_exceeded; return out;
        }
        const auto integral = numerics::integrate(
            momentum, &state, 4*panel, 4*(panel+1),
            {absolute, relative, available, p.maximum_depth});
        out.callbacks += integral.evaluations;
        if (integral.status != S::ok) { out.status = integral.status; return out; }
        sums[kind] += integral.value;
        errors[kind] += integral.error_estimate;
      }
      errors[kind] += tails[kind] + arithmetic_relative*sums[kind];
      const long double budget = p.absolute_tolerance +
                                p.relative_tolerance*sums[kind];
      if (errors[kind] > budget) {
        out.status = S::conditioning_budget_exceeded; return out;
      }
    }
    out.rho = sums[0]*scale; out.rho_error = errors[0]*scale;
    out.pressure = sums[1]/scale; out.pressure_error = errors[1]/scale;
    out.status = S::ok;
  }
  const long double rho_budget = p.absolute_tolerance +
                                p.relative_tolerance*out.rho/scale;
  const long double pressure_budget = p.absolute_tolerance +
                                    p.relative_tolerance*out.pressure*scale;
  if (out.rho_error/scale > rho_budget ||
      (need_pressure && out.pressure_error*scale > pressure_budget))
    out.status = S::conditioning_budget_exceeded;
  return out;
}
WideMoments species(const ThermalSpecies &s, double a, const ThermalPolicy &p,
                    std::size_t remaining, bool need_pressure = true) {
  WideMoments out;
  out.status = species_status(s,a);
  if (out.status != S::ok) return out;
  const long double temperature = static_cast<long double>(s.temperature_today_ev)/a;
  const long double y = static_cast<long double>(s.mass_ev)/temperature;
  out = moments(y,p,remaining,need_pressure);
  if (out.status != S::ok) return out;
  const long double t2 = temperature*temperature;
  const long double factor = s.statistical_weight*t2*t2/(2*pi*pi);
  out.rho *= factor; out.pressure *= factor;
  out.rho_error *= factor; out.pressure_error *= factor;
  return out;
}
long double critical_density_ev4(double h0) {
  const long double c = speed_of_light_m_per_s;
  const long double hbar = planck_constant_joule_second/(2*pi);
  const long double natural_to_si =
      ev_joule*ev_joule*ev_joule*ev_joule/(hbar*c*hbar*c*hbar*c);
  return detail::critical_energy_density_si(h0)/natural_to_si;
}
} // namespace

ThermalMoments evaluate_thermal_moments(double y, ThermalPolicy p) {
  ThermalMoments out;
  if (!arithmetic_supported() || !valid_policy(p)) return out;
  if (!std::isfinite(y)) { out.status = S::nonfinite_input; return out; }
  auto m = moments(y,p,p.maximum_total_callbacks);
  out.status = m.status; out.callbacks = m.callbacks;
  if (m.status != S::ok) return out;
  std::optional<double> rho, pressure; double re=0, pe=0;
  if (!store(m.rho,m.rho_error,rho,re) ||
      !store(m.pressure,m.pressure_error,pressure,pe)) {
    out.status = S::outside_domain; return out;
  }
  const long double scale=std::hypot(1.L,static_cast<long double>(y));
  if (re/scale>p.absolute_tolerance+p.relative_tolerance*(*rho/scale) ||
      pe*scale>p.absolute_tolerance+p.relative_tolerance*(*pressure*scale)) {
    out.status=S::conditioning_budget_exceeded; return out;
  }
  out.rho_moment=rho; out.pressure_moment=pressure;
  out.rho_error_estimate=re; out.pressure_error_estimate=pe;
  return out;
}
ThermalSpeciesState evaluate_thermal_species(const ThermalSpecies &s, double a,
                                            ThermalPolicy p) {
  ThermalSpeciesState out;
  if (!arithmetic_supported() || !valid_policy(p)) return out;
  auto m = species(s,a,p,p.maximum_total_callbacks);
  out.status=m.status; out.callbacks=m.callbacks;
  if (m.status != S::ok) return out;
  std::optional<double> rho, pressure; double re=0, pe=0;
  if (!store(m.rho,m.rho_error,rho,re) ||
      !store(m.pressure,m.pressure_error,pressure,pe)) {
    out.status=S::outside_domain; return out;
  }
  out.rho_ev4=rho; out.pressure_ev4=pressure;
  out.rho_error_estimate_ev4=re; out.pressure_error_estimate_ev4=pe;
  return out;
}
std::optional<std::size_t> thermal_background_payload_bound(
    std::size_t points, std::size_t count) noexcept {
  irred::detail::PayloadAccounting bytes(sizeof(ThermalBackground)+sizeof(ThermalBackgroundBatch));
  bytes.add(count,2*sizeof(ThermalSpecies));
  bytes.add(points,2*sizeof(ThermalBackgroundRow));
  return bytes.result();
}
ThermalPhysicalMapping map_thermal_physical_model(const ThermalPhysicalModel &m,
                                                 ThermalPolicy p) {
  ThermalPhysicalMapping out;
  if (!arithmetic_supported() || !valid_policy(p) || m.species.size()>16) return out;
  const auto bytes=thermal_background_payload_bound(0,m.species.size());
  if (!bytes) return out;
  if (m.species.size()>p.maximum_species || *bytes>p.maximum_native_bytes) {
    out.status=S::work_limit; return out;
  }
  for (double x:{m.h0_km_s_mpc,m.physical_baryon_density,m.physical_cdm_density,
                 m.tcmb_kelvin,m.physical_massless_nonphoton_density})
    if (!std::isfinite(x)) { out.status=S::nonfinite_input; return out; }
  if (!(m.h0_km_s_mpc>0) || !(m.tcmb_kelvin>0) || m.physical_baryon_density<0 ||
      m.physical_cdm_density<0 || m.physical_massless_nonphoton_density<0) {
    out.status=S::outside_domain; return out;
  }
  ThermalFlatModel mapped{}; mapped.h0_km_s_mpc=m.h0_km_s_mpc;
  const long double h=static_cast<long double>(m.h0_km_s_mpc)/100, h2=h*h;
  const long double t=boltzmann_ev_kelvin*m.tcmb_kelvin, t2=t*t;
  const long double photon=pi*pi/15*t2*t2/critical_density_ev4(m.h0_km_s_mpc);
  auto cast=[](long double x,double &dest) {
    if (x!=0 && !normal_positive(x)) return false;
    dest=static_cast<double>(x); return std::isfinite(dest);
  };
  if (!cast(photon,mapped.omega_gamma) ||
      !cast(m.physical_baryon_density/h2,mapped.omega_b) ||
      !cast(m.physical_cdm_density/h2,mapped.omega_cdm) ||
      !cast(m.physical_massless_nonphoton_density/h2,mapped.omega_massless_nonphoton)) {
    out.status=S::outside_domain; return out;
  }
  mapped.species.reserve(m.species.size());
  for (const auto &v:m.species) {
    if (!std::isfinite(v.temperature_today_kelvin)) { out.status=S::nonfinite_input; return out; }
    ThermalSpecies species{v.mass_ev,0,v.statistical_weight};
    if (!(v.temperature_today_kelvin>0) ||
        !cast(boltzmann_ev_kelvin*v.temperature_today_kelvin,species.temperature_today_ev)) {
      out.status=S::outside_domain; return out;
    }
    out.status=species_status(species,1);
    if (out.status!=S::ok) return out;
    mapped.species.push_back(species);
  }
  out.status=S::ok; out.model=std::move(mapped); return out;
}
ThermalBackground prepare_thermal_background(const ThermalFlatModel &m,
                                             ThermalPolicy p) {
  ThermalBackground out;
  if (!arithmetic_supported() || !valid_policy(p)) return out;
  const auto bytes=thermal_background_payload_bound(0,m.species.size());
  if (!bytes || m.species.size()>16) return out;
  if (m.species.size()>p.maximum_species || *bytes>p.maximum_native_bytes) {
    out.status_=S::work_limit; return out;
  }
  for (double x : {m.h0_km_s_mpc,m.omega_gamma,m.omega_massless_nonphoton,m.omega_b,m.omega_cdm})
    if (!std::isfinite(x)) { out.status_=S::nonfinite_input; return out; }
  if (!(m.h0_km_s_mpc>0) || m.omega_gamma<0 || m.omega_massless_nonphoton<0 ||
      m.omega_b<0 || m.omega_cdm<0) { out.status_=S::outside_domain; return out; }
  std::array<long double,20> fractions{};
  fractions[0]=m.omega_gamma; fractions[1]=m.omega_massless_nonphoton;
  fractions[2]=m.omega_b; fractions[3]=m.omega_cdm;
  for (const auto &s:m.species) {
    const auto status=species_status(s,1);
    if (status!=S::ok) { out.status_=status; return out; }
  }
  // Refuse a known component excess before acquiring any momentum integrals.
  auto supplied=fractions;
  std::sort(supplied.begin(),supplied.end(),std::greater<>());
  long double available=1;
  for (const auto x:supplied) {
    if (x>available) { out.status_=S::outside_domain; return out; }
    available-=x;
  }
  out.critical_ev4_=critical_density_ev4(m.h0_km_s_mpc);
  if (!(out.critical_ev4_>0) || !std::isfinite(out.critical_ev4_)) {
    out.status_=S::outside_domain; return out;
  }
  // Acquire storage before any FD callbacks: a source-copy allocation failure
  // must not discard already performed normalization work in an outer owner.
  // Keep failed/default source semantics; publish this owner only on success.
  ThermalFlatModel owned_source=m;
  static_assert(std::is_nothrow_move_assignable_v<ThermalFlatModel>);
  for (std::size_t i=0; i<m.species.size(); ++i) {
    auto v=species(m.species[i],1,p,p.maximum_total_callbacks-out.callbacks_,false);
    out.callbacks_+=v.callbacks;
    if (v.status!=S::ok) { out.status_=v.status; return out; }
    fractions[4+i]=v.rho/out.critical_ev4_;
    out.normalization_error_+=v.rho_error/out.critical_ev4_;
    out.omega_species_+=fractions[4+i];
  }
  // Subtract largest first: in particular 1+positive-tiny never rounds to an
  // admitted flat closure. No clipping or negative-Lambda fallback is used.
  std::sort(fractions.begin(),fractions.end(),std::greater<>());
  long double remainder=1;
  for (const auto x:fractions) {
    if (x>remainder) { out.status_=S::outside_domain; return out; }
    remainder-=x;
  }
  if (remainder<out.normalization_error_) {
    out.status_=S::conditioning_budget_exceeded; return out;
  }
  out.lambda_=remainder;
  out.source_=std::move(owned_source);
  out.method_=p.momentum_method;
  out.status_=S::ok;
  return out;
}
ThermalBackground::ThermalBackground(ThermalBackground &&other) noexcept {
  *this=std::move(other);
}
ThermalBackground &ThermalBackground::operator=(ThermalBackground &&other) noexcept {
  if (this==&other) return *this;
  status_=other.status_;
  source_=std::move(other.source_);
  method_=other.method_;
  critical_ev4_=other.critical_ev4_;
  omega_species_=other.omega_species_;
  lambda_=other.lambda_;
  normalization_error_=other.normalization_error_;
  callbacks_=other.callbacks_;
  other.status_=S::invalid_input;
  other.source_={};
  other.method_=ThermalMomentumMethod::direct_adaptive;
  other.critical_ev4_=other.omega_species_=other.lambda_=0;
  other.normalization_error_=0;
  other.callbacks_=0;
  return *this;
}
std::optional<double> ThermalBackground::omega_species_today() const noexcept {
  if (status_!=S::ok || (omega_species_!=0 && !normal_positive(omega_species_))) return {};
  return static_cast<double>(omega_species_);
}
std::optional<double> ThermalBackground::omega_lambda() const noexcept {
  if (status_!=S::ok || (lambda_!=0 && !normal_positive(lambda_))) return {};
  return static_cast<double>(lambda_);
}
ThermalScaledExpansion ThermalBackground::scaled_expansion(
    long double a, ThermalPolicy p) const {
  ThermalScaledExpansion out;
  if (!arithmetic_supported() || !valid_policy(p)) return out;
  if (status_!=S::ok) { out.status=status_; return out; }
  if (p.momentum_method!=method_) return out;
  if (!std::isfinite(a)) { out.status=S::nonfinite_input; return out; }
  if (a<0 || a>1) { out.status=S::outside_domain; return out; }
  const auto bytes=thermal_background_payload_bound(0,source_.species.size());
  if (!bytes) return out;
  if (source_.species.size()>p.maximum_species || *bytes>p.maximum_native_bytes) {
    out.status=S::work_limit; return out;
  }
  if (a==1) { out.status=S::ok; out.a4_e2=1; return out; }
  const long double a2=a*a, a4=a2*a2;
  out.a4_e2=source_.omega_gamma+static_cast<long double>(source_.omega_massless_nonphoton)+
      (source_.omega_b+static_cast<long double>(source_.omega_cdm))*a+lambda_*a4;
  out.error_estimate=normalization_error_*a4;
  for (const auto &s:source_.species) {
    const long double t=s.temperature_today_ev, t2=t*t;
    auto v=moments(static_cast<long double>(s.mass_ev)*a/t,p,
                   p.maximum_total_callbacks-out.callbacks,false);
    out.callbacks+=v.callbacks;
    if (v.status!=S::ok) { out.status=v.status; return out; }
    const long double factor=s.statistical_weight*t2*t2/(2*pi*pi*critical_ev4_);
    out.a4_e2+=factor*v.rho; out.error_estimate+=factor*v.rho_error;
  }
  out.error_estimate+=arithmetic_relative*out.a4_e2;
  if (a==0 && out.a4_e2==0 && out.error_estimate==0) out.status=S::ok;
  else if (out.a4_e2<std::numeric_limits<long double>::min() ||
           out.error_estimate<std::numeric_limits<long double>::min())
    out.status=S::outside_domain;
  else if (!(out.a4_e2>out.error_estimate) || !std::isfinite(out.a4_e2))
    out.status=S::conditioning_budget_exceeded;
  else out.status=S::ok;
  return out;
}
ThermalBackgroundBatch ThermalBackground::evaluate(
    std::span<const double> factors, unsigned requested, ThermalPolicy p) const {
  ThermalBackgroundBatch out;
  if (!arithmetic_supported() || !valid_policy(p) || (requested&~(thermal_e|thermal_h)) ||
      !requested || factors.size()>65536) return out;
  if (status_!=S::ok) { out.status=status_; return out; }
  if (p.momentum_method!=method_) return out;
  const auto bytes=thermal_background_payload_bound(factors.size(),source_.species.size());
  if (!bytes) return out;
  if (factors.size()>p.maximum_points || source_.species.size()>p.maximum_species ||
      *bytes>p.maximum_native_bytes) { out.status=S::work_limit; return out; }
  out.rows.reserve(factors.size()); out.requested_outputs=requested; out.status=S::ok;
  for (double a:factors) {
    out.rows.emplace_back(); auto &row=out.rows.back(); row.scale_factor=a;
    S status=S::ok;
    if (!std::isfinite(a)) status=S::nonfinite_input;
    else if (!(a>0) || a>1) status=S::outside_domain;
    long double e=0, e_error=0;
    if (status==S::ok && a==1) { e=1; }
    else if (status==S::ok) {
      auto local=p;
      local.maximum_total_callbacks=p.maximum_total_callbacks-out.callbacks;
      const auto scaled=scaled_expansion(a,local);
      row.callbacks+=scaled.callbacks; out.callbacks+=scaled.callbacks;
      status=scaled.status;
      if (status==S::ok) {
        const long double a2=static_cast<long double>(a)*a;
        e=std::sqrt(scaled.a4_e2)/a2;
        e_error=scaled.error_estimate/(2*std::sqrt(scaled.a4_e2-scaled.error_estimate)*a2)+
                arithmetic_relative*e;
      }
    }
    auto output=[&](unsigned mask,ThermalBackgroundValue &v,long double x,long double err) {
      if (!(requested&mask)) return;
      v.status=status;
      if (status!=S::ok) return;
      if (a==1 && normal_positive(x)) { v.value=static_cast<double>(x); v.status=S::ok; return; }
      if (!store(x,err,v.value,v.error_estimate)) v.status=S::outside_domain;
    };
    output(thermal_e,row.e,e,e_error);
    output(thermal_h,row.h_km_s_mpc,e*source_.h0_km_s_mpc,e_error*source_.h0_km_s_mpc);
  }
  return out;
}
} // namespace irred::cosmology
