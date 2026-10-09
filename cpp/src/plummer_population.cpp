#include "irred/plummer_population.hpp"
#include "payload_accounting.hpp"
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <numbers>
#include <utility>
namespace irred::gravity {
namespace {
using W=long double;
using S=numerics::Status;
constexpr W wide_epsilon=std::numeric_limits<W>::epsilon();
bool positive(double x) {return std::isfinite(x)&&x>0;}
bool arithmetic() {
  return std::fegetround()==FE_TONEAREST && std::numeric_limits<W>::digits>=64 &&
         std::numeric_limits<W>::max_exponent>=16384;
}
numerics::ScalarResult project(W value,W inherited,double tolerance) {
  double v=static_cast<double>(value);
  if(!std::isfinite(value)||!std::isfinite(v)||(value!=0&&!std::isnormal(v)))
    return {S::outside_domain};
  W error=inherited+128*wide_epsilon*std::abs(value)+std::abs(value-v)+
          8*std::numeric_limits<double>::epsilon()*std::abs(value);
  if(!std::isfinite(error)||error<0||!std::isfinite(static_cast<double>(error))||
      error>tolerance*std::abs(value)) return {S::conditioning_budget_exceeded};
  return {S::ok,v,static_cast<double>(error),0};
}
numerics::ScalarResult exponential(W logarithm,W relative_error,double tolerance) {
  W v=std::exp(logarithm);
  if(v==0) return {S::outside_domain}; // Positive support cannot underflow to zero.
  W rel=relative_error+128*wide_epsilon*(1+std::abs(logarithm));
  return project(v,std::abs(v)*rel,tolerance);
}
void fail(PlummerVelocityRow& row,S status) {
  row.status=status;
  row.distribution_function_kg_s3_m6.status=status;
  row.vector_velocity_density_s3_m3.status=status;
  row.speed_density_s_m.status=status;
}
}
std::optional<std::size_t> plummer_population_payload_bound(std::size_t n) noexcept {
  irred::detail::PayloadAccounting p(sizeof(PlummerPopulation));
  p.add(n,sizeof(PlummerPopulationRadius)); return p.result();
}
std::optional<std::size_t> plummer_population_batch_payload_bound(std::size_t n) noexcept {
  irred::detail::PayloadAccounting p(sizeof(PlummerPopulationBatch));
  p.add(n,sizeof(PlummerVelocityRow)); return p.result();
}
PlummerPopulation::PlummerPopulation(PlummerPopulation&& other) noexcept {
  *this=std::move(other);
}
PlummerPopulation& PlummerPopulation::operator=(PlummerPopulation&& other) noexcept {
  if(this!=&other) {
    status_=other.status_; source_=std::move(other.source_); radii_=std::move(other.radii_);
    log_coefficient_=other.log_coefficient_; log_coefficient_error_=other.log_coefficient_error_;
    requested_radii_=other.requested_radii_; native_evaluations_=other.native_evaluations_;
    other.status_=S::invalid_input; other.source_.reset(); other.radii_.clear();
    other.log_coefficient_=other.log_coefficient_error_=0;
    other.requested_radii_=other.native_evaluations_=0;
  }
  return *this;
}
PlummerPopulation prepare_plummer_population(PlummerSphere source,std::span<const double> radii,
                                             PlummerPopulationPreparationPolicy policy) {
  PlummerPopulation out;
  out.requested_radii_=radii.size();
  if(!arithmetic()||!positive(source.total_mass_kg)||!positive(source.scale_radius_metres)||
      !positive(source.gravitational_coupling_m3_kg_s2)||!positive(policy.relative_tolerance)||
      !policy.maximum_radii||!policy.maximum_native_evaluations||!policy.maximum_payload_bytes||
      radii.empty()) return out;
  auto bytes=plummer_population_payload_bound(radii.size());
  if(!bytes||*bytes>policy.maximum_payload_bytes||radii.size()>policy.maximum_radii||
      radii.size()>policy.maximum_native_evaluations/2) {out.status_=S::work_limit;return out;}
  out.source_=source;
  W lb=std::log(W(source.scale_radius_metres)),lg=std::log(W(source.gravitational_coupling_m3_kg_s2)),
    lm=std::log(W(source.total_mass_kg));
  W constant=std::log(24*std::sqrt(2.L)/(7*std::pow(std::numbers::pi_v<W>,3)));
  out.log_coefficient_=constant+2*lb-5*lg-4*lm;
  out.log_coefficient_error_=128*wide_epsilon*(1+std::abs(constant)+2*std::abs(lb)+
                                             5*std::abs(lg)+4*std::abs(lm));
  try {
    out.radii_.reserve(radii.size());
    auto actual=plummer_population_payload_bound(out.radii_.capacity());
    if(!actual||*actual>policy.maximum_payload_bytes) {out.status_=S::work_limit;return out;}
    out.status_=S::ok;
    for(double r:radii) {
      PlummerPopulationRadius state; state.radius_metres=r;
      auto native=evaluate_plummer_sphere(source,r,policy.relative_tolerance/16);
      auto projection=project_plummer_sphere(source,r,policy.relative_tolerance/16);
      out.native_evaluations_+=2;
      state.density_kg_m3=native.density_kg_m3;
      state.relative_potential_m2_s2=native.potential_m2_s2;
      state.relative_potential_m2_s2.value=-state.relative_potential_m2_s2.value;
      state.enclosed_mass_kg=native.enclosed_mass_kg;
      state.projected_surface_density_kg_m2=projection.surface_density_kg_m2;
      state.status=S::ok;
      // An unrelated native tidal eigenvalue can fail near its zero. Only
      // prerequisites for this population gate the retained radius state.
      for(auto field:{state.density_kg_m3,state.relative_potential_m2_s2,
                      state.projected_surface_density_kg_m2})
        if(field.status!=S::ok) {state.status=field.status;break;}
      if(state.status==S::ok) {
        W psi=state.relative_potential_m2_s2.value,error=state.relative_potential_m2_s2.error_estimate;
        W escape=std::sqrt(2*psi),factor=3*std::numbers::pi_v<W>/64;
        state.escape_speed_m_s=project(escape,error/escape,policy.relative_tolerance);
        state.one_axis_variance_m2_s2=project(psi/6,error/6,policy.relative_tolerance);
        state.projected_los_variance_m2_s2=project(factor*psi,factor*error,policy.relative_tolerance);
        for(auto field:{state.escape_speed_m_s,state.one_axis_variance_m2_s2,
                        state.projected_los_variance_m2_s2})
          if(field.status!=S::ok) {state.status=field.status;break;}
      }
      if(state.status!=S::ok && out.status_==S::ok) out.status_=state.status;
      out.radii_.push_back(state);
    }
  } catch(const std::bad_alloc&) {out.status_=S::work_limit;}
  return out;
}
PlummerPopulationBatch PlummerPopulation::evaluate(std::span<const PlummerVelocityRequest> requests,
                                                   PlummerPopulationEvaluationPolicy policy) const {
  PlummerPopulationBatch out; out.requested_rows=requests.size();
  if(status_!=S::ok) {out.status=status_;return out;}
  if(!arithmetic()||!positive(policy.relative_tolerance)||!policy.maximum_rows||
      !policy.maximum_velocity_evaluations||!policy.maximum_payload_bytes) return out;
  auto bytes=plummer_population_batch_payload_bound(requests.size());
  if(!bytes||*bytes>policy.maximum_payload_bytes||requests.size()>policy.maximum_rows||
      requests.size()>policy.maximum_velocity_evaluations) {out.status=S::work_limit;return out;}
  out.source=source_;
  try {
    out.rows.reserve(requests.size());
    auto actual=plummer_population_batch_payload_bound(out.rows.capacity());
    if(!actual||*actual>policy.maximum_payload_bytes) {out.status=S::work_limit;return out;}
    out.status=S::ok;
    for(auto request:requests) {
      PlummerVelocityRow row;row.request=request;
      ++out.velocity_evaluations;
      if(request.radius_index>=radii_.size()||!std::isfinite(request.speed_m_s)||request.speed_m_s<0) {
        fail(row,S::invalid_input);
      } else {
        const auto& state=radii_[request.radius_index];row.radius_metres=state.radius_metres;
        W speed=request.speed_m_s,kinetic=speed*speed/2,psi=state.relative_potential_m2_s2.value;
        W energy=psi-kinetic;
        W error=state.relative_potential_m2_s2.error_estimate+
                8*wide_epsilon*(std::abs(psi)+std::abs(kinetic));
        row.binding_energy_m2_s2=energy;row.binding_energy_absolute_error_m2_s2=error;
        if(energy+error<=0) {
          row.support=PlummerEnergySupport::outside;row.status=S::ok;
          row.distribution_function_kg_s3_m6={S::ok,0,0,0};
          row.vector_velocity_density_s3_m3={S::ok,0,0,0};
          row.speed_density_s_m={S::ok,0,0,0};
        } else if(energy-error<=0) {
          row.support=PlummerEnergySupport::ambiguous;fail(row,S::conditioning_budget_exceeded);
        } else {
          row.support=PlummerEnergySupport::bound;
          W log_energy=std::log(energy),log_f=log_coefficient_+3.5L*log_energy;
          W relative=log_coefficient_error_+3.5L*error/energy+
                     128*wide_epsilon*(1+std::abs(log_f)+std::abs(log_energy));
          W rho=state.density_kg_m3.value,log_rho=std::log(rho);
          W normalized_error=relative+state.density_kg_m3.error_estimate/rho+
                              128*wide_epsilon*(1+std::abs(log_rho));
          row.distribution_function_kg_s3_m6=exponential(log_f,relative,policy.relative_tolerance);
          row.vector_velocity_density_s3_m3=exponential(log_f-log_rho,normalized_error,policy.relative_tolerance);
          if(speed==0) row.speed_density_s_m={S::ok,0,0,0};
          else {
            W logarithm=std::log(4*std::numbers::pi_v<W>)+2*std::log(speed)+log_f-log_rho;
            row.speed_density_s_m=exponential(logarithm,normalized_error,policy.relative_tolerance);
          }
          row.status=S::ok;
          for(auto field:{row.distribution_function_kg_s3_m6,row.vector_velocity_density_s3_m3,
                          row.speed_density_s_m})
            if(field.status!=S::ok) {row.status=field.status;break;}
        }
      }
      if(row.status!=S::ok && out.status==S::ok) out.status=row.status;
      out.rows.push_back(row);
    }
  } catch(const std::bad_alloc&) {out.status=S::work_limit;}
  return out;
}
}
