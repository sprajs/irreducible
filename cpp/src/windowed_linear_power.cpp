#include "irred/windowed_linear_power.hpp"
#include "irred/numerics.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <span>
#include <utility>

namespace irred::windowed_linear_power {
namespace {
using W = long double;
constexpr std::size_t payload_ceiling = 16 * 1024 * 1024;
constexpr std::size_t source_ceiling = 2000000, action_ceiling = 8000000;
constexpr std::size_t callback_ceiling = 2000000;
constexpr W u = std::numeric_limits<W>::epsilon();
constexpr W ud = std::numeric_limits<double>::epsilon();
struct FloatingScope {
  std::fenv_t saved{};
  bool active = false;
  FloatingScope() {
    if (std::fegetround() == FE_TONEAREST &&
        std::numeric_limits<W>::radix == 2 &&
        std::numeric_limits<W>::digits >= 64 &&
        std::numeric_limits<double>::is_iec559 &&
        std::numeric_limits<double>::digits == 53)
      active = std::feholdexcept(&saved) == 0;
  }
  ~FloatingScope() { if (active) std::fesetenv(&saved); }
  bool range_ok() const {
    return !std::fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID |
                             FE_DIVBYZERO);
  }
};
bool normal_or_zero(W x) {
  return std::isfinite(x) && (x == 0 || std::isnormal(x));
}
bool positive(W x) { return std::isfinite(x) && x > 0 && std::isnormal(x); }
bool emit_estimate(W x, double& result) {
  if (!normal_or_zero(x) || x < 0 || x > std::numeric_limits<double>::max())
    return false;
  result = static_cast<double>(x);
  if (x == 0) return result == 0;
  if (!(result > 0) || !std::isnormal(result)) return false;
  if (static_cast<W>(result) < x)
    result = std::nextafter(result, std::numeric_limits<double>::infinity());
  return std::isfinite(result);
}
void source_payload(detail::PayloadAccounting& b, const Source& s) {
  b.string(s.identity.spectrum); b.string(s.identity.window);
  b.string(s.identity.matter_component); b.string(s.identity.growth_response);
  b.string(s.identity.grid_order); b.string(s.identity.window_calibration);
  b.vector(s.spectrum.k_per_mpc); b.vector(s.spectrum.power_mpc3);
  b.vector(s.window.theory_k); b.vector(s.window.output0_k);
  b.vector(s.window.output2_k); b.vector(s.window.output0_ids);
  b.vector(s.window.output2_ids);
  b.vector(s.window.w00); b.vector(s.window.w02);
  b.vector(s.window.w20); b.vector(s.window.w22);
}
void identity_payload(detail::PayloadAccounting& b, const SourceIdentity& s) {
  b.string(s.spectrum); b.string(s.window); b.string(s.matter_component);
  b.string(s.growth_response); b.string(s.grid_order);
  b.string(s.window_calibration);
}
struct PreparationLedger {
  std::size_t limit, total = 0;
  Status status = Status::ok;
  bool charge(std::size_t& category) {
    if (status != Status::ok) return false;
    if (total >= limit || category == std::numeric_limits<std::size_t>::max()) {
      status = Status::work_limit; return false;
    }
    ++total; ++category; return true;
  }
};
struct Ledger {
  const EvaluationPolicy& policy;
  EvaluationWork& work;
  std::size_t total = 0;
  Status status = Status::ok;
  void fail(Status s) { if (status == Status::ok) status = s; }
  bool charge(std::size_t& category) {
    if (status != Status::ok) return false;
    if (total >= policy.maximum_actions ||
        category == std::numeric_limits<std::size_t>::max() ||
        (&category == &work.angular_callbacks &&
         category >= policy.maximum_angular_callbacks)) {
      fail(Status::work_limit); return false;
    }
    ++total; ++category; return true;
  }
  std::size_t remaining() const { return policy.maximum_actions - total; }
};
bool valid_axis(const std::vector<double>& a, PreparationLedger& l,
                PreparationWork& w) {
  for (std::size_t j = 0; j < a.size(); ++j) {
    if (!l.charge(w.axis_nodes)) return false;
    if (!std::isfinite(a[j])) { l.status = Status::nonfinite_input; return false; }
    if (a[j] <= 0 || (j && a[j] <= a[j-1])) {
      l.status = Status::invalid_input; return false;
    }
  }
  return true;
}
bool valid_ids(const std::vector<std::uint64_t>& ids, PreparationLedger& l,
               PreparationWork& w) {
  for (std::size_t j = 0; j < ids.size(); ++j) {
    if (!l.charge(w.axis_nodes)) return false;
    for (std::size_t k = 0; k < j; ++k)
      if (ids[j] == ids[k]) { l.status = Status::invalid_input; return false; }
  }
  return true;
}
struct Geometry {
  W transverse_inverse = 0, parallel_inverse = 0;
  W low_q = 0, high_q = 0, volume = 0, condition = 0;
};
bool geometry(const ModelPoint& m, Geometry& g) {
  g.transverse_inverse = 1.L / m.alpha_perp;
  g.parallel_inverse = 1.L / m.alpha_parallel;
  g.low_q = std::min(g.transverse_inverse, g.parallel_inverse);
  g.high_q = std::max(g.transverse_inverse, g.parallel_inverse);
  g.volume = static_cast<W>(m.alpha_perp) * m.alpha_perp * m.alpha_parallel;
  g.condition = 1 + g.high_q / g.low_q;
  return positive(g.low_q) && positive(g.high_q) && positive(g.volume) &&
         positive(g.condition);
}
struct Angular {
  const Spectrum& spectrum;
  const ModelPoint& model;
  const Geometry& g;
  W k;
  unsigned ell;
  Ledger& ledger;
  double operator()(double mu) {
    if (!ledger.charge(ledger.work.angular_callbacks))
      return std::numeric_limits<double>::quiet_NaN();
    const W t = mu, t2 = t*t;
    const W q = std::sqrt((1-t2)*g.transverse_inverse*g.transverse_inverse +
                         t2*g.parallel_inverse*g.parallel_inverse);
    const W kt = k*q, mt = t*g.parallel_inverse/q;
    if (!positive(q) || !positive(kt) || !std::isfinite(mt) || mt < 0 || mt > 1) {
      ledger.fail(Status::unresolved_projection);
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double query = static_cast<double>(kt);
    if (!std::isfinite(query) || query <= 0 ||
        query < spectrum.k_per_mpc.front() || query > spectrum.k_per_mpc.back()) {
      ledger.fail(Status::outside_support);
      return std::numeric_limits<double>::quiet_NaN();
    }
    std::size_t lower=0, upper=spectrum.k_per_mpc.size();
    while (lower<upper) {
      if (!ledger.charge(ledger.work.bracket_comparisons))
        return std::numeric_limits<double>::quiet_NaN();
      const auto middle=lower+(upper-lower)/2;
      if (spectrum.k_per_mpc[middle]<query) lower=middle+1;
      else upper=middle;
    }
    const auto j=lower;
    if (j == spectrum.k_per_mpc.size()) {
      ledger.fail(Status::outside_support);
      return std::numeric_limits<double>::quiet_NaN();
    }
    const auto first = j == 0 ? 0 : j-1;
    if (!ledger.charge(ledger.work.interpolation_calls))
      return std::numeric_limits<double>::quiet_NaN();
    const auto p = numerics::interpolate_linear(
        std::span(spectrum.k_per_mpc).subspan(first,2),
        std::span(spectrum.power_mpc3).subspan(first,2),query);
    if (p.status != numerics::Status::ok) {
      ledger.fail(Status::unresolved_projection);
      return std::numeric_limits<double>::quiet_NaN();
    }
    const bool positive_p =
        (query == spectrum.k_per_mpc[first] ? spectrum.power_mpc3[first] > 0 :
         query == spectrum.k_per_mpc[first+1] ? spectrum.power_mpc3[first+1] > 0 :
         spectrum.power_mpc3[first] > 0 || spectrum.power_mpc3[first+1] > 0);
    if (positive_p && !(p.value > 0)) {
      ledger.fail(Status::unresolved_projection);
      return std::numeric_limits<double>::quiet_NaN();
    }
    const W response = model.b1 + static_cast<W>(model.f)*mt*mt;
    const W power = response*response*p.value/g.volume;
    const W weighted = ell == 0 ? power : 5*power*(3*t2-1)/2;
    const double out = static_cast<double>(weighted);
    if (!normal_or_zero(power) || !normal_or_zero(weighted) ||
        !std::isfinite(out) ||
        (weighted != 0 && (out == 0 || !std::isnormal(out)))) {
      ledger.fail(Status::unresolved_projection);
      return std::numeric_limits<double>::quiet_NaN();
    }
    return out;
  }
};
double angular_callback(double mu, const void* c) {
  return (*const_cast<Angular*>(static_cast<const Angular*>(c)))(mu);
}
// Uniform standard-operation diagnostic, NOT a proven libm/arithmetic bound.
// interpolate_linear reports no internal estimate: this owner supplies it.
W callback_diagnostic(const Spectrum& s, const ModelPoint& m,
                      const Geometry& g, W slope, W max_p) {
  const W amplitude = std::abs(static_cast<W>(m.b1)) + std::abs(static_cast<W>(m.f));
  const W k_error = (64*u*g.condition+ud)*s.k_per_mpc.back();
  const W p_error = slope*k_error + (64*u+ud)*max_p;
  const W response_error = 128*u*amplitude*g.condition;
  return ((2*amplitude*response_error+response_error*response_error)*max_p +
          amplitude*amplitude*p_error)/g.volume +
         (64*u*g.condition+ud)*amplitude*amplitude*max_p/g.volume;
}
bool fill_power(W value, W quadrature, W arithmetic, W grid,
                PowerValue& out) {
  if (!normal_or_zero(value) || std::abs(value)>std::numeric_limits<double>::max())
    return false;
  out.value = static_cast<double>(value);
  if (value != 0 && (out.value == 0 || !std::isnormal(out.value))) return false;
  arithmetic += std::abs(value-static_cast<W>(out.value));
  if (!emit_estimate(quadrature,out.quadrature_estimate) ||
      !emit_estimate(arithmetic,out.arithmetic_estimate) ||
      !emit_estimate(grid,out.grid_projection_estimate)) return false;
  return emit_estimate(static_cast<W>(out.quadrature_estimate)+
                       out.arithmetic_estimate+out.grid_projection_estimate,
                       out.combined_estimate);
}
bool project_power(const PowerValue& physical, W factor, W factor_error,
                   PowerValue& out) {
  const W value = static_cast<W>(physical.value)*factor;
  if (!normal_or_zero(value) || std::abs(value)>std::numeric_limits<double>::max())
    return false;
  out.value = static_cast<double>(value);
  if (value != 0 && (out.value == 0 || !std::isnormal(out.value))) return false;
  const W volume = (std::abs(static_cast<W>(physical.value))+
                    physical.combined_estimate)*factor_error +
                    (factor_error == 0 ? 0 : 4*u*std::abs(value)) +
                    std::abs(value-static_cast<W>(out.value));
  if (!emit_estimate(factor*physical.quadrature_estimate,out.quadrature_estimate) ||
      !emit_estimate(factor*physical.arithmetic_estimate,out.arithmetic_estimate) ||
      !emit_estimate(factor*physical.grid_projection_estimate,out.grid_projection_estimate) ||
      !emit_estimate(volume,out.volume_projection_estimate)) return false;
  return emit_estimate(static_cast<W>(out.quadrature_estimate)+out.arithmetic_estimate+
                       out.grid_projection_estimate+out.volume_projection_estimate,
                       out.combined_estimate);
}
} // namespace

Prepared::Prepared(Prepared&& o) noexcept
    : source_(std::move(o.source_)),
      mapped_theory_k_per_mpc_(std::move(o.mapped_theory_k_per_mpc_)),
      grid_projection_estimate_per_mpc_(std::move(o.grid_projection_estimate_per_mpc_)),
      maximum_abs_spectrum_slope_(std::exchange(o.maximum_abs_spectrum_slope_,0)),
      status_(std::exchange(o.status_,Status::invalid_owner)),
      work_(std::exchange(o.work_,{})), owns_source_(std::exchange(o.owns_source_,false)) {
  o.source_ = {};
}
Prepared& Prepared::operator=(Prepared&& o) noexcept {
  if (this != &o) {
    source_ = std::move(o.source_);
    mapped_theory_k_per_mpc_ = std::move(o.mapped_theory_k_per_mpc_);
    grid_projection_estimate_per_mpc_ = std::move(o.grid_projection_estimate_per_mpc_);
    maximum_abs_spectrum_slope_ = std::exchange(o.maximum_abs_spectrum_slope_,0);
    status_ = std::exchange(o.status_,Status::invalid_owner);
    work_ = std::exchange(o.work_,{});
    owns_source_ = std::exchange(o.owns_source_,false);
    o.source_ = {};
  }
  return *this;
}
Status Prepared::status() const noexcept { return status_; }
bool Prepared::owns_source() const noexcept { return owns_source_; }
const PreparationWork& Prepared::preparation_work() const noexcept { return work_; }
std::optional<std::size_t> Prepared::retained_payload_bytes() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  source_payload(b,source_);
  b.vector(mapped_theory_k_per_mpc_); b.vector(grid_projection_estimate_per_mpc_);
  return b.result();
}

Prepared prepare(Source&& s, const PreparationPolicy& policy) {
  Prepared out;
  out.status_ = Status::invalid_input;
  if (policy.maximum_source_actions == 0 || policy.maximum_source_actions>source_ceiling ||
      policy.maximum_live_payload_bytes == 0 || policy.maximum_live_payload_bytes>payload_ceiling)
    return out;
  FloatingScope floating;
  if (!floating.active) { out.status_ = Status::unsupported_arithmetic; return out; }
  PreparationLedger ledger{policy.maximum_source_actions};
  auto fail = [&](Status status) {
    out.status_ = status;
    std::vector<double>().swap(out.mapped_theory_k_per_mpc_);
    std::vector<W>().swap(out.grid_projection_estimate_per_mpc_);
  };
  const auto nk=s.spectrum.k_per_mpc.size(), nt=s.window.theory_k.size();
  const auto n0=s.window.output0_k.size(), n2=s.window.output2_k.size();
  if (nk<2 || nk>4096 || nt==0 || nt>2048 || n0==0 || n2==0 ||
      n0>256 || n2>256-n0 || s.spectrum.power_mpc3.size()!=nk ||
      s.window.output0_ids.size()!=n0 || s.window.output2_ids.size()!=n2 ||
      s.window.w00.size()!=n0*nt || s.window.w02.size()!=n0*nt ||
      s.window.w20.size()!=n2*nt || s.window.w22.size()!=n2*nt) return out;
  detail::PayloadAccounting requested(sizeof(Source)+2*sizeof(Prepared));
  source_payload(requested,s);
  requested.add(nt,sizeof(double)+sizeof(W));
  if (!requested.result() || *requested.result()>policy.maximum_live_payload_bytes) {
    fail(Status::payload_limit); return out;
  }
  const std::string* labels[]{&s.identity.spectrum,&s.identity.window,
      &s.identity.matter_component,&s.identity.growth_response,
      &s.identity.grid_order,&s.identity.window_calibration};
  std::size_t label_bytes=0;
  for (const auto* label:labels) {
    if (!ledger.charge(out.work_.fields)) { fail(ledger.status); return out; }
    if (label->empty() || label->size()>32768-label_bytes) return out;
    label_bytes+=label->size();
  }
  if (!ledger.charge(out.work_.fields)) { fail(ledger.status); return out; }
  if (!std::isfinite(s.spectrum.redshift)) { fail(Status::nonfinite_input); return out; }
  if (s.spectrum.redshift<0) return out;
  if (!ledger.charge(out.work_.fields)) { fail(ledger.status); return out; }
  W h=1;
  if (s.window.coordinates==Coordinates::fixed_h_reference) {
    if (!s.window.h_reference || !std::isfinite(*s.window.h_reference) ||
        *s.window.h_reference<=0) return out;
    h=*s.window.h_reference;
  } else if (s.window.coordinates!=Coordinates::physical_mpc || s.window.h_reference)
    return out;
  for (std::size_t j=0;j<nk;++j) {
    if (!ledger.charge(out.work_.spectrum_nodes)) { fail(ledger.status); return out; }
    const auto k=s.spectrum.k_per_mpc[j], p=s.spectrum.power_mpc3[j];
    if (!std::isfinite(k) || !std::isfinite(p)) { fail(Status::nonfinite_input); return out; }
    if (k<=0 || p<0 || (j && k<=s.spectrum.k_per_mpc[j-1])) return out;
    if (j) {
      const W slope=std::abs((static_cast<W>(p)-s.spectrum.power_mpc3[j-1])/
                            (static_cast<W>(k)-s.spectrum.k_per_mpc[j-1]));
      if (!normal_or_zero(slope)) { fail(Status::unresolved_projection); return out; }
      out.maximum_abs_spectrum_slope_=std::max(out.maximum_abs_spectrum_slope_,slope);
    }
  }
  if (!valid_axis(s.window.theory_k,ledger,out.work_) ||
      !valid_axis(s.window.output0_k,ledger,out.work_) ||
      !valid_axis(s.window.output2_k,ledger,out.work_) ||
      !valid_ids(s.window.output0_ids,ledger,out.work_) ||
      !valid_ids(s.window.output2_ids,ledger,out.work_)) { fail(ledger.status); return out; }
  const std::vector<double>* blocks[]{&s.window.w00,&s.window.w02,&s.window.w20,&s.window.w22};
  for (const auto* block:blocks)
    for (double coefficient:*block) {
      if (!ledger.charge(out.work_.window_coefficients)) { fail(ledger.status); return out; }
      if (!std::isfinite(coefficient)) { fail(Status::nonfinite_input); return out; }
    }
  try {
    out.mapped_theory_k_per_mpc_.resize(nt);
    detail::PayloadAccounting pending(sizeof(Source)+2*sizeof(Prepared));
    source_payload(pending,s);
    pending.vector(out.mapped_theory_k_per_mpc_);
    pending.add(nt,sizeof(W));
    if (!pending.result() || *pending.result()>policy.maximum_live_payload_bytes) {
      fail(Status::payload_limit); return out;
    }
    out.grid_projection_estimate_per_mpc_.resize(nt);
    detail::PayloadAccounting actual(sizeof(Source)+2*sizeof(Prepared));
    source_payload(actual,s);
    actual.vector(out.mapped_theory_k_per_mpc_);
    actual.vector(out.grid_projection_estimate_per_mpc_);
    if (!actual.result() || *actual.result()>policy.maximum_live_payload_bytes) {
      fail(Status::payload_limit); return out;
    }
    for (std::size_t j=0;j<nt;++j) {
      if (!ledger.charge(out.work_.mapped_axis_writes)) { fail(ledger.status); return out; }
      const W wide=static_cast<W>(s.window.theory_k[j])*h;
      const double mapped=static_cast<double>(wide);
      if (!positive(wide) || !std::isnormal(mapped)) {
        fail(Status::unresolved_projection); return out;
      }
      out.mapped_theory_k_per_mpc_[j]=mapped;
      out.grid_projection_estimate_per_mpc_[j]=
          s.window.coordinates==Coordinates::physical_mpc ? 0 :
          std::abs(wide-static_cast<W>(mapped))+2*u*std::abs(wide);
      if (j && mapped<=out.mapped_theory_k_per_mpc_[j-1]) {
        fail(Status::unresolved_projection); return out;
      }
    }
  } catch (const std::bad_alloc&) { fail(Status::payload_limit); return out; }
  if (!floating.range_ok()) { fail(Status::unresolved_projection); return out; }
  out.source_=std::move(s);
  out.owns_source_=true;
  out.status_=Status::ok;
  return out;
}

Result Prepared::evaluate(const ModelPoint& model, const EvaluationPolicy& policy) const {
  Result result;
  result.model=model;
  result.preparation_work=work_;
  if (status_!=Status::ok || !owns_source_) return result;
  auto fail=[&](Status status) {
    result.status=status;
    result.mean.clear(); result.input_multipoles.clear();
    return std::move(result);
  };
  if (!std::isfinite(policy.maximum_multipole_absolute_error_estimate_mpc3) ||
      policy.maximum_multipole_absolute_error_estimate_mpc3<=0 ||
      !std::isfinite(policy.maximum_mean_absolute_error_estimate_output_unit) ||
      policy.maximum_mean_absolute_error_estimate_output_unit<=0 ||
      policy.maximum_actions==0 || policy.maximum_actions>action_ceiling ||
      policy.maximum_angular_callbacks<3 || policy.maximum_angular_callbacks>callback_ceiling ||
      policy.maximum_callbacks_per_integral<3 || policy.maximum_callbacks_per_integral>4096 ||
      policy.maximum_depth>20 || policy.maximum_live_payload_bytes==0 ||
      policy.maximum_live_payload_bytes>payload_ceiling) return fail(Status::invalid_input);
  FloatingScope floating;
  if (!floating.active) return fail(Status::unsupported_arithmetic);
  const auto nt=source_.window.theory_k.size();
  const auto nr=source_.window.output0_k.size()+source_.window.output2_k.size();
  const auto ns=source_.spectrum.k_per_mpc.size();
  detail::PayloadAccounting requested(2*sizeof(Result)+
      2*sizeof(std::vector<PowerValue>)+2*sizeof(std::vector<double>));
  requested.embedded(retained_payload_bytes(),0);
  // Result strings are deliberate copies; allow their source capacities.
  identity_payload(requested,source_.identity);
  requested.add(nt,2*sizeof(PowerValue));
  requested.add(ns+2,sizeof(double));
  requested.add(2*nt,sizeof(double));
  requested.add(nr,sizeof(MeanRow));
  if (policy.retain_input_multipoles) requested.add(nt,sizeof(MultipoleRow));
  if (!requested.result() || *requested.result()>policy.maximum_live_payload_bytes)
    return fail(Status::payload_limit);
  std::vector<PowerValue> p0,p2;
  std::vector<double> splits,terms;
  auto actual_payload=[&]() {
    detail::PayloadAccounting b(2*sizeof(Result)+sizeof(p0)+sizeof(p2)+sizeof(splits)+sizeof(terms));
    b.embedded(retained_payload_bytes(),0);
    identity_payload(b,result.identity);
    b.vector(p0); b.vector(p2); b.vector(splits); b.vector(terms);
    b.vector(result.mean); b.vector(result.input_multipoles);
    return b.result();
  };
  try {
    result.identity=source_.identity;
    result.coordinates=source_.window.coordinates;
    result.h_reference=source_.window.h_reference;
    // Check actual old plus all still-pending requested storage at each stage.
    auto pending_fits=[&](std::size_t next_stage) {
      const auto current=actual_payload();
      if (current)
        result.known_live_payload_estimate_bytes=
            std::max(result.known_live_payload_estimate_bytes,*current);
      detail::PayloadAccounting b(2*sizeof(Result)+sizeof(p0)+sizeof(p2)+sizeof(splits)+sizeof(terms));
      b.embedded(retained_payload_bytes(),0); identity_payload(b,result.identity);
      b.vector(p0); b.vector(p2); b.vector(splits); b.vector(terms);
      b.vector(result.mean); b.vector(result.input_multipoles);
      if (next_stage==0) b.add(nt,sizeof(PowerValue));
      if (next_stage<=1) b.add(nt,sizeof(PowerValue));
      if (next_stage<=2) b.add(ns+2,sizeof(double));
      if (next_stage<=3) b.add(2*nt,sizeof(double));
      if (next_stage<=4) b.add(nr,sizeof(MeanRow));
      if (next_stage<=5 && policy.retain_input_multipoles) b.add(nt,sizeof(MultipoleRow));
      return b.result() && *b.result()<=policy.maximum_live_payload_bytes;
    };
    if (!pending_fits(0)) return fail(Status::payload_limit);
    p0.resize(nt);
    if (!pending_fits(1)) return fail(Status::payload_limit);
    p2.resize(nt);
    if (!pending_fits(2)) return fail(Status::payload_limit);
    splits.reserve(ns+2);
    if (!pending_fits(3)) return fail(Status::payload_limit);
    terms.resize(2*nt);
    if (!pending_fits(4)) return fail(Status::payload_limit);
    result.mean.reserve(nr);
    if (!pending_fits(5)) return fail(Status::payload_limit);
    if (policy.retain_input_multipoles) result.input_multipoles.reserve(nt);
    const auto actual=actual_payload();
    if (!actual || *actual>policy.maximum_live_payload_bytes) return fail(Status::payload_limit);
    result.known_live_payload_estimate_bytes=*actual;
  } catch (const std::bad_alloc&) { return fail(Status::payload_limit); }
  Ledger ledger{policy,result.work};
  const double fields[]{model.redshift,model.b1,model.f,model.alpha_perp,model.alpha_parallel};
  for (double field:fields) {
    if (!ledger.charge(result.work.model_fields)) return fail(ledger.status);
    if (!std::isfinite(field)) return fail(Status::nonfinite_input);
  }
  if (model.redshift!=source_.spectrum.redshift) return fail(Status::incompatible_epoch);
  if (model.alpha_perp<=0 || model.alpha_parallel<=0) return fail(Status::invalid_input);
  Geometry g;
  if (!geometry(model,g)) return fail(Status::unresolved_projection);
  W output_factor=1,output_factor_error=0;
  if (result.h_reference) {
    const W h=*result.h_reference;
    output_factor=h*h*h;
    output_factor_error=8*u*std::abs(output_factor);
    if (!positive(output_factor) || !positive(output_factor_error))
      return fail(Status::unresolved_projection);
  }
  W max_p=0;
  // Source extrema are inspected here under the model-field action category.
  for (double p:source_.spectrum.power_mpc3) {
    if (!ledger.charge(result.work.model_fields)) return fail(ledger.status);
    max_p=std::max(max_p,static_cast<W>(p));
  }
  const W callback_error=callback_diagnostic(source_.spectrum,model,g,
                                             maximum_abs_spectrum_slope_,max_p);
  if (!normal_or_zero(callback_error) || callback_error<0)
    return fail(Status::unresolved_projection);
  for (std::size_t j=0;j<nt;++j) {
    if (!ledger.charge(result.work.support_columns)) return fail(ledger.status);
    const W k=mapped_theory_k_per_mpc_[j], e=grid_projection_estimate_per_mpc_[j];
    const W lo=(k-e)*g.low_q,hi=(k+e)*g.high_q;
    if (!positive(lo) || !positive(hi)) return fail(Status::unresolved_projection);
    if (lo<source_.spectrum.k_per_mpc.front() || hi>source_.spectrum.k_per_mpc.back())
      return fail(Status::outside_support);
  }
  for (std::size_t j=0;j<nt;++j) {
    const W k=mapped_theory_k_per_mpc_[j];
    splits.clear(); splits.push_back(0);
    if (model.alpha_perp!=model.alpha_parallel) {
      const W first=g.transverse_inverse*g.transverse_inverse;
      const W difference=g.parallel_inverse*g.parallel_inverse-first;
      if (!normal_or_zero(difference) || difference==0)
        return fail(Status::unresolved_projection);
      for (double knot:source_.spectrum.k_per_mpc) {
        if (!ledger.charge(result.work.split_candidates)) return fail(ledger.status);
        if (knot<=k*g.low_q || knot>=k*g.high_q) continue;
        const W ratio=static_cast<W>(knot)/k;
        const W mu2=(ratio*ratio-first)/difference;
        const W wide=std::sqrt(mu2);
        const double split=static_cast<double>(wide);
        if (!(mu2>0 && mu2<1) || !positive(wide) || split<=0 || split>=1)
          return fail(Status::unresolved_projection);
        splits.push_back(split);
      }
      std::sort(splits.begin()+1,splits.end());
      for (std::size_t i=1;i<splits.size();++i)
        if (splits[i]<=splits[i-1]) return fail(Status::unresolved_projection);
    }
    splits.push_back(1);
    for (unsigned ell: {0u,2u}) {
      W value=0, quadrature=0, absolute_parts=0;
      Angular angular{source_.spectrum,model,g,k,ell,ledger};
      for (std::size_t cell=1;cell<splits.size();++cell) {
        // Reserve the worst charged callback/search/interpolation cost before
        // handing callbacks to the shared integrator. A denied callback is
        // never secretly invoked after the owner's action budget is exhausted.
        const std::size_t callback_actions=std::bit_width(ns)+2;
        const auto max_calls=std::min({policy.maximum_callbacks_per_integral,
                                      ledger.remaining()/callback_actions,
                                      policy.maximum_angular_callbacks-result.work.angular_callbacks});
        if (max_calls<3) return fail(Status::work_limit);
        const double width=splits[cell]-splits[cell-1];
        const double allowance=policy.maximum_multipole_absolute_error_estimate_mpc3*.5*width;
        if (!std::isnormal(allowance) || width<=0) return fail(Status::unresolved_projection);
        const auto part=numerics::integrate(angular_callback,&angular,
            splits[cell-1],splits[cell],{allowance,0,max_calls,policy.maximum_depth});
        if (ledger.status!=Status::ok) return fail(ledger.status);
        if (part.status!=numerics::Status::ok)
          return fail(part.status==numerics::Status::work_limit ? Status::work_limit : Status::unresolved_projection);
        value+=part.value; quadrature+=part.error_estimate;
        absolute_parts+=std::abs(static_cast<W>(part.value));
      }
      const W factor=ell==0 ? 1 : 5;
      const W amp=std::abs(static_cast<W>(model.b1))+std::abs(static_cast<W>(model.f));
      const W grid=factor*amp*amp*g.high_q*maximum_abs_spectrum_slope_/
                      g.volume*grid_projection_estimate_per_mpc_[j];
      const W arithmetic=factor*callback_error+64*u*(absolute_parts+quadrature);
      auto& p=ell==0 ? p0[j] : p2[j];
      if (!fill_power(value,quadrature,arithmetic,grid,p))
        return fail(Status::unresolved_projection);
      if (p.combined_estimate>policy.maximum_multipole_absolute_error_estimate_mpc3)
        return fail(Status::numerical_budget_exceeded);
    }
    if (policy.retain_input_multipoles) {
      if (!ledger.charge(result.work.output_writes)) return fail(ledger.status);
      MultipoleRow row;
      row.source_k=source_.window.theory_k[j];
      if (!project_power(p0[j],output_factor,output_factor_error,row.p0) ||
          !project_power(p2[j],output_factor,output_factor_error,row.p2))
        return fail(Status::unresolved_projection);
      if (static_cast<W>(row.p0.combined_estimate)/output_factor>
              policy.maximum_multipole_absolute_error_estimate_mpc3 ||
          static_cast<W>(row.p2.combined_estimate)/output_factor>
              policy.maximum_multipole_absolute_error_estimate_mpc3)
        return fail(Status::numerical_budget_exceeded);
      result.input_multipoles.push_back(row);
    }
  }
  for (unsigned ell: {0u,2u}) {
    const auto& axis=ell==0 ? source_.window.output0_k : source_.window.output2_k;
    const auto& ids=ell==0 ? source_.window.output0_ids : source_.window.output2_ids;
    const auto& block0=ell==0 ? source_.window.w00 : source_.window.w20;
    const auto& block2=ell==0 ? source_.window.w02 : source_.window.w22;
    for (std::size_t i=0;i<axis.size();++i) {
      W eq=0,ea=0,eg=0;
      for (std::size_t j=0;j<nt;++j) {
        for (unsigned l: {0u,2u}) {
          if (!ledger.charge(result.work.window_products)) return fail(ledger.status);
          const auto& p=l==0 ? p0[j] : p2[j];
          const W coefficient=l==0 ? block0[i*nt+j] : block2[i*nt+j];
          const W wide=coefficient*p.value;
          const double term=static_cast<double>(wide);
          if (!normal_or_zero(wide) || !std::isfinite(term) ||
              (wide!=0 && (term==0 || !std::isnormal(term))))
            return fail(Status::unresolved_projection);
          terms[2*j+(l==2)]=term;
          eq+=std::abs(coefficient)*p.quadrature_estimate;
          ea+=std::abs(coefficient)*p.arithmetic_estimate+
              std::abs(wide-static_cast<W>(term))+4*u*std::abs(wide);
          eg+=std::abs(coefficient)*p.grid_projection_estimate;
        }
      }
      const auto sum=numerics::compensated_sum(terms);
      if (sum.status!=numerics::Status::ok) return fail(Status::overflow);
      ea+=sum.error_estimate;
      PowerValue physical;
      MeanRow row{ell,ids[i],axis[i],{}};
      if (!fill_power(sum.value,eq,ea,eg,physical) ||
          !project_power(physical,output_factor,output_factor_error,row.power))
        return fail(Status::unresolved_projection);
      if (row.power.combined_estimate>policy.maximum_mean_absolute_error_estimate_output_unit)
        return fail(Status::numerical_budget_exceeded);
      if (!ledger.charge(result.work.output_writes)) return fail(ledger.status);
      result.mean.push_back(row);
    }
  }
  if (!floating.range_ok()) return fail(Status::unresolved_projection);
  result.status=Status::ok;
  return result;
}
} // namespace irred::windowed_linear_power
