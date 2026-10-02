#include "irred/baryon_abundance_law.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
constexpr W eps = std::numeric_limits<W>::epsilon();
constexpr std::size_t hard_states = 1024, hard_rows = 64,
    hard_state_rows = 65536, hard_axes = 384, hard_products = 1ull << 27,
    hard_solves = 65536, hard_iterations = 256, hard_evaluations = 4000000,
    hard_bytes = 1ull << 30;
bool profile() {
  return std::numeric_limits<W>::digits >= 64 &&
      std::numeric_limits<W>::max_exponent >= 16384 &&
      std::fegetround() == FE_TONEAREST;
}
bool policy(const BaryonAbundanceLawPolicy &p) {
  return p.maximum_states <= hard_states && p.maximum_rows <= hard_rows &&
      p.maximum_state_rows <= hard_state_rows && p.maximum_axes <= hard_axes &&
      p.maximum_covariance_products <= hard_products &&
      p.maximum_solves <= hard_solves && p.maximum_root_iterations <= hard_iterations &&
      p.maximum_charge_evaluations <= hard_evaluations &&
      p.maximum_native_bytes <= hard_bytes &&
      std::isfinite(p.mean_relative_sensitivity) && p.mean_relative_sensitivity > 0 &&
      p.mean_relative_sensitivity <= 1e-12 &&
      std::isfinite(p.covariance_relative_sensitivity) &&
      p.covariance_relative_sensitivity > 0 && p.covariance_relative_sensitivity <= 1e-8;
}
bool multiply(std::size_t a, std::size_t b, std::size_t &v) {
  if (b && a > std::numeric_limits<std::size_t>::max()/b) return false;
  v = a*b; return true;
}
void string_bound(irred::detail::PayloadAccounting &bytes, std::string_view s) {
  // Supported libstdc++ empty-string assignment can grow capacity to30 for
  // lengths16..29. Charge this copied-string envelope BEFORE any allocation;
  // the actual retained capacity check alone would be too late.
  if (s.size()==std::numeric_limits<std::size_t>::max()) {
    bytes.add(s.size(),2); return;
  }
  bytes.add(std::max<std::size_t>(32,s.size()+1),1);
}
std::optional<std::size_t> retained_payload(const BaryonAbundanceLawSource &s) {
  irred::detail::PayloadAccounting b(sizeof(BaryonAbundanceLaw)+sizeof(s));
  b.string(s.distribution_origin); b.string(s.dependence_origin);
  b.vector(s.rows); b.vector(s.states);
  for (const auto &r : s.rows) b.string(r.id);
  for (const auto &t : s.states) {
    b.string(t.id); b.string(t.temperature_origin); b.vector(t.temperature_kelvin);
    if (t.abundance.source()) {
      b.string(t.abundance.source()->source_origin);
      b.string(t.abundance.source()->mass_origin);
    }
  }
  return b.result();
}
// Neumaier wide summation. The empirical operation diagnostic is deliberately
// conservative and charged by the number and absolute size of actual terms.
struct Sum {
  W sum = 0, correction = 0, absolute = 0;
  std::size_t terms = 0;
  void add(W x) {
    const W next = sum+x;
    correction += std::abs(sum) >= std::abs(x) ? (sum-next)+x : (x-next)+sum;
    sum = next; absolute += std::abs(x); ++terms;
  }
  W value() const { return sum+correction; }
  W estimate() const { return (8*W(terms)+8)*eps*absolute; }
};
bool same(double x, double y) {
  return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y);
}
// A sufficient raw-input witness. No equality test on rounded output values is
// used, and source-origin text never changes the physical dependence.
bool constant(const BaryonAbundanceLawSource &s, const BaryonAbundanceLawAxis &axis) {
  if (s.states.size() == 1) return true;
  const bool h = axis.output == BaryonAbundanceLawOutput::hydrogen_nuclei;
  const bool he = axis.output == BaryonAbundanceLawOutput::helium_nuclei;
  const auto &first = *s.states.front().abundance.source();
  if (h || he) {
    bool absent = true;
    for (const auto &state : s.states) {
      const auto &x = *state.abundance.source();
      absent = absent && (x.physical_baryon_density == 0 ||
          (h ? x.helium4_mass_fraction == 1 : x.helium4_mass_fraction == 0));
    }
    if (absent) return true;
  }
  for (std::size_t i=1; i<s.states.size(); ++i) {
    const auto &x = *s.states[i].abundance.source();
    if (!same(x.physical_baryon_density,first.physical_baryon_density) ||
        !same(x.helium4_mass_fraction,first.helium4_mass_fraction)) return false;
    if (!he && !same(x.hydrogen1_effective_mass_kg,first.hydrogen1_effective_mass_kg)) return false;
    if (!h && !same(x.helium4_effective_mass_kg,first.helium4_effective_mass_kg)) return false;
    if (!h && !he && !same(s.states[i].temperature_kelvin[axis.row_index],
                         s.states[0].temperature_kelvin[axis.row_index])) return false;
  }
  return true;
}
// All nonzero casts retain their sign, including positive subnormal outputs.
// Round the reported absolute diagnostic upward, including its own cast loss.
// A nonzero quantity is never admitted as zero because of an absolute floor.
bool project(W value, W error, W scale, W budget, bool positive,
             double &out, double &diagnostic) {
  if (!std::isfinite(value) || !std::isfinite(error) || error < 0 ||
      !std::isfinite(scale) || scale < 0 || (positive && !(value > 0))) return false;
  out = static_cast<double>(value);
  if (!std::isfinite(out) || (value != 0 && (out == 0 || std::signbit(value) != std::signbit(out)))) return false;
  error += std::abs(W(out)-value);
  if (error > budget*scale) return false;
  diagnostic = static_cast<double>(error);
  if (W(diagnostic) < error) diagnostic = std::nextafter(diagnostic, INFINITY);
  return std::isfinite(diagnostic) && (error == 0 || diagnostic > 0) &&
      W(diagnostic) <= budget*scale;
}
struct Dimensions { std::size_t states=0, rows=0, axes=0, cells=0, values=0, products=0, bytes=0; };
template<class Result, class Batch, class Row>
bool dimensions(const BaryonAbundanceLawSource &s, const BaryonAbundanceLawPolicy &p,
                std::size_t outputs, bool lte, Dimensions &d) {
  d.states=s.states.size(); d.rows=s.rows.size();
  std::size_t state_rows=0, triangle=0;
  if (d.states>p.maximum_states || d.rows>p.maximum_rows ||
      !multiply(d.states,d.rows,state_rows) || state_rows>p.maximum_state_rows ||
      !multiply(d.rows,outputs,d.axes) || d.axes>p.maximum_axes ||
      !multiply(d.axes,d.axes,d.cells) || !multiply(d.states,d.axes,d.values) ||
      !multiply(d.axes,d.axes+1,triangle) ||
      !multiply(d.states,triangle/2,d.products) || d.products>p.maximum_covariance_products)
    return false;
  irred::detail::PayloadAccounting b(s.retained_native_payload_bytes);
  b.add(1,sizeof(Result)); b.add(d.axes,sizeof(BaryonAbundanceLawAxis));
  b.add(d.states,sizeof(Batch)); b.add(state_rows,sizeof(Row));
  b.add(d.axes,2*sizeof(double)); b.add(d.cells,2*sizeof(double));
  b.add(d.values,2*sizeof(W)); b.add(d.axes,2*sizeof(W)+sizeof(unsigned char));
  b.add(d.cells,2*sizeof(W));
  if (lte) {
    // One reused query batch and the simultaneous native LTE inputs/result.
    b.add(d.rows,sizeof(BaryonAbundanceLteQuery)+sizeof(atomic::HydrogenHeliumState)+
          sizeof(atomic::HydrogenHeliumRow));
    b.add(1,sizeof(atomic::HydrogenHeliumBatch));
    for (const auto &t:s.states) string_bound(b,t.temperature_origin);
  } else b.add(d.rows,sizeof(double));
  auto size=b.result();
  if (!size || *size>p.maximum_native_bytes) return false;
  d.bytes=*size; return true;
}
S moments(const BaryonAbundanceLawSource &source,
          const std::vector<BaryonAbundanceLawAxis> &axes, const Dimensions &d,
          const std::vector<W> &f, const std::vector<W> &error,
          const BaryonAbundanceLawPolicy &policy, BaryonAbundanceLawWork &work,
          std::optional<BaryonAbundanceLawMoments> &out) {
  if (!profile()) return S::invalid_input;
  std::vector<W> mean(d.axes), me(d.axes), covariance(d.cells), ce(d.cells);
  std::vector<unsigned char> fixed(d.axes);
  BaryonAbundanceLawMoments result;
  result.mean.resize(d.axes); result.mean_empirical_sensitivity.resize(d.axes);
  result.covariance.resize(d.cells); result.covariance_empirical_sensitivity.resize(d.cells);
  for (std::size_t i=0;i<d.axes;++i) {
    fixed[i]=constant(source,axes[i]);
    if (fixed[i]) { mean[i]=f[i]; me[i]=error[i]; }
    else {
      Sum value, sensitivity;
      for (std::size_t s=0;s<d.states;++s) {
        const auto &state=source.states[s];
        const W x=f[s*d.axes+i], e=error[s*d.axes+i];
        value.add(state.probability*x);
        sensitivity.add(state.probability*(e+state.probability_relative_arithmetic_estimate*(std::abs(x)+e)));
      }
      mean[i]=value.value(); me[i]=sensitivity.value()+sensitivity.estimate()+value.estimate();
    }
    if (!project(mean[i],me[i],std::abs(mean[i]),policy.mean_relative_sensitivity,
                 mean[i]!=0,result.mean[i],result.mean_empirical_sensitivity[i]))
      return S::conditioning_budget_exceeded;
    me[i]=result.mean_empirical_sensitivity[i];
  }
  for (std::size_t i=0;i<d.axes;++i) for (std::size_t j=0;j<=i;++j) {
    Sum value, sensitivity;
    if (!fixed[i] && !fixed[j]) for (std::size_t s=0;s<d.states;++s) {
      ++work.covariance_products;
      const auto &state=source.states[s];
      const W x=f[s*d.axes+i], y=f[s*d.axes+j],
          di=x-mean[i], dj=y-mean[j],
          ei=error[s*d.axes+i]+me[i]+4*eps*(std::abs(x)+std::abs(mean[i])),
          ej=error[s*d.axes+j]+me[j]+4*eps*(std::abs(y)+std::abs(mean[j])),
          product=di*dj,
          centered=std::abs(di)*ej+std::abs(dj)*ei+ei*ej;
      value.add(state.probability*product);
      sensitivity.add(state.probability*(centered+
          state.probability_relative_arithmetic_estimate*(std::abs(product)+centered)+
          12*eps*(std::abs(product)+centered)));
    }
    const W v=value.value(), e=sensitivity.value()+sensitivity.estimate()+value.estimate();
    covariance[i*d.axes+j]=covariance[j*d.axes+i]=v;
    ce[i*d.axes+j]=ce[j*d.axes+i]=e;
  }
  for (std::size_t i=0;i<d.axes;++i) {
    if (!fixed[i] && !(covariance[i*d.axes+i]>0)) return S::conditioning_budget_exceeded;
    for (std::size_t j=0;j<=i;++j) {
      const W scale=std::sqrt(covariance[i*d.axes+i])*std::sqrt(covariance[j*d.axes+j]);
      double v=0,e=0;
      if (!project(covariance[i*d.axes+j],ce[i*d.axes+j],scale,
                   policy.covariance_relative_sensitivity,i==j && !fixed[i],v,e))
        return S::conditioning_budget_exceeded;
      result.covariance[i*d.axes+j]=result.covariance[j*d.axes+i]=v;
      result.covariance_empirical_sensitivity[i*d.axes+j]=
          result.covariance_empirical_sensitivity[j*d.axes+i]=e;
    }
  }
  out=std::move(result); return S::ok;
}
template<class Value>
void collect(const Value &v, W &value, W &error, S &failure) {
  if (v.status!=S::ok || !v.value) {
    if (failure==S::ok) failure=v.status==S::ok ? S::invalid_input : v.status;
  } else { value=*v.value; error=std::abs(value)*W(v.relative_arithmetic_estimate); }
}
constexpr std::array<BaryonAbundanceLawOutput,6> lte_outputs{
    BaryonAbundanceLawOutput::hydrogen_neutral,BaryonAbundanceLawOutput::hydrogen_ionized,
    BaryonAbundanceLawOutput::helium_neutral,BaryonAbundanceLawOutput::helium_singly_ionized,
    BaryonAbundanceLawOutput::helium_doubly_ionized,BaryonAbundanceLawOutput::electron_density};
}
BaryonAbundanceLaw::BaryonAbundanceLaw(BaryonAbundanceLaw &&o) noexcept
    :status_(o.status_),source_(std::move(o.source_)) { o.status_=S::invalid_input; }
BaryonAbundanceLaw &BaryonAbundanceLaw::operator=(BaryonAbundanceLaw &&o) noexcept {
  if (this!=&o) { status_=o.status_; source_=std::move(o.source_); o.status_=S::invalid_input; }
  return *this;
}
BaryonAbundanceLaw prepare_baryon_abundance_law(
    const BaryonAbundanceLawInput &in, BaryonAbundanceLawPolicy p) {
  BaryonAbundanceLaw out;
  if (!policy(p)) return out;
  const auto ns=in.states.size(), nr=in.rows.size();
  if (!ns || !nr || in.distribution_origin.empty() || in.dependence_origin.empty()) return out;
  std::size_t state_rows=0;
  if (ns>p.maximum_states || nr>p.maximum_rows ||
      !multiply(ns,nr,state_rows) || state_rows>p.maximum_state_rows) { out.status_=S::work_limit; return out; }
  irred::detail::PayloadAccounting bytes(sizeof(BaryonAbundanceLaw)+sizeof(BaryonAbundanceLawSource));
  bytes.add(ns,sizeof(BaryonAbundanceLawRetainedState)); bytes.add(nr,sizeof(BaryonAbundanceLawOwnedRow));
  bytes.add(state_rows,sizeof(double));
  string_bound(bytes,in.distribution_origin); string_bound(bytes,in.dependence_origin);
  for (const auto &r:in.rows) { if (r.id.empty()) return out; string_bound(bytes,r.id); }
  for (const auto &s:in.states) {
    if (s.id.empty() || s.temperature_origin.empty() || s.abundance.source_origin.empty() ||
        s.abundance.mass_origin.empty() || s.temperature_kelvin.size()!=nr) return out;
    string_bound(bytes,s.id); string_bound(bytes,s.temperature_origin);
    string_bound(bytes,s.abundance.source_origin); string_bound(bytes,s.abundance.mass_origin);
  }
  auto size=bytes.result();
  if (!size || *size>p.maximum_native_bytes) { out.status_=S::work_limit; return out; }
  for (std::size_t i=0;i<nr;++i) for (std::size_t j=0;j<i;++j)
    if (in.rows[i].id==in.rows[j].id) return out;
  Sum normalization;
  for (std::size_t i=0;i<ns;++i) {
    for (std::size_t j=0;j<i;++j) if (in.states[i].id==in.states[j].id) return out;
    const double w=in.states[i].relative_mass;
    if (!std::isfinite(w)) { out.status_=S::nonfinite_input; return out; }
    if (!std::isnormal(w) || !(w>0)) { out.status_=S::outside_domain; return out; }
    normalization.add(W(w));
  }
  if (!profile()) return out;
  const W total=normalization.value(), weight_error=normalization.estimate()/total+4*eps;
  if (!std::isnormal(total) || !(total>0) || !std::isfinite(weight_error) || !(weight_error<1)) {
    out.status_=S::conditioning_budget_exceeded; return out;
  }
  auto source=std::make_shared<BaryonAbundanceLawSource>();
  source->preparation_native_payload_bound=*size;
  source->distribution_origin=in.distribution_origin; source->dependence_origin=in.dependence_origin;
  source->rows.reserve(nr); source->states.reserve(ns); source->normalization_terms=ns;
  out.status_=S::ok; out.source_=source;
  auto preserve=[&](S status) { if (out.status_==S::ok && status!=S::ok) out.status_=status; };
  for (const auto &r:in.rows) {
    source->rows.push_back({std::string(r.id),r.scale_factor});
    if (!std::isfinite(r.scale_factor)) preserve(S::nonfinite_input);
    else if (!(r.scale_factor>0) || r.scale_factor>1) preserve(S::outside_domain);
  }
  for (const auto &s:in.states) {
    BaryonAbundanceLawRetainedState state;
    state.id=s.id; state.relative_mass=s.relative_mass;
    state.probability=W(s.relative_mass)/total;
    state.probability_relative_arithmetic_estimate=weight_error/(1-weight_error)+4*eps;
    if (!std::isnormal(state.probability) || !(state.probability>0)) preserve(S::conditioning_budget_exceeded);
    state.temperature_origin=s.temperature_origin;
    state.temperature_kelvin.assign(s.temperature_kelvin.begin(),s.temperature_kelvin.end());
    for (double t:state.temperature_kelvin) {
      if (!std::isfinite(t)) preserve(S::nonfinite_input);
      else if (!std::isnormal(t) || !(t>0)) preserve(S::outside_domain);
    }
    ++source->state_preparations;
    state.abundance=prepare_baryon_abundance(s.abundance,{65536,p.maximum_native_bytes});
    preserve(state.abundance.status()); source->states.push_back(std::move(state));
  }
  const auto actual=retained_payload(*source);
  if (!actual || *actual>p.maximum_native_bytes) preserve(S::work_limit);
  else source->retained_native_payload_bytes=*actual;
  return out;
}
BaryonAbundanceDensityLawResult BaryonAbundanceLaw::evaluate_density(BaryonAbundanceLawPolicy p) const {
  BaryonAbundanceDensityLawResult out; out.source=source_;
  if (!policy(p)) return out;
  if (status_!=S::ok) { out.status=status_; return out; }
  Dimensions d;
  if (!dimensions<BaryonAbundanceDensityLawResult,BaryonAbundanceBatch,BaryonAbundanceRow>(*source_,p,2,false,d)) {
    out.status=S::work_limit; return out;
  }
  out.work.combined_native_payload_bound=d.bytes;
  out.axes.reserve(d.axes); out.attempts.reserve(d.states);
  std::vector<double> scales; scales.reserve(d.rows);
  for (std::size_t r=0;r<d.rows;++r) {
    scales.push_back(source_->rows[r].scale_factor);
    out.axes.push_back({r,BaryonAbundanceLawOutput::hydrogen_nuclei});
    out.axes.push_back({r,BaryonAbundanceLawOutput::helium_nuclei});
  }
  std::vector<W> f(d.values), error(d.values);
  S failure=S::ok;
  for (std::size_t s=0;s<d.states;++s) {
    auto batch=source_->states[s].abundance.evaluate(scales,{d.rows,p.maximum_native_bytes});
    out.work.mapping_rows+=batch.rows.size();
    if (batch.status!=S::ok || batch.rows.size()!=d.rows) {
      if (failure==S::ok) failure=batch.status==S::ok ? S::invalid_input : batch.status;
    } else for (std::size_t r=0;r<d.rows;++r) {
      const auto &row=batch.rows[r];
      if (row.admission_status!=S::ok && failure==S::ok) failure=row.admission_status;
      collect(row.hydrogen_nuclei,f[s*d.axes+2*r],error[s*d.axes+2*r],failure);
      collect(row.helium_nuclei,f[s*d.axes+2*r+1],error[s*d.axes+2*r+1],failure);
    }
    out.attempts.push_back(std::move(batch));
  }
  out.status=failure==S::ok ? moments(*source_,out.axes,d,f,error,p,out.work,out.moments) : failure;
  return out;
}
BaryonAbundanceLteLawResult BaryonAbundanceLaw::evaluate_equilibrium(unsigned mask,BaryonAbundanceLawPolicy p) const {
  BaryonAbundanceLteLawResult out; out.source=source_;
  if (!policy(p) || !mask || (mask&~63u)) return out;
  if (status_!=S::ok) { out.status=status_; return out; }
  Dimensions d;
  if (!dimensions<BaryonAbundanceLteLawResult,BaryonAbundanceLteBatch,BaryonAbundanceLteRow>(
      *source_,p,std::popcount(mask),true,d)) { out.status=S::work_limit; return out; }
  out.work.combined_native_payload_bound=d.bytes;
  out.axes.reserve(d.axes); out.attempts.reserve(d.states);
  for (std::size_t r=0;r<d.rows;++r) for (unsigned j=0;j<6;++j)
    if (mask&(1u<<j)) out.axes.push_back({r,lte_outputs[j]});
  std::vector<BaryonAbundanceLteQuery> queries(d.rows);
  std::vector<W> f(d.values),error(d.values);
  S failure=S::ok;
  for (std::size_t s=0;s<d.states;++s) {
    const auto &state=source_->states[s];
    for (std::size_t r=0;r<d.rows;++r) queries[r]={source_->rows[r].scale_factor,state.temperature_kelvin[r]};
    atomic::HydrogenHeliumPolicy child;
    child.requested_outputs=mask; child.maximum_rows=d.rows;
    child.maximum_solves=p.maximum_solves-out.work.solves;
    child.maximum_root_iterations=p.maximum_root_iterations;
    child.maximum_charge_evaluations=p.maximum_charge_evaluations-out.work.charge_evaluations;
    child.maximum_native_bytes=p.maximum_native_bytes;
    auto batch=state.abundance.evaluate_equilibrium(queries,state.temperature_origin,child,{d.rows,p.maximum_native_bytes});
    out.work.mapping_rows+=batch.rows.size(); out.work.solves+=batch.solves;
    out.work.root_iterations+=batch.root_iterations; out.work.charge_evaluations+=batch.charge_evaluations;
    if (batch.status!=S::ok || batch.rows.size()!=d.rows) {
      if (failure==S::ok) failure=batch.status==S::ok ? S::invalid_input : batch.status;
    } else {
      std::size_t a=0;
      for (const auto &row:batch.rows) {
        if (row.nuclei.admission_status!=S::ok && failure==S::ok) failure=row.nuclei.admission_status;
        for (const auto *v:{&row.nuclei.hydrogen_nuclei,&row.nuclei.helium_nuclei})
          if ((v->status!=S::ok || !v->value) && failure==S::ok) failure=v->status==S::ok ? S::invalid_input : v->status;
        if (row.equilibrium.admission_status!=S::ok && failure==S::ok) failure=row.equilibrium.admission_status;
        const auto &e=row.equilibrium;
        const std::array<const atomic::HydrogenHeliumValue*,6> values{
            &e.hydrogen_neutral,&e.hydrogen_ionized,&e.helium_neutral,
            &e.helium_singly_ionized,&e.helium_doubly_ionized,&e.electron_density};
        for (unsigned j=0;j<6;++j) if (mask&(1u<<j)) {
          collect(*values[j],f[s*d.axes+a],error[s*d.axes+a],failure); ++a;
        }
      }
    }
    out.attempts.push_back(std::move(batch));
  }
  out.status=failure==S::ok ? moments(*source_,out.axes,d,f,error,p,out.work,out.moments) : failure;
  return out;
}
} // namespace irred::cosmology
