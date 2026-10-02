#include "irred/gaussian_simulation.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace irred::statistics {
namespace {
using W = long double;
using N = numerics::Status;
void copied_string(detail::PayloadAccounting &b, const std::string &s) noexcept {
  // Fresh SSO assignment may grow to twice its inline capacity even when the
  // source's capacity is exactly its length. Charge a conservative copy envelope.
  const auto capacity = std::max(s.capacity(), std::string{}.capacity());
  if (capacity > (SIZE_MAX - 1) / 2) { b.add(SIZE_MAX, 2); return; }
  b.add(2 * capacity + 1, 1);
}
void copied_strings(detail::PayloadAccounting &b,
                    const std::vector<std::string> &v) noexcept {
  b.vector(v);
  for (const auto &s : v) copied_string(b, s);
}
void mean_bytes(detail::PayloadAccounting &b, const GeneratingMean &m) noexcept {
  b.vector(m.value); copied_strings(b,m.ordered_ids); copied_strings(b,m.coordinate_units);
  copied_string(b,m.identity); copied_string(b,m.coordinate_measure); copied_string(b,m.generating_law_identity);
}
void metadata_bytes(detail::PayloadAccounting &b, const Metadata &m) noexcept {
  copied_strings(b,m.ordered_ids);
  for (const auto *s : {&m.arithmetic_id, &m.measure, &m.table_identity,
       &m.uncertainty_identity, &m.ordering_provenance, &m.calibration_provenance,
       &m.dependence_provenance, &m.source_semantics, &m.input_matrix_convention,
       &m.treatment}) copied_string(b,*s);
}
bool normal(double x) { return std::isfinite(x) && (x == 0 || std::isnormal(x)); }
bool environment() { return std::numeric_limits<W>::digits >= 64 &&
  std::numeric_limits<W>::max_exponent >= 16384 && std::fegetround() == FE_TONEAREST; }
}
std::optional<size_t> GaussianSimulation::work_bound(size_t d, size_t n) noexcept {
  if (!d || d > SIZE_MAX / d) return {};
  detail::PayloadAccounting b(0);
  if (n && d*d > SIZE_MAX/n) return {};
  b.add(n*d*d, 32); b.add(n, 128); // sorting and row scope
  if (n && d > SIZE_MAX/n) return {};
  b.add(n*d, 256); // words, transcendentals, casting and diagnostics
  return b.result();
}
std::optional<size_t> GaussianSimulation::payload_bound(
    const Gaussian &g, const GeneratingMean &m, size_t count, unsigned outputs) noexcept {
  const auto d = g.metadata().ordered_ids.size();
  if (!d || (count && d > SIZE_MAX/count) || !outputs || (outputs & ~7u)) return {};
  detail::PayloadAccounting b(sizeof(GaussianSimulationBatch));
  b.add(g.retained_payload_bound().value_or(SIZE_MAX), 1);
  mean_bytes(b,m); metadata_bytes(b,g.metadata());
  b.add(count,sizeof(GaussianSimulationRow)+sizeof(random::Address));
  if(outputs&1u) b.add(count*d,sizeof(double));
  if(outputs&2u) b.add(count*d,sizeof(double));
  if(outputs&4u) b.add(count*d,sizeof(std::array<std::uint32_t,4>));
  b.add(d,sizeof(W));
  b.add(numerics::colouring_payload_bound(d).value_or(SIZE_MAX),1);
  return b.result();
}
GaussianSimulationBatch GaussianSimulation::simulate(
    const Gaussian &g, const GeneratingMean &m,
    std::span<const random::Address> addresses, std::uint64_t seed,
    GaussianSimulationPolicy policy) {
  GaussianSimulationBatch out;
  out.seed=seed; out.requested_vectors=addresses.size();
  out.dimension=g.metadata().ordered_ids.size(); out.outputs=policy.outputs;
  const auto d=out.dimension;
  if (!environment()) {out.status=N::outside_domain; return out;}
  if (g.status()!=DensityStatus::finite || !d || addresses.empty() ||
      !policy.outputs || (policy.outputs&~7u) || !policy.maximum_vectors ||
      !policy.maximum_elements || !policy.maximum_payload_bytes ||
      !policy.maximum_work_units || !std::isfinite(policy.maximum_scaled_arithmetic_error) ||
      policy.maximum_scaled_arithmetic_error<=0) return out;
  const auto bytes=payload_bound(g,m,addresses.size(),policy.outputs),
             work=work_bound(d,addresses.size());
  if (!bytes || !work || d>SIZE_MAX/d || d*d>policy.maximum_elements ||
      addresses.size()>policy.maximum_vectors || addresses.size()>65536 ||
      d>policy.maximum_elements/addresses.size() || *bytes>policy.maximum_payload_bytes ||
      *work>policy.maximum_work_units) {out.status=N::work_limit; return out;}
  if (g.metadata().source_semantics!="synthetic controls" ||
      !g.priors().empty() || !g.mean_shift().empty() ||
      m.ordered_ids!=g.metadata().ordered_ids || m.value.size()!=d ||
      m.coordinate_units.size()!=d || m.identity.empty() || m.generating_law_identity.empty() ||
      m.coordinate_measure.empty() || m.coordinate_measure!=g.metadata().measure)
    return out;
  for (size_t j=0;j<d;++j) if(!normal(m.value[j]) || m.coordinate_units[j].empty()) return out;
  for (auto a:addresses) if (d-1>UINT64_MAX-a.stream) return out;
  try {
    std::vector<random::Address> sorted(addresses.begin(),addresses.end());
    std::sort(sorted.begin(),sorted.end(),[](auto a,auto b){return a.sample<b.sample ||
      (a.sample==b.sample && a.stream<b.stream);});
    for(size_t i=1;i<sorted.size();++i)
      if(sorted[i-1].sample==sorted[i].sample && sorted[i].stream-sorted[i-1].stream<d) return out;
    out.mean=m; out.source_metadata=g.metadata(); out.work_units=*work;
    out.rows.resize(addresses.size());
    const auto elements=addresses.size()*d;
    if(policy.outputs&1u) out.values.resize(elements);
    if(policy.outputs&2u) out.absolute_error_estimates.resize(elements);
    if(policy.outputs&4u) out.words.resize(elements);
    std::vector<W> z(d);
    for(size_t i=0;i<addresses.size();++i) {
      auto &row=out.rows[i]; row.address=addresses[i];
      for(size_t j=0;j<d;++j) {
        const auto w=random::words(seed,{row.address.stream+j,row.address.sample});
        if(policy.outputs&4u) out.words[i*d+j]=w;
        z[j]=random::normal_cosine(w);
      }
      const auto v=numerics::colour(g.factor_,z,policy.maximum_elements,
        policy.maximum_payload_bytes,policy.maximum_scaled_arithmetic_error);
      row.status=v.status;
      if(row.status!=N::ok) continue;
      for(size_t j=0;j<d;++j) {
        const W exact=W(m.value[j])+v.value[j];
        const double value=static_cast<double>(exact);
        W error=v.absolute_error_estimates[j]+std::abs(exact-value)+
          4*std::numeric_limits<W>::epsilon()*(std::abs(W(m.value[j]))+std::abs(v.value[j]));
        if(!std::isfinite(exact)||!std::isfinite(value)||!std::isfinite(error)) {row.status=N::overflow;break;}
        if(!normal(value)||(exact!=0&&value==0)) {row.status=N::outside_domain;break;}
        double bound=static_cast<double>(error);
        if(W(bound)<error) bound=std::nextafter(bound,std::numeric_limits<double>::infinity());
        if(!normal(bound)||(error>0&&bound==0)) {row.status=N::outside_domain;break;}
        if(error/(1+std::abs(exact))>policy.maximum_scaled_arithmetic_error) {
          row.status=N::conditioning_budget_exceeded;break;}
        if(policy.outputs&1u) out.values[i*d+j]=value;
        if(policy.outputs&2u) out.absolute_error_estimates[i*d+j]=bound;
      }
    }
    out.status=N::ok;
    return out;
  } catch(const std::bad_alloc &) {out.status=N::work_limit;}
    catch(const std::length_error &) {out.status=N::work_limit;}
  // Batch allocation failure has no partial usable output; borrowed g survives.
  out.rows.clear(); out.values.clear();out.absolute_error_estimates.clear();out.words.clear();
  return out;
}
} // namespace irred::statistics
