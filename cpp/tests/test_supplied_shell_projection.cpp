#include "irred/supplied_shell_projection.hpp"
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <string>
#include <type_traits>
#include <vector>

namespace allocation_observation {
struct Record { void *address = nullptr; std::size_t bytes = 0; };
std::array<Record, 4096> records;
bool enabled = false;
std::size_t live = 0, peak = 0, calls = 0;
std::size_t fail_after = static_cast<std::size_t>(-1);
void *allocate(std::size_t n) {
  if (enabled && calls++ == fail_after) throw std::bad_alloc{};
  auto *p = std::malloc(n ? n : 1);
  if (!p) throw std::bad_alloc{};
  if (enabled) {
    for (auto &r : records) if (!r.address) {
      r = {p, n}; live += n; peak = std::max(peak, live); return p;
    }
    std::abort();
  }
  return p;
}
void release(void *p) noexcept {
  for (auto &r : records) if (p && r.address == p) {
    live -= r.bytes; r = {}; break;
  }
  std::free(p);
}
void begin(std::size_t failure = static_cast<std::size_t>(-1)) {
  if (live) std::abort();
  enabled = true; peak = calls = 0; fail_after = failure;
}
void end() { enabled = false; }
} // namespace allocation_observation
void *operator new(std::size_t n) { return allocation_observation::allocate(n); }
void *operator new[](std::size_t n) { return allocation_observation::allocate(n); }
void operator delete(void *p) noexcept { allocation_observation::release(p); }
void operator delete[](void *p) noexcept { allocation_observation::release(p); }
void operator delete(void *p, std::size_t) noexcept { allocation_observation::release(p); }
void operator delete[](void *p, std::size_t) noexcept { allocation_observation::release(p); }

using namespace irred::projection;
using S = irred::numerics::Status;
void need(bool good, const char *message) {
  if (!good) { std::cerr << "FAIL " << message << '\n'; std::exit(1); }
}
SuppliedShellView view(std::span<const double> k, std::span<const double> chi,
                       std::span<const double> amplitudes,
                       std::span<const std::string_view> ki,
                       std::span<const std::string_view> si,
                       std::span<const AtomRole> roles) {
  SuppliedShellView v;
  v.k_mpc_inverse = k; v.chi_mpc = chi; v.amplitudes = amplitudes;
  v.k_ids = ki; v.shell_ids = si; v.roles = roles;
  v.source_identity = "synthetic/supplied-shell-control";
  v.mode_origin = "one declared scalar unit-zeta control; no physical hierarchy";
  v.ordering_provenance = "literal original k-major/shell-minor order";
  return v;
}
std::size_t named_total(const ProjectionWork &w) {
  return w.source_scalar_inspections + w.source_tag_inspections +
      w.shell_visits + w.phase_products + w.Bessel_evaluations + w.series_terms +
      w.shell_products + w.signed_additions + w.radius_operations + w.output_casts;
}
void quote(std::string_view s) {
  std::cout << '"';
  for (unsigned char c : s) {
    if (c == '"' || c == '\\') std::cout << '\\' << static_cast<char>(c);
    else if (c < 32) std::cout << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<unsigned>(c) << std::dec;
    else std::cout << static_cast<char>(c);
  }
  std::cout << '"';
}
void bits(double v) {
  std::cout << '"' << std::hex << std::setw(16) << std::setfill('0')
      << std::bit_cast<std::uint64_t>(v) << std::dec << '"';
}
void emit_case(std::string_view id, SuppliedShellView v, std::span<const unsigned> ell, bool comma) {
  auto owner = acquire_supplied_shell_projection(v);
  auto result = owner.evaluate({ell});
  need(owner.status() == S::ok && result.status() == S::ok, "reference corpus admitted");
  if (comma) std::cout << ',';
  std::cout << "{\"id\":"; quote(id);
  std::cout << ",\"source\":{\"identity\":"; quote(v.source_identity);
  std::cout << ",\"mode_origin\":"; quote(v.mode_origin);
  std::cout << ",\"ordering_provenance\":"; quote(v.ordering_provenance);
  std::cout << ",\"mode\":\"scalar-unit-zeta\",\"units\":\"comoving-Mpc-no-h\",\"convention\":\"outward-exp-plus-i\",\"k_ids\":[";
  for (std::size_t i = 0; i < v.k_ids.size(); ++i) { if (i) std::cout << ','; quote(v.k_ids[i]); }
  std::cout << "],\"shell_ids\":[";
  for (std::size_t i = 0; i < v.shell_ids.size(); ++i) { if (i) std::cout << ','; quote(v.shell_ids[i]); }
  auto array_bits = [](std::string_view name, std::span<const double> data) {
    std::cout << "],\"" << name << "\":[";
    for (std::size_t i = 0; i < data.size(); ++i) { if (i) std::cout << ','; bits(data[i]); }
  };
  array_bits("k_bits", v.k_mpc_inverse); array_bits("chi_bits", v.chi_mpc);
  array_bits("amplitude_bits", v.amplitudes);
  std::cout << "],\"roles\":[";
  for (std::size_t i = 0; i < v.roles.size(); ++i) { if (i) std::cout << ','; quote(v.roles[i] == AtomRole::radial_atom ? "radial-atom" : "boundary-response"); }
  std::cout << ']';
  if (v.endpoint) {
    const auto &e = *v.endpoint;
    std::cout << ",\"endpoint\":{\"producer_identity\":"; quote(e.producer_identity);
    std::cout << ",\"coordinate_identity\":"; quote(e.coordinate_identity);
    std::cout << ",\"unit_identity\":"; quote(e.unit_identity);
    std::cout << ",\"support_lower_bits\":"; bits(e.support_lower);
    std::cout << ",\"support_upper_bits\":"; bits(e.support_upper);
    std::cout << ",\"survival_status\":" << static_cast<int>(e.survival_status);
    std::cout << ",\"survival_value_bits\":"; if (e.survival) bits(*e.survival); else std::cout << "null";
    std::cout << ",\"survival_radius_bits\":"; if (e.survival_absolute_diagnostic) bits(*e.survival_absolute_diagnostic); else std::cout << "null";
    std::cout << '}';
  }
  std::cout << "},\"requested_ell\":[";
  for (std::size_t i = 0; i < ell.size(); ++i) { if (i) std::cout << ','; std::cout << ell[i]; }
  std::cout << "],\"native\":{\"prepare_status\":" << static_cast<int>(owner.status())
      << ",\"batch_status\":" << static_cast<int>(result.status()) << ",\"rows\":[";
  for (std::size_t i = 0; i < result.rows().size(); ++i) {
    const auto &r = result.rows()[i]; if (i) std::cout << ',';
    std::cout << "{\"k_index\":" << r.k_index << ",\"ell_index\":" << r.ell_index
        << ",\"ell\":" << r.ell << ",\"status\":" << static_cast<int>(r.status) << ",\"value_bits\":";
    if (r.signed_value) bits(*r.signed_value); else std::cout << "null";
    std::cout << ",\"radius_bits\":"; if (r.absolute_numerical_radius) bits(*r.absolute_numerical_radius); else std::cout << "null";
    std::cout << '}';
  }
  const auto &w = result.work_after();
  std::cout << "],\"work\":{\"total\":" << w.total << ",\"series\":" << w.series_terms
      << ",\"radius\":" << w.radius_operations
      << ",\"source_scalar_inspections\":" << w.source_scalar_inspections
      << ",\"source_tag_inspections\":" << w.source_tag_inspections
      << ",\"shell_visits\":" << w.shell_visits
      << ",\"phase_products\":" << w.phase_products
      << ",\"Bessel_evaluations\":" << w.Bessel_evaluations
      << ",\"shell_products\":" << w.shell_products
      << ",\"signed_additions\":" << w.signed_additions
      << ",\"output_casts\":" << w.output_casts
      << ",\"source_bytes_inspected\":" << w.source_bytes_inspected
      << ",\"rejected_work_requests\":" << w.rejected_work_requests
      << "},\"payload_bytes\":" << result.peak_payload_bytes() << "}}";
}
void corpus() {
  const unsigned ell[]{0, 1, 2, 3, 4, 5, 6, 7, 8};
  std::cout << "{\"schema\":\"irred/supplied-shell-projection-reference/v1\",\"cases\":[";
  {
    const double k[]{0, 1, 2}, chi[]{0, .5, 1}, a[]{.25, -.5, .75, 1, -.25, .5, -.5, .75, -.25};
    const std::string_view ki[]{"zero", "unit", "two"}, si[]{"observer", "half", "unit"};
    const AtomRole roles[]{AtomRole::radial_atom, AtomRole::radial_atom, AtomRole::radial_atom};
    emit_case("endpoints-order", view(k, chi, a, ki, si, roles), ell, false);
  }
  {
    const double k[]{1}, chi[]{3, 7.5, 8}, a[]{1, .25, -.5};
    const std::string_view ki[]{"unit"}, si[]{"three", "seven-half", "eight"};
    const AtomRole roles[]{AtomRole::radial_atom, AtomRole::radial_atom, AtomRole::boundary_response};
    auto finite = view(k, chi, a, ki, si, roles);
    finite.endpoint = EndpointView{"synthetic finite producer", "producer redshift", "redshift", 0, 10, S::work_limit, .2, std::nullopt};
    emit_case("signed-oscillations-edge", finite, ell, true);
  }
  {
    const double k[]{1, 2}, chi[]{2, 2}, a[]{1, -1, .5, -.5};
    const std::string_view ki[]{"unit", "two"}, si[]{"duplicate", "duplicate"};
    const AtomRole roles[]{AtomRole::radial_atom, AtomRole::radial_atom};
    emit_case("exact-source-cancellation", view(k, chi, a, ki, si, roles), ell, true);
  }
  {
    const double k[]{1}, chi[]{0x1p-32}, a[]{.75};
    const std::string_view ki[]{"unit"}, si[]{"lower"};
    const AtomRole roles[]{AtomRole::radial_atom};
    const unsigned requested[]{0, 1, 8};
    emit_case("positive-lower-phase", view(k, chi, a, ki, si, roles), requested, true);
  }
  {
    const double k[]{.75}, chi[]{1.1, 2}, a[]{.5, -.75};
    const std::string_view ki[]{"three-quarters"}, si[]{"emitted-one-point-one", "two"};
    const AtomRole roles[]{AtomRole::radial_atom, AtomRole::radial_atom};
    emit_case("actual-product-radius", view(k, chi, a, ki, si, roles), ell, true);
  }
  {
    const double k[]{1, 1}, chi[]{1, 2.5}, a[]{.25, -.5, .25, -.5};
    const std::string_view ki[]{"duplicate", "duplicate"}, si[]{"one", "two-half"};
    const AtomRole roles[]{AtomRole::radial_atom, AtomRole::radial_atom};
    const unsigned requested[]{8, 0, 1, 4, 8};
    emit_case("duplicate-axis-and-ell", view(k, chi, a, ki, si, roles), requested, true);
  }
  {
    const double k[]{1}, chi[]{8}, a[]{0};
    const std::string_view ki[]{"unit"}, si[]{"eight"};
    const AtomRole roles[]{AtomRole::radial_atom};
    emit_case("raw-zero-source", view(k, chi, a, ki, si, roles), ell, true);
  }
  std::cout << "]}\n";
}

int main(int argc, char **argv) {
  static_assert(!std::is_copy_constructible_v<PreparedProjection> && !std::is_copy_constructible_v<ProjectionBatch>);
  if (argc == 2 && std::strcmp(argv[1], "--reference-cases") == 0) { corpus(); return 0; }
  need(argc == 1, "test arguments");
  need(supplied_shell_arithmetic_profile(), "strict p64/RN arithmetic profile");
  using W = long double;
  volatile W wide_one = 1, half_ulp = 0x1p-64L;
  const W tie_even = wide_one + half_ulp;
  need(tie_even == 1 && std::nextafter(W{1}, std::numeric_limits<W>::infinity()) == 1 + 0x1p-63L && std::nextafter(1.0, std::numeric_limits<double>::infinity()) == 1 + 0x1p-52, "actual RN tie and immediate representation neighbors");
  volatile W smallest = std::numeric_limits<W>::denorm_min();
  const W doubled_smallest = smallest * W{2};
  need(smallest != 0 && doubled_smallest != 0 && std::fpclassify(doubled_smallest) == FP_SUBNORMAL && doubled_smallest > smallest, "actual gradual wide underflow control");
  std::cout << "PROFILE double_digits=" << std::numeric_limits<double>::digits
      << " wide_digits=" << std::numeric_limits<W>::digits
      << " wide_min_exponent=" << std::numeric_limits<W>::min_exponent
      << " wide_max_exponent=" << std::numeric_limits<W>::max_exponent
      << " rounding=" << std::fegetround() << " sizeof_source=" << sizeof(SuppliedShellSource)
      << " sizeof_row=" << sizeof(ProjectionRow) << " sizeof_prepared_handle=" << sizeof(PreparedProjection)
      << " sizeof_result_handle=" << sizeof(ProjectionBatch) << '\n';
  double k[]{0, 1, 1}, chi[]{0, 2}, a[]{.25, .5, 1, -.25, 1, -.25};
  const std::string_view ki[]{"zero", "duplicate", "duplicate"}, si[]{"observer", "shell"};
  const AtomRole roles[]{AtomRole::radial_atom, AtomRole::boundary_response};
  auto v = view(k, chi, a, ki, si, roles);
  const unsigned ell[]{1, 0, 8, 0};
  allocation_observation::begin();
  {
    auto owner = acquire_supplied_shell_projection(v);
    need(owner.status() == S::ok && owner.source(), "acquired immutable source");
    const auto *source = owner.source();
    need(source->live_payload_bytes() == allocation_observation::live, "actual requested retained acquisition bytes");
    const auto prep = owner.cumulative_work();
    k[1] = -99; chi[1] = std::numeric_limits<double>::quiet_NaN(); a[2] = 0;
    auto r = owner.evaluate({ell});
    need(r.status() == S::ok && r.rows().size() == 12 && r.source() == source, "coarse output same source");
    need(source->k_mpc_inverse()[1] == 1 && source->chi_mpc()[1] == 2 && source->amplitudes()[2] == 1, "original input copied once");
    need(*r.rows()[0].signed_value == 0 && *r.rows()[1].signed_value == .75 && *r.rows()[2].signed_value == 0, "zero argument exact kernels and signed monopole");
    need(r.rows()[1].absolute_numerical_radius && r.rows()[1].diagnostics, "mandatory separated diagnostics");
    need(r.rows()[4].signed_value == r.rows()[8].signed_value && r.rows()[1].signed_value == r.rows()[3].signed_value, "duplicate k and ell order");
    need(r.work_before().total == prep.total && r.work_after().total == prep.total + r.work_delta().total, "complete PREP plus EVAL work");
    need(named_total(r.work_after()) == r.work_after().total && r.work_delta().phase_products == 6, "all named units and one phase per pair");
    need(source->live_payload_bytes() == allocation_observation::live && source->peak_payload_bytes() == allocation_observation::peak, "actual requested evaluation peak and scratch release");
    auto previous = std::move(r);
    auto *self = &previous; previous = std::move(*self);
    need(previous.source() == source && r.source() == nullptr && previous.rows().size() == 12, "result move self move");
    auto moved = std::move(owner); auto *owner_self = &moved; moved = std::move(*owner_self);
    need(owner.status() == S::invalid_input && !owner.source() && moved.source() == source, "prepared move self move");
    auto next = moved.evaluate({ell});
    need(next.source() == source && previous.rows()[4].signed_value == next.rows()[4].signed_value, "prior next result coexistence");
    need(source->live_payload_bytes() == allocation_observation::live && source->peak_payload_bytes() == allocation_observation::peak, "whole live multiple results exact request bytes");
    moved = PreparedProjection{};
    need(previous.source() == source && next.source() == source && *previous.rows()[1].signed_value == .75, "results outlive prepared source owner");
  }
  need(allocation_observation::live == 0, "all source-family allocation destruction");
  allocation_observation::end();
  k[1] = 1; chi[1] = 2; a[2] = 1;
  auto owner = acquire_supplied_shell_projection(v);
  const unsigned one[]{0};
  need(owner.evaluate({one, 0}).status() == S::invalid_input, "empty mask");
  need(owner.evaluate({one, 2}).status() == S::invalid_input, "unknown mask");
  need(owner.evaluate({one, 1, 33}).status() == S::outside_domain, "accuracy tightening domain");
  const unsigned bad_ell[]{9};
  auto bad_request = owner.evaluate({bad_ell});
  need(bad_request.status() == S::outside_domain && bad_request.failure().ell_index == 0 && bad_request.failure().reason == Refusal::request, "unsupported multipole retains requested index");
  auto wrong = v; wrong.units = static_cast<GeometryUnits>(42);
  need(acquire_supplied_shell_projection(wrong).status() == S::outside_domain, "wrong source units");
  wrong = v; wrong.mode = static_cast<SourceMode>(42);
  need(acquire_supplied_shell_projection(wrong).status() == S::outside_domain, "wrong source mode");
  wrong = v; wrong.convention = static_cast<FourierConvention>(42);
  need(acquire_supplied_shell_projection(wrong).status() == S::outside_domain, "wrong Fourier convention");
  const AtomRole wrong_roles[]{AtomRole::radial_atom, static_cast<AtomRole>(42)};
  wrong = v; wrong.roles = wrong_roles;
  auto wrong_role = acquire_supplied_shell_projection(wrong);
  need(wrong_role.status() == S::outside_domain && wrong_role.failure().shell_index == 1, "unsupported role retains original index");
  wrong = v; wrong.amplitudes = std::span(a, 1);
  need(acquire_supplied_shell_projection(wrong).status() == S::invalid_input, "table shape");
  wrong = v; wrong.source_identity = std::string_view("bad\xff", 4);
  auto invalid_utf8 = acquire_supplied_shell_projection(wrong);
  need(invalid_utf8.status() == S::invalid_input && invalid_utf8.cumulative_work().source_bytes_inspected == 4, "failed UTF8 byte prefix retained");
  wrong = v; wrong.source_identity = std::string_view("x\0y", 3);
  auto nul = acquire_supplied_shell_projection(wrong);
  need(nul.status() == S::invalid_input && nul.cumulative_work().source_bytes_inspected == 2, "NUL refusal retains inspected byte prefix");
  a[0] = 1.01;
  need(acquire_supplied_shell_projection(v).status() == S::outside_domain && a[0] == 1.01, "amplitude domain original preserved"); a[0] = .25;
  chi[1] = std::numeric_limits<double>::quiet_NaN();
  need(acquire_supplied_shell_projection(v).status() == S::nonfinite_input && std::isnan(chi[1]), "failed acquisition original input unchanged"); chi[1] = 2;
  auto metadata = v;
  metadata.endpoint = EndpointView{"synthetic finite history", "producer redshift", "redshift", 0, 10, S::work_limit, .2, std::nullopt};
  auto finite = acquire_supplied_shell_projection(metadata);
  auto f = finite.evaluate({one}); auto central = owner.evaluate({one});
  need(f.status() == S::ok && f.rows()[1].signed_value == central.rows()[1].signed_value && finite.source()->endpoint()->survival_status == S::work_limit && !finite.source()->endpoint()->survival_absolute_diagnostic, "literal failed upstream metadata no renormalization");
  metadata.endpoint->support_lower = 11;
  need(acquire_supplied_shell_projection(metadata).status() == S::outside_domain, "producer support order refusal");
  metadata.endpoint->support_lower = 0; metadata.endpoint->survival = 1.01;
  need(acquire_supplied_shell_projection(metadata).status() == S::outside_domain, "survival mass domain refusal");
  const int rounding = std::fegetround();
  need(std::fesetround(FE_DOWNWARD) == 0, "set hostile rounding");
  auto hostile = owner.evaluate({one});
  need(hostile.status() == S::outside_domain && hostile.work_delta().phase_products == 0 && hostile.rows().empty(), "rounding refusal before nominal arithmetic");
  need(std::fesetround(rounding) == 0, "restore original rounding");
  {
    const double tk[]{1}, tc[]{1}, ta[]{std::numeric_limits<double>::denorm_min()};
    const std::string_view tid[]{"tiny"};
    const AtomRole tr[]{AtomRole::radial_atom};
    const unsigned high_ell[]{8};
    auto tiny_owner = acquire_supplied_shell_projection(view(tk, tc, ta, tid, tid, tr));
    need(tiny_owner.status() == S::ok, "exact subnormal source imports into normal wide arithmetic");
    auto tiny_result = tiny_owner.evaluate({high_ell});
    need(tiny_result.status() == S::overflow && tiny_result.rows()[0].refusal == Refusal::nonnormal_arithmetic && !tiny_result.rows()[0].signed_value && !tiny_result.rows()[0].diagnostics, "positive unresolved public storage refuses instead of zero");
  }
  const double small_k[]{0x1p-33, 1, 2}, small_chi[]{1}, unit_a[]{1, 1, 1};
  const std::string_view small_ki[]{"below", "unit", "two"}, small_si[]{"unit"};
  const AtomRole atom[]{AtomRole::radial_atom};
  {
    const double upper_k[]{1}, upper_chi[]{8.01}, zero_a[]{0};
    const std::string_view upper_ids[]{"outside"};
    auto upper_owner = acquire_supplied_shell_projection(view(upper_k, upper_chi, zero_a, upper_ids, upper_ids, atom));
    auto upper = upper_owner.evaluate({one});
    need(upper.status() == S::outside_domain && upper.work_delta().phase_products == 1 && upper.work_delta().Bessel_evaluations == 0 && !upper.rows()[0].signed_value, "zero amplitude still guards upper phase");
  }
  auto phases = acquire_supplied_shell_projection(view(small_k, small_chi, unit_a, small_ki, small_si, atom));
  auto phase_rows = phases.evaluate({ell});
  need(phase_rows.rows().size() == 12 && phase_rows.rows()[0].status == S::outside_domain && !phase_rows.rows()[0].diagnostics && phase_rows.rows()[4].status == S::ok, "failed phase retains every ell and other k");
  {
    const double mk[]{1}, mc[]{0x1p-33, 1}, ma[]{0, 0};
    const std::string_view mid[]{"k"}, msid[]{"c0", "c1"};
    const AtomRole mr[]{AtomRole::radial_atom, AtomRole::radial_atom};
    auto mixed_view = view(mk, mc, ma, mid, msid, mr);
    mixed_view.source_identity = "s"; mixed_view.mode_origin = "m"; mixed_view.ordering_provenance = "o";
    ProjectionResources mixed_cap; mixed_cap.maximum_total_work = 26;
    auto mixed_owner = acquire_supplied_shell_projection(mixed_view, mixed_cap);
    need(mixed_owner.status() == S::ok && mixed_owner.cumulative_work().total == 21, "minimal source exact preparation work");
    auto mixed = mixed_owner.evaluate({one});
    need(mixed.status() == S::outside_domain && mixed.rows()[0].status == S::outside_domain && mixed.rows()[0].refusal == Refusal::phase && mixed.failure().shell_index == 0 && mixed.terminal_work_status() == S::work_limit && mixed.work_after().total == 26 && mixed.work_after().rejected_work_requests == 1, "first domain cause and later terminal quota both retained");
  }
  const double oscillatory_k[]{1}, oscillatory_chi[]{8}, oscillatory_a[]{1};
  const std::string_view single_id[]{"one"};
  auto oscillatory_view = view(oscillatory_k, oscillatory_chi, oscillatory_a, single_id, single_id, atom);
  auto tight = acquire_supplied_shell_projection(oscillatory_view);
  auto tighter = tight.evaluate({one, 1, 32});
  need(tighter.status() == S::conditioning_budget_exceeded && !tighter.rows()[0].signed_value, "unattainable tightened allocation refused");
  ProjectionResources cap; cap.maximum_series_terms = 1;
  auto term_owner = acquire_supplied_shell_projection(oscillatory_view, cap);
  auto term_failure = term_owner.evaluate({one});
  need(term_failure.status() == S::work_limit && term_failure.rows()[0].refusal == Refusal::work_limit && term_failure.failure().reason == Refusal::work_limit && term_failure.terminal_work_status() == S::work_limit && term_failure.work_after().series_terms == 1 && term_failure.work_after().rejected_work_requests > 0, "failed term attempt charged with actual quota cause");
  auto repeated = term_owner.evaluate({one});
  need(repeated.work_before().series_terms == 1 && repeated.work_after().series_terms == 1, "repeated failure retains original exhausted term cap");
  cap = {}; cap.maximum_source_bytes = 1;
  need(acquire_supplied_shell_projection(v, cap).status() == S::work_limit, "source byte cap preflight");
  cap = {}; cap.maximum_total_work = 1;
  auto prep_failure = acquire_supplied_shell_projection(v, cap);
  need(prep_failure.status() == S::work_limit && prep_failure.cumulative_work().total == 1, "preparation total cap includes refused work prefix");
  cap.maximum_total_work = 0;
  auto no_inspection = acquire_supplied_shell_projection(v, cap);
  need(no_inspection.status() == S::work_limit && no_inspection.cumulative_work().total == 0 && no_inspection.cumulative_work().rejected_work_requests == 1, "zero quota refuses before policy inspection");
  cap = {}; cap.maximum_total_work = 8000001;
  auto excessive = acquire_supplied_shell_projection(v, cap);
  need(excessive.status() == S::invalid_input && excessive.cumulative_work().total == 1, "above hard cap records actual first field inspection");
  cap = {}; cap.maximum_series_terms = 0;
  need(acquire_supplied_shell_projection(v, cap).cumulative_work().total == 2, "invalid second policy field actual prefix");
  cap = {}; cap.maximum_source_bytes = 0;
  need(acquire_supplied_shell_projection(v, cap).cumulative_work().total == 3, "invalid third policy field actual prefix");
  cap = {}; cap.maximum_live_bytes = 0;
  need(acquire_supplied_shell_projection(v, cap).cumulative_work().total == 4, "invalid fourth policy field actual prefix");
  allocation_observation::begin(2);
  auto allocation_failure = acquire_supplied_shell_projection(v);
  allocation_observation::end();
  need(allocation_failure.status() == S::work_limit && allocation_failure.cumulative_work().source_scalar_inspections > 0 && allocation_observation::live == 0 && allocation_failure.peak_payload_bytes() == allocation_observation::peak, "failed acquisition allocations prefix preserved and unwound");
  auto payload_probe = acquire_supplied_shell_projection(oscillatory_view);
  auto first_probe = payload_probe.evaluate({ell});
  cap = {}; cap.maximum_live_bytes = first_probe.peak_payload_bytes();
  auto bounded = acquire_supplied_shell_projection(oscillatory_view, cap);
  auto first = bounded.evaluate({ell});
  need(first.status() == S::ok, "first result fits exact observed live envelope");
  auto over = bounded.evaluate({ell});
  need(over.status() == S::work_limit && over.rows().empty() && first.source() == over.source() && *first.rows()[0].signed_value == *first_probe.rows()[0].signed_value && bounded.source()->peak_payload_bytes() <= cap.maximum_live_bytes, "prior result plus next scratch whole live refusal");
  first = ProjectionBatch{};
  auto released = bounded.evaluate({ell});
  need(released.status() == S::ok, "destruction releases actual whole-live reservation");
  std::cout << "PASS supplied-shell endpoints/order/source/lifetime/refusal/work/allocation/rounding controls\n";
}
