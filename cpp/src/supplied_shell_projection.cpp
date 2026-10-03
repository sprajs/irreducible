#include "irred/supplied_shell_projection.hpp"
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <utility>

namespace irred::projection {
using S = numerics::Status;
using W = long double;
namespace detail {
template <class T> void TypedDelete<T>::operator()(T *p) const noexcept {
  if (p) {
    std::destroy_n(p, count);
    std::allocator<T>{}.deallocate(p, count);
  }
}
struct PreparedState {
  ProjectionResources resources;
  ProjectionWork work;
};
struct ResultState {
  std::array<unsigned, 9> ell{};
  std::size_t ell_count = 0, row_count = 0;
  TypedArray<ProjectionRow> rows;
};
struct SourceStorage {
  std::size_t K = 0, S = 0, B = 0, maximum_live = 0;
  TypedArray<double> k, chi, amplitudes;
  TypedArray<IdentifierOffset> ids;
  TypedArray<AtomRole> roles;
  TypedArray<char> bytes;
  std::array<IdentifierOffset, 6> labels{};
  std::optional<EndpointView> endpoint;
};
struct SourceMaker final : SuppliedShellSource {
  explicit SourceMaker(SourceStorage &&s) : SuppliedShellSource(std::move(s)) {}
};
struct SourceAccess {
  static void seal(SourceMaker &s, std::size_t bytes) noexcept {
    s.payload_ = bytes;
    s.live_.store(bytes);
    s.peak_.store(bytes);
  }
  static bool reserve(const SuppliedShellSource &s, std::size_t n) noexcept {
    return s.reserve(n);
  }
  static void release(const SuppliedShellSource &s, std::size_t n) noexcept {
    s.release(n);
  }
};
} // namespace detail
namespace {
bool add_size(std::size_t &a, std::size_t b) noexcept {
  if (b > std::numeric_limits<std::size_t>::max() - a) return false;
  a += b;
  return true;
}
bool mul_size(std::size_t a, std::size_t b, std::size_t &v) noexcept {
  if (b && a > std::numeric_limits<std::size_t>::max() / b) return false;
  v = a * b;
  return true;
}
struct Ledger {
  ProjectionWork &w;
  const ProjectionResources &p;
  S status = S::ok;
  bool charge(std::size_t ProjectionWork::*field) noexcept {
    if (status != S::ok) return false;
    if (w.total >= p.maximum_total_work ||
        w.*field == std::numeric_limits<std::size_t>::max() ||
        (field == &ProjectionWork::series_terms &&
         w.series_terms >= p.maximum_series_terms)) {
      if (w.rejected_work_requests != std::numeric_limits<std::size_t>::max())
        ++w.rejected_work_requests;
      status = S::work_limit;
      return false;
    }
    ++(w.*field);
    ++w.total;
    return true;
  }
  bool bytes(std::size_t n) noexcept {
    if (status != S::ok) return false;
    if (n > p.maximum_source_bytes - w.source_bytes_inspected) {
      if (w.rejected_work_requests != std::numeric_limits<std::size_t>::max())
        ++w.rejected_work_requests;
      status = S::work_limit;
      return false;
    }
    w.source_bytes_inspected += n;
    return true;
  }
};
ProjectionWork difference(ProjectionWork a, const ProjectionWork &b) noexcept {
#define IRRED_SHELL_DELTA(name) a.name -= b.name
  IRRED_SHELL_DELTA(source_scalar_inspections);
  IRRED_SHELL_DELTA(source_tag_inspections);
  IRRED_SHELL_DELTA(shell_visits); IRRED_SHELL_DELTA(phase_products);
  IRRED_SHELL_DELTA(Bessel_evaluations); IRRED_SHELL_DELTA(series_terms);
  IRRED_SHELL_DELTA(shell_products); IRRED_SHELL_DELTA(signed_additions);
  IRRED_SHELL_DELTA(radius_operations); IRRED_SHELL_DELTA(output_casts);
  IRRED_SHELL_DELTA(total); IRRED_SHELL_DELTA(source_bytes_inspected);
  IRRED_SHELL_DELTA(rejected_work_requests);
#undef IRRED_SHELL_DELTA
  return a;
}
bool utf8(std::string_view s, Ledger &ledger) noexcept {
  for (std::size_t i = 0; i < s.size();) {
    if (!ledger.bytes(1)) return false;
    const auto c = static_cast<unsigned char>(s[i++]);
    if (!c) return false;
    if (c < 0x80) continue;
    unsigned more = 0, value = 0, minimum = 0;
    if (c >= 0xc2 && c <= 0xdf) { more = 1; value = c & 31; minimum = 0x80; }
    else if (c >= 0xe0 && c <= 0xef) { more = 2; value = c & 15; minimum = 0x800; }
    else if (c >= 0xf0 && c <= 0xf4) { more = 3; value = c & 7; minimum = 0x10000; }
    else return false;
    if (more > s.size() - i) return false;
    while (more--) {
      if (!ledger.bytes(1)) return false;
      const auto d = static_cast<unsigned char>(s[i++]);
      if ((d & 0xc0) != 0x80) return false;
      value = (value << 6) | (d & 63);
    }
    if (value < minimum || value > 0x10ffff ||
        (value >= 0xd800 && value <= 0xdfff)) return false;
  }
  return true;
}
struct PayloadLimit {};
struct AllocationCapture {
  std::size_t limit = 0, arrays = 0, control = 0, peak = 0;
};
template <class T> struct CaptureAllocator {
  using value_type = T;
  AllocationCapture *capture = nullptr;
  CaptureAllocator() noexcept = default;
  explicit CaptureAllocator(AllocationCapture *c) noexcept : capture(c) {}
  template <class U> CaptureAllocator(const CaptureAllocator<U> &a) noexcept
      : capture(a.capture) {}
  T *allocate(std::size_t n) {
    if (!capture) throw PayloadLimit{};
    std::size_t bytes = 0, total = capture->arrays;
    if (!mul_size(n, sizeof(T), bytes) || !add_size(total, capture->control) ||
        !add_size(total, bytes) || total > capture->limit) throw PayloadLimit{};
    T *p = std::allocator<T>{}.allocate(n);
    capture->control += bytes;
    capture->peak = std::max(capture->peak, total);
    return p;
  }
  // The construction-only capture is never dereferenced during destruction.
  void deallocate(T *p, std::size_t n) noexcept {
    std::allocator<T>{}.deallocate(p, n);
  }
  template <class U> bool operator==(const CaptureAllocator<U> &a) const noexcept {
    return capture == a.capture;
  }
};
template <class T> detail::TypedArray<T> allocate_array(std::size_t n) {
  if (!n) return {};
  auto *p = std::allocator<T>{}.allocate(n);
  try { std::uninitialized_value_construct_n(p, n); }
  catch (...) { std::allocator<T>{}.deallocate(p, n); throw; }
  return detail::TypedArray<T>(p, detail::TypedDelete<T>{n});
}
template <class T> detail::TypedArray<T> source_array(
    std::size_t n, AllocationCapture &c) {
  std::size_t bytes = 0, total = c.arrays;
  if (!mul_size(n, sizeof(T), bytes) || !add_size(total, bytes) ||
      total > c.limit) throw PayloadLimit{};
  auto p = allocate_array<T>(n);
  c.arrays = total;
  c.peak = std::max(c.peak, total);
  return p;
}
} // namespace

SuppliedShellSource::SuppliedShellSource(detail::SourceStorage &&s)
    : k_count_(s.K), shell_count_(s.S), byte_count_(s.B),
      maximum_live_(s.maximum_live), k_(std::move(s.k)), chi_(std::move(s.chi)),
      amplitudes_(std::move(s.amplitudes)), ids_(std::move(s.ids)),
      roles_(std::move(s.roles)), bytes_(std::move(s.bytes)), labels_(s.labels),
      endpoint_(s.endpoint) {
  if (endpoint_) {
    endpoint_->producer_identity = text(labels_[3]);
    endpoint_->coordinate_identity = text(labels_[4]);
    endpoint_->unit_identity = text(labels_[5]);
  }
}
SuppliedShellSource::~SuppliedShellSource() = default;
std::string_view SuppliedShellSource::text(detail::IdentifierOffset i) const noexcept {
  return i.length ? std::string_view(bytes_.get() + i.offset, i.length) : std::string_view{};
}
std::span<const double> SuppliedShellSource::k_mpc_inverse() const noexcept { return {k_.get(), k_count_}; }
std::span<const double> SuppliedShellSource::chi_mpc() const noexcept { return {chi_.get(), shell_count_}; }
std::span<const double> SuppliedShellSource::amplitudes() const noexcept { return {amplitudes_.get(), k_count_ * shell_count_}; }
std::span<const AtomRole> SuppliedShellSource::roles() const noexcept { return {roles_.get(), shell_count_}; }
std::string_view SuppliedShellSource::k_id(std::size_t i) const noexcept { return i < k_count_ ? text(ids_[i]) : std::string_view{}; }
std::string_view SuppliedShellSource::shell_id(std::size_t i) const noexcept { return i < shell_count_ ? text(ids_[k_count_ + i]) : std::string_view{}; }
std::string_view SuppliedShellSource::source_identity() const noexcept { return text(labels_[0]); }
std::string_view SuppliedShellSource::mode_origin() const noexcept { return text(labels_[1]); }
std::string_view SuppliedShellSource::ordering_provenance() const noexcept { return text(labels_[2]); }
std::optional<EndpointView> SuppliedShellSource::endpoint() const noexcept { return endpoint_; }
bool SuppliedShellSource::has_boundary_response() const noexcept {
  return std::find(roles_.get(), roles_.get() + shell_count_, AtomRole::boundary_response) != roles_.get() + shell_count_;
}
std::size_t SuppliedShellSource::source_payload_bytes() const noexcept { return payload_; }
std::size_t SuppliedShellSource::live_payload_bytes() const noexcept { return live_.load(); }
std::size_t SuppliedShellSource::peak_payload_bytes() const noexcept { return peak_.load(); }
bool SuppliedShellSource::reserve(std::size_t n) const noexcept {
  auto old = live_.load();
  for (;;) {
    if (old > maximum_live_ || n > maximum_live_ - old) return false;
    if (live_.compare_exchange_weak(old, old + n)) break;
  }
  auto peak = peak_.load();
  while (peak < old + n && !peak_.compare_exchange_weak(peak, old + n)) {}
  return true;
}
void SuppliedShellSource::release(std::size_t n) const noexcept {
  auto old = live_.load();
  for (;;) {
    if (n > old || old - n < payload_) std::terminate();
    if (live_.compare_exchange_weak(old, old - n)) return;
  }
}
bool supplied_shell_arithmetic_profile() noexcept {
#ifdef __FAST_MATH__
  return false;
#else
  return std::numeric_limits<double>::radix == 2 &&
      std::numeric_limits<double>::digits == 53 &&
      std::numeric_limits<double>::min_exponent == -1021 &&
      std::numeric_limits<double>::max_exponent == 1024 &&
      std::numeric_limits<W>::radix == 2 && std::numeric_limits<W>::digits == 64 &&
      std::numeric_limits<W>::min_exponent == -16381 &&
      std::numeric_limits<W>::max_exponent == 16384 &&
      std::fegetround() == FE_TONEAREST;
#endif
}

PreparedProjection::PreparedProjection() = default;
PreparedProjection::~PreparedProjection() { reset(); }
void PreparedProjection::reset() noexcept {
  state_.reset();
  if (source_ && charge_) detail::SourceAccess::release(*source_, charge_);
  charge_ = 0; source_.reset(); fallback_ = S::invalid_input;
}
PreparedProjection::PreparedProjection(PreparedProjection &&o) noexcept { *this = std::move(o); }
PreparedProjection &PreparedProjection::operator=(PreparedProjection &&o) noexcept {
  if (this != &o) {
    reset(); source_ = std::move(o.source_); state_ = std::move(o.state_);
    fallback_ = o.fallback_; fallback_work_ = o.fallback_work_; failure_ = o.failure_;
    charge_ = std::exchange(o.charge_, 0); fallback_peak_ = o.fallback_peak_;
    o.fallback_ = S::invalid_input; o.fallback_work_ = {}; o.failure_ = {}; o.fallback_peak_ = 0;
  }
  return *this;
}
S PreparedProjection::status() const noexcept { return fallback_; }
const SuppliedShellSource *PreparedProjection::source() const noexcept { return source_.get(); }
const ProjectionWork &PreparedProjection::cumulative_work() const noexcept { return state_ ? state_->work : fallback_work_; }
Failure PreparedProjection::failure() const noexcept { return failure_; }
std::size_t PreparedProjection::peak_payload_bytes() const noexcept { return source_ ? source_->peak_payload_bytes() : fallback_peak_; }

ProjectionBatch::ProjectionBatch() = default;
ProjectionBatch::~ProjectionBatch() { reset(); }
void ProjectionBatch::reset() noexcept {
  state_.reset();
  if (source_ && charge_) detail::SourceAccess::release(*source_, charge_);
  charge_ = 0; source_.reset(); fallback_ = S::invalid_input;
}
ProjectionBatch::ProjectionBatch(ProjectionBatch &&o) noexcept { *this = std::move(o); }
ProjectionBatch &ProjectionBatch::operator=(ProjectionBatch &&o) noexcept {
  if (this != &o) {
    reset(); source_ = std::move(o.source_); state_ = std::move(o.state_);
    fallback_ = o.fallback_; before_ = o.before_; after_ = o.after_; delta_ = o.delta_;
    failure_ = o.failure_; charge_ = std::exchange(o.charge_, 0); peak_ = o.peak_;
    o.fallback_ = S::invalid_input; o.before_ = {}; o.after_ = {}; o.delta_ = {}; o.failure_ = {}; o.peak_ = 0;
  }
  return *this;
}
S ProjectionBatch::status() const noexcept { return fallback_; }
const SuppliedShellSource *ProjectionBatch::source() const noexcept { return source_.get(); }
std::span<const ProjectionRow> ProjectionBatch::rows() const noexcept { return state_ ? std::span<const ProjectionRow>(state_->rows.get(), state_->row_count) : std::span<const ProjectionRow>{}; }
std::span<const unsigned> ProjectionBatch::requested_ell() const noexcept { return state_ ? std::span<const unsigned>(state_->ell.data(), state_->ell_count) : std::span<const unsigned>{}; }
const ProjectionWork &ProjectionBatch::work_before() const noexcept { return before_; }
const ProjectionWork &ProjectionBatch::work_after() const noexcept { return after_; }
const ProjectionWork &ProjectionBatch::work_delta() const noexcept { return delta_; }
Failure ProjectionBatch::failure() const noexcept { return failure_; }
std::size_t ProjectionBatch::retained_payload_bytes() const noexcept { return charge_; }
std::size_t ProjectionBatch::peak_payload_bytes() const noexcept { return peak_; }

PreparedProjection acquire_supplied_shell_projection(SuppliedShellView v, ProjectionResources p) {
  PreparedProjection out;
  Ledger ledger{out.fallback_work_, p};
  auto fail = [&](S s, Refusal r, std::size_t k = static_cast<std::size_t>(-1), std::size_t sh = static_cast<std::size_t>(-1)) {
    out.fallback_ = ledger.status == S::ok ? s : ledger.status;
    out.failure_ = {r, k, sh, static_cast<std::size_t>(-1)};
  };
  if (!p.maximum_total_work || p.maximum_total_work > 8000000 ||
      !p.maximum_series_terms || p.maximum_series_terms > 6000000 ||
      !p.maximum_source_bytes || p.maximum_source_bytes > 65536 ||
      p.maximum_live_bytes < sizeof(detail::PreparedState) || p.maximum_live_bytes > 4 * 1024 * 1024) {
    fail(S::invalid_input, Refusal::resource_policy); return out;
  }
  for (unsigned field = 0; field < 4; ++field)
    if (!ledger.charge(&ProjectionWork::source_tag_inspections)) { fail(S::work_limit, Refusal::resource_policy); return out; }
  if (!ledger.charge(&ProjectionWork::source_tag_inspections)) { fail(S::work_limit, Refusal::source_shape); return out; }
  const auto K = v.k_mpc_inverse.size(), N = v.chi_mpc.size();
  std::size_t KS = 0;
  if (!K || K > 256 || !N || N > 64 || !mul_size(K, N, KS) ||
      v.amplitudes.size() != KS || v.k_ids.size() != K || v.shell_ids.size() != N || v.roles.size() != N) {
    fail(S::invalid_input, Refusal::source_shape); return out;
  }
  for (unsigned field = 0; field < 3; ++field) {
    if (!ledger.charge(&ProjectionWork::source_tag_inspections)) { fail(S::work_limit, Refusal::source_tags); return out; }
    const bool valid = field == 0 ? v.mode == SourceMode::scalar_unit_zeta :
        field == 1 ? v.units == GeometryUnits::comoving_mpc_no_h :
        v.convention == FourierConvention::outward_exp_plus_i;
    if (!valid) { fail(S::outside_domain, Refusal::source_tags); return out; }
  }
  for (std::size_t j = 0; j < N; ++j) {
    if (!ledger.charge(&ProjectionWork::source_tag_inspections)) { fail(S::work_limit, Refusal::source_tags, -1, j); return out; }
    if (v.roles[j] != AtomRole::radial_atom && v.roles[j] != AtomRole::boundary_response) { fail(S::outside_domain, Refusal::source_tags, -1, j); return out; }
  }
  std::size_t B = 0;
  auto inspect_text = [&](std::string_view s) {
    if (!ledger.charge(&ProjectionWork::source_tag_inspections)) return false;
    if (s.empty()) return false;
    if (s.size() > p.maximum_source_bytes - ledger.w.source_bytes_inspected) {
      if (ledger.w.rejected_work_requests != std::numeric_limits<std::size_t>::max()) ++ledger.w.rejected_work_requests;
      ledger.status = S::work_limit; return false;
    }
    if (!utf8(s, ledger) || !add_size(B, s.size())) return false;
    return true;
  };
  if (!inspect_text(v.source_identity) || !inspect_text(v.mode_origin) || !inspect_text(v.ordering_provenance)) { fail(S::invalid_input, Refusal::source_strings); return out; }
  for (auto s : v.k_ids) if (!inspect_text(s)) { fail(S::invalid_input, Refusal::source_strings); return out; }
  for (auto s : v.shell_ids) if (!inspect_text(s)) { fail(S::invalid_input, Refusal::source_strings); return out; }
  auto number = [&](double x, bool nonnegative) {
    if (!ledger.charge(&ProjectionWork::source_scalar_inspections)) return false;
    if (!std::isfinite(x)) { ledger.status = S::nonfinite_input; return false; }
    if (nonnegative && x < 0) { ledger.status = S::outside_domain; return false; }
    return true;
  };
  if (v.endpoint) {
    auto &e = *v.endpoint;
    if (!inspect_text(e.producer_identity) || !inspect_text(e.coordinate_identity) || !inspect_text(e.unit_identity)) { fail(S::invalid_input, Refusal::source_strings); return out; }
    if (!number(e.support_lower, false) || !number(e.support_upper, false) || e.support_lower > e.support_upper) { fail(S::outside_domain, Refusal::source_number); return out; }
    if (!ledger.charge(&ProjectionWork::source_tag_inspections) || static_cast<unsigned>(e.survival_status) > static_cast<unsigned>(S::conditioning_budget_exceeded)) { fail(S::invalid_input, Refusal::source_tags); return out; }
    if (e.survival && (!number(*e.survival, true) || *e.survival > 1)) { fail(S::outside_domain, Refusal::source_number); return out; }
    if (e.survival_absolute_diagnostic && !number(*e.survival_absolute_diagnostic, true)) { fail(S::outside_domain, Refusal::source_number); return out; }
  }
  for (std::size_t i = 0; i < K; ++i) if (!number(v.k_mpc_inverse[i], true)) { fail(S::outside_domain, Refusal::source_number, i); return out; }
  for (std::size_t j = 0; j < N; ++j) if (!number(v.chi_mpc[j], true)) { fail(S::outside_domain, Refusal::source_number, -1, j); return out; }
  for (std::size_t t = 0; t < KS; ++t) if (!number(v.amplitudes[t], false) || std::abs(v.amplitudes[t]) > 1) { fail(S::outside_domain, Refusal::source_number, t / N, t % N); return out; }
  if (!supplied_shell_arithmetic_profile()) { fail(S::outside_domain, Refusal::arithmetic_profile); return out; }
  AllocationCapture capture{p.maximum_live_bytes - sizeof(detail::PreparedState)};
  try {
    detail::SourceStorage d;
    d.K = K; d.S = N; d.B = B; d.maximum_live = p.maximum_live_bytes;
    d.k = source_array<double>(K, capture); d.chi = source_array<double>(N, capture);
    d.amplitudes = source_array<double>(KS, capture);
    d.ids = source_array<detail::IdentifierOffset>(K + N, capture);
    d.roles = source_array<AtomRole>(N, capture); d.bytes = source_array<char>(B, capture);
    std::copy(v.k_mpc_inverse.begin(), v.k_mpc_inverse.end(), d.k.get());
    std::copy(v.chi_mpc.begin(), v.chi_mpc.end(), d.chi.get());
    std::copy(v.amplitudes.begin(), v.amplitudes.end(), d.amplitudes.get());
    std::copy(v.roles.begin(), v.roles.end(), d.roles.get());
    std::size_t offset = 0;
    auto store = [&](std::string_view s) { detail::IdentifierOffset id{offset, s.size()}; std::copy(s.begin(), s.end(), d.bytes.get() + offset); offset += s.size(); return id; };
    d.labels[0] = store(v.source_identity); d.labels[1] = store(v.mode_origin); d.labels[2] = store(v.ordering_provenance);
    for (std::size_t i = 0; i < K; ++i) d.ids[i] = store(v.k_ids[i]);
    for (std::size_t j = 0; j < N; ++j) d.ids[K + j] = store(v.shell_ids[j]);
    if (v.endpoint) { d.endpoint = v.endpoint; d.labels[3] = store(v.endpoint->producer_identity); d.labels[4] = store(v.endpoint->coordinate_identity); d.labels[5] = store(v.endpoint->unit_identity); }
    auto source = std::allocate_shared<detail::SourceMaker>(CaptureAllocator<detail::SourceMaker>(&capture), std::move(d));
    detail::SourceAccess::seal(*source, capture.arrays + capture.control);
    out.source_ = source;
    if (!detail::SourceAccess::reserve(*source, sizeof(detail::PreparedState))) throw PayloadLimit{};
    out.charge_ = sizeof(detail::PreparedState);
    out.state_ = std::make_unique<detail::PreparedState>(detail::PreparedState{p, out.fallback_work_});
    out.fallback_ = S::ok;
  } catch (const PayloadLimit &) {
    fail(S::work_limit, Refusal::allocation);
  } catch (const std::bad_alloc &) {
    fail(S::work_limit, Refusal::allocation);
  }
  if (!out.state_ && out.source_ && out.charge_) { detail::SourceAccess::release(*out.source_, out.charge_); out.charge_ = 0; }
  out.fallback_peak_ = capture.peak;
  return out;
}

namespace {
constexpr W unit_roundoff = 0x1p-64L;
bool normal_or_zero(W x) noexcept { return x == 0 || std::isnormal(x); }
struct Bounds {
  Ledger &ledger;
  S status = S::ok;
  bool begin() noexcept {
    if (status != S::ok) return false;
    if (!ledger.charge(&ProjectionWork::radius_operations)) {
      status = ledger.status; return false;
    }
    return true;
  }
  W finish(W r, bool upward, bool zero_witness) noexcept {
    if (r == 0 && zero_witness) return 0;
    if (!std::isnormal(r) || r < 0) { status = S::overflow; return 0; }
    W bound = std::nextafter(r, upward ? std::numeric_limits<W>::infinity() : -std::numeric_limits<W>::infinity());
    if (!std::isnormal(bound) || bound < 0) { status = S::overflow; return 0; }
    return bound;
  }
  W add(W a, W b, bool up = true) noexcept {
    if (!begin()) return 0;
    if (!normal_or_zero(a) || !normal_or_zero(b) || a < 0 || b < 0) { status = S::overflow; return 0; }
    return finish(a + b, up, a == 0 && b == 0);
  }
  W sub(W a, W b, bool up = true) noexcept {
    if (!begin()) return 0;
    if (!normal_or_zero(a) || !normal_or_zero(b) || a < b || b < 0) { status = S::overflow; return 0; }
    return finish(a - b, up, a == b);
  }
  W mul(W a, W b, bool up = true) noexcept {
    if (!begin()) return 0;
    if (!normal_or_zero(a) || !normal_or_zero(b) || a < 0 || b < 0) { status = S::overflow; return 0; }
    return finish(a * b, up, a == 0 || b == 0);
  }
  W div(W a, W b, bool up = true) noexcept {
    if (!begin()) return 0;
    if (!normal_or_zero(a) || !std::isnormal(b) || a < 0 || b <= 0) { status = S::overflow; return 0; }
    return finish(a / b, up, a == 0);
  }
  W gamma(unsigned d) noexcept {
    if (!d) return 0; // The operation count is an exact zero witness.
    W numerator = mul(static_cast<W>(d), unit_roundoff);
    W denominator = sub(1, numerator, false);
    if (status != S::ok || numerator >= 1 || denominator <= 0) { if (status == S::ok) status = S::outside_domain; return 0; }
    return div(numerator, denominator);
  }
  W eta(unsigned d) noexcept {
    W g = gamma(d);
    if (!d || status != S::ok) return 0;
    W denominator = sub(1, g, false);
    if (status != S::ok || g >= 1 || denominator <= 0) { if (status == S::ok) status = S::outside_domain; return 0; }
    return div(g, denominator);
  }
  double emit(W value) noexcept {
    if (status != S::ok) return 0;
    if (!ledger.charge(&ProjectionWork::output_casts)) { status = ledger.status; return 0; }
    if (!normal_or_zero(value) || value < 0) { status = S::overflow; return 0; }
    double stored = static_cast<double>(value);
    if (value == 0) return 0;
    if (!std::isnormal(stored) || stored <= 0) { status = S::overflow; return 0; }
    if (static_cast<W>(stored) < value) stored = std::nextafter(stored, std::numeric_limits<double>::infinity());
    if (!std::isnormal(stored)) { status = S::overflow; return 0; }
    return stored;
  }
};
struct Kernel {
  S status = S::invalid_input;
  W value = 0, tail = 0, arithmetic = 0;
};
Kernel bessel(W x, unsigned ell, Ledger &ledger) noexcept {
  Kernel result;
  if (!ledger.charge(&ProjectionWork::Bessel_evaluations)) { result.status = ledger.status; return result; }
  if (x == 0) { result.status = S::ok; result.value = ell == 0 ? 1 : 0; return result; }
  Bounds b{ledger};
  if (!ledger.charge(&ProjectionWork::series_terms)) { result.status = ledger.status; return result; }
  W term = 1;
  for (unsigned r = 1; r <= ell; ++r) {
    W factor = x / static_cast<W>(2 * r + 1);
    if (!std::isnormal(factor)) { result.status = S::overflow; return result; }
    term = term * factor;
    if (!std::isnormal(term)) { result.status = S::overflow; return result; }
  }
  W square = x * x;
  if (!std::isnormal(square)) { result.status = S::overflow; return result; }
  W sum = 0, absolute_sum = 0, term_error_sum = 0;
  W term_error = b.mul(b.eta(2 * ell), std::abs(term));
  const W tail_goal = b.div(1, 100000000000000.L, false);
  unsigned n = 0, additions = 0;
  for (;;) {
    if (b.status != S::ok) { result.status = b.status; return result; }
    if (!ledger.charge(&ProjectionWork::signed_additions)) { result.status = ledger.status; return result; }
    sum = sum + term;
    if (!normal_or_zero(sum)) { result.status = S::overflow; return result; }
    ++additions;
    absolute_sum = b.add(absolute_sum, std::abs(term));
    term_error_sum = b.add(term_error_sum, term_error);
    if (!ledger.charge(&ProjectionWork::series_terms)) { result.status = ledger.status; return result; }
    const unsigned next = n + 1;
    const unsigned denominator = 2 * next * (2 * ell + 2 * n + 3);
    const W product = term * square;
    if (!std::isnormal(product)) { result.status = S::overflow; return result; }
    const W next_term = -product / static_cast<W>(denominator);
    if (!std::isnormal(next_term)) { result.status = S::overflow; return result; }
    const W next_error = b.mul(b.eta(2 * ell + 3 * next), std::abs(next_term));
    if (n >= 3) {
      const W omitted = b.add(std::abs(next_term), next_error);
      if (b.status != S::ok) { result.status = b.status; return result; }
      if (omitted <= tail_goal) {
        result.value = sum; result.tail = omitted;
        result.arithmetic = b.add(term_error_sum, b.mul(b.gamma(additions), absolute_sum));
        result.status = b.status; return result;
      }
    }
    if (n == 32) { result.status = S::conditioning_budget_exceeded; return result; }
    term = next_term; term_error = next_error; n = next;
  }
}
bool exact_power_two(double value) noexcept {
  const auto bits = std::bit_cast<std::uint64_t>(value);
  const auto exponent = (bits >> 52) & 2047;
  return !(bits >> 63) && exponent && exponent != 2047 && !(bits & ((std::uint64_t{1} << 52) - 1));
}
struct Accumulator {
  S status = S::ok;
  W sum = 0, absolute_products = 0;
  W argument = 0, tail = 0, arithmetic = 0, products = 0;
  unsigned additions = 0;
};
struct EvaluationScratch { std::array<Accumulator, 9> accumulators; };
void finish_row(ProjectionRow &row, Accumulator &a, unsigned q, Ledger &ledger) noexcept {
  row.status = a.status; row.check = ProjectionCheck::refused;
  if (a.status != S::ok) return;
  Bounds b{ledger};
  if (!ledger.charge(&ProjectionWork::output_casts)) { row.status = ledger.status; return; }
  const double value = static_cast<double>(a.sum);
  if (!std::isfinite(value) || (a.sum != 0 && !std::isnormal(value))) { row.status = S::overflow; return; }
  const W readback = static_cast<W>(value);
  W cast_loss = 0;
  if (a.sum != readback) {
    const W high = std::max(std::abs(a.sum), std::abs(readback));
    const W low = std::min(std::abs(a.sum), std::abs(readback));
    if (low && high <= 2 * low) cast_loss = high - low; // Sterbenz, same sign.
    else cast_loss = b.sub(high, low);
    if (!normal_or_zero(cast_loss)) { row.status = S::overflow; return; }
  }
  const W signed_error = b.mul(b.gamma(a.additions), a.absolute_products);
  const std::array<W, 6> components{a.argument, a.tail, a.arithmetic, a.products, signed_error, cast_loss};
  std::array<double, 6> stored{};
  W upper = 0, lower = 0;
  for (std::size_t j = 0; j < components.size(); ++j) {
    stored[j] = b.emit(components[j]);
    upper = b.add(upper, static_cast<W>(stored[j]));
    lower = b.add(lower, static_cast<W>(stored[j]), false);
  }
  const double radius = b.emit(upper);
  const W radius_readback = static_cast<W>(radius);
  const double padding = b.emit(b.sub(radius_readback, lower));
  const W scale = static_cast<W>(std::uint64_t{1} << q);
  W allocation = b.div(1, 10000000000.L, false);
  allocation = b.mul(allocation, b.add(1, std::abs(readback), false), false);
  allocation = b.div(allocation, scale, false);
  if (b.status != S::ok) { row.status = b.status; return; }
  if (radius_readback > allocation) { row.status = S::conditioning_budget_exceeded; return; }
  row.signed_value = value; row.absolute_numerical_radius = radius;
  row.diagnostics = ProjectionDiagnostics{stored[0], stored[1], stored[2], stored[3], stored[4], stored[5], padding};
  row.status = S::ok; row.check = ProjectionCheck::numerical_radius_within_allocation;
}
} // namespace

ProjectionBatch PreparedProjection::evaluate(ProjectionRequest request) {
  ProjectionBatch out;
  out.source_ = source_;
  out.before_ = cumulative_work(); out.after_ = out.before_;
  auto finish = [&] {
    out.after_ = cumulative_work(); out.delta_ = difference(out.after_, out.before_);
    out.peak_ = source_ ? source_->peak_payload_bytes() : fallback_peak_;
  };
  if (fallback_ != S::ok || !state_ || !source_) {
    out.fallback_ = fallback_; out.failure_ = failure_; finish(); return out;
  }
  auto &work = state_->work;
  Ledger ledger{work, state_->resources};
  auto fail = [&](S s, Refusal r) { out.fallback_ = ledger.status == S::ok ? s : ledger.status; out.failure_.reason = r; finish(); };
  if (!ledger.charge(&ProjectionWork::source_tag_inspections)) { fail(S::work_limit, Refusal::request); return out; }
  if (request.outputs != supplied_shell_signed_projection || request.ell.empty() || request.ell.size() > 9) { fail(S::invalid_input, Refusal::request); return out; }
  if (!ledger.charge(&ProjectionWork::source_tag_inspections)) { fail(S::work_limit, Refusal::request); return out; }
  if (request.accuracy_reduction_power > 32) { fail(S::outside_domain, Refusal::request); return out; }
  for (auto ell : request.ell) {
    if (!ledger.charge(&ProjectionWork::source_tag_inspections)) { fail(S::work_limit, Refusal::request); return out; }
    if (ell > 8) { fail(S::outside_domain, Refusal::request); return out; }
  }
  if (!supplied_shell_arithmetic_profile()) { fail(S::outside_domain, Refusal::arithmetic_profile); return out; }
  const auto K = source_->k_mpc_inverse().size(), N = source_->chi_mpc().size(), L = request.ell.size();
  std::size_t count = 0, row_bytes = 0, result_bytes = sizeof(detail::ResultState), reservation = 0;
  if (!mul_size(K, L, count) || !mul_size(count, sizeof(ProjectionRow), row_bytes) ||
      !add_size(result_bytes, row_bytes)) { fail(S::overflow, Refusal::allocation); return out; }
  reservation = result_bytes;
  if (!add_size(reservation, sizeof(EvaluationScratch)) || !detail::SourceAccess::reserve(*source_, reservation)) { fail(S::work_limit, Refusal::allocation); return out; }
  out.charge_ = result_bytes;
  std::unique_ptr<EvaluationScratch> scratch;
  bool scratch_reserved = true;
  try {
    out.state_ = std::make_unique<detail::ResultState>();
    out.state_->ell_count = L; out.state_->row_count = count;
    std::copy(request.ell.begin(), request.ell.end(), out.state_->ell.begin());
    out.state_->rows = allocate_array<ProjectionRow>(count);
    scratch = std::make_unique<EvaluationScratch>();
    for (std::size_t i = 0; i < K; ++i) for (std::size_t j = 0; j < L; ++j) {
      auto &row = out.state_->rows[i * L + j];
      row.k_index = i; row.ell_index = j; row.ell = request.ell[j];
      row.check = ProjectionCheck::refused;
    }
    out.fallback_ = S::ok;
    for (std::size_t i = 0; i < K; ++i) {
      scratch->accumulators = {};
      S phase_status = S::ok;
      for (std::size_t sh = 0; sh < N; ++sh) {
        if (!ledger.charge(&ProjectionWork::shell_visits) || !ledger.charge(&ProjectionWork::phase_products)) { phase_status = ledger.status; break; }
        const double k = source_->k_mpc_inverse()[i], chi = source_->chi_mpc()[sh];
        const W x = static_cast<W>(k) * static_cast<W>(chi);
        W argument_radius = 0;
        Bounds b{ledger};
        if (x != 0 && !std::isnormal(x)) phase_status = S::overflow;
        else if (k == 0 || chi == 0) { /* Exact source zero witness. */ }
        else if (exact_power_two(k) || exact_power_two(chi)) {
          if (x < 0x1p-32L || x > 8) phase_status = S::outside_domain;
        } else {
          argument_radius = b.mul(b.div(unit_roundoff, b.sub(1, unit_roundoff, false)), x);
          const W lo = b.sub(x, argument_radius, false), hi = b.add(x, argument_radius);
          if (b.status != S::ok) phase_status = b.status;
          else if (lo < 0x1p-32L || hi > 8) phase_status = S::outside_domain;
        }
        if (phase_status != S::ok) { if (out.failure_.reason == Refusal::none) out.failure_ = {Refusal::phase, i, sh, static_cast<std::size_t>(-1)}; continue; }
        const W amplitude = static_cast<W>(source_->amplitudes()[i * N + sh]);
        if (amplitude == 0) continue;
        for (std::size_t j = 0; j < L; ++j) {
          auto &a = scratch->accumulators[j];
          if (a.status != S::ok) continue;
          const Kernel kernel = bessel(x, request.ell[j], ledger);
          if (kernel.status != S::ok) { a.status = kernel.status; continue; }
          Bounds coefficient{ledger};
          const W weight = std::abs(amplitude);
          a.argument = coefficient.add(a.argument, coefficient.mul(weight, coefficient.div(argument_radius, 2)));
          a.tail = coefficient.add(a.tail, coefficient.mul(weight, kernel.tail));
          a.arithmetic = coefficient.add(a.arithmetic, coefficient.mul(weight, kernel.arithmetic));
          if (!ledger.charge(&ProjectionWork::shell_products)) { a.status = ledger.status; continue; }
          const W product = amplitude * kernel.value;
          if ((kernel.value != 0 && !std::isnormal(product)) || !normal_or_zero(product)) { a.status = S::overflow; continue; }
          a.products = coefficient.add(a.products, coefficient.mul(coefficient.eta(1), std::abs(product)));
          a.absolute_products = coefficient.add(a.absolute_products, std::abs(product));
          if (!ledger.charge(&ProjectionWork::signed_additions)) { a.status = ledger.status; continue; }
          a.sum = a.sum + product; ++a.additions;
          if (!normal_or_zero(a.sum)) a.status = S::overflow;
          else a.status = coefficient.status;
        }
      }
      for (std::size_t j = 0; j < L; ++j) {
        auto &row = out.state_->rows[i * L + j];
        auto &a = scratch->accumulators[j];
        if (phase_status != S::ok) a.status = phase_status;
        finish_row(row, a, request.accuracy_reduction_power, ledger);
        if (row.status != S::ok && out.fallback_ == S::ok) {
          out.fallback_ = row.status;
          if (out.failure_.reason == Refusal::none) out.failure_ = {row.status == S::conditioning_budget_exceeded ? Refusal::numerical_allocation : Refusal::nonnormal_arithmetic, i, static_cast<std::size_t>(-1), j};
        }
      }
    }
  } catch (const std::bad_alloc &) {
    out.state_.reset();
    detail::SourceAccess::release(*source_, out.charge_); out.charge_ = 0;
    fail(S::work_limit, Refusal::allocation);
  }
  scratch.reset();
  if (scratch_reserved) detail::SourceAccess::release(*source_, sizeof(EvaluationScratch));
  finish();
  return out;
}
} // namespace irred::projection
