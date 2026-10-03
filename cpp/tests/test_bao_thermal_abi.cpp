// Exact original3-row input; standalone/ABI parity shares native ancestry.
// GL8 direct-z/sqrt(a) polynomial control reuses the existing peer's original
// source algorithm/constants, with its own512/256 refinement gate (not a new
// independent external reference). No runtime/reference receipt earned yet.
#include "irred/abi.h"
#include "irred/bao_thermal.hpp"
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <numbers>
#include <stdexcept>
#include <string_view>
using namespace irred;
// Track requested heap bytes, excluding this header/allocator overhead/RSS.
// Atomic observable state and noinline operators survive Release optimization.
namespace allocation {
struct alignas(std::max_align_t) Header {
  size_t bytes;
};
std::atomic<size_t> live = 0, calls = 0, peak = 0, baseline = 0, fail_on = 0;
std::atomic<bool> armed = false;
void start(size_t fail = 0) {
  baseline = live.load();
  calls = 0;
  peak = 0;
  fail_on = fail;
  armed = true;
}
} // namespace allocation
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(size_t n) {
  if (allocation::armed && ++allocation::calls == allocation::fail_on)
    throw std::bad_alloc();
  if (n > SIZE_MAX - sizeof(allocation::Header))
    throw std::bad_alloc();
  auto *h = static_cast<allocation::Header *>(
      std::malloc(sizeof(allocation::Header) + (n ? n : 1)));
  if (!h)
    throw std::bad_alloc();
  h->bytes = n;
  const auto current = allocation::live.fetch_add(n) + n;
  if (allocation::armed && current >= allocation::baseline &&
      current - allocation::baseline > allocation::peak)
    allocation::peak = current - allocation::baseline;
  return h + 1;
}
NOINLINE void *operator new[](size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept {
  if (p) {
    auto *h = static_cast<allocation::Header *>(p) - 1;
    allocation::live -= h->bytes;
    std::free(h);
  }
}
NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
NOINLINE void operator delete(void *p, size_t) noexcept {
  ::operator delete(p);
}
NOINLINE void operator delete[](void *p, size_t) noexcept {
  ::operator delete(p);
}
namespace {
using W = long double;
void need(bool x, const char *message) {
  if (!x)
    throw std::runtime_error(message);
}
void near(W a, W b, W allowance, const char *message) {
  need(std::abs(a - b) <= allowance, message);
}
irred_bytes text(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
template <size_t N> irred_f64_buffer values(const std::array<double, N> &v) {
  return {sizeof(irred_f64_buffer),
          IRRED_ABI_VERSION,
          2,
          0,
          v.data(),
          N,
          sizeof(v)};
}
template <size_t N> irred_strings ids(const std::array<irred_bytes, N> &v) {
  return {v.data(), N, sizeof(v)};
}
struct Result {
  irred_bao_thermal_result *p = nullptr;
  ~Result() { irred_bao_thermal_result_destroy(p); }
  Result() = default;
  Result(const Result &) = delete;
  Result &operator=(const Result &) = delete;
  irred_bao_thermal_view view() {
    irred_bao_thermal_view v{};
    need(irred_bao_thermal_result_view(p, &v) == IRRED_OK, "view");
    return v;
  }
};
struct Fixture {
  std::array<double, 3> y{24, 17, 22};
  std::array<double, 9> C{64, 16, 0, 16, 80, 16, 0, 16, 96};
  std::array<irred_bao_query, 3> queries{{{1, 0, 0}, {1, 1, 0}, {1, 2, 0}}};
  std::array<irred_bytes, 3> row{text("synthetic-DM-z1"),
                                 text("synthetic-DH-z1"),
                                 text("synthetic-DV-z1")};
  std::array<irred_bao_thermal_species, 2> sp{};
  std::array<irred_bao_thermal_model, 2> models{};
  irred_bao_thermal_source source{};
  irred_bao_thermal_batch batch{};
  irred_bao_thermal_policy policy{};
  Fixture() {
    for (auto &x : sp)
      x = {sizeof(x), IRRED_ABI_VERSION, 0, 0, 0, 1.95, 2};
    for (size_t i = 0; i < 2; ++i) {
      auto &m = models[i];
      m.struct_size = sizeof(m);
      m.abi_version = IRRED_ABI_VERSION;
      m.id = text(i ? "synthetic-H0-72-fixed-physical-inputs"
                    : "synthetic-H0-70-fixed-physical-inputs");
      m.h0_km_s_mpc = i ? 72 : 70;
      m.physical_baryon_density = .0245;
      m.physical_cdm_density = .1225;
      m.tcmb_kelvin = 2.7255;
      m.physical_massless_nonphoton_density = 1e-5;
      m.species = &sp[i];
      m.species_count = 1;
      m.species_byte_length = sizeof(sp[i]);
      m.z_drag = 1059.95;
      m.drag_origin = text(
          "explicit synthetic supplied drag; no physical endpoint prediction");
      m.source_origin = text("same explicit omega_b/omega_cdm/other density, "
                             "Tcmb and mass0 FD Kelvin state; H0 varies");
    }
    source.struct_size = sizeof(source);
    source.abi_version = IRRED_ABI_VERSION;
    source.role = 1;
    source.covariance_unit = 0;
    source.queries = queries.data();
    source.query_count = 3;
    source.query_byte_length = sizeof(queries);
    source.observed = values(y);
    source.covariance = values(C);
    source.ordered_ids = source.covariance_axis_ids = ids(row);
    source.table_identity = text("synthetic-three-ratio-fixed-observation/v1");
    source.covariance_identity = text("16-times-integer-SPD3/v1");
    source.ordering_provenance = text(
        "DM,DH,DV at exact source z=1; full covariance in that same order");
    source.calibration_provenance = text(
        "explicit synthetic ratios; no measured calibration qualification");
    source.dependence_provenance =
        text("one full supplied three-row covariance; no cross-probe or "
             "cross-model independence asserted");
    source.redshift_convention = text(
        "synthetic shared FLRW model redshift; no observed frame conversion");
    batch = {sizeof(batch), IRRED_ABI_VERSION, models.data(), 2,
             sizeof(models)};
    auto &p = policy;
    p.struct_size = sizeof(p);
    p.abi_version = IRRED_ABI_VERSION;
    p.arithmetic = 1;
    p.requested = 7;
    p.maximum_rows = 64;
    p.maximum_matrix_elements = 4096;
    p.maximum_models = 16;
    p.maximum_species_per_model = 16;
    p.maximum_string_bytes = 65536;
    p.maximum_output_array_elements = 2048;
    p.maximum_native_bytes = 16 << 20;
    p.maximum_preparation_native_bytes = p.maximum_evaluation_native_bytes =
        8 << 20;
    p.maximum_total_callbacks = 200000000;
    p.maximum_forward_sensitivity = p.maximum_projection_log_density_error =
        1e-8;
    p.predictions = {sizeof(p.predictions),
                     IRRED_ABI_VERSION,
                     0,
                     30,
                     30,
                     0,
                     1e-8,
                     2e-10,
                     1e-10,
                     5e-10,
                     1e-12,
                     2e-12,
                     20000000,
                     100000000,
                     16 << 20,
                     200000,
                     100000000,
                     16 << 20};
  }
  cosmology::ThermalObservableRequest model(size_t i) const {
    auto &m = models[i];
    return {{m.h0_km_s_mpc,
             m.physical_baryon_density,
             m.physical_cdm_density,
             m.tcmb_kelvin,
             m.physical_massless_nonphoton_density,
             {{sp[i].mass_ev, sp[i].temperature_today_kelvin,
               sp[i].statistical_weight}}},
            m.z_drag,
            std::string(reinterpret_cast<const char *>(m.drag_origin.data),
                        m.drag_origin.length),
            std::string(reinterpret_cast<const char *>(m.source_origin.data),
                        m.source_origin.length)};
  }
  bao::PreparedDensity observation() const {
    bao::DensityInput a;
    a.queries = {{1, bao::Observable::transverse_over_ruler},
                 {1, bao::Observable::hubble_over_ruler},
                 {1, bao::Observable::volume_over_ruler}};
    a.observed.assign(y.begin(), y.end());
    a.covariance.assign(C.begin(), C.end());
    for (auto &x : row)
      a.ordered_ids.emplace_back(reinterpret_cast<const char *>(x.data),
                                 x.length);
    a.role = bao::RowRole::synthetic_control;
    a.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
    auto str = [](irred_bytes v) {
      return std::string(reinterpret_cast<const char *>(v.data), v.length);
    };
    a.table_identity = str(source.table_identity);
    a.covariance_identity = str(source.covariance_identity);
    a.ordering_provenance = str(source.ordering_provenance);
    a.calibration_provenance = str(source.calibration_provenance);
    a.dependence_provenance = str(source.dependence_provenance);
    return bao::prepare_density(std::move(a),
                                {64, 4096, 65536, 8 << 20, 1e-8,
                                 numerics::Arithmetic::longdouble_cpu_v1});
  }
};
template <class F> W gl(F f, W end, unsigned n) {
  constexpr W x[]{
      .183434642495649804939476142360184L, .525532409916328985817739049189246L,
      .796666477413626739591553936475831L, .960289856497536231683560868569473L};
  constexpr W w[]{
      .362683783378361982965150449277195L, .313706645877887287337962201986601L,
      .222381034453374470544355994426240L, .101228536290376259152531354309962L};
  W sum = 0, h = end / (2 * n);
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < 4; ++j) {
      W mid = (2 * i + 1) * h;
      sum += h * w[j] * (f(mid - h * x[j]) + f(mid + h * x[j]));
    }
  return sum;
}
W critical(W H) {
  const W pi = std::numbers::pi_v<W>, c = 299792458.L, ev = 1.602176634e-19L,
          hbar = 6.62607015e-34L / (2 * pi),
          mpc = 648000000000.L * 149597870700.L / pi;
  const W hs = H * 1000 / mpc;
  return 3 * hs * hs * c * c / (8 * pi * 6.67430e-11L) / ev *
         std::pow(hbar * c / ev, 3);
}
std::array<W, 3> polynomial(cosmology::ThermalObservableRequest q, W z,
                            unsigned n) {
  auto m = q.model;
  const W pi = std::numbers::pi_v<W>, kb = 1.380649e-23L / 1.602176634e-19L,
          h2 = std::pow(W(m.h0_km_s_mpc) / 100, 2),
          b = m.physical_baryon_density / h2,
          g = pi * pi / 15 * std::pow(kb * m.tcmb_kelvin, 4) /
              critical(m.h0_km_s_mpc),
          matter = b + m.physical_cdm_density / h2;
  W radiation = g + m.physical_massless_nonphoton_density / h2;
  for (auto sp : m.species) {
    need(sp.mass_ev == 0, "polynomial only massless FD limit");
    radiation += sp.statistical_weight * 7 * pi * pi / 240 *
                 std::pow(kb * sp.temperature_today_kelvin, 4) /
                 critical(m.h0_km_s_mpc);
  }
  W lambda = 1 - radiation - matter, scale = 299792.458L / m.h0_km_s_mpc;
  auto P = [&](W a) { return radiation + matter * a + lambda * a * a * a * a; };
  auto E = [&](W z0) {
    W a = 1 / (1 + z0);
    return std::sqrt(P(a)) / (a * a);
  };
  W rs = scale / std::sqrt(3.L) *
         gl(
             [&](W u) {
               W a = u * u;
               return 2 * u / std::sqrt(P(a)) /
                      std::sqrt(1 + 3 * b * a / (4 * g));
             },
             1 / std::sqrt(1 + W(q.z_drag)), n),
    dm = scale * gl([&](W t) { return 1 / E(t); }, z, n), dh = scale / E(z);
  return {dm / rs, dh / rs, std::cbrt(dm * dm * z * dh) / rs};
}

void numeric_controls(Fixture &f, Result &r) {
  need(irred_bao_thermal_evaluate(&f.source, &f.batch, &f.policy, &r.p) ==
               IRRED_OK &&
           r.p,
       "original admitted ABI");
  auto v = r.view();
  need(v.source_factor_completed && v.thermal_batch_call_attempted &&
           v.model_count == 2,
       "native phases");
  auto obs = f.observation();
  need(obs.status() == statistics::DensityStatus::finite, "standalone factor");
  std::array<cosmology::ThermalObservableRequest, 2> requests{f.model(0),
                                                              f.model(1)};
  bao::ThermalDensityPolicy p;
  p.maximum_models = 16;
  p.maximum_queries = 64;
  p.maximum_string_bytes = 65536;
  p.maximum_native_bytes = 8 << 20;
  p.maximum_total_callbacks = 200000000;
  p.maximum_forward_sensitivity = 1e-8;
  p.maximum_projection_log_density_error = 1e-8;
  p.requested = 7;
  auto native = obs.evaluate_thermal(requests, p);
  need(native.slots.size() == 2, "standalone batch");
  constexpr W inverse[]{29, -6, 1, -6, 24, -4, 1, -4, 19};
  for (size_t i = 0; i < 2; ++i) {
    auto &x = v.rows[i];
    auto &s = native.slots[i];
    need(x.density.state.availability == 1 &&
             x.prediction_state.availability == 1 &&
             x.residual_state.availability == 1 && s.result,
         "original requested groups");
    auto fine = polynomial(requests[i], 1, 512),
         coarse = polynomial(requests[i], 1, 256);
    W residual[3], quad = 0, oldq = 0;
    for (size_t j = 0; j < 3; ++j) {
      const W allowance = 1e-10L + 5e-10L * std::abs(fine[j]);
      near(fine[j], coarse[j], .05L * allowance, "reference ratio refinement");
      near(x.predictions.data[j], fine[j], allowance, "direct-z sqrt(a) ratio");
      need(x.predictions.data[j] == s.predictions[j] &&
               x.residuals.data[j] == s.residuals[j],
           "standalone ABI bit parity");
      residual[j] = W(f.y[j]) - fine[j];
    }
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k) {
        quad += residual[j] * inverse[j * 3 + k] * residual[k] / 1760;
        oldq += (W(f.y[j]) - coarse[j]) * inverse[j * 3 + k] *
                (W(f.y[k]) - coarse[k]) / 1760;
      }
    near(quad, oldq, .05L * 1e-8L, "reference cofactor refinement");
    const W ld = std::log(450560.L),
            norm = 3 * std::log(2 * std::numbers::pi_v<W>);
    near(x.density.quadratic, quad, 1e-8L, "named cofactor quadratic");
    near(x.density.log_determinant, ld, 1e-8L, "named determinant");
    near(x.density.log_normalization, norm, 1e-8L, "named normalization");
    near(x.density.log_density, -(quad + ld + norm) / 2, 1e-8L,
         "named log density");
    need(x.density.log_density == s.result->density.log_value,
         "standalone ABI density bits");
    need(x.callbacks == x.outer_callbacks + x.momentum_callbacks &&
             x.preparation_callbacks <= x.momentum_callbacks,
         "work no double count");
  }
  need(v.callbacks == v.outer_callbacks + v.momentum_callbacks &&
           v.callbacks <= 200000000,
       "original whole work cap");
  need(std::abs(v.rows[0].predictions.data[0] - v.rows[1].predictions.data[0]) >
           1e-8,
       "fixed physical H0 has no cancellation");
}
} // namespace
int main() {
  try {
    Fixture f;
    std::string caller_id = "synthetic-H0-70-fixed-physical-inputs",
                caller_redshift = "synthetic shared FLRW model redshift; no "
                                  "observed frame conversion";
    f.models[0].id = text(caller_id);
    f.source.redshift_convention = text(caller_redshift);
    Result original;
    numeric_controls(f, original);
    const auto baseline = allocation::live.load();
    {
      Result r;
      allocation::start();
      auto q = f.policy;
      q.maximum_preparation_native_bytes = 0;
      auto status = irred_bao_thermal_evaluate(&f.source, &f.batch, &q, &r.p);
      allocation::armed = false;
      need(status == IRRED_OK && r.p, "source subcap owned refusal");
      auto v = r.view();
      need(v.phase == 0 && !v.source_available &&
               !v.source_prepare_call_attempted && v.numerical_status == 4 &&
               v.model_count == 0,
           "source preflight withheld");
      const size_t minimum = v.retained_payload_bytes;
      need(allocation::calls == 1 && allocation::peak == minimum,
           "allocation-free minimum constructor");
      for (size_t cap : {size_t(0), minimum - 1, minimum}) {
        Result small;
        auto z = f.policy;
        z.maximum_native_bytes = cap;
        z.maximum_preparation_native_bytes = z.maximum_evaluation_native_bytes =
            0;
        allocation::start();
        auto result =
            irred_bao_thermal_evaluate(&f.source, &f.batch, &z, &small.p);
        allocation::armed = false;
        need(result ==
                 (cap < minimum ? IRRED_BAO_THERMAL_QUOTA_REFUSED : IRRED_OK),
             "minimum exact boundary namespace");
        need(bool(small.p) == (cap >= minimum) &&
                 allocation::calls == (cap < minimum ? 0 : 1),
             "no hidden quota allocation");
      }
    }
    need(allocation::live == baseline, "minimum lifetime released");
    {
      Result r;
      allocation::start();
      auto status =
          irred_bao_thermal_evaluate(&f.source, &f.batch, &f.policy, &r.p);
      allocation::armed = false;
      need(status == IRRED_OK && r.p, "measured original");
      auto v = r.view();
      const auto sites = allocation::calls.load();
      need(allocation::peak <= v.required_global_peak_bytes &&
               allocation::live - baseline <= v.retained_payload_bytes,
           "whole native peak and retained charge");
      std::cout << "allocation sites=" << sites << " peak=" << allocation::peak
                << " retained=" << allocation::live - baseline
                << " declared_peak=" << v.required_global_peak_bytes
                << " declared_retained=" << v.retained_payload_bytes << '\n';
      for (size_t site = 1; site <= sites; ++site) {
        Result failed;
        allocation::start(site);
        auto code = irred_bao_thermal_evaluate(&f.source, &f.batch, &f.policy,
                                               &failed.p);
        allocation::armed = false;
        need(code == IRRED_ALLOCATION_FAILURE && !failed.p,
             "every allocation fault releases owner");
      }
    }
    need(allocation::live == baseline, "whole call lifetime released");
    {
      auto invalid = f;
      invalid.C[0] = -64;
      invalid.source.covariance = values(invalid.C);
      Result r;
      need(irred_bao_thermal_evaluate(&invalid.source, &invalid.batch,
                                      &invalid.policy, &r.p) == IRRED_OK,
           "invalid covariance owned");
      auto v = r.view();
      need(v.source_prepare_call_attempted && !v.source_factor_completed &&
               v.source_available && !v.actual_source_arithmetic_id.length &&
               !v.thermal_batch_call_attempted && v.model_count == 0,
           "failed factor raw source only");
    }
    {
      auto q = f.policy;
      q.maximum_evaluation_native_bytes = 0;
      Result r;
      need(irred_bao_thermal_evaluate(&f.source, &f.batch, &q, &r.p) ==
               IRRED_OK,
           "eval subcap owned");
      auto v = r.view();
      need(v.source_factor_completed && v.source_available &&
               !v.thermal_batch_call_attempted && v.model_count == 0 &&
               v.numerical_status == 4,
           "retained source eval refusal");
    }
    for (unsigned mask : {1u, 2u, 4u}) {
      auto q = f.policy;
      q.requested = mask;
      Result r;
      need(irred_bao_thermal_evaluate(&f.source, &f.batch, &q, &r.p) ==
               IRRED_OK,
           "output mask");
      auto v = r.view();
      for (size_t i = 0; i < v.model_count; ++i) {
        auto &x = v.rows[i];
        need((x.predictions.length != 0) == bool(mask & 2) &&
                 (x.residuals.length != 0) == bool(mask & 4),
             "output omission arrays");
      }
    }
    {
      auto q = f.policy;
      q.maximum_projection_log_density_error = 1e-30;
      Result r;
      need(irred_bao_thermal_evaluate(&f.source, &f.batch, &q, &r.p) ==
               IRRED_OK,
           "projection refusal");
      auto v = r.view();
      need(v.model_count == 2 &&
               v.rows[0].density.state.availability ==
                   static_cast<uint32_t>(cosmology::Availability::failed) &&
               v.rows[0].prediction_state.availability == 1 &&
               v.rows[0].residual_state.availability == 1,
           "projection preserves successful groups");
      need(v.rows[0].density.quadratic == 0 &&
               v.rows[0].density.log_density == 0 &&
               v.rows[0].density.has_projection_estimate,
           "failed density suppression retains diagnostic");
    }
    {
      auto q = f.policy;
      q.maximum_total_callbacks = 0;
      Result r;
      need(irred_bao_thermal_evaluate(&f.source, &f.batch, &q, &r.p) ==
               IRRED_OK,
           "zero work owned");
      auto v = r.view();
      need(v.callbacks == 0 && v.model_count == 2 &&
               v.rows[0].prediction_state.availability ==
                   static_cast<uint32_t>(cosmology::Availability::failed),
           "zero cap no work or dropped models");
    }
    // Hostile descriptors refuse before allocation/dereference, never
    // manufacture owner.
    auto hostile = [&](const irred_bao_thermal_source *a,
                       const irred_bao_thermal_batch *b,
                       const irred_bao_thermal_policy *q, uint32_t expected) {
      irred_bao_thermal_result *out = nullptr;
      allocation::start();
      auto status = irred_bao_thermal_evaluate(a, b, q, &out);
      allocation::armed = false;
      need(status == expected && !out && !allocation::calls,
           "hostile no allocation/owner");
    };
    hostile(reinterpret_cast<const irred_bao_thermal_source *>(uintptr_t(1)),
            &f.batch, &f.policy, IRRED_INVALID_INPUT);
    {
      auto x = f.source;
      x.abi_version = 1;
      hostile(&x, &f.batch, &f.policy, IRRED_ABI_MISMATCH);
    }
    {
      auto x = f.source;
      x.query_count = UINT64_MAX;
      hostile(&x, &f.batch, &f.policy, IRRED_INVALID_INPUT);
    }
    {
      auto x = f.source;
      x.covariance.byte_length--;
      hostile(&x, &f.batch, &f.policy, IRRED_INVALID_INPUT);
    }
    {
      auto x = f.source;
      x.observed.reserved = 1;
      hostile(&x, &f.batch, &f.policy, IRRED_INVALID_INPUT);
    }
    {
      auto x = f.batch;
      x.model_byte_length--;
      hostile(&f.source, &x, &f.policy, IRRED_INVALID_INPUT);
    }
    {
      auto x = f.policy;
      x.predictions.reserved = 1;
      hostile(&f.source, &f.batch, &x, IRRED_INVALID_INPUT);
    }
    {
      auto x = f.policy;
      x.predictions.momentum_method = 1;
      hostile(&f.source, &f.batch, &x, IRRED_INVALID_INPUT);
    }
    // Original caller arrays/text may disappear without invalidating retained
    // views.
    auto lifetime = original.view();
    std::string().swap(caller_id);
    std::string().swap(caller_redshift);
    need(std::string_view(
             reinterpret_cast<const char *>(lifetime.rows[0].source.id.data),
             lifetime.rows[0].source.id.length) ==
             "synthetic-H0-70-fixed-physical-inputs",
         "model ID owned after caller text destruction");
    need(std::string_view(reinterpret_cast<const char *>(
                              lifetime.source.redshift_convention.data),
                          lifetime.source.redshift_convention.length) ==
             "synthetic shared FLRW model redshift; no observed frame "
             "conversion",
         "redshift convention owned after caller text destruction");
    f.y.fill(-999);
    f.C.fill(-999);
    f.queries.fill({5, 2, 0});
    need(lifetime.source.observed.data[0] == 24 &&
             lifetime.source.covariance.data[0] == 64 &&
             lifetime.source.queries[0].z == 1,
         "owned source independent input lifetime");
    std::cout << "thermal BAO ABI original controls passed\n";
    return 0;
  } catch (const std::exception &e) {
    allocation::armed = false;
    std::cerr << e.what() << '\n';
    return 1;
  }
}
