// Independent boundary review. Native parity is transport evidence; independent
// piecewise equations/refinements are in test_piecewise_background_hostile.
#include "irred/abi.h"
#include "irred/piecewise_background.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
static long long fail_after = -1, live = 0;
static unsigned checks = 0;
void *operator new(std::size_t n) {
  if (fail_after == 0)
    throw std::bad_alloc();
  if (fail_after > 0)
    --fail_after;
  void *p = std::malloc(n ? n : 1);
  if (!p)
    throw std::bad_alloc();
  ++live;
  return p;
}
void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
void operator delete(void *p, std::size_t) noexcept { ::operator delete(p); }
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }
#define CHECK(x)                                                               \
  do {                                                                         \
    ++checks;                                                                  \
    if (!(x))                                                                  \
      throw std::runtime_error("piecewise ABI: " #x);                          \
  } while (false)
using namespace irred::cosmology;
struct Fixture {
  std::array<cosmo_piecewise_parameters, 3> parameters{
      {{70, {-0.4, -0.4, -0.2, 0.1, 0.3}},
       {100, {-1, -1, -1, -1, -1}},
       {70, {3, 0, 0, 0, 0}}}};
  std::array<cosmo_background_query, 7> queries{{{0, 0, 0, 0},
                                                 {.1, .1, 0, 0},
                                                 {.3, .3, 0, 0},
                                                 {.6, .6, 0, 0},
                                                 {1, 1, 0, 0},
                                                 {2.5, 2.5, 0, 0},
                                                 {.7, .71, 1, 0}}};
  cosmo_piecewise_batch batch{sizeof(cosmo_piecewise_batch),
                              COSMO_ABI_VERSION,
                              parameters.data(),
                              parameters.size(),
                              sizeof(parameters),
                              queries.data(),
                              queries.size(),
                              sizeof(queries)};
  cosmo_piecewise_policy policy{sizeof(cosmo_piecewise_policy),
                                COSMO_ABI_VERSION,
                                0,
                                8,
                                16,
                                128,
                                1000000,
                                1000};
};
struct View {
  const cosmo_piecewise_slot *rows = nullptr;
  uint64_t n = 0, segments = 0;
  uint32_t status = 999, numerical = 999;
};
View view(cosmo_piecewise_result *r) {
  View v;
  CHECK(cosmo_piecewise_result_view(r, &v.rows, &v.n, &v.status, &v.numerical,
                                    &v.segments) == COSMO_OK);
  return v;
}
cosmo_piecewise_result *run(const cosmo_piecewise_batch &b,
                            const cosmo_piecewise_policy &p) {
  cosmo_piecewise_result *r = nullptr;
  CHECK(cosmo_piecewise_evaluate(&b, &p, &r) == COSMO_OK);
  CHECK(r);
  return r;
}
void failure(const cosmo_piecewise_slot &r) {
  CHECK(!r.has_geometry && !r.has_q && !r.has_jerk && !r.has_q0);
  CHECK(r.q_convention == COSMO_PIECEWISE_Q_CONVENTION_NOT_ASSESSED);
  CHECK(r.jerk_availability == COSMO_PIECEWISE_JERK_AVAILABILITY_NOT_ASSESSED);
  CHECK(r.expansion_E == 0 && r.radial_integral == 0 && r.jerk == 0 &&
        r.assigned_q == 0);
}
int main() {
  try {
    static_assert(sizeof(cosmo_piecewise_parameters) == 48 &&
                  sizeof(cosmo_piecewise_batch) == 56 &&
                  sizeof(cosmo_piecewise_policy) == 56 &&
                  sizeof(cosmo_piecewise_slot) == 328);
    CHECK(offsetof(cosmo_piecewise_parameters, q) == 8);
    CHECK(offsetof(cosmo_piecewise_slot, parameters) == 24);
    Fixture f;
    auto *r = run(f.batch, f.policy);
    auto v = view(r);
    CHECK(v.n == 21 && v.status == COSMO_BACKGROUND_STATUS_OK &&
          v.numerical == COSMO_NUMERICAL_STATUS_OK);
    std::array<Query, 7> q;
    for (std::size_t j = 0; j < q.size(); ++j)
      q[j] = {f.queries[j].z_expansion, f.queries[j].z_observer,
              static_cast<Convention>(f.queries[j].convention)};
    std::size_t remaining = f.policy.maximum_total_segment_visits, total = 0;
    for (std::size_t i = 0; i < f.parameters.size(); ++i) {
      std::array<double, 5> coeff;
      std::copy_n(f.parameters[i].q, 5, coeff.begin());
      auto native =
          prepare_piecewise_q({f.parameters[i].h0_km_s_mpc, coeff})
              .evaluate_batch(q, {f.policy.maximum_queries, remaining});
      remaining -= native.segments_processed;
      total += native.segments_processed;
      for (std::size_t j = 0; j < q.size(); ++j) {
        auto &row = v.rows[i * q.size() + j];
        CHECK(row.parameter_index == i && row.query_index == j);
        CHECK(row.parameters.h0_km_s_mpc == f.parameters[i].h0_km_s_mpc);
        for (int k = 0; k < 5; ++k)
          CHECK(row.parameters.q[k] == coeff[k]);
        CHECK(row.query.z_expansion == q[j].z_expansion &&
              row.query.z_observer == q[j].z_observer);
        CHECK(std::string_view(
                  reinterpret_cast<const char *>(row.model_id.data),
                  row.model_id.length) == PiecewiseBackground::model_id);
        if (native.status != Status::ok) {
          CHECK(row.status == static_cast<uint32_t>(native.status));
          failure(row);
          continue;
        }
        auto &s = native.slots[j];
        CHECK(row.status == static_cast<uint32_t>(s.status));
        CHECK(row.numerical_status ==
              static_cast<uint32_t>(s.numerical_status));
        CHECK(row.segments_processed == s.segments_processed);
        if (s.status != Status::ok) {
          failure(row);
          continue;
        }
        CHECK(row.has_geometry && bool(row.has_q) == s.assigned_q.has_value() &&
              bool(row.has_jerk) == s.jerk.has_value() &&
              bool(row.has_q0) == s.q0_within_piecewise_model.has_value());
        CHECK(row.q_convention == static_cast<uint32_t>(s.q_convention) &&
              row.jerk_availability ==
                  static_cast<uint32_t>(s.jerk_availability));
        if (s.assigned_q)
          CHECK(row.assigned_q == *s.assigned_q);
        if (s.jerk)
          CHECK(row.jerk == *s.jerk);
        if (s.q0_within_piecewise_model)
          CHECK(row.q0_within_piecewise_model == *s.q0_within_piecewise_model);
        auto &g = s.geometry;
        CHECK(row.expansion_E == g.expansion_E &&
              row.h_km_s_mpc == g.h_km_s_mpc &&
              row.radial_integral == g.radial_integral);
        CHECK(row.radial_mpc == g.radial_mpc &&
              row.transverse_mpc == g.transverse_mpc &&
              row.angular_diameter_mpc == g.angular_diameter_mpc &&
              row.luminosity_mpc == g.luminosity_mpc);
        CHECK(row.dimensionless_luminosity_shape ==
                  g.dimensionless_luminosity_shape &&
              row.lookback_seconds == g.lookback_seconds &&
              row.volume_mpc3_per_sr_per_redshift ==
                  g.volume_mpc3_per_sr_per_redshift);
      }
    }
    CHECK(v.segments == total);
    CHECK(v.rows[0].has_q0 && v.rows[0].radial_mpc == 0);
    CHECK(v.rows[1].has_jerk);
    CHECK(!v.rows[2].has_jerk &&
          v.rows[2].q_convention ==
              COSMO_PIECEWISE_Q_CONVENTION_RIGHT_LIMIT_AT_INTERNAL_JUMP);
    CHECK(v.rows[5].jerk_availability ==
          COSMO_PIECEWISE_JERK_AVAILABILITY_ONE_SIDED_ENDPOINT);
    const double retained = v.rows[0].parameters.q[0];
    f.parameters[0].q[0] = 1;
    f.queries[0].z_expansion = 1;
    CHECK(v.rows[0].parameters.q[0] == retained &&
          v.rows[0].query.z_expansion == 0);
    cosmo_piecewise_result_destroy(r);
    f = Fixture{};
    f.batch.parameters = f.parameters.data();
    f.batch.queries = f.queries.data();
    // Global work admits whole queries: first model5, next z1 rejects0 then .1
    // consumes1.
    std::array<cosmo_background_query, 2> work{{{1, 1, 0, 0}, {.1, .1, 0, 0}}};
    auto b = f.batch;
    b.parameter_count = 2;
    b.parameter_bytes = 2 * sizeof(cosmo_piecewise_parameters);
    b.queries = work.data();
    b.query_count = 2;
    b.query_bytes = sizeof(work);
    auto p = f.policy;
    p.maximum_total_segment_visits = 6;
    r = run(b, p);
    v = view(r);
    CHECK(v.n == 4 && v.segments == 6);
    CHECK(v.rows[2].status == COSMO_BACKGROUND_STATUS_WORK_LIMIT &&
          v.rows[2].numerical_status == COSMO_NUMERICAL_STATUS_WORK_LIMIT);
    CHECK(v.rows[2].segments_processed == 0);
    failure(v.rows[2]);
    CHECK(v.rows[3].status == COSMO_BACKGROUND_STATUS_OK &&
          v.rows[3].segments_processed == 1);
    CHECK(v.rows[0].segments_processed + v.rows[1].segments_processed +
              v.rows[2].segments_processed + v.rows[3].segments_processed ==
          v.segments);
    cosmo_piecewise_result_destroy(r);
    // A geometrically unrepresentable query still consumes admitted analytic
    // work.
    work[0] = {1e-200, 1e-200, 0, 0};
    b.parameter_count = 1;
    b.parameter_bytes = sizeof(cosmo_piecewise_parameters);
    p.maximum_total_segment_visits = 1;
    r = run(b, p);
    v = view(r);
    CHECK(v.segments == 1 && v.rows[0].segments_processed == 1);
    CHECK(v.rows[0].status != COSMO_BACKGROUND_STATUS_OK);
    failure(v.rows[0]);
    CHECK(v.rows[1].status == COSMO_BACKGROUND_STATUS_WORK_LIMIT);
    failure(v.rows[1]);
    cosmo_piecewise_result_destroy(r);
    // Equal/interior jump flags around both adjacent floating-point edges.
    for (double edge : piecewise_q_edges) {
      for (double z : {std::nextafter(edge, -INFINITY), edge,
                       std::nextafter(edge, INFINITY)}) {
        cosmo_background_query query{z, z, 0, 0};
        b = f.batch;
        b.parameter_count = 1;
        b.parameter_bytes = sizeof(cosmo_piecewise_parameters);
        b.queries = &query;
        b.query_count = 1;
        b.query_bytes = sizeof(query);
        r = run(b, f.policy);
        v = view(r);
        auto direct =
            prepare_piecewise_q({70, {-0.4, -0.4, -0.2, 0.1, 0.3}})
                .evaluate_batch(
                    std::span<const Query>(std::array<Query, 1>{
                        {{z, z, Convention::geometric_same_redshift}}}),
                    {1, 100});
        CHECK(v.rows[0].status ==
              static_cast<uint32_t>(direct.slots[0].status));
        if (direct.slots[0].status == Status::ok) {
          CHECK(v.rows[0].q_convention ==
                static_cast<uint32_t>(direct.slots[0].q_convention));
          CHECK(bool(v.rows[0].has_jerk) == direct.slots[0].jerk.has_value());
        } else
          failure(v.rows[0]);
        cosmo_piecewise_result_destroy(r);
      }
    }
    b = f.batch;
    p = f.policy;
    p.maximum_slots = 1048576;
    r = run(b, p);
    cosmo_piecewise_result_destroy(r);
    p.maximum_slots = 1048577;
    r = nullptr;
    CHECK(cosmo_piecewise_evaluate(&b, &p, &r) == COSMO_INVALID_INPUT && !r);
    // Unknown convention is a tagged row failure, not invented geometry.
    auto queries = f.queries;
    queries[0].convention = UINT32_MAX;
    b = f.batch;
    b.queries = queries.data();
    r = run(b, f.policy);
    v = view(r);
    CHECK(v.rows[0].status != COSMO_BACKGROUND_STATUS_OK);
    failure(v.rows[0]);
    cosmo_piecewise_result_destroy(r);
    // Counts/caps reject before dereferencing falsely oversized caller storage.
    for (int mode = 0; mode < 5; ++mode) {
      b = f.batch;
      p = f.policy;
      if (mode == 0) {
        b.parameter_count = 4097;
        b.parameter_bytes = 4097 * sizeof(cosmo_piecewise_parameters);
      }
      if (mode == 1) {
        b.query_count = 4097;
        b.query_bytes = 4097 * sizeof(cosmo_background_query);
      }
      if (mode == 2)
        p.maximum_slots = 20;
      if (mode == 3)
        p.maximum_native_output_bytes = 1;
      if (mode == 4)
        p.maximum_parameters = 0;
      r = run(b, p);
      v = view(r);
      CHECK(v.n == 0 && v.status == COSMO_BACKGROUND_STATUS_WORK_LIMIT &&
            v.numerical == COSMO_NUMERICAL_STATUS_WORK_LIMIT);
      cosmo_piecewise_result_destroy(r);
    }
    // Failed source identities retain raw bits, including signed zero and a
    // quiet-NaN payload. Absence flags, rather than finite sentinels, gate use.
    auto params = f.parameters;
    params[0].q[0] = -0.0;
    params[1].h0_km_s_mpc = std::bit_cast<double>(uint64_t{0x7ff8000000000123});
    params[2].q[0] = std::bit_cast<double>(uint64_t{0x7ff8000000000456});
    b = f.batch;
    b.parameters = params.data();
    r = run(b, f.policy);
    v = view(r);
    CHECK(std::bit_cast<uint64_t>(v.rows[0].parameters.q[0]) ==
          std::bit_cast<uint64_t>(-0.0));
    CHECK(std::bit_cast<uint64_t>(v.rows[7].parameters.h0_km_s_mpc) ==
          uint64_t{0x7ff8000000000123});
    CHECK(std::bit_cast<uint64_t>(v.rows[14].parameters.q[0]) ==
          uint64_t{0x7ff8000000000456});
    failure(v.rows[7]);
    failure(v.rows[14]);
    cosmo_piecewise_result_destroy(r);
    // Structural errors return transport cause/null owned result.
    for (int mode = 0; mode < 6; ++mode) {
      b = f.batch;
      p = f.policy;
      if (mode == 0)
        b.abi_version++;
      if (mode == 1)
        p.struct_size--;
      if (mode == 2)
        p.reserved = 1;
      if (mode == 3)
        b.parameter_bytes--;
      if (mode == 4)
        b.parameters = reinterpret_cast<const cosmo_piecewise_parameters *>(
            reinterpret_cast<const char *>(f.parameters.data()) + 1);
      if (mode == 5) {
        b.parameter_count = UINT64_MAX;
        b.parameter_bytes = UINT64_MAX;
      }
      r = reinterpret_cast<cosmo_piecewise_result *>(uintptr_t(1));
      auto status = cosmo_piecewise_evaluate(&b, &p, &r);
      CHECK(status == (mode == 0 ? COSMO_ABI_MISMATCH : COSMO_INVALID_INPUT));
      CHECK(!r);
    }
    b = f.batch;
    b.parameter_count = 0;
    b.parameter_bytes = 0;
    b.parameters = nullptr;
    r = run(b, f.policy);
    v = view(r);
    CHECK(v.n == 0 && v.status == COSMO_BACKGROUND_STATUS_OK &&
          v.segments == 0);
    cosmo_piecewise_result_destroy(r);
    CHECK(cosmo_piecewise_result_destroy(nullptr) == COSMO_OK);
    View bad;
    bad.rows = reinterpret_cast<const cosmo_piecewise_slot *>(uintptr_t(1));
    bad.n = 17;
    CHECK(cosmo_piecewise_result_view(nullptr, &bad.rows, &bad.n, &bad.status,
                                      &bad.numerical,
                                      &bad.segments) == COSMO_INVALID_INPUT);
    CHECK(!bad.rows && bad.n == 0 && bad.segments == 0);
    // Exhaust every actual allocation; exception boundary and ownership remain
    // valid after each partial construction. No fixed allocation count assumed.
    const auto baseline = live;
    unsigned failures = 0;
    bool succeeded = false;
    for (int n = 0; n < 100; ++n) {
      r = nullptr;
      fail_after = n;
      auto status = cosmo_piecewise_evaluate(&f.batch, &f.policy, &r);
      fail_after = -1;
      if (status == COSMO_ALLOCATION_FAILURE) {
        ++failures;
        CHECK(!r);
      } else {
        CHECK(status == COSMO_OK && r);
        cosmo_piecewise_result_destroy(r);
        succeeded = true;
      }
      CHECK(live == baseline);
      if (succeeded)
        break;
    }
    CHECK(succeeded && failures > 0);
    std::printf("piecewise ABI hostile checks %u PASS allocation_failures %u\n",
                checks, failures);
    return 0;
  } catch (const std::exception &e) {
    fail_after = -1;
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
