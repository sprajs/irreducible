#include "irred/abi.h"
#include "irred/background.hpp"
#include <bit>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
namespace c = irred::cosmology;
namespace {
unsigned checks = 0;
long fail_after = -1, live = 0;
void check(bool x, const char *why) {
  ++checks;
  if (!x)
    throw std::runtime_error(why);
}
irred_f64_buffer f64(const double *p, size_t n) {
  return {sizeof(irred_f64_buffer), IRRED_ABI_VERSION, 2, 0, p, n,
          n * sizeof(double)};
}
struct Fixture {
  double q[2]{0, 3};
  irred_expansion_spec models[2]{};
  irred_expansion_request queries[3]{};
  irred_expansion_batch b{};
  irred_expansion_policy p{};
  Fixture() {
    for (int i = 0; i < 2; ++i)
      models[i] = {sizeof(models[i]), IRRED_ABI_VERSION, 1, 0, f64(q + i, 1)};
    for (int i = 0; i < 3; ++i) {
      queries[i].struct_size = sizeof(queries[i]);
      queries[i].abi_version = IRRED_ABI_VERSION;
      queries[i].requested = 32;
      queries[i].z_expansion = i == 2 ? .5 : 1;
    }
    b.struct_size = sizeof(b);
    b.abi_version = IRRED_ABI_VERSION;
    b.models = models;
    b.model_count = 1;
    b.model_byte_length = sizeof(models[0]);
    b.queries = queries;
    b.query_count = 3;
    b.query_byte_length = sizeof(queries);
    p.struct_size = sizeof(p);
    p.abi_version = IRRED_ABI_VERSION;
    p.maximum_models = 2;
    p.maximum_queries = 3;
    p.maximum_slots = 6;
    p.maximum_native_bytes = 1000000;
  }
};
struct View {
  const irred_expansion_row *rows = nullptr;
  uint64_t count = 0, callbacks = 0, segments = 0;
  uint32_t status = 0, num = 0;
  void read(irred_expansion_result *r) {
    check(irred_expansion_result_view(r, &rows, &count, &status, &num,
                                      &callbacks, &segments) == IRRED_OK,
          "result view");
  }
};
} // namespace
void *operator new(size_t n) {
  if (fail_after == 0)
    throw std::bad_alloc();
  if (fail_after > 0)
    --fail_after;
  auto p = std::malloc(n ? n : 1);
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
void operator delete(void *p, size_t) noexcept { operator delete(p); }
int main() {
  try {
    Fixture f;
    irred_expansion_result *r = nullptr;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_OK && r,
          "E-only owner");
    View v;
    v.read(r);
    check(v.count == 3 && v.callbacks == 0 && v.segments == 0,
          "E only zero work");
    auto native = c::prepare(c::ConstantQ(0), c::FlatFLRW{});
    c::EvaluationPolicy np;
    np.maximum_queries = 3;
    np.maximum_native_bytes = 1000000;
    np.maximum_callbacks = 0;
    np.maximum_segment_visits = 0;
    c::Request nr[3]{{1, 32}, {1, 32}, {.5, 32}};
    auto expected = native.evaluate(nr, np);
    for (size_t i = 0; i < 3; ++i) {
      check(v.rows[i].expansion.state.availability == 1, "available E");
      check(std::bit_cast<uint64_t>(v.rows[i].expansion.E) ==
                std::bit_cast<uint64_t>(
                    expected.slots[i].expansion.value->expansion_E),
            "native E bits");
      check(v.rows[i].radial.state.availability == 0, "unrequested radial");
    }
    const irred_expansion_node *nodes = nullptr;
    uint64_t nn = 0;
    check(irred_expansion_result_nodes(r, &nodes, &nn) == IRRED_OK && nn == 2,
          "exact duplicate nodes");
    check(v.rows[0].node_index == v.rows[1].node_index,
          "duplicate map identity");
    const irred_expansion_model_view *ms = nullptr;
    uint64_t nm = 0;
    uint32_t available = 0;
    check(irred_expansion_result_models(r, &ms, &nm, &available) == IRRED_OK &&
              nm == 1 && available,
          "model source");
    f.q[0] = -1;
    check(ms[0].source.parameters.data[0] == 0, "owned model parameter");
    const irred_expansion_request *qs = nullptr;
    uint64_t nq = 0;
    check(irred_expansion_result_queries(r, &qs, &nq, &available) == IRRED_OK &&
              nq == 3,
          "query source");
    f.queries[0].z_expansion = 2;
    check(qs[0].z_expansion == 1, "owned query");
    irred_expansion_result_destroy(r);
    f = Fixture{}; // Rebind copied internal pointers explicitly.
    f.models[0].parameters = f64(f.q, 1);
    f.models[1].parameters = f64(f.q + 1, 1);
    f.b.models = f.models;
    f.b.queries = f.queries;
    f.b.model_count = 2;
    f.b.model_byte_length = sizeof(f.models);
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_OK,
          "mixed model owner");
    v.read(r);
    check(v.count == 3, "invalid model no fabricated rows");
    irred_expansion_result_models(r, &ms, &nm, &available);
    check(nm == 2 && ms[1].preparation_status != 0 && ms[1].row_count == 0 &&
              ms[1].source.parameters.data[0] == 3,
          "invalid attempted model retained");
    irred_expansion_result_destroy(r);
    f.b.model_count = 1;
    f.b.model_byte_length = sizeof(f.models[0]);
    f.p.has_integration = 1;
    f.p.max_depth = 24;
    f.p.absolute_tolerance = 1e-12;
    f.p.relative_tolerance = 1e-12;
    f.p.integration_max_evaluations = 100000;
    f.p.maximum_callbacks = 100000;
    f.queries[0].requested = 36;
    f.queries[0].presence_flags = 2;
    f.queries[0].h0_km_s_mpc = 0;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_OK,
          "nested H/time owner");
    v.read(r);
    check(v.rows[0].expansion.state.availability == 1 &&
              v.rows[0].expansion.H_km_s_mpc.state.availability == 3,
          "invalid H preserves E");
    check(v.rows[0].clock.state.availability == 1 &&
              v.rows[0].clock.lookback_seconds.state.availability == 3,
          "invalid seconds preserves clock");
    irred_expansion_result_destroy(r);
    f.queries[0].requested = 64;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_INVALID_INPUT && !r,
          "unknown mask structural");
    f.queries[0].requested = 32;
    f.p.maximum_native_bytes = 0;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_OK && r,
          "zero quota diagnostic owner");
    v.read(r);
    check(v.count == 0 && v.status == 5 && v.callbacks == 0 && v.segments == 0,
          "quota empty worklimit");
    irred_expansion_result_models(r, &ms, &nm, &available);
    check(!available && nm == 0, "quota source unavailable");
    irred_expansion_result_destroy(r);
    f.p.maximum_native_bytes = 1000000;
    f.b.model_count = 65537;
    f.b.model_byte_length = 65537 * sizeof(irred_expansion_spec);
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_INVALID_INPUT && !r,
          "huge count before dereference");
    f.b.model_count = 1;
    f.b.model_byte_length = sizeof(f.models[0]);
    f.p.abi_version = 0;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_ABI_MISMATCH && !r,
          "version resets owner");
    f.p.abi_version = IRRED_ABI_VERSION;
    f.queries[0].presence_flags = 0;
    f.queries[0].h0_km_s_mpc = 0;
    // Two analytic models share one global segment allowance, never reset.
    double bins[5]{-1, 0, -.5, .2, 1};
    auto oldmodel = f.models[0];
    for (auto &m : f.models)
      m = {sizeof(m), IRRED_ABI_VERSION, 3, 0, f64(bins, 5)};
    f.b.model_count = 2;
    f.b.model_byte_length = sizeof(f.models);
    for (auto &query : f.queries)
      query.requested = 1;
    f.p.maximum_segment_visits = 5;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_OK,
          "global analytic owner");
    v.read(r);
    check(v.count == 6 && v.segments <= 5 && v.callbacks == 0,
          "global segments capped");
    check(v.rows[0].radial.state.availability == 1 &&
              v.rows[3].radial.state.availability == 3 &&
              v.rows[3].radial.state.numerical_status ==
                  (uint32_t)irred::numerics::Status::work_limit,
          "second model retains precise exhausted work cause");
    irred_expansion_result_destroy(r);
    f.models[0] = oldmodel;
    f.b.model_count = 1;
    f.b.model_byte_length = sizeof(f.models[0]);
    for (auto &query : f.queries)
      query.requested = 32;
    f.p.maximum_segment_visits = 0;
    auto count = f.b.query_count;
    auto bytes = f.b.query_byte_length;
    f.b.query_count = f.b.query_byte_length = 0;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_OK,
          "zero query owner");
    v.read(r);
    check(v.count == 0 && v.status == 0 && v.callbacks == 0 && v.segments == 0,
          "zero query success");
    irred_expansion_result_destroy(r);
    f.b.query_count = count;
    f.b.query_byte_length = bytes;
    auto mc = f.b.model_count;
    auto mb = f.b.model_byte_length;
    f.b.model_count = f.b.model_byte_length = 0;
    check(irred_expansion_evaluate(&f.b, &f.p, &r) == IRRED_OK,
          "zero model owner");
    v.read(r);
    check(v.count == 0 && v.status == 0, "zero model success");
    irred_expansion_result_destroy(r);
    f.b.model_count = mc;
    f.b.model_byte_length = mb;
    long baseline = live;
    unsigned failures = 0;
    for (long i = 0; i < 100; ++i) {
      fail_after = i;
      auto rc = irred_expansion_evaluate(&f.b, &f.p, &r);
      fail_after = -1;
      if (rc == IRRED_OK) {
        check(r != nullptr, "allocation sweep success owner");
        irred_expansion_result_destroy(r);
        check(live == baseline, "success cleanup");
        break;
      }
      check(rc == IRRED_ALLOCATION_FAILURE && !r, "allocation failure null");
      check(live == baseline, "allocation cleanup");
      ++failures;
    }
    check(failures > 5, "multiple allocation sites");
    check(irred_expansion_result_destroy(nullptr) == IRRED_OK, "null destroy");
    std::cout << "current expansion ABI peer PASS " << checks
              << " allocation failures " << failures << '\n';
  } catch (const std::exception &e) {
    fail_after = -1;
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
