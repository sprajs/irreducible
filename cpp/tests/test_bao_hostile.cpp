// Independent peer: scale-factor composite Simpson, direct antiderivatives,
// and 2x2 adjugate Gaussian. Shared libm/equations, distinct quadrature route.
// Fixed observable budget 2e-12+2e-10|reference|; reference refinement <=10%.
// Synthetic normalized density budget 1e-8; no original asset/engine
// dependency.
#include "irred/bao.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
int checks = 0;
long double worst = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
void close(double x, long double y) {
  auto b = 2e-12L + 2e-10L * std::abs(y);
  worst = std::max(worst, std::abs(x - y) / b);
  check(std::abs(x - y) <= b, "observable budget");
}
long double radial(double z, double om, double w, double wa, int n) {
  long double lo = 1 / (1 + (long double)z), h = (1 - lo) / n, s = 0;
  for (int i = 0; i <= n; ++i) {
    long double a = lo + i * h;
    auto E =
        std::sqrt(om / (a * a * a) + (1 - om) * std::pow(a, -3 * (1 + w + wa)) *
                                         std::exp(-3 * wa * (1 - a)));
    s += (i == 0 || i == n ? 1 : i % 2 ? 4 : 2) / (a * a * E);
  }
  return s * h / 3;
}
} // namespace
int main() {
  try {
    bao::Policy p;
    const long double scale = 299792.458L / 10000;
    for (double q : {-1., -.5, 0., .5})
      for (double z : {0., .0001, .51, 2.33, 5.}) {
        auto bg = cosmology::prepare(
            {cosmology::Model::constant_q_flat_v1, 70, 0, q});
        std::array<bao::Query, 3> qs{
            {{z, bao::Observable::transverse_over_ruler},
             {z, bao::Observable::hubble_over_ruler},
             {z, bao::Observable::volume_over_ruler}}};
        auto out = bao::evaluate(bg, bao::Ruler(10000), qs, p);
        check(out.status == cosmology::Status::ok, "batch");
        long double I = q == 0
                            ? std::log1p((long double)z)
                            : -std::expm1(-q * std::log1p((long double)z)) / q;
        long double dm = scale * I,
                    dh = scale / std::pow(1 + (long double)z, 1 + q);
        for (auto &s : out.slots)
          check(s.status == cosmology::Status::ok, "analytic row");
        close(out.slots[0].dimensionless_value, dm);
        close(out.slots[1].dimensionless_value, dh);
        close(out.slots[2].dimensionless_value, std::cbrt(z * dm * dm * dh));
      }
    for (auto pars : {std::array<double, 3>{.3, -.9, .4},
                      std::array<double, 3>{.3, -1.1, -.4}})
      for (double z : {.1, .7, 2.33, 5.}) {
        auto bg = cosmology::prepare_cpl({70, pars[0], pars[1], pars[2]});
        std::array<bao::Query, 1> qs{
            {{z, bao::Observable::transverse_over_ruler}}};
        auto out = bao::evaluate(bg, bao::Ruler(10000), qs, p);
        auto a = scale * radial(z, pars[0], pars[1], pars[2], 8192),
             b = scale * radial(z, pars[0], pars[1], pars[2], 16384);
        check(std::abs(a - b) <= .1L * (2e-12L + 2e-10L * std::abs(b)),
              "independent refinement");
        close(out.slots[0].dimensionless_value, b);
      }
    bao::DensityInput in;
    in.queries = {{1, bao::Observable::transverse_over_ruler},
                  {2, bao::Observable::hubble_over_ruler}};
    in.observed = {(double)scale + 2, (double)scale - 3};
    in.covariance = {4, 1, 1, 9};
    in.ordered_ids = {"transverse", "hubble"};
    in.role = bao::RowRole::synthetic_control;
    in.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
    in.table_identity = "synthetic";
    in.covariance_identity = "exact integer covariance";
    in.ordering_provenance = "explicit IDs";
    in.calibration_provenance = "unknown";
    in.dependence_provenance = "unknown";
    bao::DensityPolicy dp{p,
                          10,
                          100,
                          10000,
                          1000000,
                          1e-10,
                          numerics::Arithmetic::longdouble_cpu_v1};
    auto prepared = bao::prepare_density(in, dp);
    check(prepared.status() == statistics::DensityStatus::finite,
          "density prepare");
    auto bg =
        cosmology::prepare({cosmology::Model::constant_q_flat_v1, 70, 0, -1});
    std::array<bao::ModelQuery, 1> models{{{bg, bao::Ruler(10000)}}};
    auto batch = prepared.evaluate(models, dp);
    check(batch.slots[0].result.density.status ==
              statistics::DensityStatus::finite,
          "normalized density");
    auto r0 = (long double)in.observed[0] - scale,
         r1 = (long double)in.observed[1] - scale;
    auto q = (9 * r0 * r0 - 2 * r0 * r1 + 4 * r1 * r1) / 35;
    auto ref = -.5L * (q + std::log(35.L) + 2 * std::log(2 * std::acos(-1.L)));
    check(std::abs(batch.slots[0].result.density.log_value - ref) < 1e-8L,
          "adjugate normalized reference");
    check(prepared.metadata().dependence_provenance == "unknown",
          "unknown dependence retained");
    auto permuted = in;
    std::swap(permuted.queries[0], permuted.queries[1]);
    std::swap(permuted.observed[0], permuted.observed[1]);
    std::swap(permuted.ordered_ids[0], permuted.ordered_ids[1]);
    permuted.covariance = {9, 1, 1, 4};
    auto reordered = bao::prepare_density(permuted, dp);
    auto rebatch = reordered.evaluate(models, dp);
    check(reordered.metadata().ordered_ids == permuted.ordered_ids,
          "declared source order preserved");
    check(std::abs(rebatch.slots[0].result.density.log_value -
                   batch.slots[0].result.density.log_value) < 1e-12,
          "consistent complete coordinate permutation invariant");
    auto moved = std::move(reordered);
    check(moved.evaluate(models, dp).slots[0].result.density.status ==
              statistics::DensityStatus::finite,
          "moved-to prepared target usable");
    auto abandoned = reordered.evaluate(models, dp);
    check(abandoned.slots.empty() || abandoned.slots[0].result.density.status !=
                                         statistics::DensityStatus::finite,
          "moved-from target cannot publish finite density");
    bao::PreparedDensity empty;
    check(empty.evaluate(models, dp).status !=
              statistics::DensityStatus::finite,
          "default target rejected");
    auto bad = in;
    bad.ordered_ids[1] = bad.ordered_ids[0];
    check(bao::prepare_density(bad, dp).status() !=
              statistics::DensityStatus::finite,
          "duplicate IDs");
    bad = in;
    bad.covariance[1] = 2;
    check(bao::prepare_density(bad, dp).status() !=
              statistics::DensityStatus::finite,
          "asymmetry");
    auto tight = dp;
    tight.maximum_forward_sensitivity = 1e-30;
    auto failed = prepared.evaluate(models, tight);
    check(failed.slots[0].numerical_status ==
              numerics::Status::conditioning_budget_exceeded,
          "precise conditioning cause");
    check(failed.slots[0].predictions.empty(), "failed payload cleared");
    auto cap = dp;
    cap.maximum_models = 0;
    check(prepared.evaluate(models, cap).numerical_status ==
              numerics::Status::work_limit,
          "model cap");
    std::array<bao::Query, 1> invalid{{{1, static_cast<bao::Observable>(99)}}};
    check(bao::evaluate(bg, bao::Ruler(10000), invalid, p).slots[0].status ==
              cosmology::Status::invalid_input,
          "unknown observable");
    auto huge = std::span<const bao::Query>(
        (const bao::Query *)nullptr, std::numeric_limits<std::size_t>::max());
    check(bao::evaluate(bg, bao::Ruler(10000), huge, p).status ==
              cosmology::Status::work_limit,
          "pre-dereference cap");
    check(bao::evaluate(bg, bao::Ruler(4999), invalid, p).status ==
              cosmology::Status::unsupported_domain,
          "ruler domain");
    std::printf("{\"suite\":\"bao_peer\",\"checks\":%d,\"maximum_budget_"
                "fraction\":%.17Lg,\"passed\":true}\n",
                checks, worst);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
