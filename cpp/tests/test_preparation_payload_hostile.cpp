//! Allocation-size instrumentation of preparation peak, independent of bound
//! formula.
#include "irred/supernova.hpp"
#include "irred/bao.hpp"
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
namespace {
struct alignas(std::max_align_t) Header {
  size_t size;
  bool tracked;
};
bool tracking = false;
size_t retained = 0, peak = 0;
unsigned checks = 0;
void check(bool x, const char *why) {
  ++checks;
  if (!x)
    throw std::runtime_error(why);
}
} // namespace
void *operator new(size_t n) {
  auto raw = std::malloc(sizeof(Header) + (n ? n : 1));
  if (!raw)
    throw std::bad_alloc();
  auto h = new (raw) Header{n, tracking};
  if (tracking) {
    retained += n;
    if (retained > peak)
      peak = retained;
  }
  return h + 1;
}
void operator delete(void *p) noexcept {
  if (p) {
    auto h = static_cast<Header *>(p) - 1;
    if (h->tracked)
      retained -= h->size;
    std::free(h);
  }
}
void operator delete(void *p, size_t) noexcept { operator delete(p); }
int main() {
  using namespace irred;
  try {
    observations::Input in;
    in.profile = observations::Profile::typed_magnitude_covariance;
    in.role = observations::Role::observed_measurement;
    in.unit = observations::Unit::magnitude;
    in.calibration = observations::Calibration::unknown;
    in.uncertainty = observations::Uncertainty::covariance;
    in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
    in.component = observations::Component::total;
    in.table_sha256 = std::string(64, 'a');
    in.uncertainty_sha256 = std::string(64, 'b');
    in.ordering_provenance = std::string(2500, 'o');
    in.calibration_provenance = std::string(1800, 'c');
    in.dependence_provenance = std::string(1700, 'd');
    for (size_t i = 0; i < 5; ++i) {
      in.values.push_back(20 + i);
      in.measurement_ids.push_back(std::string(500 + i * 73, 'a' + i));
      in.event_ids.push_back(in.measurement_ids.back());
      in.missing.push_back(0);
      in.quality.push_back(0);
    }
    in.uncertainty_axis_ids = in.measurement_ids;
    in.uncertainty_matrix.assign(25, .1);
    for (size_t i = 0; i < 5; ++i)
      in.uncertainty_matrix[i * 5 + i] = 2 + i;
    auto prepared = observations::prepare(std::move(in), {5, 25, 100000});
    check(prepared.status() == observations::Status::ok,
          "source prepared outside measurement");
    auto shared =
        std::make_shared<const observations::Prepared>(std::move(prepared));
    for (auto arithmetic : {numerics::Arithmetic::binary64_legacy_v1,
                            numerics::Arithmetic::longdouble_cpu_v1})
      for (bool reserve_heavy : {false, true}) {
        auto bound =
            supernova::preparation_payload_bound(*shared, 2, arithmetic);
        check(bound.has_value(), "declared peak");
        tracking = true;
        retained = peak = 0;
        supernova::SelectedMagnitudeSource selected{
            shared,
            {0, 4},
            {shared->source().measurement_ids[0],
             shared->source().measurement_ids[4]},
            {{.1, {.1, cosmology::Convention::geometric_same_redshift}},
             {.4, {.4, cosmology::Convention::geometric_same_redshift}}}};
        size_t extra = 0;
        if (reserve_heavy) {
          const auto old = selected.source_indices.capacity();
          selected.source_indices.reserve(113);
          extra = (selected.source_indices.capacity() - old) * sizeof(size_t);
        }
        // Caller construction is outside preparation; transferred storage
        // remains in scope.
        peak = retained;
        const size_t allowance = *bound + extra;
        {
          auto consumer = supernova::prepare(
              std::move(selected), {arithmetic, 2, 4, 1e-8, allowance});
          check(consumer.status() == supernova::Status::ok,
                "admitted partial selection");
          const auto observed = peak;
          tracking = false;
          check(observed <= allowance, "actual preparation heap peak bounded");
          check(consumer.selected_source().source.get() == shared.get(),
                "shared source not duplicated");
          std::cout << "arithmetic " << (unsigned)arithmetic << " reserve "
                    << reserve_heavy << " observed " << observed << " bound "
                    << allowance << '\n';
        }
        check(retained == 0, "tracked preparation allocations recovered");
      }
    for (bool reserve_heavy : {false, true}) {
      tracking = true;
      retained = peak = 0;
      {
        bao::DensityInput input;
        input.queries = {{.1, bao::Observable::transverse_over_ruler},
                         {.2, bao::Observable::hubble_over_ruler},
                         {.3, bao::Observable::volume_over_ruler}};
        input.observed = {4, 30, 15};
        input.covariance = {4, 1, 0, 1, 3, 1, 0, 1, 2};
        input.ordered_ids = {std::string(503, 'a'), std::string(601, 'b'),
                             std::string(709, 'c')};
        input.role = bao::RowRole::synthetic_control;
        input.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
        input.table_identity = std::string(64, 'a');
        input.covariance_identity = std::string(64, 'b');
        input.ordering_provenance = std::string(1700, 'o');
        input.calibration_provenance = std::string(1500, 'c');
        input.dependence_provenance = std::string(1300, 'd');
        if (reserve_heavy) {
          input.covariance.reserve(317);
          input.ordered_ids[0].reserve(1801);
          input.dependence_provenance.reserve(3001);
        }
        const auto source = bao::retained_source_payload_bound(input);
        check(source.has_value(), "BAO transferred capacity admitted to bound");
        size_t identity = 0;
        for (const auto &id : input.ordered_ids)
          identity += id.capacity() + 1;
        for (const auto *s : {&input.table_identity, &input.covariance_identity,
                              &input.ordering_provenance,
                              &input.calibration_provenance,
                              &input.dependence_provenance})
          identity += s->capacity() + 1;
        const auto bound = bao::preparation_payload_bound(3, *source, identity);
        check(bound.has_value(), "BAO preparation bound representable");
        peak = retained;
        auto density = bao::prepare_density(
            std::move(input), {3, 9, 20000, *bound, 1e-8,
                              numerics::Arithmetic::longdouble_cpu_v1});
        check(density.status() == statistics::DensityStatus::finite,
              "BAO instrumented preparation succeeds");
        const auto observed = peak;
        tracking = false;
        check(observed <= *bound, "BAO actual simultaneous heap within bound");
        std::cout << "BAO reserve " << reserve_heavy << " observed " << observed
                  << " bound " << *bound << '\n';
      }
      check(retained == 0, "BAO tracked preparation allocations recovered");
    }
    // Borrowed matrix storage is intentionally outside these Gaussian scopes.
    const std::array<double,4> matrix{4,1,1,9};
    const std::array<double,2> response{1,1};
    for (auto kind : {statistics::MatrixKind::covariance,
                      statistics::MatrixKind::precision}) {
      tracking = true;
      retained = peak = 0;
      {
        statistics::Metadata metadata;
        metadata.ordered_ids = {std::string(501,'a'),std::string(603,'b')};
        metadata.measure = "magnitude";
        metadata.ordering_provenance = std::string(1701,'o');
        metadata.ordering_provenance.reserve(3001);
        metadata.ordered_ids.reserve(113);
        auto bound = statistics::gaussian_preparation_payload_bound(
            2, kind, numerics::Arithmetic::longdouble_cpu_v1, metadata);
        check(bound.has_value(), "Gaussian preparation bound representable");
        peak = retained;
        auto gaussian = statistics::prepare_gaussian(
            matrix, kind, std::move(metadata), 4, 1e-8,
            numerics::Arithmetic::longdouble_cpu_v1);
        check(gaussian.status() == statistics::DensityStatus::finite,
              "Gaussian reserve-heavy metadata preparation succeeds");
        check(peak <= *bound, "Gaussian actual preparation within bound");
        const auto preparation_peak = peak;
        // Independent observed simultaneous heap, with the prepared owner still
        // live. The helper covers sequential native scratch, not caller buffers.
        const std::array<double,2> residual{1,2};
        for (bool profile : {false,true}) {
          const auto base=retained;
          auto scratch_bound=gaussian.evaluation_payload_bound(3,profile);
          check(scratch_bound.has_value(), "evaluation scratch bound representable");
          peak=retained;
          for (size_t row=0;row<3;++row) {
            if (profile) {
              auto value=gaussian.profile_offset(residual,response,gaussian.metadata().ordered_ids,1e-8);
              check(value.status==statistics::DensityStatus::finite,"profile evaluation succeeds");
              check(peak-base<=*scratch_bound,"profile simultaneous native heap bounded");
            } else {
              auto value=gaussian.evaluate(residual,gaussian.metadata().ordered_ids,1e-8);
              check(value.density.status==statistics::DensityStatus::finite,"density evaluation succeeds");
              check(peak-base<=*scratch_bound,"density simultaneous native heap bounded");
            }
            check(retained==base,"sequential evaluation workspace recovered");
          }
          std::cout << "Gaussian evaluation profile " << profile << " observed " << peak-base << " bound " << *scratch_bound << '\n';
        }

        // Distinct proper priors grow latent/history identity storage. This
        // measures old owner + candidate coexistence, not candidate alone.
        for (size_t prior : {0u,1u}) {
          std::string identity(701+prior*127,'x'+prior);
          auto transform_bound = gaussian.proper_offset_payload_bound(identity);
          check(transform_bound.has_value(), "proper transform bound");
          peak = retained;
          auto next = gaussian.proper_offset(response, gaussian.metadata().ordered_ids,
              .25, 2, identity, true, 4, 1e-8);
          check(next.status() == statistics::DensityStatus::finite,
                "proper prior transform succeeds");
          check(peak <= *transform_bound,
                "old and candidate simultaneous proper-prior payload bounded");
          gaussian = std::move(next);
        }
        tracking = false;
        std::cout << "Gaussian kind " << (unsigned)kind << " preparation observed "
                  << preparation_peak << " bound " << *bound << '\n';
      }
      check(retained == 0, "Gaussian tracked allocations recovered");
    }
    for (auto kind : {observations::Uncertainty::covariance,
                      observations::Uncertainty::precision}) {
      auto source_input = shared->source();
      source_input.uncertainty = kind;
      if (kind == observations::Uncertainty::precision) {
        source_input.profile = observations::Profile::gaussian_fixture_v1;
        source_input.role = observations::Role::synthetic_control;
      }
      source_input.uncertainty_unit = kind == observations::Uncertainty::covariance
          ? observations::UncertaintyUnit::magnitude_squared
          : observations::UncertaintyUnit::inverse_magnitude_squared;
      auto source = observations::prepare(std::move(source_input), {5,25,100000});
      check(source.status() == observations::Status::ok,"partial source prepared");
      auto bound = statistics::selected_gaussian_preparation_payload_bound(
          source, 2, numerics::Arithmetic::longdouble_cpu_v1);
      check(bound.has_value(), "selected Gaussian preparation bound");
      tracking = true;
      retained = peak = 0;
      {
        const std::array<size_t,2> indices{0,4};
        auto gaussian = statistics::prepare_selected_observations(source, indices,
            25, 1e-8, numerics::Arithmetic::longdouble_cpu_v1);
        check(gaussian.status() == statistics::DensityStatus::finite,
              "selected covariance/full precision marginal succeeds");
        const auto observed = peak;
        tracking = false;
        check(observed <= *bound,"selected Gaussian actual peak bounded");
        check(gaussian.metadata().matrix_validation_scope ==
            (kind == observations::Uncertainty::covariance
             ? statistics::MatrixValidationScope::selected_covariance_only
             : statistics::MatrixValidationScope::full_precision_then_marginal),
             "selection matrix semantics preserved");
        std::cout << "selected kind " << (unsigned)kind << " observed "
                  << observed << " bound " << *bound << '\n';
      }
      check(retained == 0,"selected Gaussian tracked allocations recovered");
    }
    std::cout << "Independent preparation payload PASS " << checks << '\n';
  } catch (const std::exception &e) {
    tracking = false;
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
