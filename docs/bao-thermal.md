# Thermal conditional BAO density

The standalone C++20 `irred/bao_thermal.hpp` API evaluates the retained
`bao::PreparedDensity` observation owner with a coarse batch of
`cosmology::ThermalObservableRequest` physical sources. `evaluate_thermal`
prepares one thermal state per model and uses it for all ordered supplied-drag
ratios. The observation matrix and Gaussian factor stay retained and are neither
copied nor refactored for model evaluations. This consumer has no CLI or C ABI
route. The existing massless conditional density and free-ruler hypothesis keep
their distinct identities.

Read [thermal observables](thermal-observables.md) for the explicit physical
baryon/CDM/other-massless densities, Kelvin temperatures, species state weights,
SI conventions, flat closure and supplied drag. Physical density inputs are
`omega_i=Omega_i h^2`; the mapped state holds fractions. Changing H0 with fixed
physical densities and temperatures changes the fractions and generally changes
ratios and density. H0 cancellation applies only when all relevant dimensionless
fractions and relic mass/temperature evolution are preserved explicitly. This
consumer does not implicitly complete a Planck source or predict drag.

The prepared observation owner retains its source role, full ordered SPD ratio
covariance, IDs, calibration and dependence metadata. It currently admits
redshifts 0 through 5. With ordered observed ratios y and model predictions mu,
all declared rows enter the same normalized ratio-coordinate density as the
[massless conditional consumer](bao-conditional.md):

    log p = -1/2 [(y-mu)^T C^-1 (y-mu) + log det C + n log(2 pi)].

The thermal provider supplies D_M/r_s, D_H/r_s and D_V/r_s from one thermal E/H
state and a caller-supplied finite nonnegative drag redshift. Only ratio types
appearing in the queries are requested. The existing `bao::Output` masks request
density, predictions and residuals separately; omitted values remain absent.
The prediction vector always follows the retained observation order. A failed
model or ruler cannot produce a fabricated quadratic. Density projection
refusal preserves successful requested predictions and residuals. Batch finite
status means processing admission, not universal per-model success.

Projection uses the shared conditional-density helper. For r=y-mu, a=C^-1 r
and componentwise provider/subtraction/storage diagnostics eps, it estimates

    ||a||inf ||eps||1 + 1/2 estimated||C^-1||inf ||eps||inf ||eps||1.

The inverse norm estimate comes from the retained factor's condition estimate
divided by the covariance infinity norm. The default absolute log-density
projection allowance stays **1e-8**. Thermal ratio admission is separately
1e-10 + 5e-10 relative by default and cannot guarantee that arbitrary covariance
and residuals satisfy the density allowance. A refusal reports the estimate and
cause; neither observational covariance inflation nor a tolerance relaxation is
performed. Diagnostics are empirical numerical estimates, not certified bounds
or independent random uncertainties. Factor/solve sensitivity remains a
separate caller requirement. Density arithmetic must match the retained
Gaussian; thermal work requires its supported wide arithmetic and rounding.

`ThermalDensityPolicy` bounds model/query count, copied origin bytes,
simultaneous owned payload and whole-batch callbacks before result allocation.
Each model retains an independently owned physical source in its slot; mapping
or closure failures remain inspectable. Provider resource limits apply as well.
The payload helper includes all returned source species/origins/vectors and
serial mapping/preparation/provider/solve scratch. Retained observations,
borrowed inputs, allocator overhead, stack frames and RSS are excluded.

`preparation_status` reports mapping/background preparation admission.
`preparation_callbacks` records present-day relic normalization work, including
failed normalization. `outer_callbacks` counts distance/ruler quadrature;
`momentum_callbacks` includes preparation and nested momentum work. Thus
`callbacks=outer_callbacks+momentum_callbacks`, while preparation is a subset
of momentum and must not be added twice. Preparation and evaluation share the
provider per-model total and momentum ceilings; their actual combined work also
reduces the consumer whole-batch ceiling before the next model. Mapping itself
performs no quadrature callbacks. No observation readback, runtime model plugin,
server or unbounded cache is needed.

Permanent synthetic controls retain massive direct-momentum high-precision and
matched CLASS facts from the thermal peer's explicit source; they independently
compose ratio density with 3x3 covariance cofactors. Original massless-polynomial
controls use direct-z and sqrt(a) composite Gauss–Legendre integration with
refinement at most 5% of the frozen ratio and propagated density allocations.
Density quadratic, determinant, normalization and log density retain an absolute
1e-8 comparison allocation. Tests cover fixed physical versus fixed fractional
scaling, ordered covariance permutations, partial failures, projection refusal,
quotas, nested preparation accounting and source lifetime. Shared physical
constants and equations are shared ancestry, not independent observations.

These controls supply named numerical evidence only. Released BAO fitted
compressions require a separate physical-validity review for this thermal
hypothesis. No observational qualification, posterior, cross-probe combination,
full thermal BAO compression qualification, recombination or drag prediction is
claimed.
