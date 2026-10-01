# Conditional BAO density

The native `irred/bao_conditional.hpp` consumer evaluates the retained
`bao::PreparedDensity` observation object against a coarse batch of
`cosmology::SoundHorizonRequest` points. Observation acquisition, exact row
order, full SPD ratio covariance, factorization, source role, calibration and
dependence metadata remain owned by that prepared object. The factor is reused;
no covariance is copied or refactored per point. This native consumer leaves the
free-ruler likelihood and CLI as their distinct physical hypothesis.

One pressureless-matter, massless-radiation, flat Lambda model supplies both the
late distances and the conditional tight-coupling sound horizon. H0 is in
km/s/Mpc, distances and the ruler in Mpc, baryons are a subset of matter and
photons a subset of radiation. The source supplies finite nonnegative drag
redshift and a nonempty origin. The precise domain and equations are those of
[early/late predictions](early-late.md) and [the sound horizon](sound-horizon.md).
The prepared observations currently admit redshifts 0 through 5. No independent
free H0*r_d enters this consumer.

For retained ordered ratios y, predictions mu and full covariance C, the
normalized density in the product of dimensionless ratio coordinates is

    log p = -1/2 [(y-mu)^T C^-1 (y-mu) + log det C + n log(2 pi)].

All declared covariance rows remain present. Numerical error is never added to
observational C. At fixed density fractions and supplied drag, both distances
and sound horizon scale as 1/H0: these ratios and this density are H0 invariant.
This is neither drag prediction, CMB or growth qualification, nor observational
H0 inference. Using a released fitted compression remains conditional on its
physical validity for the supplied hypothesis; numerical acceptance cannot
establish that validity.

Existing `bao::Output` bits separately request normalized density, predictions
and residuals. Only required ratio types are requested from the shared provider;
unrequested vectors/results are absent. Each point and requested output group
retains its own cause. Prediction success survives projection-budget refusal;
a failed ruler supplies neither usable ratios nor fabricated quadratic.
Batch finite status indicates admission, not success of every point.

Policies bound point count, retained query count, origin string bytes,
simultaneous owned payload and global callbacks before result allocation or
numerical scans. Global callbacks include every supplied-ruler evaluation and
all distance work, including failures. The provider's own limits apply too.
The payload helper charges returned slots, copied origins, requested vectors,
serial early/late rows and numerical scratch; retained observation objects,
borrowed input, allocator overhead, stack frames and RSS are excluded.

Projection admission is separate from factor/solve sensitivity. With residual
r, solve a=C^-1 r and componentwise prediction/subtraction/cast diagnostics eps,
the consumer estimates a log-density perturbation by

    ||a||inf ||eps||1 + 1/2 estimated||C^-1||inf ||eps||inf ||eps||1.

The inverse norm estimate uses the retained factor's condition estimate divided
by the covariance infinity norm. This and provider quadrature diagnostics are
empirical estimates, not rigorous universal error bounds. The strict caller
projection allowance defaults to 1e-8. An absent density result and failed state
report refusal without invalidating already successful prediction outputs.
Arithmetic must match the retained Gaussian; provider arithmetic additionally
requires round-to-nearest and the supported long-double host contract.

Permanent synthetic controls compare direct-z distances and a sqrt(a) ruler
with composite 8-point Gauss-Legendre refinement, distinct from production's
coordinates and adaptive rule. Ratios use the inherited 1e-11+5e-11*|reference|
comparison allocation; quadratic, log determinant, normalization and log density
use a separately named absolute 1e-8 allocation. Reference refinement occupies
at most 5% of its allocation. Explicit 3x3 covariance cofactors are independent
of production Cholesky. Radiation analytic limits, H0 cancellation, supplied
drag and baryon loading, row/covariance permutation, mixed failed points,
projection refusal, omitted outputs, quotas and owner moves are covered.
These controls qualify their named synthetic cases only; no original released
BAO data comparison or joint inference is claimed by them.
