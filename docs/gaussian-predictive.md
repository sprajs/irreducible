# Retained joint Gaussian prediction

The native C++ `irred/gaussian_predictive.hpp` consumer predicts distinct
synthetic calibration observations after conditioning the existing
[proper Gaussian parameter posterior](gaussian-posterior.md). There is no
CLI/C ABI operation. A supplied proper prior, independent noise declarations,
fixed design and exact coordinate order define this calculation; they do not
qualify observed calibration uncertainty or a released finite-box posterior.

For beta~N(m,S), r=X beta+epsilon and y*=A beta+epsilon*, the original prior
and training noise are independent. Future noise epsilon*~N(0,R*) is explicitly
independent of both that prior and training noise. C, S and R* are full SPD
matrices; future noise may be correlated across future rows. The posterior
supplies beta|r~N(mu,V), giving

    b = A mu,
    W = R* + A V A^T,
    log p(y*|r) = -1/2 [(y*-b)^T W^-1 (y*-b)
                         + log det W + k log(2 pi)].

The density uses the declared product of future-coordinate Lebesgue measures.
For every nonzero v, v^T W v = v^T R* v + (A^T v)^T V(A^T v) > 0.
An invertible Cholesky change of variables therefore gives integral one over
the entire future coordinate space. Numerical admission covers resolved
inputs of this model; refusal does not mean zero probability or define
truncated support. No box normalization or finite-support prior is supplied.

Shared beta creates training/future and cross-future dependence. Distinct
events alone do not establish noise independence. Training/future cross-noise
covariance, singular future noise, selection/censoring, nonlinear response,
improper or zero-width priors and measured-data qualification are unsupported.
Never include the same beta/calibration prior both in the response and an
already-marginalized future covariance. The original-prior correlated
calibration operator and this conditional predictive law have separate roles.

## Preparation, identities and ownership

`GaussianPredictive::prepare` borrows a `const GaussianPosterior&`, its exact
ordered training vector/row IDs, a `const Gaussian&` representing future noise,
row-major A and `const PredictiveMetadata&`. Metadata copying occurs after
resource admission inside the allocation exception handler. Every borrowed
input remains usable on success and failure. Training is conditioned once;
W is assembled and factored once. No source factor or complete source
covariance is copied into the predictive owner.

Training row order must exactly match the posterior source. A columns and
their parameter units must exactly match the prior's ordered parameter IDs
and units. Future axes come from the supplied noise Gaussian. Measurement
IDs must be disjoint across training/future sets; parallel event-ID vectors
must be nonempty, unique within each set and disjoint across sets. Event IDs
and measurement IDs are different namespaces. These are checked caller
declarations, not independently verified astronomical event linkage.

This bounded consumer requires `source_semantics="synthetic controls"` on
both sources. Future noise must be unshifted and have no applied latent prior.
Its arithmetic must match the posterior source. The caller supplies a
nonempty conditional-noise identity and explicit declarations that future
noise is conditional on beta and independent of the original beta prior and
training noise. The returned law is a derived conditional prediction, never
an additional observation or a new original generative prior.

One future unit is supported. It must match the prior's residual unit;
`future_covariance_unit` is exactly `future_unit+"^2"`, `future_measure` is
`"product d("+future_unit+")"`, matching the future Gaussian's measure.
Each response column declares `future_unit+"/"+parameter_unit`. Unit labels
perform no conversion. Changed coordinates require explicit transformations
of data, responses, noise and parameter moments; normalized densities acquire
the corresponding Jacobian.

The move-only owner retains b, W and its normalized-density factor, their
numerical diagnostics, A, exact training values, ordered row/event and parameter
identities, original source metadata and prior/design/dependence identities.
It retains neither original noise factor nor a cloned posterior. After
preparation the inputs may be destroyed. Const mean/covariance and diagnostic
spans borrow the predictive owner and are invalidated by moves/destruction.
Moves explicitly invalidate the old owner and clear its retained views;
self move retains it. Empty/failed owners and move cleanup allocate nothing.

`log_density` evaluates one ordered whole future vector with retained b/W.
It neither reconditions training nor reassembles, copies or refactors
covariance. It reads no original training or future-noise covariance; retained
W is used by the solve and its residual diagnostics. It does not multiply
independent per-row likelihoods. Changing r, A or
R* requires another preparation. Status must be finite before interpreting
numeric result fields; failed preparation returns no usable mean/covariance,
and failed density evaluation supplies no usable likelihood.

## Numerical and resource contract

At least one training and future coordinate, and at least two parameter
coordinates, are admitted, matching the existing parameter-posterior contract.
Finite normal binary64 inputs are required; zero is allowed. Wide products
require round-to-nearest and long double with at least 64 mantissa bits and
maximum exponent at least 16384. Arithmetic follows the existing source/noise
Gaussian identity without fallback. Underflow, overflow, nonfinite values,
unsupported arithmetic, unresolved conditioning or quotas cause refusal.
There is no jitter, clipping, dropped coordinate or numerical noise added to
any statistical covariance.

Mean error propagates each posterior mean error through absolute response
weights, wide summation and actual binary64 storage rounding. B=A V and
W=R*+B A^T are accumulated in wide arithmetic. Each symmetric entry is
computed once and mirrored exactly. Inherited posterior covariance error,
absolute product/accumulation diagnostics and actual storage rounding produce
separate covariance-entry numerical errors. These are distinct from W's
physical predictive variance and from uncertainty in externally supplied
covariance assets.

Let eta be the retained W inverse-norm estimate times its numerical covariance
error's infinity norm. Preparation refuses eta>=0.01 or relative estimate
eta/(1-eta) above the policy. Density admission adds mean/subtraction error
projected through that inverse norm, quadratic sensitivity q*eta/(1-eta)/2,
log-determinant sensitivity -k*log1p(-eta)/2 and retained Gaussian solve/cast
diagnostics. The scaled log-density estimate must meet the requested limit.
Positive unresolved diagnostics cannot silently round to zero. These are
empirical arithmetic/conditioning estimates, not rigorous universal bounds.

Default limits are one million matrix elements, 256 MiB payload,
100 million declared work units and 1e-10 maximum forward sensitivity.
Preparation checks the conservative dimension count
32*(n*n+n*p+p*p+k*p*p+k*k*p+8*k*k*k); evaluation checks 32*k*k.
These counts define bounded work scope, not measured elapsed time.
Evaluation uses the tighter retained/caller ceilings.

Preparation's simultaneous payload envelope includes both already-retained
borrowed owners once, candidate storage/factors, conditioning and assembly
scratch, numeric outputs and copied metadata. Callers must not charge those
owners twice. Retained and evaluation bounds use actual vector/string
capacities; evaluation reports scratch separately and admission sums retained
plus scratch. Borrowed arrays, allocator metadata, stack and RSS are excluded.
All checked products/sums precede allocations. Internal allocation failure
returns typed work-limit status; caller construction of inputs is outside the
operation. Original inputs and an already-prepared predictive owner remain
usable after failed preparation/evaluation.

## Named synthetic acceptance

The original two-parameter offset/slope control has two training rows and two
distinct future rows in mag. Its proper prior and full correlated future noise
are declared synthetic assumptions. Exact rational mean/W/determinant and
normalized-density controls use 2e-12*(1+absolute reference) for moments and
2e-11*(1+absolute reference) for log density, quadratic and log determinant.
Zero response gives the independent future-noise law; null/rank-deficient
training remains lawful under the proper prior. Shared responses retain
cross-future covariance, even when future noise is diagonal.

An independently authored joint-covariance/cofactor route and original direct
prior-times-training-times-future quadrature test this consumer. Native future
density integration uses the same retained owner, an independently integrated
future coordinate measure and explicit omitted-tail accounting. The integral
allocation is 2e-11 absolute; reference refinement consumes at most 5% of each
named allocation. Shared Gaussian kernels are ancestry, not independent
evidence. Permanent controls also challenge units/axes, duplicate events,
prior-bearing noise, full covariance, ownership/lifetime, arithmetic refusal,
all measured allocation sites and payload/work limits. These controls do not
qualify arbitrary supplied inputs, a physical prior, observational held-out
calibration, H0 posterior/coverage or the released ladder box target.
