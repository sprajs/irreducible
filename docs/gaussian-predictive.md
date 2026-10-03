# Retained joint Gaussian prediction

The native C++ `irred/gaussian_predictive.hpp` consumer predicts distinct
synthetic calibration observations after conditioning the existing
[proper Gaussian parameter posterior](gaussian-posterior.md). The strict inline
`statistics.gaussian_predictive` CLI operation and coarse ABI2 consumer use the
same native repeated-conditioning owner. A supplied proper prior, independent noise declarations,
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

## Repeated-conditioning contract

For unchanged C, X, m, S, A and R, repeated training vectors obey
`mu(r)=m+V X^T C^-1(r-Xm)` and `b(r)=A mu(r)`. The normalized future
covariance `W=R+A V A^T` is independent of r. `GaussianPredictiveConditioning::prepare` consumes a posterior only after
successful invariant preparation. `evaluate` borrows contiguous row-major
training/future pools, exact coordinate IDs and a positive count. It owns
the original posterior, response, source/parameter/event lineage and one W
factor; no posterior or input noise factor is cloned. Ownership transfers only
after successful preparation. Original vectors belong to the caller and each
coarse batch preserves their order and individual refusals. A changed design,
prior, response or noise requires another owner.

`PredictiveOutputs` selects future means (with their absolute numerical errors),
whole-vector log densities, or both. Density-only evaluation still computes the
conditional mean internally but does not retain it in output. Per-vector
refusal withholds its numeric outputs; global admission failure returns no rows.
Finite normal/zero inputs, arithmetic and sensitivity allocations match the
one-conditioned-law consumer. Named comparisons retain 2e-12*(1+abs(reference))
for moments and 2e-11*(1+abs(reference)) for density parts; existing independently
derived rational/Decimal references remain the scientific ancestry.

Batch admission charges every requested row before allocations, including rows
that subsequently refuse. Its conservative per-row scope is
`32*(n*n+n*p+p*p+k*p)`, plus `32*k*k` for requested densities.
Exact training translation and original-coordinate mean translation are each
charged once; requested densities additionally charge future translation.
The peak payload includes retained owners once, pooled requested outputs and
one row of condition/projection/solve scratch. Borrowed original arrays,
allocator bookkeeping, stack and RSS are excluded. Neither work counts nor
matched timers establish a universal performance or scientific qualification.

`GaussianPredictiveConditioning::posterior()` requires finite owner status.
The posterior remains available to other const library consumers; the batch
owner does not return a cloned parameter law. Ordinary owner controls observe
one setup Cholesky call, then zero Cholesky calls with three training whitenings
and three density solves for a three-admitted/one-refused batch. 64-bit Linux GCC/Clang
ELF forwarding hooks observe the original routines without replacing their
arithmetic. Portable named mathematical checks do not depend on those hooks.

## Strict inline predictive CLI

```sh
irred run tests/fixtures/gaussian-predictive.json STORE
```

Schema version 2 supplies the original fixed training noise/design/proper prior,
complete pooled training vectors, fixed future noise/response, and explicit
prediction declarations. `source_semantics` is exactly `synthetic_controls`.
Every request section is an object with no unknown or duplicate members;
positional object arrays, null sections and nonfinite JSON numbers are rejected.
Training/future vectors remain in their original order. Cases are alternative
conditionals under the same fixed model; their labels do not declare a joint
law across cases or independent samples.

`outputs.means` requests `predictive_means`, including absolute numerical errors.
`outputs.joint_log_densities` requests `joint_predictive_log_densities`, one
normalized density for each whole ordered future vector. At least one output is
required. `future_vectors` must be present as an object exactly when density is
requested, and absent otherwise. The conditional mean and joint covariance are
still internal dependencies of density; the CLI exports no covariance getter.
Each requested group has a required numerical check. Omitted groups are
`not_requested` and carry no numerical payload or receipt check.

One native call prepares C, the proper posterior, R and W once and evaluates
both borrowed contiguous vector pools in a coarse batch. The result moves the
native output pools and exposes borrowed views until serialization finishes.
No per-case FFI, source covariance readback or repeated factorization is needed.
The original prior is applied once. Training/future measurement IDs and parallel
event IDs are unique and disjoint in their respective namespaces; all response,
parameter, unit and vector orders match exactly. Unit labels perform no
conversion. Supplied calibration identities state that the same beta prior is
excluded from both conditional noise covariances.

Global preparation/admission failure publishes no rows. Each admitted case
retains its original index and typed refusal. When both groups are requested,
a case whose joint density fails withholds its mean and density atomically;
successful neighboring cases remain. Refusal is never zero probability and no
case is dropped or renormalized. Phase flags record attempted and completed
preparations separately, including an empty-batch refusal before all native
preparations. Method and arithmetic records identify the attempted algorithm;
unexecuted phases have no actual method/arithmetic claim.

The hard request ceiling is 16 MiB, one million aggregate numeric fields
including requested outputs, 65,536 cases, 1 MiB aggregate UTF-8 text, and
1–256 bytes per text field. Native limits are 256 MiB payload, 100 million
declared work units, one thread, and forward sensitivity at most 1e-10.
`resource_policy` may request smaller limits, including zero resource ceilings;
they remain original ceilings and produce completed resource refusals. A native
byte ceiling smaller than the allocation-free result envelope uses the scoped
ABI `IRRED_GAUSSIAN_PREDICTIVE_QUOTA_REFUSED` disposition with no result owner.
Invalid wire/version/alignment/length/UTF-8 and allocation transport failures
remain distinct execution failures.

Before metadata or scientific allocations, checked whole-request reservation is

    WC = 32 n^3,
    WP = 32 (n^2 p + n p^2 + 8 p^3),
    WR = 32 k^3,
    WA = 32 (n^2 + n p + p^2 + k p^2 + k^2 p + 8 k^3),
    WB = 32 B (n^2 + n p + p^2 + k p + D k^2).

Here n, p, k and B are training, parameter, future and case counts; D is one
exactly when density is requested. The sum charges every requested case,
including refused cases. These reservations bound declared work, rather than
reporting measured kernel calls. The stored predictive preparation policy
retains the original work allowance minus WC/WP/WR; evaluation receives that
allowance minus WA, preserving legal larger batches. Payload admission takes
the maximum simultaneous conversion/preparation/evaluation lifetime, retaining
borrowed native owners once. Records include minimum result, peak and retained
payload and all five work components. Borrowed Rust request/wire buffers, JSON
serialization, stack, allocator bookkeeping and RSS are outside native payload.

The named rational/cofactor fixture uses the native guide's unchanged moment,
density and normalization budgets. Independent direct-product quadrature,
refinement, tail and arithmetic uncertainty must fit inside the same total
comparison allocation; their reference share is at most 5%, not an added
allowance. Native/ABI/CLI and installed-consumer parity share numerical ancestry
and qualify transport and ownership only. Allocation-prefix and no-elision
controls test the selected runtime's actual payload behavior. Generic synthetic
responses do not qualify physical future designs, cross-noise calibration,
selection, a measured population, parameter-posterior coverage or observational
inference. The proper-prior ladder remains a separate native-only consumer.
