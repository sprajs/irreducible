# Proper Gaussian parameter posterior

The native `irred/gaussian_posterior.hpp` calculation conditions an explicit
proper Gaussian parameter prior on a fixed linear observation model. It is
available through the C++ SDK. The strict-inline `statistics.gaussian_posterior`
CLI/C ABI consumer returns a batch of means and one invariant covariance;
parameter-density evaluation remains SDK-only. For explicitly
independent beta~N(m,S) and epsilon~N(0,C), r=X beta+epsilon gives

    P = S^-1 + X^T C^-1 X,
    V = P^-1,
    mu = m + V X^T C^-1 (r-Xm).

The returned covariance V and conditional mean mu define a normalized Gaussian
parameter density in the caller's declared product of parameter-coordinate
Lebesgue measures. `log_density` includes the determinant and p log(2 pi).
This is a parameter posterior under the supplied proper prior. It is distinct
from the relative design profile, the conditional estimator sampling law and
the observed-residual density obtained by marginalizing a calibration prior.
It does not fit a nonlinear cosmological model or reproduce a released posterior.

The source Gaussian owns the supplied observation covariance and its retained
factor. `ParameterPrior` specifies ordered unique parameter IDs, units, shared
nuisance IDs (a unique subset), mean, full SPD covariance S, prior/design IDs,
residual unit, parameter measure and an explicit noise/prior independence
declaration. Unit/measure labels are provenance, not physical conversions.
At least two parameters are supported. Observation/design axes must match
exactly. Prior-bearing or mean-shifted source Gaussians are rejected, preventing
an accidental second application of the same prior. A null or rank-deficient
design is permitted because S is proper; no jitter or implicit regularization
is introduced. Zero-width box support and improper flat priors are unsupported.

Preparation whitens design columns with the retained C factor and identity
columns with a temporary S factor. Wide Gram sums form P, with diagnostics for
whitening, accumulation and final binary64 storage. The existing precision-input
Gaussian preparation produces V and its retained normalized-density factor.
The owner retains the design, whitened design, precision and prior. All candidate
allocations and checks finish before it consumes the source; failed preparation
leaves that source usable. Moves invalidate the old owner and self move preserves
it. Evaluation uses retained factors and matrix-vector operations without
refactorization, copying the source covariance or scanning it again.

Finite normal binary64 inputs are required; zero is allowed. Arithmetic requires
round-to-nearest and long double with at least 64 mantissa bits/exponent range
16384. Conditioning, underflow/overflow or unrepresentable positive diagnostics
cause refusal. Element/payload limits apply before allocations. The preparation
work admission uses the conservative dimension count
32*(n*n*p+n*p*p+8*p*p*p), with a default 100-million ceiling; this bounds this
implementation's declared work scope and is not a measured wall-time guarantee.
The payload estimates cover owned/transient matrix, vector and metadata storage,
excluding borrowed input, allocator overhead, stack and RSS.

The posterior covariance's relative numerical error estimate propagates the
precision diagnostic through its inverse norm. Mean outputs carry numerical
absolute-error estimates separately from physical posterior variance. Density
admission propagates the mean/subtraction estimates through P and adds covariance
quadratic/log-determinant sensitivity. Default maximum scaled forward sensitivity
is 1e-10. These are empirical arithmetic/conditioning estimates, not rigorous
continuous or universal error certificates. Numerical noise is not added to S
or C. A tighter unattainable request fails without changing its allocation.

Permanent owner controls use exact rational two-dimensional cofactor references
derived before implementation, plus a separately composed Bayes density identity,
null/rank-deficient designs, parameter-unit Jacobian, exact axis order, proper
prior rejection, quotas and owner moves. The fixed named mean/covariance and
density comparison allocation is 2e-12*(1+absolute reference). Shared factor and
Gaussian kernels in the Bayes comparison are ancestry; they are not independent
evidence. A separately derived covariance-conditioning/cofactor reference and
direct joint prior-times-likelihood quadrature now test normalization, moments,
rank-deficient/null and underdetermined designs, axis permutations and unit
Jacobians. Named mean/covariance controls use 2e-12*(1+absolute reference),
log density uses 2e-11*(1+absolute reference); reference refinement consumes at
most 5% of those allocations. Two initial reference failures remain preserved.
These controls do not qualify arbitrary inputs or physical priors. No observational calibration uncertainty,
confidence coverage or systematic-error claim follows from these controls.

## Strict-inline family of conditionals

Run the original two-parameter synthetic example with:

```sh
irred run tests/fixtures/gaussian-posterior.json STORE
```

The request supplies full row-major noise covariance C, fixed design X, proper
prior mean m and covariance S, then complete ordered training vectors. One
coarse native call prepares the source and posterior once and conditions every
vector using those retained factors. The result owns one posterior; its common
V view borrows that owner and each mean is moved into the pooled result. No
source covariance readback, per-vector FFI, reapplication of the prior or new
factorization is needed between means. The standalone C++ API remains usable.

This is a family of conditionals under one fixed model. It does not declare a
joint posterior across cases, independent generated draws or a sampling law.
The source role is restricted to synthetic controls. Noise/design/conditioning
row order and parallel event order must match exactly; design/prior parameter
IDs and units match exactly, and shared nuisance IDs name one subset of those
parameters. Event labels may repeat to declare shared events; neither event
labels nor distinct case IDs prove independence. Full C retains the declared
cross-row noise. C entries have residual-unit squared dimensions, S entries
have the product of their two parameter units, and each design column declares
residual-unit/parameter-unit dimensions. Labels perform no conversions.

The fixed required outputs are `posterior_covariance` and `posterior_means`.
V carries the native relative numerical-error estimate; each mean carries its
native absolute-error vector, backward residual and scaled sensitivity. These
empirical numerical estimates remain separate from physical posterior variance.
A case refusal has no numeric mean/error payload and remains in the original
order; successful neighbors and V survive. Any refused case fails the required
means check. Global preparation/resource refusal withholds both moments. No
jitter, dropped case, weaker sensitivity or changed prior repairs a refusal.
The operation supplies output checks, method and arithmetic to the recorder.
The immutable raw and resolved request preserve all original vectors and prior
facts; numerical assurance does not qualify their physical interpretation.

The original example uses the named rational two-dimensional control: its first
mean is (20604/25511, -7125/51022) and V=(15160,4466;4466,8592)/25511. The second
vector is exactly X*m in ideal arithmetic, giving mean m. These are inherited
control facts; interface parity is shared native ancestry, not new independent
scientific evidence. No density, evidence, nonlinear fit, released-box target or
fixed-truth estimator sampling output is supplied by this CLI operation.

All nested objects reject unknown/duplicate JSON members. Input is bounded to
16 MiB, numerical pools to one million elements, cases to 65536, total metadata
to 1 MiB and each text field to 256 UTF-8 bytes. Requested native payload is at
most 256 MiB and declared work at most 100 million units. Smaller original
ceilings remain binding. For n rows, p parameters and B cases, checked admission
charges

    32*(n^3+n^2*p+n*p^2+8*p^3+B*(n^2+n*p+p^2))

before allocation/evaluation, including cases that later refuse. This is a
conservative declared work scope, not measured FLOPs or elapsed time. The
pooled input/output element count is n*n+n*p+p+2*p*p+B*n+2*B*p. Native phase
payload includes the result owner, conversion metadata, source/posterior
preparation, retained factors once, pooled mean/error/case capacities and one
serial conditioning scratch space. Borrowed Rust input/wire descriptors, JSON
serialization, stack, allocator overhead and RSS are excluded.

An empty failure result has a disengaged optional posterior and empty vectors;
its constructor allocates only the charged result object. If the requested cap
cannot hold even that object, the new evaluator alone returns scoped transport
`GAUSSIAN_POSTERIOR_QUOTA_REFUSED` (6) with a null owner before allocation. Rust
records a completed resource refusal with no native payload, both moments
withheld, requested method/arithmetic and no actual execution claim. At an
admissible cap an owned failure envelope remains charged; global cleanup
releases all scientific/vector ownership. Failed native error-vector capacities
remain charged internally even though failed views expose no numeric payload.

Owned refusals retain native preparation attempt/completion flags and the
number of conditioning calls attempted. The reported method identifies the
attempted stage; a refused noise preparation does not establish a completed
factor or actual arithmetic evaluation. Its returned arithmetic contract is
labelled as attempted, while completed source/posterior preparation is recorded
separately. Admission-only refusals keep actual execution fields absent.

The predictive CLI reuses these bounded inline training input types, source and
proper-posterior owners. Its additional fixed future law and whole-vector density
are described in [Gaussian prediction](gaussian-predictive.md); neither route
reapplies a marginalized prior or changes this parameter-posterior calculation.
