# Synthetic proper-prior ladder and joint held-out prediction

The native `irred/calibration_predictive.hpp` consumer composes the compiled
[ladder equations](calibration-ladder.md) with the retained
[Gaussian posterior](gaussian-posterior.md) and
[joint Gaussian predictive](gaussian-predictive.md) owners. There is no CLI,
Rust or C ABI operation. This is a bounded synthetic supplied-shape calculation,
not a released-data ladder, observational H0 inference or a finite-box posterior.

For the exact existing ordered host moduli, M_Cep, b, gamma, M_SN, eta and delta,

    y = o + X beta + epsilon,
    beta ~ N(m,S), epsilon ~ N(0,C),
    y* = o* + A beta + epsilon*, epsilon* ~ N(0,R*),
    V = (S^-1 + X^T C^-1 X)^-1,
    mu = m + V X^T C^-1 (y-o-Xm),
    E[y*|y] = o* + A mu, W = R* + A V A^T.

S, C and R* must be proper full SPD matrices. The original parameter prior,
training noise and future noise are explicitly independent; future noise is
conditional on beta. Internal correlations in each matrix remain retained.
Parameter density uses the prior's declared product coordinate measure; future
density uses `product d(mag)`. Translation by an offset has unit Jacobian.
No calibration prior is applied twice: delta is one shared parameter, and
conditional C/R* exclude a contribution already marginalized over that delta.
A calibration measurement remains a data row, separate from the proper prior.

`linearize(Model, DesignPolicy)` obtains design, offsets and metadata from the
existing single compiled row equation. Its model validation is prediction-valid:
proper prior conditioning may use rank-deficient or missing-rung training rows.
The relative `Ladder` preserves its two anchored hosts/rung and full-rank gates.
The host dictionary still has at least two unique ordered entries. Future rows
share the exact ordered host dictionary, magnitude/metallicity conventions,
reference H0, supplied distance shape and calibration/dependence identities.
Future conditional covariance may have its own explicit identity.

Every original row/host/event ID remains in the retained model. Gaussian
predictive lineage is exactly `event:` followed by a Cepheid/SN event ID, or
`measurement:` followed by the anchor/calibration row ID. Empty astronomical
anchor event IDs do not become empty predictive lineage. All row IDs and actual
source-event IDs must be unique and disjoint between training/future sets.
Shared hosts and equal responses are allowed; different events alone do not
prove conditional-noise independence. Caller PredictiveMetadata must match the
compiled parameter axes/units, tagged event vectors, design/dependence identity
and conditional noise identity exactly; it also declares independent future
noise and a nonempty conditioning identity. Observation source metadata must
match synthetic-control role, exact row order, mag measure and declared
calibration/dependence/conditional-covariance provenance. Labels perform no
physical conversions and do not verify astronomical provenance.

`LadderPosterior::prepare` validates and acquires training covariance once,
consuming its Gaussian only after every candidate allocation/check succeeds.
Failure preserves the source. The `posterior()` getter requires finite status;
an empty wrapper retains no parent allocation. Static bounds reject an invalid
borrowed posterior. The retained posterior exposes normalized parameter density,
conditional mean/covariance and original source/prior metadata.
`LadderPredictive::prepare` borrows that owner, the exact ordered training vector,
future conditional-noise owner, future model and metadata. Borrowed owners remain
usable on success/failure. Future covariance is assembled and factored once.
The resulting owner retains model/design/offset/lineage, mean/W, numerical errors
and the parent metadata, without cloning the posterior or source factors.
Changing training values, future response or covariance requires preparation.
Evaluation performs whole-vector solves without copying/refactoring covariance.
Views borrow their owner and are invalidated by moves/destruction; old owners
are invalid after moves, and self move preserves the owner.

## Offset domain, numerical and resource admission

Observation values are supplied in their original magnitude coordinates.
Offset subtraction must be exactly representable as a finite normal binary64
value (zero allowed). Wide subtraction and an independent error-free TwoSum
remainder check reject a rounded translation. This deliberately excludes such
evaluations rather than claiming propagation of an omitted offset error;
original data must not be changed to pass. An underflowed, subnormal or overflowing
translation fails with no usable result. Future mean offset addition propagates
parent mean error, wide accumulation and actual storage rounding into separate
absolute error estimates and checks the requested sensitivity. Those errors are
numerical diagnostics, not physical predictive variance or covariance jitter.

Existing wide arithmetic and parent sensitivity contracts apply unchanged.
Preparation reserves 64*n*(p+1) wrapper work units before parent preparation;
prediction additionally reserves 32 times the training row count. Evaluation
reserves its translated row count before the parent solve. H0 projection
additionally reserves p*p+32 units for its full covariance norm and projection. Parent work admission
uses the remaining ceiling, so wrapper work does not bypass the combined limit.
The default wrapper policies retain one million elements, 256 MiB payload,
100 million declared combined work units and 1e-10 forward sensitivity. Static
preparation bounds conservatively include model/design/offset/lineage construction,
parent preparation storage, retained borrowed owners once, and simultaneous
centering/mean scratch. Retained bounds use actual vector/string capacities;
evaluation bounds report scratch, including the translated whole vector. Admission
adds wrapper retained storage to parent/evaluation scratch and uses tighter
retained/caller limits. Payload excludes allocator bookkeeping, borrowed arrays,
stack and RSS; work counts are declared scope, not measured elapsed-time bounds.
Allocation failures return typed refusal; caller argument construction lies outside
the operation. For exact byte-boundary checks, compute the preparation bound on
the constructed model/prior/metadata objects and move those exact objects into
preparation: a by-value copy may have different string/vector capacities.
No jitter, dropped rows, changed budgets or covariance repair occurs.

## H0 posterior projection

`h0_projection` returns the posterior eta mean/variance and the derived lognormal
median, expectation and log-standard-deviation under the shared coordinate
H0 = Href exp(a eta), a = ln(10)/5:

    median = Href exp(a mu_eta),
    expectation = Href exp(a mu_eta + a^2 V_eta_eta/2),
    log standard deviation = a sqrt(V_eta_eta).

A nonpositive, nonfinite, subnormal or unrepresentable output refuses the entire
projection. The result retains the eta-mean and eta-variance absolute arithmetic
errors,
and separate relative errors for median, expectation and log-standard-deviation.
Absolute variance error uses the full posterior covariance infinity norm times
its inherited relative error; an entrywise relative error is not assumed.
Diagnostic casts round upward and each requested projection must meet the
sensitivity ceiling, separately from physical variance. This uses the proper
prior; it is distinct
from relative recovery and the explicitly centered estimator sampling law.
It supplies no arbitrary quantile inversion or observational H0 result.

## Named checks and coverage interpretation

The owner control has two hosts/eight parameters, one anchor training row and a
proper diagonal Gaussian prior. An independently derived scalar conjugate update
fixes posterior host mean/variance. Two future Cepheid/Hubble rows with correlated
noise retain the shared delta covariance and a nonzero magnitude offset.
Named mean/covariance controls use 2e-12*(1+absolute reference); normalized density
parts use 2e-11*(1+absolute reference). Relative fit rank refusal is not a
proper-prior refusal. Axis/lineage mismatch, rounded translation, byte boundaries,
move/lifetime and every measured preparation allocation site have distinct checks.

Fixed-truth recovery differs from prior-predictive calibration. With
K = V X^T C^-1 and fixed beta0, conditional posterior-mean bias is
(I-KX)(m-beta0) and its sampling covariance is K C K^T. A posterior credible
interval therefore need not have nominal frequentist coverage at that truth.
H0 interval coverage equals eta interval coverage because its projection is
monotonic. Fixed-truth future error has bias A(I-KX)(beta0-m) and covariance
R* + A K C K^T A^T, generally different from W. Under the separately declared
prior-predictive ensemble, the conditional future law instead has covariance W;
its whitened joint ellipsoid has the exact corresponding chi-square probability.
Analytic Gaussian controls precede stochastic coverage; a seeded empirical
receipt cannot replace those distinctions or qualify observational coverage.

The independent peer freezes a separate eleven-row, two-host/eight-parameter
rational fixture with correlated full prior/noise and nonzero Hubble offsets.
Covariance conditioning and independently formed prior/noise precision agree
exactly in Fraction arithmetic before native evaluation. A prior-whitened
two-coordinate quadrature integrates prior × training × future noise directly;
16/24-order composite refinement is below five percent of the density allocation.
A noiseless fixed-truth response checks posterior bias and the resulting eta/H0
credible-interval coverage formula. Numerical integration of the conditional
two-row joint ellipsoid agrees with the chi-square2 mass 1-exp(-radius²/2),
with separate refinement. These are analytic synthetic law checks; no empirical
coverage campaign or observational held-out claim is made.
