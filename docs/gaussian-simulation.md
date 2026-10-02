# Addressed synthetic Gaussian generation and recovery

The native `irred/gaussian_simulation.hpp` batch consumes a retained covariance
factor through a borrowed, unshifted synthetic-control `Gaussian`. A supplied
ordered generating mean defines the ideal law `x ~ N(mean,C)`. There is no CLI
or C ABI operation. Mean axes and coordinate measure must match exactly; units
are declared labels and perform no conversion. Prior-bearing or shifted source
Gaussians are refused.

The emitted law is a discrete pseudo-Gaussian approximation: shared
Philox4x32-10 words, halfbin uniforms `(word+0.5)/2^32`, the original detector's
words1/2 cosine Box–Muller expression, retained triangular colouring and final
binary64 rounding. Finite normal tails, represented-factor error, platform
transcendental arithmetic and observation rounding are separate from the ideal
Gaussian law. The continuous posterior/predictive density is not the exact mass
function of those emitted vectors. Integer words replay exactly; transcendental
values require the same supported arithmetic/libm environment. Counter addresses
and empirical cross moments do not prove statistical independence.

Vector address `(stream_base,sample)` assigns coordinate `j` to
`(stream_base+j,sample)`. Stream overflow and overlapping coordinate ranges for
the same sample fail admission. Across-call replay is allowed; a consumer owns
the global seed/address ledger and distinct prior/training/future role ranges.
All requested vectors retain status and address, including refused vectors.
Pooled vector-major values, arithmetic errors and random words have an explicit
output mask; numeric slots of failed vectors are unusable. A whole-batch admission
failure retains seed/request count without allocating per-vector successes.

`numerics::colour` multiplies the retained Cholesky triangle with wide products
and sums. It copies no factor, reads no original covariance and refactors nothing.
Its error estimates concern accumulation against that triangle, not a certificate
for the supplied statistical covariance or a universal libm error bound. The
sampler adds mean/storage diagnostics and upward-casts positive estimates. It
requires round-to-nearest and long double with at least 64 mantissa bits and
exponent range 16384. Nonfinite/subnormal values, unresolved positive diagnostics,
overflow, unsupported arithmetic and unmet allocations cause refusal.

Default admission bounds are 65536 vectors, one million numeric elements,
256 MiB simultaneous payload, 100 million declared work units and `1e-10`
scaled arithmetic sensitivity. `work_bound` charges `32*N*d*d+256*N*d+128*N`;
these are conservative operation-scope counts, not wall-time guarantees.
Payload admission includes the borrowed Gaussian once, copied mean/metadata with conservative fresh-string capacity envelopes,
requested pooled outputs, row/address-sort storage and simultaneous wide colour
scratch. Borrowed input arrays, allocator overhead, stack and RSS are excluded.
Allocation failures return typed refusal and preserve the borrowed Gaussian.

## Named empirical controls

`test_gaussian_recovery` and `test_calibration_recovery` are explicit native
campaign executables, registered in ordinary CTest and required native CI. Each main cell has 32768 attempted
realizations, split between seeds `0x123456789abcdef0` and
`0xfedcba9876543210`, in chunks of 1024. Full prior/training/future covariances
and ordered disjoint row/event lineage are retained. A shared parameter draw
generates both training and future responses in the prior-predictive ensemble;
the fixed-truth ensemble instead holds its declared parameter vector constant.

The generic offset/slope fixture reuses the independent predictive peer's
two-dimensional correlated law. The ladder uses a separately identified narrow,
fully correlated dyadic law with eleven training rows, eight parameters and two
future rows. Its prior mean is `(30,31,-4,-3,1/4,-19,1/4,1/8)`; fixed truth
changes eta to `1/2`. Prior/noise factors have diagonals `1/16` and `1/8`,
respectively, with alternating first-column entries `±1/128`. Future noise factor
is `(1/8,0;1/32,1/8)`. This target is distinct from the broad rational peer
fixture, whose unchanged inputs are retained as a refusal challenge.

The finite correlated-colour extrema keep the narrow law's nonzero-offset
magnitudes inside `[16,32]`. Integer offsets 40/41/42 then subtract exactly in
binary64. The ideal Gaussian has unbounded tails and a different translation
domain; the finite-law proof does not qualify its entire support. Observations
are never snapped, clipped or regenerated to pass the exact translation gate.

Independent Fraction covariance/precision routes fix `K,V,W`, fixed-truth mean
bias and sampling covariances. Independently refined CDF facts fix coordinate
coverage. Fixed-truth posterior-mean error has bias `(I-KX)(m-beta0)` and covariance
`K C K^T`; future error has bias `A(I-KX)(beta0-m)` and covariance
`R+A K C K^T A^T`. These generally differ from the proper posterior/predictive
covariances. In the shifted narrow ladder, ideal eta/H0 95% credible-interval
coverage is approximately `0.0004294447488094294`, rather than 95%. Monotonic H0
projection preserves the eta interval event. Only prior-predictive joint
ellipsoids use a central chi-square reference; fixed-truth joint errors are
biased and anisotropic.

Empirical mean/full covariance criteria use six ideal standard errors plus
separate `1e-6*(1+abs(reference))` discretization and `1e-10` arithmetic
allocations. Coverage uses `6*sqrt(p*(1-p)/N)+6/N`, plus frozen `5e-4`
discretization and `1e-10` arithmetic allocations. An independently derived
halfbin coupling allocation is below `5e-4` for these named fixture dimensions:
it excludes radial-uniform endpoints with alpha `1e-6`, charges their union mass,
and bounds ideal whitened boundary shells. This concerns an independent-uniform
mathematical comparator; it is not a proof of actual Philox independence or
platform libm accuracy. These are named seeded discrepancy criteria, not
statistical coverage guarantees.

The campaign ledger records actual seeds, role stream ranges, ordered generating
identities and each successful attempt with its original parameter/training/future
vectors, means, quadratics and interval predicates. Admission refusals record
planned and actually attempted counts before stopping. Failed attempts retain both
observation vectors, stage and concrete numerical/density statuses; nonfinite
values retain IEEE binary64 bits. Every refusal is retained in the native receipt. Coverage bounds use the
attempted denominator: `[hits/N,(hits+refusals)/N]`. Any main-cell refusal withholds
empirical acceptance; the broad refusal probe makes no coverage claim. Numerical
acceptance and synthetic agreement establish neither observed calibration
uncertainty nor measured-data qualification. Existing predictive preparation is
charged for every changed training vector; no new reconditioning operation is
introduced. Receipts report generation work and conditional preparation/density
time separately, and total campaign time includes stage and attempt logging.
The matched coarse/one-vector generation comparison excludes equality verification
from both timed arms and separately charges retained comparison outputs; it is
not a consolidated campaign throughput benchmark.
