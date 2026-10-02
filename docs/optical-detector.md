# Supplied joint optical law through detector censoring

This bounded native C++20 composition answers one synthetic question: given a
finite joint optical-transmission law and fixed sampled spectrum, what joint
count or censored-record law follows after marginalizing its shared state?
The first consumer is a two-state/two-band SDK experiment. State masses and
instrument/source declarations are supplied synthetic controls; numerical
acceptance does not qualify a measured passband, camera or source population.
No CLI, C/Rust ABI, temporal source model or population inference is added.

## Frozen measure and physical assumptions

Reuse [sampled photometry](photometry.md) and its [finite calibration law](photometry-calibration.md)
for each ordered state/band expected transmitted photon signal N_sb. Relative
masses are positive finite normal binary64, with the calibration owner's single
wide normalization p_s=w_s/sum(w). Optical transmission precedes the fixed
wavelength-independent quantum efficiency. Background, dark current, exposure,
read noise, gain and bias use [the detector owner's law](detector-selection.md).
The calibration band's observer exposure is also the dark-current exposure;
there is no second exposure or supplied photon expectation that could disagree.

Conditional on one optical state, all declared band counts are independent
Poisson arrivals and their Gaussian read noises are independent of those counts
and each other. This assumption must be explicitly selected and its origin
supplied. It does not follow from different band labels or from a shared state.
Different joint records are evaluated separately; they are not multiplied into
a cross-record likelihood. A temporal/coherent calibration state across multiple
records would require a separately declared joint-record law.

A record contains exactly one ordered detector observation per band. Detection
means Y_b >= threshold_b. A nondetection retains its threshold and supplies no
realized value. At zero read noise a detected value uses the integer electron
count and discrete counting measure; at positive read noise it uses ADU and a
continuous ADU density. A joint likelihood can therefore have product units
ADU^(-q), where q is the number of continuous detected bands; count masses and
censoring probabilities are dimensionless. The outer record chooses joint
censoring or selected-only measure. Inner observations always use the joint
measure, so per-band selected denominators cannot be applied accidentally.

For record d, define c_sb(d_b) as the parent's mass, ADU density, or nondetection
probability and a_sb=P(Y_b >= threshold_b | s). The only selection event in this
slice is E=all declared bands detected at their record thresholds:

    L(d) = sum_s p_s product_b c_sb(d_b)
    P(E) = sum_s p_s product_b a_sb
    L(d | E) = L(d) / P(E).

Selected-only records must have every band detected. Its denominator is this
same joint mixture event, rather than a product of marginal band probabilities.
Marginalizing the optical state generally creates band dependence. Expected
signal calibration covariance is never added as independent detector noise.
Structural zero support has an explicit flag; unresolved positive values fail.
Every required state remains in input order with its conditional expectation
and detector diagnostics. A required state refusal withholds its record's
mixture without dropping states or renormalizing their masses. Required optical
state/output refusals withhold all mixtures. The explicit state-only request omits population moments and their unrelated unresolved-spread
gate; their default request retains that gate. No moment-only optical approximation
is used.

## Frozen numerical allocation

The optical parent assigns empirical forward sensitivity 3e-12 abs(N_sb).
Evaluate the detector law at the nominal N and outward binary64 endpoints
N +/- that allocation; include endpoint rounding and refuse a nonrepresentable
positive endpoint. Zero photons have an exact zero witness and need no endpoint
calls. Each state-band record carries its nominal, lower and upper attempts.
Maximum endpoint changes in conditional log values and event log probabilities
are added to the parent detector diagnostics. These are named empirical
sensitivity controls, rather than interval guarantees or physical calibration
uncertainties. An endpoint leaving the detector's bounded domain is a refusal.

Products and mixtures use wide log coordinates. Mixture log diagnostics combine
conditional errors with their actual nonnegative mixture responsibilities,
wide accumulation and final binary64 rounding. Selected-only diagnostics include
both numerator and event denominator. The default final dimensionless allocation
is 5e-12 + 1e-8*(1+abs(log L)); the relative scale allows accumulated inherited
photon sensitivity across at most 64 bands near the detector's lambda<=64 limit.
Individual parent gates remain unchanged. A failed gate does not admit a zero,
Gaussian approximation, jitter or a smaller state set.
The returned joint numerator and all-band event log each have their own final
rounding diagnostic and gate, including for joint records. The selected-only
value separately carries both contributions. A rare unresolved event can
therefore conservatively refuse this combined output request.

Named two-state/two-band controls require actual log agreement within
2e-10*(1+abs(reference log L)). Constant-spectrum polynomial photon integrals
supply independent conditional expectations; finite Poisson sums establish
joint normalization, censoring and cross-band dependence. Positive-read controls
use original characteristic-function inversion in each band and explicit finite
state sums, with independently refined quadrature consuming at most 5% of that
allocation. Shared SI constants, mathematical Poisson/Gaussian laws and system
transcendentals are ancestry, not independent instrument qualification.

## Ownership and resources

Inputs are synchronous borrowed spans; results own origin, state, band and record
IDs and scalar attempts, with no retained curve views. Bounds precede sample
scans and allocations. The hard dimensions are 1024 states, 64 bands, 1024
records and 65536 state-band-record cells; configured caps can lower them.
There are at most three detector evaluations per state-band, each including one
threshold-only nondetection probe per record. All attempts and their Poisson
work are retained, including failed endpoints. The admitted output/scratch
payload includes the conservative optical parent quota, result vectors, string
characters, detector rows and bounded mixture scratch. Borrowed caller storage,
allocator overhead, scalar stack and RSS are excluded. Allocation failure returns
a work-limit diagnostic and preserves available attempts without an aggregate.

Ordinary native owner and peer tests require no external data or reference runtime.
