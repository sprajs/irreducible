# Bounded uniform-box Gaussian target

The C++ `GaussianBox` consumer conditions the supplied full-rank fixed linear
Gaussian likelihood on an explicitly finite, positive box. For active ordered
coordinates beta, B=product[l_j,u_j], q=(y-X beta)^T C^-1(y-X beta), the prior is
1/Vol(B) on B. Any fixed coordinate has a declared point-mass provenance outside
this active product Lebesgue measure. The proper-Gaussian prior route is separate.

Completing the quadratic gives q=q_min+(beta-mu)^T H(beta-mu), H=X^T C^-1 X.
With retained AP=QR and original column scale D,
logdet(H)=2 sum log(D_j)+2 sum log(abs(R_jj)). The permutation changes no
determinant magnitude. The relative box integral is

    log I_B = -q_min/2 + p log(2 pi)/2 - logdet(H)/2 + log P_G(B).

The prior-normalized relative evidence subtracts log Vol(B). The normalized
observation-density evidence additionally subtracts
(logdet(C)+n log(2 pi))/2. They are distinct outputs; a source sampler's arbitrary
constant log prior supplies neither normalization. q_min is evaluated from the
wide QR transformed residual tail; the profile's actual post-cast q is retained
separately and is not silently substituted for the minimum.

For every coordinate, use its own mean and marginal variance to form both
standardized endpoint distances. For t>=0, Q(t)<=exp(-t^2/2)/2, and
exp(-x)<=(1+x/65536)^(-65536). This rational majorant is evaluated with outward
basic arithmetic and sixteen squarings. Negative distances receive upper bound
one. If E is the sum of all endpoint tail upper bounds, P_G(B) is enclosed by
[1-E,1] when E<1, irrespective of correlation. E>=1 or a log-width beyond the
request withholds the inadmitted normalization/marginal output. The proper finite
box still exists mathematically when this algorithm refuses it.

For a requested coordinate with Gaussian marginal CDF F, arbitrary dependence
gives max(0,(F-E)/(1-E))<=F_B<=min(1,F/(1-E)). Therefore its u-quantile is within
[F^-1(u(1-E)), F^-1(u+(1-u)E)], intersected with its supplied box interval.
CDF inversion uses outward composite Simpson quadrature of exp(-z^2/2), with
the global fourth-derivative bound 12; the error is bounded by
12*b*(b/N)^4/180 before the normalizing constant. A positive Taylor series for
exp(v), v<=73/256, with a geometric remainder, eight squarings and reciprocal
encloses exp(-x^2/2). Logarithms use the positive atanh series after binary range
reduction. Parsed decimal bounds for pi are widened before use. No certificate
is obtained by putting nextafter around an unbounded libm exp/log result.
The inverse search endpoints at±12 are witnessed by the rational tail bound,
avoiding long-interval quadrature merely to admit that search.

The elementary interval contract requires round-to-nearest binary floating
arithmetic, long double with at least 64 mantissa bits and exponent range16384,
finite normal values or exact zero, disabled fast-math and contraction. Square
root proposals are accepted only after outward squared inequalities witness
both bounds. Overflow, subnormal arithmetic, unresolved interval width or
resource exhaustion refuses the calculation. Every CDF node/attempt is charged,
including failed refinement and bisection. Default allocations are log box-mass
width1e-9, requested original-coordinate quantile width1e-8 and standard-normal
quadrature absolute radius2e-11, with8192 nodes (including endpoints) per CDF,96 bisections per
inverse and256 CDF attempts per evaluation. Thus the total CDF node budget is
at most8192*256, and every refinement attempt remains within that charge.
The reported box-conditional endpoint CDF intervals also include the excluded
mass bound and can be wider than that quadrature allocation; their actual
widths are retained. Final original-coordinate quantile width is admitted
independently, including binary64 outward reporting.

These are **conditional arithmetic and truncation enclosures for the reported
Gaussian completion**. The retained covariance/QR completion has the existing
empirical conditioning and arithmetic diagnostics; these do not certify the
original X,C,y. Source-level log integral acceptance requires a separately pinned
independent completion/reference with its own refinement allocation. The source
fixed44 active46 target freezes log integral comparison at1e-7 absolute and
reference refinement at5% of that allocation. Its requested first marginal is
the original-axis46 median, in the source's5log10 numerical coordinate. Numerical
acceptance does not identify original axis44 physically or qualify observed H0.

Preparation consumes the move-only retained DesignProfile only after support,
metadata and resource validation. Evaluation does not copy/refactor C or X.
The observation normalization reads logdet(C) from the admitted retained factor;
it does not solve an artificial zero residual to retrieve determinant metadata.
Full covariance solve conditioning is distinct from the retained QR completion
and determinant comparison. Original-input determinant accuracy still requires
the independent reference gate.
All active ordered IDs and bounds must match exactly; no dropped coordinates,
jitter, widened box or changed budget repairs a refusal. Public tests use
analytic separable truncated controls, correlated broad-box witnesses, original
unit/axis transformations and near-bound/tail refusals. External original-source
inputs and failed receipts remain in their experiment owner.

Each result carries an operation-owned stage and explicit availability flags.
The separate typed `completion_step` identifies the attempted profile,
whitening, QR, marginal-variance or covariance-determinant gate even when no
completion is available. `completion_parameter_index` is present only during
a marginal-variance attempt, in active source order. A completed Gaussian keeps
its complete completion step through a later box refusal.
An input/setup refusal is unassessed, with no completion. A later tail refusal
retains the completed Gaussian and all endpoint margins; a later inversion
refusal retains already admitted normalization intervals and actual work.
Unadmitted intervals stay unavailable even though their storage defaults to zero.
Only a complete result receives finite status. Preserved diagnostic stages do
not convert a refused requested calculation into a qualified result.

The synthetic suite also retains an isolated observation with variance2^-40 and
zero design response. Its determinant remains in the three-row observation
normalization; an observation2^-20 in that row adds exactly one to q_min.
An overstrict inherited profile sensitivity policy still refuses before any
completion output, with the attempted gate preserved.

The named synthetic owner suite includes these analytic controls. A separate original Decimal
reference uses Machin's pi identity, the integrated Gaussian alternating power
series and monotone bisection at90 and120 digits. It checks the finite[-8,8]
quantiles, asymmetric[-7,9] median and mass, and the rational2D cofactor
normalizing prefactor. Maximum reference refinement is2.67e-70; these fixtures
have ancestry distinct from native Simpson/rational-tail arithmetic. The two
affected retained-design suites pass84 owner and626 peer controls. These tests
qualify their named small cases. The integrated strict Release build, five
affected native suites and fresh installed-library normalization/median consumer
also pass. Actual released active46 source normalization and its independently
refined numerical reference remain separate acceptance work.
