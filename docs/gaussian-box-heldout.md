# One heldout row under an original finite box

`irred/gaussian_box_heldout.hpp` is a standalone C++ consumer for a supplied
fixed linear design, full correlated Gaussian observation noise and a proper
normalized finite uniform parameter box. It prepares one ordered training
factor and retained QR, then evaluates coarse pools of training vectors,
heldout candidate values and ordered request pairs. It has no CLI or C ABI.
The source implementation, build/runtime acceptance and independent released
reference remain separate gates; execution is currently unearned.

This law is conditional on the supplied design, offsets, covariance and box.
It is generally a mixture, not a scalar Gaussian parameter posterior. A
mathematical training box weight does not qualify observational inference.
The proper Gaussian prior consumers have a different target and are not used.

For the original ordered observations, partition one row `t` and all remaining
rows `T`, preserving their original order. Let `r_T=y_T-o_T`, `r_t=v-o_t`,
`k=C_Tt`, `u=C_TT^{-1} k`, `s=C_tt-k^T u`,
`a=x_t-u^T X_T`, and `b=r_t-u^T r_T`. At a fixed parameter vector beta,

```
y_t | y_T,beta ~ N(o_t+x_t beta+u^T(y_T-o_T-X_T beta), s)
q_J(beta)=q_T(beta)+(b-a beta)^2/s
 det(C)=det(C_TT)*s.
```

For the same active box `B` with volume `V`, define
`I_T=integral_B exp(-q_T/2) d beta` and the analogous `I_J`. The two normalized
observation evidences are

```
Z_T=(2*pi)^(-|T|/2)*det(C_TT)^(-1/2)*I_T/V
Z_J=(2*pi)^(-(|T|+1)/2)*det(C)^(-1/2)*I_J/V
p_B(v|y_T)=Z_J/Z_T
log p_B=log I_J-log I_T-0.5*(log det(C)-log det(C_TT)+log(2*pi)).
```

The consumer composes the two reported relative-integral enclosures directly.
It cancels `V` and the common observation constants symbolically. Its remaining
constant uses both actual retained covariance log determinant reports, including
the original full covariance factor report captured before source release.
The Schur computation is a preparation diagnostic and whitening coordinate; it
does not silently replace that reported normalization constant.

The original parameter IDs and their active original indices are explicit.
Exactly one separately named fixed coordinate is literal zero and a point mass
outside the active measure. The design contains all original rows with active
columns in the declared original-axis map; original full-column bytes remain
external immutable lineage. Its value, unknown label where applicable and provenance are
retained. Finite box bounds must have positive width in the original native
parameter units. No zero-width Gaussian prior, Jacobian change, dropped axis,
refitted box or sign/orientation repair is introduced. Offsets are supplied in
the original observation units and subtracted in wide arithmetic; actual
reported coefficients and adjusted residuals remain binary64 and their
stationarity is checked against the retained normalized design.

Preparation requires an authoritative covariance `MatrixKind` recorded by
`prepare_gaussian`; a caller-authored convention string cannot establish it.
It requires full declared matrix validation and empty applied latent priors,
mean shifts and selection history. An external owner must still pin original
bytes, dtype, adapters, row/event/calibrator lineage and scientific identities.
The source remains valid after any failed setup. Only complete successful
preparation consumes it.

One principal `C_TT` factor is retained. With `h=L_T^{-1}k`, the conditional
variance is `C_tt-sum(h_i^2)`. The empirical cancellation screen is

```
K_s=(abs(C_tt)+sum(h_i^2))/s
E_s=K_s*(2*cross_whitening_rounding_estimate+(2*|T|+4)*epsilon_long_double).
```

A positive computed variance with this screen is not an original-input
certificate. Negative, unresolved, subnormal or nonfinite intermediates refuse;
no jitter, variance floor or independence fallback is used. The appended
normalized response is `(x_t/D-h^T B)/sqrt(s)` for the retained
`B=L_T^{-1}X_T D^{-1}`. Givens rotations update `[R_T; g_P]` once, retaining the
same scaling and pivot order. Updated triangular conditions and all marginal
variance screens are re-admitted before source consumption. A deficient training
design refuses even when a finite box integral exists mathematically.

Each training pool member is whitened and transformed once. A bounded cache
retains exactly `2*p+2` wide scalars: the QR head, `p` stationarity denominator
sums, orthogonal-tail quadratic and `h^T` whitened RHS. Joint values reuse that
cache. Each value still whitens its actual rounded adjusted residual to check
post-cast stationarity. No original `C` refactor, per-row FFI or generic cache
framework is involved. Requests and all pool orders are retained.

The shared box normalizer uses the existing rational-tail omitted-mass bound,
outward elementary operations and reported Gaussian completion. It performs no
CDF or quantile work. The enclosure scope is explicitly conditional on both
reported completions and does not certify them against original `C,X,y`.
The training and joint integral widths each retain the existing `1e-9` box
allocation. The final reported conditional log-density width is at most
`3e-9`, including covariance-constant and reporting arithmetic. Numerical
refusal for a finite candidate is not zero probability or `outside_support`.

`BoxHeldoutPolicy` requires explicit nonzero preparation and evaluation work
caps. Public source-derived upper-bound functions expose their complete logical
loop/copy/helper charges, including failed prefixes and independent candidate
pool scanning; they are not FLOP or instruction measurements. Evaluation may
only tighten the frozen setup policy. Payload receipts count simultaneous
retained/setup/scratch/output ownership, vector/string capacities and no-NRVO
headers; borrowed arrays, allocator bookkeeping and RSS are separately excluded.
The released protocol additionally requires its own wall-time/address limits
and shared job lease. Default/failed/moved-from accessors are safe, failed
receipts move intact, and self moves preserve state.

Preparation metadata work depends on logical contents, so reserves and copies
of the same metadata preserve its work upper; payload still counts actual
capacities. Its complete copy/move/default/literal ledger includes an explicit
fixed 417-unit term for literal length scans and metadata headers. This changes
the source-derived required work rather than a configured cap or default: the
same cap can refuse through the existing precharge and retain the source.

A failure retains only availability-qualified diagnostics. Optional raw refusal
witnesses may be nonfinite and must be transported as tagged values, not bare
JSON `nan`/`inf` or invented finite zero. The installed consumer source includes
accepted, Schur-refused, shape-refused and serializer-only nonfinite and partial
allocation-layout controls. A batch earns `output_layout_available` only after
all public arrays finish initialization; a partial layout can have different
array lengths and missing flags withhold their scalar as `null`. Per-item
allocation failures retain their actual attempted work and earned diagnostics,
while a failed training cache withholds every dependent request. Late retained
QR allocations return the actual whitening prefix. The separate bounded
allocation-ordinal source control checks tiny setup/profile/pooled-evaluation
prefixes, requested C++ heap peaks and leak freedom; it disables injection
before serialization. Argument construction before entry can still throw and
cannot consume the original source.
Synthetic analytic and independent direct rectangle-quadrature tests are
separate from the released reference and its qualification.

For NEXT05 the preregistered source identity retains all original 3492 rows,
47 original parameters and original covariance. Only zero-based original row
3213 is withheld. Original axes 0--43,45,46 remain the 46-dimensional measure;
axis44 is literal zero outside it. All `C_TT,C_Tt,C_tT,C_tt` entries and shared
calibrators remain. Physical/H0/release-equivalent inference remains unsupported
while the NEXT04 reduction, axis/anchor sign and full-fit prior/preprocessing
issues remain open. The root-approved released comparison allocation is
`2.1e-7` with shares `(1e-7,1e-7,5e-9,5e-9)`; the complete independent-reference
allocation is `1.05e-8` with shares `(5e-9,5e-9,2.5e-10,2.5e-10)`. Acceptance
compares the farthest native interval endpoint against the independently bounded
reference center plus complete reference radius. These allocations do not alter
older box budgets or establish scientific interpretation.
