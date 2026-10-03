# Synthetic two-state retained H/He history SDK proof caller

The bounded SDK proof caller in `cpp/tests/test_installed_abundance_history_cohort.cpp`
composes existing public abundance and H/He history owners. Its helper remains
local test/consumer source: no public SDK API, CLI operation or C/Rust ABI is
added. Numerical, installed, resource and complete independent-reference
execution are unverified until their recorded gates pass. This guide states
the fixed source and arithmetic contract before runtime.

The source identity is
`synthetic-two-equal-state-H1-He4-neutral-mass-photon-temperature-retained-HHe-working-history-law/v1`.
It propagates an original **synthetic** two-state law through the bounded
[H/He singlet history](hydrogen-helium-history.md). It does not qualify NEXT-12's
physical abundance/mass/temperature probability law. The existing
[finite LTE law](baryon-abundance-law.md) accepts arbitrary supplied matter-T
vectors. Those cannot be dropped, truncated or relabeled as this caller's
supplied photon temperature. Here T_m starts at T_CMB,0*(1+z_i) and is evolved;
it is an output rather than a supplied trajectory.

## Original complete support and outputs

Every scalar, ID, label and binary64 bit assertion is frozen in
`cpp/tests/abundance_history_cohort_source.hpp` before any runtime. A then B
each have relative mass1 and exact probability1/2. Both have empty explicit
species, physical omega_b=.02237, omega_cdm=.12 and other-massless density1.7e-5.
A supplies H0=63 km/s/Mpc, photon T_CMB,0=2.7 K and Y_He4=.24; B supplies
H0=77, T_CMB,0=2.75 K and Y_He4=.26. Both supply neutral effective masses
m_H1=1.6735328383153192e-27 kg and m_He4=6.646479071583153e-27 kg, with the
existing installed synthetic abundance-fixture ancestry explicitly recorded.
These source choices supply no observed abundance or mass probability law.

The complete physical source, neutral effective mass convention, fixed G/atomic
assets, state/order/row IDs and source/mass/photon-temperature origins remain
owned. One [abundance map](baryon-abundance.md) emits nuclei at a=1 per state.
Those exact emitted binary64 nuclei and the history's exact mapped thermal
background are the working kinetic inputs. Original source-map diagnostics,
once-captured thermal witnesses and emitted bits are preserved provenance.
They are **not propagated through the ODE** or added as random covariance.
Source-to-trajectory sensitivity for a continuous original physical input law
remains open. No remapping, new physical equation or conditional drag root is
used by this caller.

Boundaries are z_i=2700 and z_L=300. The original ordered rows are 2700,1300,300
with mask24, matter temperature then per-redshift Thomson opacity at each row.
Native pooled queries add no ODE or background preparation work. The six-axis
population moments are

```
mu_i=(f_Ai+f_Bi)/2
Delta_i=f_Ai-f_Bi
C_ij=Delta_i*Delta_j/4.
```

Cross-row and cross-output covariance stays in one full symmetric matrix.
Temperature means have K units, opacity means per-redshift units, and covariance
has the corresponding product units. This exact law has rank at most one;
singularity is valid. Projected binary64 entries need not have exact rank/PSD
eigenstructure. No jitter, diagonal replacement, PSD repair or independence
assumption changes the law. Any required state/row/group or aggregate numerical
failure withholds all moments while retaining both original state masses and
actual attempts. Missing native rows remain missing. Equal rounded outputs from
different relevant raw inputs do not establish zero variance.

## Frozen finite-operation diagnostic graph

Let W be the supported strict nearest long double, u=epsilon(W)/2. Each
nonzero normal wide add/subtract/multiply result r receives the finite-
operation round allowance `r_round=abs(r)*u/(1-u)`, rounded upward in W.
Zero from exact cancellation has no round allowance; a nonzero product/input
lost to wide zero, overflow or unsupported wide subnormal refuses. These
per-operation allowances are empirical floating arithmetic diagnostics. The
complete reference gate is separate. Each diagnostic sum/product is rounded
upward with `nextafter` at the operation, with exact-zero terms kept zero.

For a mean, compute s=x_A+x_B in W, then halve by an exact binary scaling.
`a_mu` is half the s round allowance plus the measured binary64 output-cast
loss. For a difference, compute d=x_A-x_B in W; `a_Delta` is that subtraction's
round allowance. The difference is retained in W and has no binary64 cast.
For covariance compute p=d_i*d_j in W, then quarter by exact binary scaling.
`a_C` is one quarter of p's round allowance plus measured covariance output-
cast loss. Diagnostic projection is outward; its actual reported value is
tested against the allocation. Binary scaling must retain a nonzero normal W
coordinate or refuse. No cast loss or diagnostic rounding is hidden in a method
name or a count-times-epsilon slogan.

With each parent's reported absolute error e_si:

```
E_mu_i=(e_Ai+e_Bi)/2+a_mu_i
E_Delta_i=e_Ai+e_Bi+a_Delta_i
E_C_ij=(abs(d_i)*E_Delta_j+abs(d_j)*E_Delta_i+
        E_Delta_i*E_Delta_j)/4+a_C_ij.
```

The mean and covariance errors are deterministic empirical numerical sensitivity,
separate from physical C. Positive diagnostics must project upward and remain
representable. Means retain T_m's `2e-4 K+2e-6*abs(mu)` and opacity's
`2e-9+3e-6*abs(mu)` allocations. Covariance sensitivity must fit one percent
of the geometric variance scale. The reducer uses the exact two-state identity
`sqrt(C_ii*C_jj)=abs(d_i*d_j)/4` to avoid an unnecessary libm square root;
the scale is conservatively reduced by the product's wide rounding allowance.
This new one-percent spread resolution was explicitly admitted for this caller;
LTE moment budgets do not transfer. Unresolved small positive spread refuses.

Bit-identical complete relevant physical raw inputs at common boundaries/query/
policy witness a constant axis. At the initial row alone, identical raw photon
T_CMB,0 at common z_i is a sufficient narrower T_m witness. A witnessed constant
axis has exact zero covariance rows/columns and zero covariance sensitivity;
its mean retains forward numerical error. Nonconstant axes need strictly
positive resolved variance. Positive mean/variance or signed nonzero covariance
lost on projection refuses, even within an absolute allocation. Resolvable
binary64 subnormals charge measured cast loss. Off-diagonal zero is not a
declaration of independence.

## Work, complete live storage and ownership

This scope explicitly requests preparation cap8,000,006, logical evaluation
cap1024, total two-evaluation cap8,002,054 and whole-live cap128 MiB. These new
scope limits do not change native child defaults N8192/fine65536/work4M/32MiB or
their numerical allocations. The N16384 comparison is a separate recorded
request at unchanged caps/allocations. Preparation counts two complete source
intakes, four abundance prepare/map calls and each of the native four work
categories once, including failure. Every child receives the remaining global
work/byte ceiling; requested and served values are separately owned.

Evaluation logical counts are two dispatches, six queried rows, twelve required
scalar collections, six bounded raw-input witnesses, six mean coordinates,
six difference coordinates, twenty-one upper-triangle products and twenty-seven
result projections. Failed operations still count. These are bounded logical
categories, not CPU instruction counts. Physical preparation work is reported
separately; output work sums never supply ledger freshness.

The whole-live expression reserves two complete public
`hydrogen_helium_history_payload_bound` envelopes using the admitted maximum
fine-interval bound, three outputs and supported origin-copy capacity. They
remain reserved throughout the cohort lifetime, even after native scratch is
released. They bound opaque history storage; they do not measure its private
node capacities. Caller-known actual capacities are reported separately.
Add two local owner headers for unelided return/move, the immutable source
header and all string/species copy envelopes, abundance source/owner/temporary
headers/strings, a history request header/origin, map/emitted snapshots,
both prior and new complete query receipt headers/native row capacities,
fixed six-axis/full-36-cell moments and diagnostics, and wide reduction scratch.
All fixed arrays are charged at their full sizeof storage; their stack location
does not exempt them. Shared source is charged once. Checked products/additions
and supported copied-string `max(32,length+1)` envelopes precede allocation;
actual capacities are checked afterward. Allocator/control-block bookkeeping,
borrowed input, scalar stack frames and RSS are excluded.

The local owning wrapper deletes its copy route, moves each history once and
queries through const references. Moves clear source ownership and invalidate
the moved-from wrapper; self move preserves it. Receipts share immutable source
and own actual native query attempts/moments without copying histories. At most
prior and new query receipts coexist in the installed two-evaluation scope;
the final receipt can outlive caller release. Allocation failures use RAII and
preserve acquired source/work/actual child evidence without fabricating rows.
The byte request is not an earned measured resource claim before threshold,
allocation-site and complete library/caller no-elision checks execute.

## Verification status and physical limits

Permanent exact-rational half-law, signed covariance, refusal, ordering,
resource, allocation, lifetime, strict-policy and no-elision controls accompany
the caller. Their actual execution must be recorded before they are described
as passing. A fresh installed-header/static-library caller gate is independent
of simply linking tests inside the build tree.

An independent complete trajectory/moment comparison must use exact emitted
working inputs and restore original row order and duplicates. All affecting
arithmetic, nonlinear-solve, interpolation/projection and refinement error must
fit <=5% of each unchanged group/moment allocation. The existing physical-
remapping/sorting reference cannot be run unchanged for this law. Complete
five-percent reference qualification remains withheld; inherited Decimal rate
facts and long-double/system-libm stiff agreement are not a full high-precision
trajectory certificate. Fixed support/allocations cannot be tuned after failure.

This synthetic pushforward supplies neither a physical joint abundance law,
BBN, atomic/model uncertainty, arbitrary matter-T history, full helium stages,
physical drag endpoint, reionization, present-day CMB visibility nor a posterior.
