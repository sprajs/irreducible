# Finite supplied radial-shell projection

The C++20 `irred/supplied_shell_projection.hpp` source defines a bounded
synthetic projection of already integrated, real isotropic radial atoms:

\[
 \Theta_\ell(k_i)=\sum_{s=0}^{S-1}A_{is}j_\ell(k_i\chi_s).
\]

The caller supplies the dimensionless signed amplitudes `A`, including every
shell integration weight. The calculation does not infer a source from an
expansion history, recombination history or perturbation solution. It has no
CLI or C ABI route. Numerical qualification requires the native, installed and
independent controls for the exact source/compiler candidate; the source
implementation and its proposed controls alone do not establish that result.

## Source and conventions

`SuppliedShellView` borrows coarse spans of `k_mpc_inverse`, `chi_mpc` and
`amplitudes`, ordered wavenumber and shell identifiers, and shell roles. The
amplitude array is wavenumber-major and shell-minor, with exactly `K*S` values.
Acquisition validates the original input once, then retains one immutable
source. Duplicate coordinates, identifiers and requested multipoles preserve
their supplied order. Nonempty UTF-8 `source_identity`, `mode_origin` and
`ordering_provenance` describe the caller's source; labels do not certify their
contents.

Coordinates use comoving Mpc and inverse Mpc, with no factor of `h`. The fixed
source mode is one supplied scalar unit-\(\zeta\) response. Its regular-mode,
gauge and species origin remain the caller's responsibility. The outward
line-of-sight convention and Fourier factor are \(e^{+i\mathbf k\cdot
\mathbf x}\). With standard orthonormal complex spherical harmonics and the
Condon–Shortley phase,

\[
 e^{ix\mu}=\sum_{\ell=0}^{\infty}(2\ell+1)i^\ell j_\ell(x)P_\ell(\mu),
 \qquad
 j_\ell(x)=\frac{1}{2i^\ell}\int_{-1}^{1}e^{ix\mu}P_\ell(\mu)\,d\mu.
\]

The corresponding scalar harmonic response has the `4*pi*i^ell` convention.
The API returns only individually requested signed `Theta` coefficients for
`ell=0..8`. It supplies no angular reconstruction or remainder for the
infinite multipole expansion. `radial_atom` and `boundary_response` roles
distinguish supplied terms; they share the stated isotropic kernel. A boundary
role does not create an omitted boundary amplitude.

An optional `EndpointView` retains the producer's identity, coordinate/unit
labels, finite ordered support endpoints, survival status, optional survival
and optional absolute diagnostic. Its coordinate is the producer axis; it is
not converted into `chi`. A present survival must lie in `[0,1]`. Failed or
missing upstream diagnostics remain failed or absent. The projection does not
renormalize survival, supply missing boundary mass, or include upstream source
or shell-discretization error in its numerical radius.

Raw Doppler `mu` factors, Bessel derivatives, complex or multipole-dependent
amplitude tables, polarization, vector, tensor and isocurvature sources are
outside this input law. Kelvin conversion and a primordial spectrum are
separate operations.

## Domain and ownership

Acquisition admits `1<=K<=256`, `1<=S<=64`, finite nonnegative coordinates and
finite `abs(A)<=1`. Evaluation checks every ordered `(k,chi)` pair, including
zero amplitudes: its exact product must be zero or conservatively enclosed in
`[2^-32,8]`. An unresolved boundary enclosure refuses; no phase is clipped or
replaced by zero. Each request contains one to nine ordered multipoles and the
exact `supplied_shell_signed_projection` output mask. Unknown bits and
unsupported multipoles refuse.

`PreparedProjection` and `ProjectionBatch` are move-only owners. Each result
retains the same immutable source object, without copying its tables or
preparing another background. Source and result spans are borrowed views whose
owner must remain alive. Results may outlive the prepared owner. Prior and next
results may coexist within the shared source's whole-live payload cap. Moves
transfer storage and receipts; self moves preserve them. Calls to one prepared
owner are serial, since they share its cumulative work ledger.

For example, after filling the source view and ordered identifiers:

```cpp
auto prepared = irred::projection::acquire_supplied_shell_projection(source);
const unsigned ell[] = {0, 1, 8};
auto batch = prepared.evaluate({ell});
// Inspect preparation, batch and row statuses before optional values.
// batch.source() is the same retained source as prepared.source().
```

Every row carries its original `k_index`, `ell_index` and `ell`, a numerical
status and an operation-owned check. Accepted rows contain a signed binary64
value, a complete outward numerical radius and separated diagnostics. Refused
quantities remain absent. A failed pair refuses all requested multipoles for
that wavenumber; independent later wavenumbers remain visible when work permits.
Batch and work receipts preserve the first failure and the charged prefix.
Each refused row retains its operation cause. A separate optional
`terminal_work_status()` records later exhausted work without replacing an
earlier phase-domain failure.

## Numerical method and arithmetic

The native kernel uses the direct spherical-Bessel series

\[
 t_0=\frac{x^\ell}{(2\ell+1)!!},\qquad
 t_{n+1}=-\frac{t_nx^2}{2(n+1)(2\ell+2n+3)},\qquad j_\ell=\sum_n t_n.
\]

The initial factors and one reused square determine the operation graph.
Computed-term rounding is bounded using `d_n=2*ell+3*n`,
`gamma_d=d*u/(1-d*u)` and `eta_d=gamma_d/(1-gamma_d)`. The ordered signed sum
has its separate absolute-sum rounding bound. Only after included index three
may the decreasing alternating first-omitted term stop the series, at a
downward-enclosed exact `1/10^14`. There are at most 33 included and 34 attempted
terms per kernel. An unattained tail refuses.

The proposed arithmetic profile is radix-two binary64 with precision 53 and
wide `long double` with precision 64, wide normal minimum exponent `-16381`
and maximum exponent `16384` in `numeric_limits` notation. It requires
`FE_TONEAREST`, no fast math and disabled implicit FMA contraction. Other
profiles refuse. Nonzero nominal intermediates and public numerical bounds
must remain finite normal. Algebraic input zero is distinguished from rounded
cancellation. A cancellation-zero value can retain a positive radius.

Elementary radius operations use one rounded wide operation and an immediate
representable outward neighbor. Positive denominators and the acceptance
allocation are enclosed downward. No library Bessel, sine/cosine, square root
or power routine enters the kernel. The bound includes argument error through
`abs(j_l')<=1/2`, the Bessel tail and arithmetic, amplitude multiplication,
ordered signed accumulation and the measured single value-cast loss.
Diagnostics retain final component/radius-storage assembly padding separately.
It is already represented in the complete radius and must not be added again.

The complete emitted radius must fit the conservative exact allocation

\[
 A_\mathrm{native}=\frac{1+|\widehat\Theta_\ell|}{10^{10}2^q},
 \qquad 0\le q\le32.
\]

`accuracy_reduction_power=q` only tightens acceptance. It does not change the
tail goal or silently expand a work cap. These bounds concern the exact
emitted finite-source law, excluding physical input uncertainty, missing
source closure and continuous line-of-sight discretization error.

## Work and payload limits

`ProjectionResources` may tighten the hard limits of 6,000,000 attempted series
terms, 8,000,000 total named nonbyte work units, 65,536 inspected source bytes
and 4 MiB simultaneous requested engine payload. Structural maxima can refuse
these limits. Preparation, failed requests and repeated evaluations share one
ledger; a move or refusal does not reset it. Named work includes each source
scalar/tag inspection, shell visit, phase product, kernel, attempted term,
shell product, signed addition, elementary radius operation and output cast.
An unavailable quota refuses before execution and records a rejected request.

The source owns exact typed arrays and one provenance byte block. Let `C` be
the actual shared-source/control allocation request, `I` the size of its
identifier-offset pair, `G` the shell-role size and `B` the stored string bytes.
Its retained payload is

\[
 C+\operatorname{sizeof}(double)(K+S+KS)+I(K+S)+GS+B.
\]

The shared live ledger also includes the prepared state, every live result
state and its `K*L` rows, and the active scratch allocation. Integer additions
and products are checked before allocation; prior results remain charged.
Blocks are destroyed before their reservation is released. Payload describes
actual requested engine storage, excluding caller spans/handles/containers,
ordinary call frames, allocator and exception metadata and process RSS.

## Independent controls and remaining gates

The independent Python reference consumes a preserved native test corpus;
its checker never launches the native binary. It imports exact float bits and
uses separately directed Decimal precision-100 interval arithmetic, original
exact-rational Sturm isolation of all 16 Legendre roots, positive GL16 weight
enclosures, its own signed trigonometric Taylor series, and one-, two- and
four-panel angular refinement. Analytic GL remainder and root, weight,
arithmetic, trigonometric and refinement diagnostics remain separate. Exact
zero and `ell=0,1` analytic controls complement the angular algorithm.

For its retained center `C_ref`, radius `E_ref`, native radius `E_native`,
outward difference `Delta` and downward allocation
`A_ref=1e-10*(1+abs(C_ref))`, all four external gates are required:

1. `E_ref <= .05*A_ref`.
2. `E_native <= .90*A_ref`.
3. `Delta+E_native+E_ref <= A_ref`.
4. The independently enclosed intervals overlap.

The reference hard caps are `K,S<=16`, nine requested multipoles, 262,144
angular node visits, 33,554,432 trigonometric terms, 128 term attempts per
trigonometric evaluation, 8,192 root checks, 16,384 temporary integer bits,
40,000,000 complete named work units and 8 MiB whole-live payload. Root/node
preparation and discarded refinements count. Actual Python/libmpdec profile,
allocation controls, compiler flags and matched-quality costs require saved
execution evidence; proposed formulas alone do not qualify either route.

The [primary CMB contract](primary-cmb-projection.md) remains open. This finite
projection supplies no photon/baryon/species/metric hierarchy, gauge closure,
present-day visibility or reionization history, continuous source integration,
primordial spectrum, temperature `C_l`, measured likelihood or posterior.
[Finite-endpoint histories](recombination-drag.md) retain survival mass and
their own source/error gates. No background, thermal, visibility, variance or
likelihood owner is duplicated here.
