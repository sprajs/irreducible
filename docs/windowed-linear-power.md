# Supplied linear power and window mean

The native C++ owner in `irred/windowed_linear_power.hpp` retains a supplied
piecewise-linear spectrum and four named window blocks. One evaluation returns
the coupled monopole/quadrupole on the required theory grid and the ordered
window means. The model is
`supplied-linear-Kaiser-zero-dispersion-geometric-AP-window02/v1`.
It is a conditional forward response. It supplies no observed likelihood,
nonlinear spectrum, tracer calibration or fitted cosmological input.

The source implementation and permanent controls are present; compilation,
independent numerical/refinement passage, resource fit and installed-consumer
validation are separate gates. Executable CLI/ABI discovery has no such route.

## Physical identity and equations

At the supplied redshift, fix one emitted linear matter spectrum in physical
units, k in Mpc^-1 and P in Mpc^3. Node powers must be nonnegative. The declared
working law is linear interpolation in k between the exact supplied binary64
nodes; it is not a bound on an underlying continuous CLASS spectrum. Retain the
matter component, spectrum/growth ancestry, epoch and grid order. A single
scale-independent supplied f is a model choice; radiation, photons and massive
relic perturbations are not evolved here. [GR growth](gr-growth.md) and
[perfect-fluid transfer](linear-transfer.md) retain their separate assumptions.

The linear redshift-space law is delta_s=(b1+f*mu_true^2)*delta_m, with fixed
finite real b1/f, zero velocity dispersion and zero extra stochastic/shot noise.
Negative response coefficients are mathematically allowed; P2 and window means
remain signed. Negative physical spectrum nodes are refused without clipping.

Supply positive geometric alpha_perp=D_A/D_A,fid and
alpha_parallel=H_fid/H at this same epoch. No distance, H, sound ruler or prior
is computed. For mu=mu_fid in [0,1],

```
Q = sqrt((1-mu^2)/alpha_perp^2 + mu^2/alpha_parallel^2)
k_true = k_fid*Q
mu_true = (mu/alpha_parallel)/Q
P_fid = (b1+f*mu_true^2)^2*P_L(k_true)/(alpha_perp^2*alpha_parallel)
P0 = integral_0^1 P_fid dmu
P2 = 5*integral_0^1 P_fid*(3*mu^2-1)/2 dmu
```

The inverse-volume factor follows from x_true=A*x_fid and k_true=A^-T*k_fid.
Historical ruler-conditioned fitted alphas include r_s,fid/r_s; they cannot be
inserted as geometric ratios when changing the full spectrum.

Use the literal supplied dimensionless, row-major four-block action:

```
m0_i = sum_j (W00_ij*P0_j + W02_ij*P2_j)
m2_i = sum_j (W20_ij*P0_j + W22_ij*P2_j)
```

All output-0 rows precede all output-2 rows, with their own retained numeric IDs
and axes. Columns refer to one increasing theory grid. No extra dk, transpose,
row normalization or integral-constraint correction is applied. Unsorted public
column permutations refuse, even when all four blocks were permuted together.

## Units, support and ownership

The window can use physical coordinates, or a fixed supplied h_ref convention:
k_phys=h_ref*k_h and P_h=h_ref^3*P_phys, in (Mpc/h_ref)^3. Preparation retains
original axes and maps its theory grid once. W is unchanged. Window sums use
physical powers; each mean projects once afterward. Optional serialized input
multipoles project separately. h_ref is retained across model evaluations and
is distinct from CLASS h_model and varied H0. CLASS source conversion belongs
to the pinned upstream adapter, not this model.

Every retained theory column requires the whole transformed support
[k/max(alpha),k/min(alpha)] inside supplied spectrum support. Grid-conversion
diagnostics additionally require the displaced interval
[(k-e_k)/max(alpha),(k+e_k)/min(alpha)] inside that support before slope
forwarding. Exact physical-axis endpoint copies have e_k=0; a fixed-h boundary
with positive mapping displacement may refuse. There is no extrapolation,
clamping, dropped column or repaired split. Angular integration splits the
piecewise-linear knot crossings; unresolved numerical projection refuses.

Prepare once by moving a `Source` into `prepare`. Validation failure leaves
the input unconsumed. A valid `Prepared` is move-only and evaluates through a
const coarse call; moves invalidate the source owner, self move is a no-op.
No covariance or observation object is accepted. Results own bounded source
metadata and requested numerical vectors. Required late failure withholds the
whole numerical batch and preserves acquired identities and all attempted
work. Early policy/payload/arithmetic rejection before result-metadata
acquisition reports empty metadata; the prepared source remains retained.

## Numerical and resource scope

The compiled interpolation equation is the existing two-knot
`numerics::interpolate_linear`. The new owner brackets a validated retained
grid rather than validating the whole spectrum on each callback. The existing
adaptive Simpson integrator is applied to the knot-resolved angular cells.
Half the explicit Mpc^3 multipole allowance goes to width-weighted absolute
quadrature; relative tolerance is zero. Callback, interpolation, geometry,
cast and accumulation diagnostics are owned here, because the interpolation
helper itself reports zero internal error estimate.

Error components are quadrature, arithmetic, grid mapping and output volume
projection. Absolute window coefficients transport input diagnostics; signed
signal cancellation does not cancel them. Final means must pass the requested
output-unit allowance. Optional serialized multipoles must independently pass
their total, including h^3/cast diagnostics, converted back against the original
physical allowance, even when W annihilates that multipole. These are empirical
estimators/operation diagnostics, not complete physical or arithmetic enclosures.

There are no inferred physical defaults. Explicit policies may select smaller
limits than the maximum supported slice: spectrum4096 knots, theory2048
columns, total256 rows; preparation2M actions; evaluation8M actions and2M
angular callbacks;4096 callbacks per smooth integral and depth20. Callback
admission reserves its worst bracket/interpolation action cost before entering
the shared integrator. Failures count; there is no automatic cap growth/retry.
Each evaluation is a transaction, with preparation cost reported separately;
no global model-campaign budget is claimed.

The chosen known-payload admission ceiling is16MiB. Requested old+pending
storage is checked before allocation; actual vector/string capacities are
checked afterward. Fixed reservations avoid growth. Inventory includes owner,
scratch, result metadata, original capacities and simultaneous move/return
headers. An unexpected allocation profile refuses. A live-peak qualification
requires the actual selected allocator/library and requested-new observer in
named cases. It excludes allocator bookkeeping, C malloc, opaque library
frames and whole RSS; no universal reserve(n)==capacity(n) or whole-core
no-elision qualification follows from this implementation.

The [synthetic SDK example](../cpp/examples/windowed_linear_power.cpp) supplies
the exact emitted nearly affine table, signed cross-block window and fixed
geometry. Its reference must use that table's interpolation, rather than
replace the decimal nodes by an exact untabled affine function. Permanent
analytic affine controls use dyadic nodes/powers. The selected caller budgets
are1e-7Mpc^3 multipole diagnostics and1e-6 output-unit mean diagnostics, a
synthetic regression resolution rather than an observational error requirement.
The independent Legendre16/32/64 peer also compares each measured native
discrepancy with the emitted native combined diagnostic plus that peer's own
32-to64 change. This is a required empirical consistency control. Neither the
last refinement change nor its agreement with the native estimator is a proven
reference enclosure.

## Source and unassessed model errors

[Gil-Marin et al., 1509.06386v2](https://arxiv.org/abs/1509.06386v2) supplies the
half-range multipole convention and literal four-block discrete action. Its
historical nonlinear/nonlocal-bias/TNS/FoG fit is a different model: the exact
Lorentzian FoG kernel is delegated and its fitted alphas include a ruler factor.
The preceding window-construction index/normalization ambiguities are not
repaired here. Public availability/citation instructions establish no new
asset redistribution license; no external numerical assets are vendored.

Even no-AP Kaiser has P4=8*f^2*P_L/35; anisotropic AP creates further even
multipoles. Their omission from the supplied 0/2 operator is UNASSESSED model
discrepancy, not zero numerical error. Released window calibration/support,
mode effects and omitted integral constraint remain unassessed. Only6/106
released covariance diagonals agree with the candidate order. There is no
inverse covariance, released likelihood, row repair/cut or inference claim.
