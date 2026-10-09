# Diagonal Bianchi I expansion and directional photon redshift

The native C++ SDK supplies an anisotropic homogeneous GR model, with isotropic
perfect-fluid dust/radiation/Lambda stress and trace-free shear. It returns a
retained batch of mean and principal-axis expansion, anisotropy, shear fraction
and photon redshift for supplied observer-frame directions. This is a Bianchi I
model, with three directional scale factors; it is not an isotropic FLRW model
with an extra density label.

Use `irred/bianchi_i.hpp` and `evolve_bianchi_i`; the
[SDK example](../cpp/examples/bianchi_i.cpp) is a complete native consumer.
The identity is
`diagonal-Bianchi-I-GR-isotropic-perfect-fluid-dust-radiation-Lambda-shear/v1`.
The ABI and CLI do not expose it. It supplies no angular-diameter/luminosity
distance, beam focusing, light-cone optical propagation or observational fit.
In particular, perfect-fluid radiation with isotropic pressure is an **assumed
stress closure**. A freely streaming collisionless radiation distribution may
have anisotropic stress and does not obey this shear evolution.

## Metric, normalization and owned inputs

The metric and comoving observer anchor are

```
ds^2 = -dt^2 + a^2 sum_i exp(2 beta_i) (dx_i)^2
sum_i beta_i = 0
a(observer) = 1, beta_i(observer) = 0
```

Supply `H_anchor>0` in inverse seconds. It is a normalization scale at the chosen
finite observer epoch, not an implicitly inferred present-day `H0`. Fluid
fractions `Omega_m`, `Omega_r`, `Omega_Lambda` are nonnegative and refer to
`rho_c,anchor=3 H_anchor^2/(8 pi G)` at that same epoch. No numerical value of `G`
is needed for these normalized outputs.

Supply two canonical shear coordinates `s_x`, `s_y`. The engine owns
`s_z=-(long double(s_x)+s_y)` and `Omega_shear=sum_i s_i^2/6`. Here
`s_i=dot(beta_i)(anchor)/H_anchor`. The third coordinate is derived once, not
accepted as an inconsistent redundant input and repaired. The represented
fluid-plus-shear sum must differ from one by at most eight binary64 epsilons.
Fractions are never renormalized: actual `E(1)` equals the square root of that
represented sum, so the represented mean anchor H differs from `H_anchor` by
its admitted arithmetic closure discrepancy. Output trace/volume identities
retain their ordinary floating-point residuals, without forced cancellation.

A query contains `a` and a supplied unit direction `n_i` in the observer's
orthonormal principal-axis frame. Its represented squared norm must differ from
one by at most eight binary64 epsilons. It is never normalized. Consequently
`1+z(1,n)=sqrt(sum_i n_i^2)` retains the admitted norm-roundoff discrepancy from
one. Direction signs are preserved; only their squares enter this homogeneous
redshift operator. These are observed directions, not coordinate-direction
cosines at emission.

## Independent physical derivation

For diagonal trace-free shear, `sigma^2=(1/2) sum_i dot(beta_i)^2`.
The Einstein Hamiltonian constraint and trace-free spatial equations, with no
anisotropic matter stress, give

```
3 H^2 = 8 pi G rho_total + sigma^2
dot(dot(beta_i)) + 3 H dot(beta_i) = 0
rho_m proportional to a^-3; rho_r proportional to a^-4; rho_Lambda constant
```

Thus `dot(beta_i)=H_anchor s_i a^-3` and

```
E(a)^2 = Omega_m a^-3 + Omega_r a^-4 + Omega_Lambda + Omega_shear a^-6
beta_i(a) = s_i integral_1^a db / (b^4 E(b))
H_i(a) = H_anchor [ E(a) + s_i a^-3 ]
f_shear(a) = Omega_shear a^-6 / E(a)^2
```

The spatial translation Killing fields conserve each photon's covariant spatial
momentum `k_i`. At the anchor `a=1,beta_i=0`, those momenta are proportional to
the supplied orthonormal direction `n_i`. Comoving observers measure
`omega^2=sum_i k_i^2/[a^2 exp(2 beta_i)]`; therefore

```
1+z(a,n) = sqrt(sum_i n_i^2 exp(-2 beta_i(a))) / a
```

The mean expansion is positive. Individual `H_i` may be positive, zero or
negative. A contracting axis can produce `z<0` for an earlier emitter; the engine
retains this blueshift and never enforces positive redshift or positive axis H.

Jacobs,
[Astrophysical Journal 153, 661 (1968)](https://doi.org/10.1086/149694), is historical
context. Original-paper access was unavailable during this implementation
(HTTP 503), so no source-equation locator or paper-reproduction claim is made.
The implementation and tests use the explicit Einstein/shear/null-momentum
derivation above; no third-party code or paper text is copied.

## Numerical contract and retained batches

Queries proceed from observer toward the past: nonincreasing scales in
`[1e-4,1]`, including duplicates. Output preserves query order and every supplied
direction. An anisotropy integral retained across segments is reused for repeated
epochs and different directions. The complete initial model and derived anchor
shear are retained in the batch, with no surviving borrowed span or callback.
Empty batches validate the initial state and do no quadrature.

Set `s_max=max_i |s_i|` and
`D(a)=Omega_m a^3+Omega_r a^2+Omega_Lambda a^6+Omega_shear`.
For `x=ln a`, the positive quadrature integrates `s_max/sqrt(D(exp(x)))` toward
the observer, accumulating `J=-s_max integral_1^a db/(b^4 E(b))`.
Each beta is `-(s_i/s_max) J`, so the absolute quadrature allocation directly
controls the largest beta coordinate. This scaled integrand avoids integrating
an unnecessarily enormous unscaled shear integral for very small initial shear.
It uses the shared `numerics::integrate` adaptive Simpson kernel, with an
absolute allocation proportional to each log-scale segment's length and a
relative allocation per positive segment. Pure isotropy has beta exactly zero
and needs no quadrature. Pure shear uses the analytic logarithmic integral.

Defaults are `1e-10` absolute beta tolerance and `1e-9` relative tolerance,
maximum depth 30, one million callbacks, 4096 rows and 4 MiB native payload.
Nonzero tolerances must be at least `128 epsilon(binary64)`; the absolute one
must be positive. Round-to-nearest is required. The retained beta estimate sums
quadrature estimates and arithmetic diagnostics; it must fit the requested
absolute-plus-relative cumulative allocation. The directional `1+z` relative
estimate propagates the largest beta estimate with `expm1(error)` and includes
arithmetic diagnostics. These are empirical estimates, **not proven bounds**.
Errors in supplied parameters, directions or stress closure are excluded.
Consumers should check refinement at their actual inputs and error allocation.

All callback attempts, including refused quadrature attempts, count against the
single batch cap. No hidden clamp, composition repair or automatic policy
increase occurs. Before numerical work the payload gate counts the result
header, requested retained rows and a fixed 32768-byte operation/recursive-kernel
allowance. Borrowed input, caller-held copies, allocator metadata and RSS are
excluded. Mean and nonzero axis H, directional scale factors and `1+z` must be
normal in platform `long double` arithmetic; a mathematically zero axis H is
valid. The binary64 quadrature callback must remain normal and positive.
Unrepresentable nonzero quantities refuse rather than becoming zero. The
explicit additional `|beta_i|<=64` gate prevents unsafe exponential projection.
Within the normalized trace-free model, the stated scale interval keeps the
physical anisotropy far below that exponent gate.

Any failure withholds all rows, retaining the model and attempted callback
counter. Nonfinite inputs yield `nonfinite_input`; unsupported physical,
ordering, direction or representability inputs yield `outside_domain`;
resource/depth exhaustion yields `work_limit`; invalid numerical policy yields
`invalid_input`; an unresolved cumulative diagnostic allocation yields
`conditioning_budget_exceeded`. Refusal is not a partial accepted trajectory.

## Permanent independent controls

[Native tests](../cpp/tests/test_bianchi_i.cpp) check exact isotropic
matter/radiation/Lambda limits and vacuum Kasner. For `s=(2,-1,-1)`,
`E=a^-3`, `beta=s ln a`, and directional scale factors are `(a^3,1,1)`.
The other tested Kasner state `s=(-2,1,1)` has a contracting principal axis and
an earlier directional blueshift. Axis and off-axis photon energies are checked
against the conserved-momentum geometry, including `a=1e-4`.

An independent composite four-node Gauss-Legendre rule uses analytically derived
roots/weights of the fourth Legendre polynomial, rather than adaptive Simpson.
It is refined from 128 to 256 cells. A separate proper-time RK4 evolves the scale
factor, all three beta coordinates and all three shear rates, including their
Hamiltonian contribution, and refines from 8192 to 16384 steps. Its time endpoint
is fixed by independently refined Gauss-Legendre lookback integration. Each
reference refinement must consume less than 5% of the named `2e-8`
absolute/relative comparison allocation. Native tolerance refinement is also
checked independently. Named mixtures span initial shear amplitudes `0.01`,
`0.2`, `0.8` with positive dust/radiation/Lambda; these do not qualify every
admitted state. Tests cover signed axes, trace-free volume, duplicates, retained
ownership, invalid normalization/directions/ranges, arithmetic admission and
callback/row/payload/depth caps. They do not qualify observational data or optics.
