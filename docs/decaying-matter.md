# Pressureless dark matter decaying to massless radiation

The native C++ SDK evolves a bounded, homogeneous, flat-GR initial-value
problem with pressureless parent dark matter, massless daughter dark radiation,
optional stable pressureless matter and a cosmological constant. It provides
retained ordered expansion/composition/energy-transfer trajectories for synthetic
background comparisons. It does not provide perturbations, a transfer function,
BBN, recombination, CMB observables or a qualified inference model.

Use `irred/decaying_matter.hpp` and `evolve_decaying_matter`. The physical identity
is `flat-GR-finite-anchor-pressureless-decay-massless-daughters/v1`; the numerical
method is `positive-comoving-daughter-log-scale-RK4-step-doubling/v1`. This is
separate from the thermal species background and the independent effective fluid.
The [SDK example](../cpp/examples/decaying_matter.cpp) is a complete native
consumer. The ABI/CLI does not expose this model.

## Initial state, equations and source

Supply a finite initial scale factor `a_i>0`, an expanding initial Hubble rate
`H_i>0` **in inverse seconds**, and stable/parent/daughter/Lambda fractions at that
same epoch. These refer to `rho_ci=3 H_i^2/(8 pi G)`, not today's critical density.
`H_i` is an initial anchor, not `H0`; the scale factor need not equal one at any
present epoch. The fraction sum must equal one to within eight binary64 epsilons;
values are never renormalized. Consequently the actual represented initial
`E=H/H_i` equals the square root of that represented sum, within its admitted
roundoff discrepancy from one. This API reports densities divided by `rho_ci`;
conversion to a physical density requires a separately chosen `G` and unit system.

The supplied constant `g=Gamma/H_i>=0` defines a proper-cosmic-time decay rate
`Gamma=g H_i` in inverse seconds. Positive rates decay; zero gives the analytic
noninteracting matter/radiation/Lambda limit. With a dot denoting proper time,

```
dot(rho_s) + 3 H rho_s = 0
dot(rho_d) + 3 H rho_d = -Gamma rho_d
dot(rho_r) + 4 H rho_r = +Gamma rho_d
dot(rho_Lambda) = 0
H^2 = (8 pi G/3) (rho_s + rho_d + rho_r + rho_Lambda)
p_total = rho_r/3 - rho_Lambda
```

These continuity equations are independently implemented from the equations in
Audren, Lesgourgues, Mangano, Serpico and Tram, *Strongest model-independent bound
on the lifetime of Dark Matter*, [arXiv:1407.2418v1](https://arxiv.org/abs/1407.2418v1),
Section 2.1, equations (2.1) and (2.2). There their derivative is conformal time and
the transfer is `a Gamma rho_d`; `dt=a d_eta` gives the proper-time equations
above. No CLASS code is copied. Unlike that source's present-fraction shooting
and asymptotic initial daughter boundary, this consumer takes **both initial
parent and daughter densities explicitly at a finite anchor** and only evolves
forward. Source perturbation and observational claims do not transfer to this
background implementation.

## Positive forward coordinates and outputs

For `x=ln(a/a_i)`, `tau=H_i (t-t_i)` and comoving daughter density
`R=(rho_r/rho_ci) exp(4x)`, the retained evolution is

```
rho_s/rho_ci = f_s exp(-3x)
rho_d/rho_ci = f_d exp(-3x-g tau)
rho_r/rho_ci = R exp(-4x)
rho_Lambda/rho_ci = f_Lambda
E = sqrt(sum of these four densities)
dtau/dx = 1/E
dR/dx = g f_d exp(x-g tau)/E
tau(0)=0, R(0)=f_r
```

Here `f_d` is the **initial** parent fraction. Parent survival is analytic;
daughter production is nonnegative. Classical RK4 stages in these two coordinates
have positive derivatives and cannot produce negative retained coordinates.
There are no density clamps, composition renormalizations or inferred past
boundary conditions.

Each row reports `E`, `H` in inverse seconds, elapsed seconds and `tau`, all four
density fractions and scaled densities, total `w=f_r/3-f_Lambda`, deceleration
`q=(1+3w)/2`, and `Q/(H rho_total)=g f_d(a)/E` with `Q=Gamma rho_d`. The continuity
residual checks the algebraic cancellation of transfer in the total derivative;
it is not an independent numerical conservation certificate. The survival
identity and independent time-coordinate reference supply stronger tests.

## Domain, error and resources

Requested scale factors must be finite, nondecreasing and lie in
`a_i<=a<=a_i exp(16)`; duplicates are retained in order. Empty batches are valid
and perform no evolution. The model supports `0<=g<=10^4`, `0<=tau<=10^6`, and,
when an initial parent exists, at most 600 elapsed parent lifetimes (`g tau<=600`).
These are explicit computational/representation bounds, not cosmological
restrictions. A trial stage outside these bounds also refuses the calculation.
Positive parent/daughter values that cannot be represented remain refusals,
never substituted zeros. Physical H and positive elapsed seconds must be normal
in the platform's `long double` arithmetic. Heterogeneous lifetime/rate ranges
are not covered by a uniform observational qualification.

The solver uses RK4 full-step/two-half-step differences divided by 15, retains
the two-half-step state and adjusts steps. Per-step budgets allocate the supplied
absolute plus relative state tolerance in proportion to log-scale step length
across the full requested span (at least one in the denominator). Defaults are
`1e-11` absolute and `1e-9` relative. Absolute tolerance must be positive; both
nonzero tolerances must exceed the platform arithmetic floor
`128 epsilon(long double)`. Round-to-nearest arithmetic is required. Public
estimates accumulate local `tau` and `R` differences; the expansion estimate
propagates those through the exact density projection and includes arithmetic
roundoff. These are empirical diagnostics, **not proven bounds**; accumulated
local differences do not bound nonlinear global errors. The tolerance controls
the two integration coordinates, not every arbitrarily tiny parent fraction.
Consumers should check refinement at their actual state and error allocation.

Default caps are one million RHS evaluations, 4096 rows and 4 MiB native payload.
Every attempted stage counts, including rejected steps. Before allocation or
work, the payload gate counts retained rows plus the result owner and an 8192-byte
operation allowance. It excludes borrowed input, caller-held copies, allocator
metadata and process RSS; copy/move lifetime is owned by the caller. Failed batches
withhold all rows and retain attempted-work counters and the supplied initial
state. `work_limit` reports a cap, `outside_domain` reports unsupported physical
or representation ranges, `nonfinite_input` preserves nonfinite refusals, and
`conditioning_budget_exceeded` reports an unresolved step. Invalid numerical
policy is `invalid_input`. No partial trajectory is usable as an accepted result.

## Permanent numerical evidence

[The native test](../cpp/tests/test_decaying_matter.cpp) compares zero-decay
expansion against analytic matter/radiation/Lambda densities and pure-species
proper-time integrals. It checks exact particle survival, positive daughter
emergence, monotone comoving daughter energy, late radiation domination, closure
and continuous transfer cancellation. Independent fixed-step RK4 evolves the
**physical densities and scale factor in proper time**, rather than the
production log-scale coordinates. Its 16384/32768-step refinement must consume
less than 5% of the named `2e-8` absolute/relative comparison allocation. Native
refinement separately changes both tolerances by a factor 100. Named source-free
synthetic cases cover decay rates `0.01`, `1`, `10`, mixed species and pure decay;
they do not qualify the entire admitted computational domain. Permanent refusal
and ownership checks cover invalid closure/rates/epochs, lifetime/scale/work/
payload/row caps, nonfinite inputs, duplicates, and retained copy/move lifetime.
