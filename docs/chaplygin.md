# Generalized Chaplygin unified barotropic fluid

The C++ SDK supplies a flat-GR homogeneous generalized Chaplygin fluid, with
separately conserved ordinary pressureless matter and isotropic perfect-fluid
radiation. It retains background composition, equation of state, a formal
barotropic derivative and optional flat distances from the same expansion law.
Ordinary matter is a separately supplied component, not a renamed total dark
matter density. The unified fluid has its own physical identity and is not the
thermal background, quintessence or the nominal effective-fluid model.

Use `irred/chaplygin.hpp` and `evolve_chaplygin`; the complete
[SDK example](../cpp/examples/chaplygin.cpp) requests expansion and distances.
The identity is
`flat-GR-generalized-Chaplygin-unified-barotropic-fluid-ordinary-dust-radiation/v1`.
There is no ABI/CLI route, perturbation evolution, ruler, CMB closure or inference
qualification. Its formal `dp/de` output is not a collisionless signal speed or
a qualified rest-frame propagation law, especially on the constant-density
vacuum endpoint.

## Physical law, inputs and source

With energy density `e` and pressure `p` in the same units, the barotropic law is
`p=-A e^(-alpha)`. Thus `A` has units `[energy density]^(1+alpha)`. Normalize at the
chosen observer epoch `a=1`, with
`A_s=A/e_anchor^(1+alpha)`. Supply `0<=A_s<=1`, `0<=alpha<=1`, a positive anchor
Hubble normalization in **km/s/Mpc**, and nonnegative ordinary-matter/radiation
fractions at that same epoch. The computational H range is `[1,1000] km/s/Mpc`.
It is the supplied finite observer anchor, not an inferred cosmological epoch.

The engine owns the remaining density fraction canonically:
`Omega_cg=1-(long double(Omega_ordinary)+Omega_radiation)>0`. It rejects a zero or
negative component rather than reporting an absent fluid as if its EOS were
measured. No input fraction is renormalized or clamped. This rule uses the actual
stored binary64 inputs: decimal `.95` and `.05` have a wide sum slightly below
one and therefore leave a tiny positive unified-fluid fraction; exact binary
`.75` and `.25` leave zero and are refused. Actual `E(1)` is obtained from the
represented density sum, not forced to one after a closure calculation.
Normalized energy densities refer to `3 H_anchor^2 c^2/(8 pi G)`; no numerical
`G` is needed for the emitted ratios.

Continuity `de/dln(a)=-3(e+p)` gives

```
d[e^(1+alpha)]/dln(a) = -3(1+alpha) [e^(1+alpha)-A]
F(a) = e(a)/e_anchor
     = [A_s + (1-A_s) a^(-3(1+alpha))]^(1/(1+alpha))
w_cg = -A_s / [A_s + (1-A_s) a^(-3(1+alpha))]
dp/de = -alpha w_cg
E(a)^2 = Omega_ordinary a^-3 + Omega_radiation a^-4 + Omega_cg F(a)
```

These are independently implemented from Bento, Bertolami and Sen,
*Generalized Chaplygin Gas, Accelerated Expansion and Dark Energy-Matter
Unification*, [gr-qc/0202064v1](https://arxiv.org/abs/gr-qc/0202064v1), equations
(1), (12) and (13). Original PDF equation access was verified; the source file
SHA256 is `0e246430a70015a450f4994cecef6f748272295f9e7a77523358f9419362aa04`.
The paper states `A>0` and `0<alpha<=1`. Native `A_s=0` and `alpha=0` are
separately declared mathematical endpoints extending that source domain.
`alpha=1` gives the original Chaplygin law; `A_s=1` is constant density.
Born-Infeld/microphysics and perturbation claims from the source are not imported
into this background-only consumer. No paper text or third-party code is copied.

## State and flat-distance outputs

Each row reports `E`, H in km/s/Mpc, `F`, all three instantaneous component
fractions, `w_cg`, an independently positive expression for `1+w_cg`, the formal
barotropic slope, total `w=Omega_cg(a) w_cg+Omega_radiation(a)/3`, and
`q=(1+3w)/2`. The enthalpy fraction is computed from the positive evolving term,
not by subtracting a nearly vacuum `w` from one.

Exact `A_s=0` is dust, and exact `A_s=1` is vacuum, handled by their analytic
expressions. General inputs use a log-sum-exp of the positive vacuum and evolving
terms, followed by their normalized fractions. `F(1)=1` follows the defining
anchor normalization and does not repair ordinary/radiation input closure.

Request `chaplygin_distances` to obtain

```
chi(a) = (c/H_anchor) integral_a^1 db/(b^2 E(b))
D_A(a) = a chi(a), D_L(a) = chi(a)/a
```

These use the same compiled E state, the shared speed-of-light definition and
flat homogeneous geometry. They are synthetic model distances, with no source,
instrument or observational calibration. Unrequested distances are absent,
not zero-valued placeholders. At the observer endpoint all three distances are
exactly zero by their integration interval.

## Numerical, domain and ownership contract

Queries proceed from observer toward the past: finite nonincreasing `a` in
`[1e-4,1]`. Duplicates are retained in their input order and reuse the accumulated
distance. This is a bounded past-anchor calculation; future queries are explicit
refusals. The scale range bounds the largest logarithmic exponent at about55.3
for `alpha=1`; log-sum arithmetic prevents direct-term overflow or cancellation.
Very small positive EOS/enthalpy/components that cannot be represented remain
refusals, never replaced with endpoint zeros. States and physical projections
must remain normal in platform `long double`; distance callbacks must remain
normal in binary64. Round-to-nearest is required.

Distances accumulate positive segment integrals in `x=ln(a)`, with callback
`exp(-x)/E(exp(x))`, using shared adaptive Simpson quadrature. Defaults are
`1e-7 Mpc` absolute distance tolerance and `1e-9` relative tolerance, maximum
quadrature depth30, one million **quadrature callbacks**, 4096 rows and 4 MiB
native payload. Analytic background row evaluations do not count as quadrature
callbacks; they are bounded by the admitted row count. Background-only requests
perform no quadrature and may use a zero callback cap. No callback budget is
silently enlarged. All attempted quadrature callbacks, including refused
attempts, count against the single batch cap.

A positive absolute distance tolerance is required; relative tolerance must be
at least `128 epsilon(binary64)`. Each segment receives a fraction of the absolute
allocation by log-scale length, divided by eight and scaled by the smallest
requested `a` so even `D_L=chi/a` fits its own Mpc allocation. Its relative
allocation is divided by eight. Retained distance estimates add quadrature,
callback arithmetic and background arithmetic diagnostics and must pass the
requested absolute-plus-relative gate for each emitted distance. The background
relative diagnostic must also fit the requested relative tolerance. These are
empirical estimates, **not proven error bounds**. They exclude input-parameter,
stress-closure and observational uncertainties. Consumers should check
refinement at their actual inputs and output error allocation.

The result retains the complete supplied model, derived anchor fluid fraction,
requested-output mask, rows and callback count. No borrowed span or callback
survives the call. Before numerical work the payload gate counts the result
header, retained rows and a 32768-byte operation/recursive-kernel allowance.
Borrowed inputs, caller-held copies, allocator metadata and RSS are excluded.
Empty batches validate the initial state and perform no numerical work.

Any refusal withholds **all rows** while retaining model and attempted-work
metadata. Nonfinite inputs yield `nonfinite_input`; physical/order/representation
ranges yield `outside_domain`; row/payload/callback/depth exhaustion yields
`work_limit`; invalid numerical/output policy yields `invalid_input`; an
unresolved diagnostic allocation yields `conditioning_budget_exceeded`. A failed
batch is not a partially accepted trajectory.

## Independent permanent controls

[Native tests](../cpp/tests/test_chaplygin.cpp) check exact `A_s=0` dust and
`A_s=1` vacuum expansion/distances for several alpha values, and `alpha=0`
equivalence to a separately conserved dust-plus-Lambda density decomposition.
An independent fixed-step RK4 evolves `F_x=-3(F-A_s F^-alpha)` from the anchor,
instead of using the production closed solution. An independently derived
four-node Gauss-Legendre rule supplies distances instead of adaptive Simpson.
Density steps refine8192→16384 and quadrature cells128→256; each reference
refinement must occupy less than5% of its named `2e-8` comparison allocation.
Native distance tolerances separately refine by a factor100.

Four critical-density fixtures near `A_s=0`, near `A_s=1` at both sides of its
transition, and fractional alpha use direct-algebra mpmath1.3.0 at100/160 decimal
digits with exact binary64 integer-ratio inputs. The independent precision
refinement is below `1e-90` scaled; permanent comparisons use `3e-14` relative
for `F`, small negative `w`, positive `1+w` and E. Fixture ancestry and literals
are retained in the test. These are numerical controls, not external
cosmological truth. Named mixture cases span `A_s=.01,.7,.999` and
`alpha=0,.5,1`; they do not qualify the entire admitted computational domain.
Other controls retain duplicates, absent outputs, copy/move lifetime, tiny
positive canonical fluid fractions and invalid model/query/resource policies.
