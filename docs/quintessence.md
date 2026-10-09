# Canonical exponential quintessence background

`irred/quintessence.hpp` supplies a retained, flat GR background with a minimally
coupled canonical scalar, separately conserved pressureless dust and massless
radiation. This is a field initial-value model, with identity
`GR/flat-canonical-exponential-quintessence-dust-radiation-IVP/v1`, rather than a
prescribed CPL equation of state. It supplies no perturbations, scalar sound
horizon, growth, recombination, likelihood, parameter posterior or observational
qualification. The consumer is C++ SDK only.

With reduced Planck mass `Mpl`, the potential is
`V(phi)=Vstar exp(-lambda phi/Mpl)`, and the source equations are

```text
3 Mpl² H² = rho_m + rho_r + phi_dot²/2 + V
phi_ddot + 3 H phi_dot + dV/dphi = 0
rho_m ∝ a^-3; rho_r ∝ a^-4
```

The scalar is canonical in the GR cosmological proper-time frame. Define
`N=ln a`, `x=phi_dot/(sqrt(6) Mpl H)`, `y=sqrt(V)/(sqrt(3) Mpl H)` and
`Q=3x²+(3/2)Omega_m+2Omega_r`. The one compiled evolution owner uses

```text
x'       = -3x + sqrt(3/2) lambda y² + x Q
y'       = y [-sqrt(3/2) lambda x + Q]
Omega_m' = Omega_m [-3 + 2Q]
Omega_r' = Omega_r [-4 + 2Q]
(ln E)'  = -Q
```

The prime means `d/dN`. `Omega_phi=x²+y²` and
`w_phi=(x²-y²)/(x²+y²)`. Separate dust/radiation evolution preserves their exact
zero boundaries; the code monitors `|x²+y²+Omega_m+Omega_r-1|`. It never clamps
negative densities, renormalizes closure, or modifies supplied fractions.

The original [Copeland, Liddle and Wands scaling analysis](https://arxiv.org/abs/gr-qc/9711068)
provides the exponential-potential autonomous variables and single-barotropic
fluid fixed points (equations 6–10 and Table I). The simultaneously retained
dust/radiation equations above follow their separate GR conservation laws;
they are a distinct mixture, not a source claim about one fixed fluid. The
focused source inspection used version 2, PDF SHA256
`520595adc32264c9c0c2c9837929b25c340aa53eafd6ccaf0507a064bdf5df8c`.

## Initial state and supported domain

The supplied anchor is exactly `a=1`. `h_anchor_km_s_mpc` is its positive Hubble
rate, at most `1e6 km/s/Mpc`; `E=H/H_anchor`. The signed kinetic root is `x`,
`potential_fraction=y²` and `radiation_fraction=Omega_r`. These finite supplied
fractions must obey `x²+y²+Omega_r<=1`, with nonnegative potential/radiation.
The remaining fraction is dust. `lambda` is finite in `[-20,20]`. Inputs are
initial conditions, not fitted cosmological parameters or an implied attractor.
`Vstar` and a field origin need not be separately specified because only their
initial potential product affects this background. A zero potential is the
explicit kinetic/no-potential boundary; it is not a strictly positive
exponential potential at finite field value.

Requested scale factors lie in `[exp(-8),exp(4)]`. Both evolution directions
start from the same immutable anchor. Backward scalar evolution can be
sensitive to initial conditions; positivity, closure, arithmetic or work gates
may refuse an otherwise syntactically admitted trajectory. No attractor or
present-day normalization is imposed. Every input row retains its original
position, even after a refusal or exhausted batch budget.

Each row reports `E`, dimensional `H`, `w_phi`, and scalar/dust/radiation
fractions. Exactly zero or numerically unresolved scalar energy withholds
`w_phi` with `singular` status while retaining the other outputs. Optional flat
radial and luminosity distances use the same evolved `ln E` and integrate
`dJ/dN=exp(-N)/E`. For `a<=1`, `D_M=-c J/H_anchor` and `D_L=D_M/a`, with
`c=299792.458 km/s`. Future (`a>1`) distances are refused independently while
background outputs survive. Distances use the usual `1+z=1/a` relation for this
flat GR observer; no frame conversion or source effects are supplied.

## Numerical contract and evidence

The portable CPU implementation requires nearest rounding and `long double`
with at least 64 mantissa bits. It uses classical RK4 full versus two half
steps; the accepted state is the two-half-step result, without an implicit
projection or Richardson state correction. Each attempt counts all twelve RHS
callbacks, including rejected attempts. Per-row and whole-batch callback caps,
halving limits, checked row-storage accounting and fixed stack workspace bound
work/storage. Defaults are not raised on failure. The payload accounting is a
requested object-storage bound, not whole-process RSS or allocator accounting.

Policy absolute/relative tolerances allocate **local** step-doubling and
arithmetic diagnostics in the integration coordinates and separately admit the
Friedmann residual. The local allowance is
`(abs + rel max(|old_i|,|new_i|)) |step|/|target_N|`. Reported
`absolute_error_estimate` values transform accumulated accepted local
step/arithmetic diagnostics at the endpoint. They are not propagated global
state-error intervals and do not certify final observable accuracy under
amplification. A policy tolerance is not a promise of an absolute error in Mpc
or an end-to-end likelihood budget. Consumers must earn their own final-output
comparison gate and input-sensitivity contract. Singular denominators and
unsupported arithmetic refuse rather than acquiring an artificial tolerance.

The permanent native test uses a frozen `3e-8*(1+|reference|)` analytic
comparison allocation for the zero-velocity `lambda=0` Lambda/dust/radiation
limit, both signs of pure kinetic stiff matter, scalar-dominated exponential
power-law, dust scaling and radiation scaling. Exact stiff distances and
independent composite Simpson Lambda distances challenge the distance state.
A transient case instead integrates dimensional Klein–Gordon variables
`q=phi_dot/(Mpl H_anchor)` and `U=V/(3Mpl² H_anchor²)` with algebraic Friedmann
matter/radiation densities. Its final `E`, `w_phi` and distance comparisons use
`2e-7*(1+|reference|)`, with 4096→8192 reference refinement below one percent of
that allocation. This independent formulation and producer refinement check
named initial states; they do not qualify all allowed trajectories. The tests
also preserve field-reflection symmetry, zero scalar energy, invalid states,
row order, partial availability, capped work, storage overflow and impossible
arithmetic allocations.

The standalone installed consumer is [cpp/examples/quintessence.cpp](../cpp/examples/quintessence.cpp):

```cpp
#include <irred/quintessence.hpp>
using namespace irred::cosmology;
auto model = prepare_quintessence({1.2, 67.4, .08, .65, .00009});
const double a[]{1., .8, .5};
auto result = model.evaluate(a, true);
```

Inspect batch and per-output statuses before using optional values. The
example initial state is a declared numerical control, not an observational
fit or a Planck state.
