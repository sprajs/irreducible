# Conditional pure-hydrogen thermal, opacity and truncated drag history

The native C++20 SDK `irred/recombination_drag.hpp` prepares an owned, bounded
pure-hydrogen non-equilibrium history. It returns the electron fraction `x_e`,
a positive baryon-drag optical depth measured from a supplied late endpoint,
and a separately accepted unit-depth root. The compiled identity is
`pure-H-Peebles3level-F1-prescribed-Tm-equals-Tr-truncated-drag/v1`.

The default remains the prescribed effective three-level approximation with
`T_m=T_r=T_CMB,0(1+z)`. Explicitly selecting
`HydrogenTemperatureModel::evolved_compton_adiabatic` evolves matter temperature
and qualifies additional finite-endpoint Thomson opacity and visibility groups.
`model_identity()` and `method_identity()` declare the selected physical and
numerical profiles.
There is no helium, reionization, molecular physics, calibrated multilevel
correction, perturbation calculation, CMB likelihood or observational
qualification. Residual ionization and the finite late endpoint affect the
root: it is not a physical cosmological drag epoch. The existing supplied-drag
thermal distance/ruler and BAO consumers retain their separate route.

## Equations, atomic assets and source boundary

[Peebles 1968, equations 26, 30 and 31](https://articles.adsabs.harvard.edu/pdf/1968ApJ...153....1P)
and [Seager, Sasselov and Scott 1999 v2, equation 1, equation 3 and the unnumbered beta relation on page 4](https://arxiv.org/abs/astro-ph/9909275v2)
supply the effective three-level structure and hydrogen rate fit. This slice
uses `F=1`; it does not transfer their empirically calibrated `F=1.14` to a
pure-hydrogen model. The latter paper warns about the approximate history below
redshift 300, which is outside this history's domain.

For `T=T_m` (equal to `T_r` in the prescribed variant), `t=T/10000 K`, hydrogen nuclei density `n_H` and the same owned
background's physical Hubble rate `H`:

```
alpha = 1e-19 * 4.309 * t^(-0.6166)/(1+0.6703*t^0.5300) [m^3/s]
E21 = h*c/lambda_Lya, B2 = B1-E21
beta = alpha*(2*pi*m_e*kB*T/h^2)^(3/2)*exp(-B2/(kB*T)) [1/s]
K = lambda_Lya^3/(8*pi*H)
C = (1+K*Lambda*n_H*(1-x))/(1+K*(Lambda+beta)*n_H*(1-x))
dx/dz = C*(n_H*alpha*x^2-beta*(1-x)*exp(-E21/(kB*T)))/(H*(1+z))
R = 3*Omega_b/(4*Omega_gamma*(1+z))
tau_drag(z;z_late) = integral[z_late,z] c*sigma_T*n_H*x/(H*(1+z)*R) dz
```

The excited-state photoionization threshold in `beta` is `B2`, following
Peebles equation 26. The additional ground-to-excited exponential makes the
stationary balance use `B1`. The explicit `B2=B1-h*c/lambda_Lya` energy closure
is an approximation, not an unchanged reproduction of every original atomic
asset. `lambda_Lya=121.5682 nm` and `Lambda=8.22458 /s` follow the 1999 source.

The shared [hydrogen equilibrium](hydrogen-equilibrium.md) constants fix
`m_e=9.1093837139e-31 kg` and `B1=13.598434599702 eV` to CODATA 2022 and NIST
ASD 5.12. This history also fixes `m_p=1.67262192595e-27 kg` and
`sigma_T=6.6524587051e-29 m^2` to [CODATA 2022](https://physics.nist.gov/cuu/pdf/all.pdf).
Exact SI `h`, `kB`, `c` and eV definitions come from `quantities.hpp`.
Atomic uncertainties, the rate fit, energy closure and prescribed temperature
are separate physical approximations, excluded from numerical diagnostics.
Original equations are implemented here; no external solver code or table is
copied into the library or permanent reference.

The physical-density source is mapped once with the thermal background API.
Pure hydrogen means every mapped baryon belongs to a hydrogen nucleus, with
fixed effective rest mass `m_H=m_p+m_e-B1/c^2`:

```
n_H,0 = rho_critical(H100)*omega_b_physical/m_H
n_H(z) = n_H,0*(1+z)^3
```

A private electron quantum-density equation is shared with the Saha consumer,
preserving its original scalar arithmetic and phase-space prefactor. The
private critical-density helper shares the original thermal background's
fixed `G=6.67430e-11 SI` and arithmetic. Binding-energy changes do not feed back
into that pressureless background. `R`, photons, `H` and `n_H` use this same
physical model. A ground-state Saha calculation supplies only the initial
condition through the shared atomic SDK; subsequent points solve the
non-equilibrium equation. Its translational mass-ratio approximation remains
one, as declared by the atomic consumer.

## Admission and retained ownership

The prescribed slice admits `1550<=z_initial<=1650`; the evolved slice admits
`1550<=z_initial<=1600`. The latter is a bounded numerical qualification scope,
not a physical boundary. Both admit
`300<=z_late<=600`, `60<=H0<=80 km/s/Mpc`, physical baryon density
`0.015<=omega_b<=0.03`, physical CDM density `0.08<=omega_cdm<=0.15`,
`2.70<=T_CMB,0<=2.75 K` and other physical massless density from zero to
`3e-5`. At most 16 explicit thermal species are admitted. Zero-mass species
retain the thermal mapper's existing finite-positive temperature/weight and
flat-closure admission. Up to three species can instead have explicit positive
mass `0<m<=0.3 eV`, present temperature `1.8<=T_nu,0<=2.0 K` and populated-state
weight `1<=g<=2`. These are bounded numerical qualification profiles, not a
mass-sum, hierarchy or effective-species mapping. An explicit species already
contributes its complete energy density; `physical_massless_nonphoton_density`
must not count that species again. No alternative expansion function or
radiation-free sound/ruler admission is introduced.

`prepare_pure_hydrogen_history(request, policy)` owns the source, one prepared
thermal state and history mesh. Copies retain independent owned storage;
moves invalidate the source object. Preparation/query allocation failures return
`work_limit`; ordinary value copies can throw `std::bad_alloc`. Copy assignment
constructs a complete replacement first, preserving the destination if allocation
fails. Ordered `evaluate(redshifts, mask)` batches
request independently masked scalar groups listed below. Invalid
redshifts fail individually, preserving valid rows. An unrequested group has no
value. Evaluation performs interpolation, with no new background or ODE work.
The unit-depth root has its own status; its numerical refusal does not erase
accepted history coordinates.

## Numerical profile and resource accounting

Strict floating-point compilation, `FE_TONEAREST` and long double with at least
64 mantissa bits and exponent range 16384 are required. The prescribed solver uses
backward Euler with the unique monotone scalar implicit root, safeguarded
Newton/bisection and a fixed 80-iteration ceiling. Three meshes at `N`, `2N`
and `4N` support Richardson extrapolation. Nonpositive or super-unit fractions
are refused; there is no physical-state clipping. Positive trapezoid integration
sets `tau(z_late)=0`. Linear query/root interpolation includes local curvature
and propagated history/background estimates.

Default accepted numerical allocations are `1e-8+1e-6*abs(x)` for electron
fraction, `2e-7+1e-6*abs(tau)` for optical depth, and `0.002` redshift for the
unit-depth root. Refinement and arithmetic diagnostics are empirical estimates,
not universal rigorous bounds and not model-uncertainty estimates. A tighter
requested allocation can refuse one output while preserving another.

The default base mesh has 8192 intervals; the finest mesh has a hard ceiling of
65536. Combined work counts initial equilibrium solves, physical background
evaluations, nested momentum callbacks and nonlinear right-hand-side calls,
including failed attempts. Its default ceiling is four million. The default
output cap is 4096 points, hard maximum 65536. Explicit payload bounds cover
retained sources/background, coefficient/history meshes, temporary mesh buffers
and output capacities before allocation. The default cap is 16 MiB and hard
maximum 1 GiB; allocator metadata and process RSS are excluded.

## Independent controls

Permanent tests contain original direct-SI equations with an independent
resolved explicit RK4 algorithm, without engine headers or production
readback. The constants and physical closure ancestry are shared and declared;
algorithm and arithmetic implementation are independently authored. Both routes
use long double and the same system `libm`; this is an algorithm/refinement
control, not an independent high-precision arithmetic certificate. The separate
Decimal controls of the atomic consumer qualify Saha fractions, not this whole
history. Steps `0.01` and
`0.005` redshift must refine to within five percent of each frozen named
allocation, maintain stage positivity and have a resolved stiffness diagnostic.
Controls cover three distinct physical models, initial and late endpoint
changes, finite-difference ODE residuals, numerical mesh refinement, supplied
Saha departure, masks, source lifetime, quotas and invalid closures. These are
synthetic model controls, not additional measured data or a validation of the
physical closure against the real universe.


## Evolved matter temperature and finite-endpoint photon measure

The evolved identity is
`pure-H-Peebles3level-F1-SSS1999-Tm-rates-Compton-adiabatic-finite-endpoint-Thomson/v1`.
It uses the literal source's `T_m` throughout `alpha`, `beta` and the
hydrogen ground-to-excited exponential. `T_r` enters the radiation energy
and the Compton target. This convention is explicit; it is not a mixed
radiation-temperature photoionization scheme or calibrated RECFAST/HyRec.
Equation **5** of SSS1999v2 supplies the temperature closure; equation 4 is a
helium rate and equation 6 is HeII Saha, neither the temperature nor beta source.
For pure hydrogen `x_e=x_p=x`, `f_He=0`:

```
u = 1+z, T_r = T_CMB,0*u
u_gamma = Omega_gamma*rho_critical_energy(H0)*u^4 [J/m^3]
Gamma_C = 8*sigma_T*u_gamma/(3*m_e*c)*x/(1+x) [1/s]
dT_m/dz = (Gamma_C*(T_m-T_r)+2*H*T_m)/(H*u) [K per redshift]
q_T = c*sigma_T*n_H*x/(H*u) [per redshift]
tau_T(z;z_late) = integral[z_late,z] q_T dz
g_z = q_T*exp(-tau_T) [per redshift]
S(z;z_late) = exp(-tau_T)
```

The radiation energy reuses the thermal mapper's photon state. Its exact SI
identity is `a_R=8*pi^5*kB^4/(15*h^3*c^3)`, without a separately rounded
radiation constant. Mapped photon rounding and H estimates enter the temperature
diagnostic separately. The heat equation includes only adiabatic expansion and
Compton coupling; other heating, cooling and a variable-particle-number energy
term are excluded. The boundary is one shared Saha fraction with
`T_m(z_initial)=T_r(z_initial)`, rather than equilibrium imposed throughout.
The finite endpoint excludes all later ionization and reionization.

A finite-cell positive linear opacity has an exact quadratic integral for
`tau_T`. Both query `q_T` and `g_z` use that same cell definition, giving
`integral g_z dz = 1-S(z_initial)` plus surviving boundary mass
`S(z_initial)`. Visibility is not renormalized. This per-redshift measure is
not a distribution per conformal time, today's optical depth or a full CMB
last-scattering prediction. Drag retains the separate `q_T/R` equation and
endpoint-dependent unit-depth root.

| Mask | Row group | Units |
|---|---|---|
| 1 | `electron_fraction` | dimensionless |
| 2 | `drag_depth` | dimensionless |
| 4 | `matter_temperature_kelvin` | K |
| 8 | `thomson_depth` | dimensionless |
| 16 | `thomson_opacity_per_redshift` | per redshift |
| 32 | `visibility_per_redshift` | per redshift |
| 64 | `finite_endpoint_survival` | dimensionless |

The prescribed variant exposes its explicit radiation temperature under mask4;
its new photon groups are individually `outside_domain`, preserving accepted
legacy groups. Unknown mask bits and an empty mask are invalid. Each requested
group has its own optional scalar, status and absolute empirical error estimate.

The coupled solver eliminates the backward-Euler temperature equation at each
trial fraction. With decreasing-redshift step `h`, `A=Gamma_C/(H*u)`:

```
T_new = (T_previous+h*A(x_new)*T_r,new)/(1+h*A(x_new)+2*h/u_new)
R_x = x_new-x_previous+h*f(x_new,T_new)
J = 1+h*(f_x+f_T*dT_new/dx_new)
```

A safeguarded bracket and at most80 Newton/bisection trials check the effective
Jacobian and both equation residuals. No global coupled monotonicity/uniqueness
proof is asserted. Unresolved, ill-conditioned or nonphysical steps refuse;
there is no clipping or Saha substitution. Nested quadratic-redshift meshes
`z_i=z_initial-(z_initial-z_late)*(i/N)^2` resolve the initial stiff Compton
boundary. Richardson states and interpolation must retain `0<x<1` and
`0<T_m<=T_r`. The prescribed uniform mesh and scalar output bits retain their
original arithmetic identity.

Additional accepted allocations are `2e-5 K+2e-7*abs(T_m)` for temperature,
`2e-7+1e-6*abs(tau_T)` for Thomson depth, `2e-9+2e-6*abs(q_T)` per redshift
for opacity, `2e-9+3e-6*abs(g_z)` per redshift for visibility and
`1e-5*abs(S)` for survival. The finite visibility-mass allocation is `2e-7`
absolute. Opacity propagates fraction and H errors. Partial-cell depth errors integrate
nonnegative endpoint opacity errors with the same positive primitive weights;
cumulative refinement and curvature allowances are retained. Exponential diagnostics
use `e_S=S*expm1(e_tau)` and
`e_g=S*(e_q+(abs(q)+e_q)*expm1(e_tau))`, plus arithmetic. A tiny positive
value/product/diagnostic cannot silently become zero; projection refuses if it
is unrepresentable or exceeds its own allocation.

Evolved preparation counts both ionization and temperature RHS evaluations,
including failed trials. Actual enlarged row/header sizes and evolved auxiliary
storage are charged by the payload helper; tight old byte caps may now refuse.
The default finest32768 grid fits the default16MiB preparation cap; finest65536
requires an explicitly larger cap. These are payload bounds, not RSS estimates.

Independent thermal controls use original direct-SI equations and a separately
authored two-stage order3 L-stable RadauIIA algorithm in logarithmic expansion,
with stage Newton/finite-difference Jacobian, resolved initial Compton layer,
positive stages and no clipping. Reference refinement must use at most five
percent of every named downstream allocation. Constant-coupling/adiabatic
analytic limits and independent opacity quadrature challenge temperature signs,
units, finite-interval visibility mass and boundary survival. Both compiled
routes use long double/system libm with declared shared constants; this is
algorithm/refinement evidence, not an independent high-precision certificate.
Atomic, rate-fit and truncated thermal-model errors remain excluded from these
numerical allocations.


## Explicit massive thermal-relic background composition

A positive-mass source selects `pure_hydrogen_relic_history_id` or
`evolved_hydrogen_relic_history_id`. All-zero-mass sources keep the original
model IDs and scalar arithmetic; the prescribed and evolved ODE method IDs,
atomic assets, source bounds, initial Saha/temperature conditions and every
numerical allocation remain unchanged. The new identities qualify only the
composition of this bounded hydrogen closure with explicit collisionless relics.

The same retained [thermal background](thermal-neutrino.md) owns the momentum
integral, present normalization and flat Lambda closure. With
`q=p*a/T_nu,0` in the declared natural-unit convention,

```
I_rho(y) = integral[0,infinity] q^2*sqrt(q^2+y^2)/(exp(q)+1) dq
rho_nu(a) = g*T_nu,0^4*I_rho(m*a/T_nu,0)/(2*pi^2*a^4)
```

This collisionless phase-space law follows the explicit momentum formulation
in [Lesgourgues and Tram 2011, equation 2.1](https://arxiv.org/abs/1104.2935).
The production history calls the already qualified retained `scaled_expansion`
provider at each finest-mesh node. It reuses those physical coefficients across
the three hydrogen meshes; it introduces no expansion table or replacement law.
Relics change `H`, including the source-defined present flat closure. They do
not become hydrogen nuclei or photon energy, and do not change `R`, Saha,
the pure-hydrogen heat capacity or the atomic rates. Provider H estimates
continue to propagate into hydrogen, temperature and opacity diagnostics.

The default four-million combined work cap remains unchanged. An explicit
caller cap is required for the qualified positive-mass profiles. For the
synthetic `H0=67.4`, `omega_b=.0224`, `omega_cdm=.12`, `T_CMB=2.7255`, other
massless density zero, species `(m,T0,g)=(.06,1.95,2),(0,1.95,2),(0,1.95,2)`
and evolved endpoints1600/300, the default cap refused at exactly4,000,000
operations after733 background samples. The unchanged default8192 base mesh
with an explicit500,000,000 cap accepted with178,512,172 operations, including
178,135,936 momentum callbacks,32769 background evaluations,171733 ionization
and171733 temperature RHS evaluations and one Saha solve. The measured single
GCC16.2.1 Release pilot took15.25 seconds; this is a pinned case cost, not a
portable runtime guarantee. Failed provider calls contribute their actual work;
queries perform no new momentum or ODE calls. The existing payload accounting,
16MiB default,1GiB hard cap and output masks remain in force.

Separate original direct-SI momentum and Radau controls use composite GL16/32
on `[0,64]` with an exponential tail control, independent of production adaptive
momentum integration. Reference stage coefficients may be reused at their own
stage redshift during nonlinear iteration. Named reference gates include the
light/heavy transition, standard-like one-massive plus two-massless profile,
three positive species, split populated-state weights, initial boundary and
fixed-physical-density H0 behavior. The H comparison allocation is2e-10 relative;
reference refinement consumes at most five percent of each unchanged downstream
allocation. These controls share fixed SI/atomic input ancestry and long-double
system arithmetic, and qualify numerical composition rather than full
recombination, a physical drag epoch, CMB inference or observations.
