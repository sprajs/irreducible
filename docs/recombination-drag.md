# Conditional pure-hydrogen recombination and truncated drag history

The native C++20 SDK `irred/recombination_drag.hpp` prepares an owned, bounded
pure-hydrogen non-equilibrium history. It returns the electron fraction `x_e`,
a positive baryon-drag optical depth measured from a supplied late endpoint,
and a separately accepted unit-depth root. The compiled identity is
`pure-H-Peebles3level-F1-prescribed-Tm-equals-Tr-truncated-drag/v1`.

This is a conditional effective three-level approximation with prescribed
matter temperature `T_m=T_r=T_CMB,0(1+z)`. It does not evolve matter temperature.
There is no helium, reionization, molecular physics, calibrated multilevel
correction, perturbation calculation, CMB likelihood or observational
qualification. Residual ionization and the finite late endpoint affect the
root: it is not a physical cosmological drag epoch. The existing supplied-drag
thermal distance/ruler and BAO consumers retain their separate route.

## Equations, atomic assets and source boundary

[Peebles 1968, equations 26, 30 and 31](https://articles.adsabs.harvard.edu/pdf/1968ApJ...153....1P)
and [Seager, Sasselov and Scott 1999, equations 1–3 and 6](https://arxiv.org/abs/astro-ph/9909275)
supply the effective three-level structure and hydrogen rate fit. This slice
uses `F=1`; it does not transfer their empirically calibrated `F=1.14` to a
pure-hydrogen model. The latter paper warns about the approximate history below
redshift 300, which is outside this history's domain.

For `T=T_m=T_r`, `t=T/10000 K`, hydrogen nuclei density `n_H` and the same owned
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

The first bounded slice admits `1550<=z_initial<=1650`,
`300<=z_late<=600`, `60<=H0<=80 km/s/Mpc`, physical baryon density
`0.015<=omega_b<=0.03`, physical CDM density `0.08<=omega_cdm<=0.15`,
`2.70<=T_CMB,0<=2.75 K` and other physical massless density from zero to
`3e-5`. At most 16 explicit thermal species are admitted, all with zero mass;
the existing thermal mapper enforces their temperature/weight and flat closure
domains. Massive relics are explicitly refused. No alternative expansion
function or radiation-free sound/ruler admission is introduced.

`prepare_pure_hydrogen_history(request, policy)` owns the source, one prepared
thermal state and history mesh. Copies retain independent owned storage;
moves invalidate the source object. Ordered `evaluate(redshifts, mask)` batches
request `hydrogen_electron_fraction`, `hydrogen_drag_depth` or both. Invalid
redshifts fail individually, preserving valid rows. An unrequested group has no
value. Evaluation performs interpolation, with no new background or ODE work.
The unit-depth root has its own status; its numerical refusal does not erase
accepted history coordinates.

## Numerical profile and resource accounting

Strict floating-point compilation, `FE_TONEAREST` and long double with at least
64 mantissa bits and exponent range 16384 are required. The solver uses
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
algorithm and arithmetic implementation are independent. Steps `0.01` and
`0.005` redshift must refine to within five percent of each frozen named
allocation, maintain stage positivity and have a resolved stiffness diagnostic.
Controls cover three distinct physical models, initial and late endpoint
changes, finite-difference ODE residuals, numerical mesh refinement, supplied
Saha departure, masks, source lifetime, quotas and invalid closures. These are
synthetic model controls, not additional measured data or a validation of the
physical closure against the real universe.
