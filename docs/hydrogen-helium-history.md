# Bounded shared-charge H/He singlet kinetic history

The native history owner `irred/hydrogen_helium_history.hpp` is a separately
identified approximation for H I/H II and He I/He II with one free-electron
density and evolved matter temperature. It uses an independently supplied
hydrogen and helium nuclei inventory and the same retained massless thermal
background for every equation. It supplies fractions, shared electron density,
matter temperature and Thomson opacity per redshift. There is no CLI/C ABI,
physical drag epoch, full recombination or CMB calculation.

## Frozen physical source contract

The model identity is
`HII-HeII-singlet-RecfastCLASS-convention-NIST-central-FH1-TphotoTm-Compton-adiabatic-bounded/v1`.
This is an original compiled equation implementation, with shared modern atomic
central energies and explicitly chosen source-code conventions. It is not a
literal reproduction of SSS1999 or the calibrated/default RECFAST model.

The retained paper is
[Seager, Sasselov and Scott 1999 v2](https://arxiv.org/pdf/astro-ph/9909275v2),
SHA256 `fe4232be14d12ccad7ad46f795dd21c4fa1f66d6aba18c6e512c4341cc70bc07`.
The separately reviewed source is
[RecfastCLASS at commit 64bbab707faf4de4779a9e04edd180fef18d98fa](https://github.com/lesgourg/class_public/tree/64bbab707faf4de4779a9e04edd180fef18d98fa/external/RecfastCLASS),
with these immutable source identities:

| Source member | SHA256 |
| --- | --- |
| `wrap_recfast.c` | `e5ea80edd4652eb1a752cb4a2b835f4b74db91ca90248e5b5d04a8a3eada7cd5` |
| `wrap_recfast.h` | `1e5609a3e5b1ddd8f7ad99df336311a9f769867e52f203e0cb6481b617e2de73` |

The source inspection and restart receipts remain local and ignored. No source
solver, switching code, fitted correction or external table is copied into the
library or permanent reference. Scalar central facts and equations are stated
below. Copying external software would require a separate licensing review.

There are unresolved differences between the paper and the code:

- Paper equation 2 prints a negative He singlet-splitting exponent; the pinned
  code uses a positive exponent. The new model explicitly selects the code
  convention. An independent Boltzmann/Sobolev derivation supports that sign;
  this does not establish an author erratum.
- The code supplies the He photoionization factor four and excited-state
  binding threshold `chi-E_2s`; the paper's unnumbered beta relation omits
  the factor four and its `nu_2s` notation is ambiguous when read with the
  preceding transition definitions. The selected law has the existing shared
  ground-state Saha statistical convention as its stationary balance.
- Code `T_0=10^0.477121 K` is explicitly retained rather than paper `3 K`.
  Hydrogen retains the existing `F=1`, `Lambda_H=8.22458/s` and Ly-alpha
  wavelength, rather than paper `F=1.14` or code `8.2245809/s`.
- The shared NIST ionization energies replace the source's older ionization
  wavenumbers. He excitation wavenumbers retain the pinned source values.
  The resulting excited-state energy closure is its own approximation.
- The pinned header labels the helium line wavenumbers as `eV` in comments;
  the source's `h*c*L` arithmetic establishes inverse metres. That source
  comment error is preserved in this ledger, not adopted as a unit convention.
- The full wrapper switches to additional helium physics at its late switch
  and has other calibrated corrections. The isolated singlet closure omits
  those switches and does not claim full-wrapper output fidelity.

## Equations and assets

Let `p=n_HII/n_H`, `q=n_HeII/n_He`, `u=1+z`, and
`n_H=n_H0*u^3`, `n_He=n_He0*u^3`. Both supplied densities count all nuclei of
their respective species. The shared charge closure is exactly
`n_e=n_H*p+n_He*q`. The two fractions use their own nuclei denominator.
Neither `q` nor the helium/hydrogen number ratio is a helium mass fraction.

The shared quantum-density owner and atomic assets are the same as
[ground-state H/He LTE](hydrogen-helium-equilibrium.md):

Their adopted engine serializations at the original base are pinned by SHA256:
`hydrogen_equilibrium.hpp` is
`a96ec02bc050f837b857f48dfa015d5bbf6cbb6a52164dea5db62c876eaf0597`,
`hydrogen_helium_equilibrium.hpp` is
`eb8e4de5f4894b11762eb61fd2c4b97623171b5b18995c0971533f525b2258cc`,
and the shared quantum-density helper is
`3030df1e422bf77c8417f055d66fb04e74a8a1eae2ae082d72988a5113f4e102`.
These identify compiled scalar/equation bytes, not original NIST query or
CODATA download bytes. Earlier engine source-read records support adoption
of the central values; this history packet does not supply a new independent
audit or original-download serialization hash for those atomic sources.
That original atomic-source hash/serialization qualification remains separate
from the conditional IVP and its numerical acceptance.

```
Q_e = (2*pi*m_e*kB*T_m/h^2)^(3/2)
chi_H = 13.598434599702 eV; chi_HeI = 24.587389011 eV
E_H = h*c/(121.5682 nm)
E_He2s = h*c*(1.66277434e7 /m)
E_He2p = h*c*(1.71134891e7 /m)
delta = E_He2p-E_He2s
t=T_m/10000 K
alpha_H = 1e-19*4.309*t^(-0.6166)/(1+0.6703*t^0.5300) [m^3/s]
alpha_He = 10^(-16.744)/[sqrt(T_m/10^0.477121 K)
             *(1+sqrt(T_m/10^0.477121 K))^(1-0.711)
             *(1+sqrt(T_m/10^5.114 K))^(1+0.711)] [m^3/s]
beta_H = alpha_H*Q_e*exp[-(chi_H-E_H)/(kB*T_m)] [1/s]
beta_He = 4*alpha_He*Q_e*exp[-(chi_HeI-E_He2s)/(kB*T_m)] [1/s]
b_H = alpha_H*Q_e*exp[-chi_H/(kB*T_m)] [1/s]
b_He = 4*alpha_He*Q_e*exp[-chi_HeI/(kB*T_m)] [1/s]
K_H = (121.5682 nm)^3/(8*pi*H)
K_He = (1/(1.71134891e7 /m))^3/(8*pi*H)
v_H = K_H*n_H*(1-p); v_He = K_He*n_He*(1-q)
C_H = (1+8.22458*v_H)/(1+(8.22458+beta_H)*v_H)
w = exp[-delta/(kB*T_m)]
C_He = (w+51.3*v_He)/(w+(51.3+beta_He)*v_He)
dp/dz = C_H*[alpha_H*n_e*p-b_H*(1-p)]/[H*u]
dq/dz = C_He*[alpha_He*n_e*q-b_He*(1-q)]/[H*u]
```

The stable He expression is algebraically the selected positive-exponent
escape law. It requires no arbitrary exponential cap or small-ionization
zeroing. It retains all finite positive He states. `b_i` combines the two
detailed-balance exponentials before underflow; it is the same declared
`beta_i*exp(-E_i/kB/T_m)` equation.

Paper equation 5 supplies the truncated Compton/adiabatic heat structure,
with the same shared charge and the explicitly supplied helium/hydrogen nuclei
ratio. Equations 6/7 are the paper's helium Saha relations, not heat equations.
The only retained heat processes are Compton coupling and adiabatic expansion:

```
T_r = T_CMB,0*u
u_gamma = Omega_gamma*rho_critical_energy(H0)*u^4 [J/m^3]
Gamma_C = 8*sigma_T*u_gamma/(3*m_e*c)*n_e/(n_H+n_He+n_e) [1/s]
dT_m/dz = [Gamma_C*(T_m-T_r)+2*H*T_m]/[H*u]
q_T = c*sigma_T*n_e/[H*u] [per redshift]
```

The retained thermal mapper owns photon energy, H and flat closure. Exact SI
constants, CODATA2022 `m_e` and `sigma_T=6.6524587051e-29 m^2` are fixed.
H/He atomic uncertainty, fitted-rate uncertainty, the singlet/energy closure,
the absent HeIII dynamics and other heating/cooling are excluded from numerical
diagnostics. The source helium fit's quoted accuracy applies to 4000–10000 K;
its use below that interval is explicitly an extrapolation of this model.

The initial temperature is `T_m=T_r`. A restricted two-stage Saha neutrality
root uses `p=A/(n_e+A)`, `q=B/(n_e+B)`,
`A=Q_e exp(-chi_H/kT)`, `B=4 Q_e exp(-chi_HeI/kT)`, and the same shared charge
closure. It is the stationary balance of the subsequent rate equations at a
fixed state. The existing three-stage LTE result is not projected, clipped or
renormalized into this model. Saha supplies only this initial condition.

`D_HeIII/n_e=Q_e exp(-chi_HeII/kT_m)/n_e` is a retained excluded-stage activity
witness. Its maximum logarithm is evaluated over every retained history node
and must be at most `log(1e-12)`. It is an adjacent-stage LTE activity and a
physical scope gate; it is not a numerical-error term or a proof of a kinetic
HeIII bound between nodes.

## Domain and ownership

The first profile admits `2600<=z_initial<=2800`, `300<=z_late<=600`,
`60<=H0<=80 km/s/Mpc`, `0.015<=omega_b<=0.03`, `0.08<=omega_cdm<=0.15`,
`2.70<=T_CMB,0<=2.75 K` and `0<=omega_massless_nonphoton<=3e-5`.
At most sixteen explicit zero-mass thermal species are admitted under their
existing temperature/weight/flat-closure gates. Positive-mass relics are
explicitly outside this first profile; the existing pure-H relic consumer
retains its separate identity.

Independent present nuclei inputs require `0.1<=n_H0<=0.3 /m^3`,
`0.005<=n_He0<=0.035 /m^3`, and `0.04<=n_He0/n_H0<=0.12`, with a nonempty
`nuclei_origin` label. These are numerical scope bounds, not a cosmic abundance
law. Consistency with the separate background baryon mass density is caller
input responsibility. The [supplied abundance mapper](baryon-abundance.md)
can form both densities once from its explicit mass convention; this history
does not infer masses, Y, BBN or an observational abundance distribution.

The owner retains the request, background and prepared nodes. Copies own their
storage; moves invalidate the source owner; copy assignment builds a complete
replacement before changing its destination. Evaluation is an ordered coarse
batch and performs no new ODE/background work. Requested groups have separate
status, optional value and absolute error estimate. Unrequested or refused
values are never zeros. Invalid query rows preserve other valid rows.

`thermal_mapping_witnesses()` exposes the four once-captured scalar witnesses
of the original complete [physical map](thermal-observables.md), in photon,
baryon, CDM and other-massless order. A later preparation refusal retains
these already earned fields. Copying retains their values; moving clears
them from the source owner. They are arithmetic metadata, and are not
propagated into the current kinetic-history errors or interpreted as a
physical abundance law.

One private within-cell owner supplies the interpolation used by every
history query. It retains the existing linear fraction/temperature law,
neighboring-slope curvature estimates and interpolated nodal coefficient
`A=c*sigma_T/[H*(1+z)]`. Opacity is `A*n_e` with the same shared electron
density and error composition. It performs no new interior background query
and owns no second full history vector. Extracting this common cell law does
not provide a drag primitive, late opacity tail or physical drag epoch.

| Output mask | Group | Units |
| --- | --- | --- |
| 1 | HII fraction per H nucleus | dimensionless |
| 2 | HeII fraction per He nucleus | dimensionless |
| 4 | shared electron density | electrons/m^3 |
| 8 | matter temperature | K |
| 16 | Thomson opacity | per redshift |

## Numerical contract and independent evidence

Production uses decreasing-redshift backward Euler. It eliminates the linear
temperature equation and solves the two coupled fractions with an analytic
2x2 Jacobian, determinant/conditioning gates and a physical-box line search.
The convergence check uses both relative resolvable Newton root corrections;
their accumulated component/temperature estimates enter the retained errors. Raw stiff
residuals alone are not accurate root-distance measures. There is no asserted
global uniqueness theorem for the temperature-eliminated joint root. Failure
refuses the history; no fraction clipping or equilibrium substitution occurs.

Nested quadratic-redshift meshes at N, 2N and 4N form order-two Richardson
states. Their separate old/new extrapolants, interpolation curvature, inherited
H diagnostics and wide/cast arithmetic provide empirical error estimates.
Intermediate middle/coarse interpolation can cancel in the mesh difference.
Writing their signed errors as `e_m,e_c`, the leading Richardson error is
`-Delta/3-2*e_m+e_c/3`; the estimate therefore separately charges
`2*E_m+E_c/3`, plus the final query curvature. The preserved first independent
comparison found an underestimated diagnostic before these terms were retained;
the consumer allocations were unchanged.
Strict floating point, nearest rounding and long double with at least 64
mantissa bits and exponent range 16384 are required.

The named default numerical allocations are `1e-8+2e-6*abs(x)` for each
fraction, `2e-4 K+2e-6*abs(T_m)` for temperature, and
`2e-9+3e-6*abs(q_T)` per redshift for opacity. The electron-density allocation
is the positive weighted sum of the two fraction allocations, with shared
nuclei arithmetic. These are numerical consumer allocations, not physical
rate/atomic uncertainty. Independent reference refinement must consume at
most five percent of each named allocation.

Default N is 8192 (minimum two intervals for a curvature stencil), finest nodes
are bounded by 65536 intervals, work is
bounded by four million charged background evaluations, momentum callbacks,
initial neutrality evaluations and coupled RHS/Jacobian evaluations. Failed
trials consume the original work cap. The default explicit payload cap is
32 MiB and the hard cap is 1 GiB. Bounds include simultaneous preparation
buffers, retained request/background and returned rows, excluding allocator
metadata and process RSS. Allocation failure returns `work_limit`.

The named full-group comparison consumer explicitly requests N=16384 with
the original work, byte and numerical allocations. The default N=8192 remains
unchanged and can refuse a trace-He fraction while preserving H, electron,
temperature and opacity groups. The strengthened interpolation diagnostic
first exposed that refusal at redshift 1300; both the original diagnostic
failure and the subsequent default full-group refusal are preserved. A larger
mesh is an explicit caller request, not an automatic fallback or budget change.

The installed standalone consumer in `cpp/tests/test_installed_consumer.cpp`
uses the public `irred/hydrogen_helium_history.hpp` header and linked static
library. It prepares synthetic supplied nuclei with the explicit N=16384
policy, retains a copy after replacing the original owner and changing caller
inputs, moves that copy, and evaluates an ordered three-redshift batch. It
checks all five output groups and the same shared-charge equation without
additional physical work. `tools/check_install.py` compiles and runs this
consumer against a fresh install. This exercises SDK composition and lifetime;
it does not qualify cosmic nuclei inputs or the absent physical drag tail.

The permanent independent comparison uses an original three-variable order3
L-stable Radau IIA algorithm in log expansion, a direct scaled-cubic Saha
initialization and a separately authored escape-rate form. It imports no
engine helpers or native results. Three named nominal/boundary states include
an explicit massless thermal species. Two reference resolutions must refine
within five percent of each unchanged allocation; every actual discrepancy
also must fit the reported native diagnostic plus reference refinement and
wide reference arithmetic. A default/finer native comparison checks the
accepted groups separately from the trace-He refusal.

Frozen rate/escape facts at 800, 5000 and 7600 K come from original
Decimal110/150 direct-SI calculations with Machin pi and decimal exp/ln/sqrt.
Their maximum relative refinement is `3.837e-108`; the permanent native
allocation is `2e-15` relative. They use no engine outputs or system libm and
are not regenerated during builds. The complete stiff histories still use
long double/system libm on both routes; these rate facts are not a
high-precision certificate for the entire trajectory.

Shared constants/equations are declared physical ancestry. Numerical agreement
does not qualify the physical model or provide independent atomic/rate
uncertainty. Owner tests retain charge, positivity, temperature, source lifetime,
masked outputs, strict rounding, explicit resource refusal and invalid inputs.

The absent late history remains consequential. A root integrated only to
redshift 300–600 cannot be presented as physical `z_drag`. A source-defined
endpoint needs a separately qualified late opacity/residual/reionization law
and its positive tail/boundary, plus a downstream ruler budget. This API
deliberately supplies no drag root or present-day visibility.
