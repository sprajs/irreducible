# Conditional supplied-depth H/He drag and ruler

The native C++ request in `irred/conditional_hydrogen_helium_drag.hpp` retains one
supplied abundance map, one [bounded H/He history](hydrogen-helium-history.md),
and that history's thermal background. It asks how a supplied **synthetic**
late-depth interval changes a conditional unit-depth root and tight-coupling
comoving ruler. Its identity is
`exact-emitted-H1-He4-nuclei-and-retained-thermal-state-supplied-D-conditional-drag-ruler/v1`.
It supplies no reionization law, full HeIII history, measured late-depth prior,
probability on an interval, or source-qualified physical drag prediction.

The source uses one physical thermal model, a supplied He4 neutral-effective
mass fraction and supplied H1/He4 effective masses in kg. The
[abundance owner](baryon-abundance.md) emits nuclei densities at a=1. Those exact
binary64 values and the complete original mapped thermal state are then fixed
working inputs. The original physical densities, temperature, masses, source
origins, emitted bits, abundance diagnostics and four once-captured thermal-map
witnesses remain available. Their operation/cast errors are provenance under
this working law; they are **unpropagated**, and are not a physical-input error
channel or measured uncertainty. Nonzero mass/photon reconstruction residuals
remain visible. H0 is copied exactly. The four mapper records are ordered
photon, baryon, CDM, other massless fractional density.

The first conditional profile has empty explicit thermal species. The existing
standalone thermal and H/He domains retain their own gates. H/He starts at the
restricted, artificial ground-state LTE initialization in 2600<=z_i<=2800 and
ends in 300<=z_L<=600. It retains shared electron charge, the singlet-only
code convention, temperature evolution, HeIII activity witness and low-T fit
extrapolation of the H/He guide. The literal SSS1999 equations and RecfastCLASS
binding/statistical/escape conventions remain distinct; this consumer does not
repair or promote their discrepancy. Original atomic source qualifications and
fixed central-value assets are inherited from the H/He guide.

For u=1+z, a=1/u and independently emitted nuclei nH0,nHe0, the shared cell law
is ne=u³(nH0 p+nHe0 q). Its opacity coefficient A is the affine interpolation
of nodal c sigma_T/(H u); a direct interior background quotient would be a
different numerical law. One private cell owner serves history queries and
the conditional positive primitive. The same retained background supplies H,
P=a⁴E², loading R0=3 Omega_b/(4 Omega_gamma), R=R0/u and sound speed
c/sqrt[3(1+R)]. H is in s⁻¹ for optical depth, ne in m⁻³, sigma_T in m² and c
in m/s. The adopted native sigma_T is 6.6524587051e-29 m² (CODATA2022); fixed
G is 6.67430e-11 SI (CODATA2018). They have distinct source provenance and are
not random numerical errors or newly measured atomic masses.

The physical definition uses depth from today:

```
wD(z) = qT(z)/R(z) = c sigma_T ne/[H u R]
D = tau_D(z_L),  K(z) = integral[z_L,z] wD dz,  tau_D(z)=D+K(z)
t = 1-D
[z_min,z_max] = [K^-1(1-D_upper), K^-1(1-D_lower)]
rs(z) = c/(H0 sqrt3) integral[0,1/(1+z)] da/sqrt[P(a)(1+R0 a)]  [Mpc]
drs/dz [Mpc] = -cs/(H * megaparsec_in_metres)
```

The selected CLASS drag definition has source pin
`e85808324f51fc694d12e3ed7439552a3c3f9540`; the redshift change of variables is
reviewer algebra. The H/He convention pin
`64bbab707faf4de4779a9e04edd180fef18d98fa` is separate ancestry and supplies no
late-history closure. **Physical D remains unqualified.** D=0 is an explicitly
truncated mathematical control, rather than evidence of negligible late depth.
D=1 has the exact lower-boundary root z_L. No supplied D means no requested
conditional endpoints.

Assuming all electrons arise from the fixed H1/He4 inventory, with no other
charged sector or pair production, positivity gives
0<=ne<=(nH0+2nHe0)u³ and

```
0 <= D <= U
U = c sigma_T (nH0+2nHe0)/(R0 H0_SI)
    integral[a_L,1] da/[a³ sqrt(P(a))].
```

The factor two permits unknown late HeIII despite the early singlet history.
U can exceed one. This ceiling supplies no lower bound, tail shape, useful late
prior or probability. A supplied interval must pass both the numerical capacity
gate and retained K support. Unsupported or unresolved endpoints refuse; the
interval is never trimmed to available support. Other ordered rows remain
visible. A ruler refusal retains its accepted root when both roots are supported.

Affine u,p,q,A give degree-six positive drag coefficients. Endpoint/curvature
and loading diagnostics give degree at most eight. The implementation restricts
nonnegative Bernstein coefficients before direct partial-cell integration and
retains one compensated positive prefix value/error pair. It does not subtract
nearby global antiderivatives or fit projected history queries.

Each root uses a locator bracket and at most eight expanded proposals within
128 endpoint trials. The **same expanded bracket** must contain the locator and
all depth/target shifts, have uniform EK, a strictly positive rate lower bound m,
and a common-background upper bound L for cs/H. Robust endpoint signs and
rho=(EK+ET)/m close the bracket. A central derivative alone is insufficient.
Nonpositive bounds, support crossing, exhaustion and unresolved positive
quantities refuse. The rationalized loading correction is
xi=ER a_max/[sqrt(B)(sqrt(A)+sqrt(B))], A=1+R0 a_max,
B=1+(R0-ER)a_max; a tiny positive correction must resolve. The shared loading
ratio uses 32 wide eps times |R0| plus measured reporting-cast loss. Mapper
errors remain separate from that propagated ratio diagnostic.

Every endpoint has a combined absolute allocation of 1e-3 Mpc, relative zero:
5e-4 inherited depth/history/loading, 1e-4 locator, 3e-4 complete shared ruler
(quadrature, background, loading, arithmetic), and 1e-4 target/additional
projection. Each slot and the sum must pass. The same projected redshift serves
the root and ruler. Physical interval width is separate from numerical error.
Refinement and arithmetic diagnostics are empirical; Richardson differences
do not certify an original physical enclosure. The independent reference must
earn its **combined refinement plus all affecting arithmetic <=5e-5 Mpc for each
endpoint**. A work-count/epsilon or Newton proxy alone does not earn that final
qualification; the source-only reference records that gate as withheld.

Hard caps are four million aggregate work and 32 MiB owned live payload, four
intervals/eight endpoints, 200000 outer callbacks per integral and depth30.
The meaningful installed synthetic caller explicitly chooses N=16384/fine65536,
requested history work3999998 and bytes33546240. These are requests, not proof
of an 8 KiB parent envelope. The owner checks actual remaining work/live bytes
and retains requested/served history ceilings and phase-attempt diagnostics.
Defaults remain unchanged. Each scoped batch includes preparation once and all
failed attempts; child fields transfer individually through checked additions.
Public work checksums can refuse overflow and cannot supply ledger freshness.

Payload bounds include actual string/vector capacities, owners, mapping records,
history preparation scratch, one prefix, bounded coefficient scratch, ordered
outputs and simultaneous return headers without NRVO. Borrowed
unrelated owners, allocator bookkeeping, recursion stack and RSS are excluded.
The source-only caller flag disables elision in its test translation unit. It
does not disable return elision in the compiled library; complete no-NRVO
validation of the library and caller remains a separate qualification gate.
This new conditional owner is move-only: its actual consumer prepares once and
borrows immutable retained history through const coarse evaluation. It supplies
no whole-owner copy constructor or copy assignment; the existing standalone
H/He history copy interface is unchanged. Source-clearing moves are noexcept,
and self-move retains the owner. Ordinary allocation failures use RAII. The
installed owning receipt outlives its input request and intervals. Source and
budget diagnostics move with it; a moved-from owner is invalid. Source admission
requires a valid policy, bounded nonempty origin strings, the empty explicit
species profile and initial payload admission. A request refused before that
boundary has no retained source or budget record and has zero attempted work.
The admitted preparation phase retains requested caps; failure to allocate the
source copy can still leave the source record absent. Once the request has been
acquired, later refusals preserve it, requested caps, earned work and whichever
phase diagnostics have actually been reached. Source admission
and committed controls do not establish numerical, resource, installed-library
or observational qualification. No production CLI/ABI extension is supplied.
