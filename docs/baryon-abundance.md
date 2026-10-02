# Supplied H1/He4 abundance and a shared-electron LTE consumer

The native C++20 `irred/baryon_abundance.hpp` API maps an explicitly supplied
baryon density and helium mass fraction into ordered physical nuclei densities.
Its minimal consumer supplies those densities to the existing shared-electron
[H/He ground-state LTE operator](hydrogen-helium-equilibrium.md), with an
explicitly supplied matter temperature. There is no CLI or C ABI route.

## Physical identity

`physical_baryon_density` means `omega_b=Omega_b*h^2`, with `h=H0/100`.
The source supplies neutral ground-state H1 and He4 effective masses in kg,
with `source_origin` and `mass_origin` labels. Those labels record caller
provenance; they do not certify a mass measurement. No helium mass is inferred,
no number-fraction interpretation is selected, and `m_He=4*m_H` is not assumed.
`helium4_mass_fraction=Y` is the He4 share of this specified fixed effective
rest-mass density. The source contains only these two species.

The unchanged private thermal critical-density helper owns the SI equation,
with the existing IAU Mpc, exact SI c and fixed CODATA2018 G convention:

```
rho0 = rho_critical(H100) * omega_b
rho(a) = rho0/a^3
n_H(a) = rho(a)*(1-Y)/m_H
n_He(a) = rho(a)*Y/m_He
n_He/n_H = [Y/(1-Y)] * m_H/m_He
Y = m_He*n_He/(m_H*n_H + m_He*n_He)
```

Nuclei densities count every ionization stage in physical nuclei per cubic
metre. Masses include the caller's declared neutral electron/binding convention,
remain fixed across ionization, and do not feed binding-energy changes into the
pressureless thermal background. The current pure-H history's effective mass,
arithmetic and scientific identity are untouched. This API supplies no BBN
prediction, abundance measurement, isotope mixture, helium kinetic rates,
full recombination, expansion feedback, physical drag epoch or CMB closure.

## Admission, ownership and arithmetic

Prepare an immutable source once, then evaluate ordered coarse scale-factor
batches with `0<a<=1`. Source admission requires finite `omega_b>=0`, finite
`0<=Y<=1`, normal positive finite masses and nonempty provenance labels. A
structurally admitted source is retained even when subsequent physical admission
fails. Copies own their strings; moves transfer ownership and invalidate the
moved source; self moves preserve the owner.

The mapping supports exact zero densities at `omega_b=0`, `Y=0` and `Y=1`.
A zero density carries a zero diagnostic. Positive outputs use wide intermediate
density/mass arithmetic and independently admitted H and He projections: an
unrepresentable H output does not destroy a representable He output. Positive
underflow is refused, never silently set to zero. Subnormal binary64 outputs
are accepted only if their measured cast loss meets the same relative budget.

The arithmetic profile requires nearest rounding and long double with at least
64 mantissa bits and exponent range 16384. The empirical operation diagnostic
charges the unchanged critical-density sequence, mapping products/divisions,
and actual binary64 cast loss, conservatively rounding diagnostics upward.
Mapping admission requires a relative estimate at most `1e-15`. This is an
arithmetic diagnostic, not a universally certified bound; mass, G, abundance
uncertainty and model assumptions are separately excluded.

The LTE consumer retains supplied `{scale_factor, temperature_kelvin}` rows and
a nonempty `temperature_origin`. Temperature is supplied matter temperature;
there is no implicit photon-temperature scaling or thermal evolution. It
requires positive baryon density and `0<Y<1`, because the existing mixture
operator excludes fractions of absent species. Endpoint requests return an
explicit `outside_domain` LTE refusal, while their mapping remains available.
No pure-H or pure-He route is substituted.

One existing LTE batch owns charge neutrality, Saha coefficients and solver
work. Its temperature, dilute-gas, requested-output and numerical gates remain
in force. With maximum relative mapping diagnostic `d`, the additional absolute
log-density allowance is `d/(1-d)`. Monotone charge neutrality limits the induced
log-electron error to this allowance; stage fractions receive at most twice it.
The consumer adds these allowances to the native output diagnostics and
requires each combined relative estimate to meet `1e-13`. Failed or omitted
outputs have no value; they must never be interpreted as zero.

Batch success means processing admission. Mapping-row admission and individual
H/He groups, then LTE-row and individual LTE groups, remain authoritative.
Invalid mappings consume no charge solves. Global LTE solve/iteration/evaluation
counters are retained, including partially failed batches.

Hard limits are 65536 rows and 1 GiB explicit native payload. Checked accounting
charges retained source strings, returned rows and temperature provenance,
and simultaneous temporary mapped LTE inputs and native LTE results. The LTE
callee's separate byte/work limits also apply. These bounds exclude borrowed
inputs, allocator bookkeeping and RSS. `std::bad_alloc` propagates with ordinary
RAII cleanup; no failed allocation is converted into a scientific value.

## Evidence

Owner tests challenge dimensional inversion, unequal-mass number ratios,
`a^-3` scaling, endpoints and adjacent binary64 values, independently refused
representations, masks, global quotas, moves/copies and each measured allocation
failure. Independent peer tests retain exact-binary-input high-precision
critical-density/mapping facts and a separately solved charge-neutrality
reference. Shared physical constants or atomic ancestry remain explicit;
transport agreement with the native LTE operator alone is not independent
scientific evidence. Named comparisons establish bounded numerical evidence,
not measured abundance or cosmological qualification.
