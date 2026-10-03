# Homogeneous ground-state hydrogen–helium equilibrium

The native C++20 SDK `irred/hydrogen_helium_equilibrium.hpp` evaluates ordered
batches of independently supplied temperature and physical total hydrogen and
helium nuclei densities. It supplies separately represented H I/H II and
He I/He II/He III fractions, and their shared charge-neutral electron density.
It is a possible initial-state prerequisite for a later helium consumer. It
does not advance a kinetic ionization history or infer abundance, cosmological
density, expansion, drag or CMB observables. There is no CLI/C ABI route.

## Physical identity and equations

The fixed model is homogeneous ground-state H/4He LTE, with nonrelativistic
ideal Maxwell–Boltzmann species, a common temperature, zero photon chemical
potential and adjacent ion/atom translational mass ratios approximated by one.
The two supplied densities `n_H,n_He` count all ionization stages of each
species, in physical nuclei per cubic metre. Their ratio is caller input.
Neither a mass fraction nor a cosmological density is silently substituted.

The unchanged [hydrogen equilibrium](hydrogen-equilibrium.md) quantum-density
owner supplies

```
Q_e = (2*pi*m_e*kB*T/h^2)^(3/2)
A = Q_e*exp(-chi_H/(kB*T))
B = 4*Q_e*exp(-chi_HeI/(kB*T))
D = Q_e*exp(-chi_HeII/(kB*T)).
```

Ground electronic degeneracies H/H+/HeI/HeII/HeIII/e are `2/1/1/2/1/2`.
The helium states are `1s^2 1S_0`, `1s 2S_(1/2)` and a bare nucleus. The
neutral and singly ionized electronic counts follow `g=2J+1` from the selected
ASD ground states. Free-electron and bare-ion counts are supplied conventions;
common nuclear factors are omitted consistently in adjacent ratios. This
reader-derived count algebra does not audit nuclear-spin values or qualify
LTE/partition applicability. The selected helium level notes state 4He ancestry.
With shared free-electron density `e>0`, the hydrogen probabilities are
`h0=e/(e+A)`, `h1=A/(e+A)`. Helium probabilities `f0,f1,f2` normalize the
weights `1,B/e,B*D/e^2`. Charge neutrality closes the mixture:

```
C(e) = n_H*h1 + n_He*(f1+2*f2)
e = C(e), h0+h1=1, f0+f1+f2=1.
```

Independently solving hydrogen and helium with different electron densities
does not satisfy this closure. Every stage fraction is computed directly, not
as the complement of a rounded partner. A dominant binary64 fraction may round
to one while its separately represented positive partner remains available.
Independent rounded fractions are not an exact binary64 conservation identity.

Admission requires finite `1<=T<=100000 K` and strictly positive finite
`n_H,n_He`. Zero-density species are outside this mixture interface because
fractions of an absent species are undefined. The separate pure-H consumer
remains unchanged. A conservative dilute electron-gas gate requires
`(n_H+2*n_He)/Q_e<=0.001`, with a wide-arithmetic margin. Definite violations
return `outside_domain`; unresolved thresholds return conditioning refusal.
This sufficient total-density bound does not qualify LTE, ideal-gas accuracy,
ground-state partition truncation or translational-mass approximations.

## Fixed atomic assets and provenance

The model reuses the existing SI h/kB/eV definitions, CODATA2022 electron mass
and ASD hydrogen ionization energy described in [hydrogen equilibrium](hydrogen-equilibrium.md).
It adds these fixed central helium facts from the ASD 5.12 query:

- He I ionization: `24.587389011 eV`, quoted uncertainty `2.5e-8 eV`.
- He II ionization: `54.4177655282 eV`, quoted uncertainty `1e-9 eV`;
  ASD encloses this energy in parentheses, as it also does the hydrogen energy.

The accepted binding is the [seven-scalar Prospector serialization at commit
669093656632a7600fa2a4d883455ccd9a9c3da8](https://github.com/sprajs/prospector/blob/669093656632a7600fa2a4d883455ccd9a9c3da8/register/source-data/hhe-atomic-central-asd512-codata2022-v1.json),
SHA256 `d5c189e4623b5eadc2f44f5035dc959e95ff85df69fc199cf8fc0f57ca9c1839`,
with its [bounded source review](https://github.com/sprajs/prospector/blob/669093656632a7600fa2a4d883455ccd9a9c3da8/register/reviews/2026-10-03-atomic-source-claims.json).
Original snapshot hashes and selected locators preserve the H/He
query markup, official parenthesis legend, ground-state notes and SI/CODATA
facts. All seven central decimal values and four nonzero uncertainty decimal
values match the current `long double` declarations across the two equilibrium
headers and `quantities.hpp`. Decimal agreement does not qualify platform bits,
binary64 conversion or finite representation. SI definitional exactness remains
distinct from finite binary storage and from nonexact atomic central values.

ASD's energy parentheses identify ab-initio values or values otherwise not
derived from evaluated experimental data. He I has no such marker; its absence
does not certify an independently measured likelihood. The checked ASD help
does not specify standard uncertainty, a sigma multiplier, coverage probability
or a probability law for the displayed ionization uncertainties. CODATA
explicitly labels the electron-mass uncertainty as standard; its `(28)` notation
has a different meaning. No joint scalar covariance is established, and repeated
query rows/page/PDF representations are not independent evidence. The displayed
ASD eV uncertainties retain the provider's conversion contribution without
adding a second conversion uncertainty.

The selected He I/He II level notes identify ground angular momentum and 4He
ancestry. The query cites ASD 5.12, while the level datasets were prepared for
5.11 and 5.10 respectively; the He II note mentions CODATA2018. This is not
evidence that all seven scalars share a fresh CODATA2022 adjustment. Underlying
atomic papers and excited-level/partition data were not reviewed in this bounded
audit. Older rounded handbook energies are not substituted. Asset identity is
`SI2019-CODATA2022-me-NIST-ASD5.12-H1-He4-fixed`. The header exposes the new
central energies and uncertainties for native consumers.

Atomic uncertainty is excluded from numerical diagnostics. LTE, partition
truncation, ideal-gas and translational-mass assumptions are also separate;
they are not independent noise terms to add to arithmetic error. The compiled
model retains the original algebra; the new source audit supports selected
scalar facts and conditional ground electronic counts. It does not reconstruct
full source-derived Saha/LTE, partition or kinetic closure. No external solver,
partition table or database asset is bundled. NIST distinguishes SRD compilation
copyright from other government works; see its [licensing statement](https://www.nist.gov/open/license).

## Log root and separately admitted outputs

Production uses `ell=log(e)` and solves `F(ell)=ell-log(C(exp(ell)))=0`.
The Saha coefficients, probabilities and charge sum are formed in log space.
The deterministic initial bracket is

```
ell_upper = log(n_H+2*n_He)
ell_lower = log(C(exp(ell_upper))).
```

Since charge decreases as electron density increases, this brackets the unique
positive root. Safeguarded Newton steps stay in the bracket, with bisection
fallback. The derivative is `1+[n_H*h0*h1+n_He*Var(J)]/C`, where
`Var(J)=f0*f1+4*f0*f2+f1*f2`. Production forms these ratios from logarithms;
it does not divide two underflowed charge values. Tiny charge or electron
density need not be materialized to calculate available neutral fractions.

The arithmetic profile requires nearest rounding, strict floating-point flags,
and long double with at least 64 mantissa bits and exponent range 16384.
Unresolved nonfinite arithmetic or iteration/evaluation exhaustion returns an
explicit failure. The finite root residual/bracket distance and a baseline
64-wide-epsilon allowance accompany Saha logs, density logs and evaluation.
Implicit neutrality sensitivities propagate coefficient/density allowances to
log electron density. Fraction sensitivities propagate them to each requested
stage; helium's log-electron sensitivity magnitude is `abs(mean_charge-stage)`.
Each output also receives wide arithmetic and binary64 cast loss.

Runtime admission requires a relative numerical diagnostic at most `1e-13`,
with no absolute floor. This is an empirical arithmetic/libm/root/cast
diagnostic, not a universally certified bound or physical error propagation.
A requested value must be strictly positive. Underflow, nonfinite results or
insufficient cast precision return `conditioning_budget_exceeded` without a
value. Positive subnormals qualify only when their cast loss meets the same
allowance. At low temperature, neutral fractions can remain available while
charged stages and electron density fail. Unrequested or failed partners must
never be interpreted as zero.

The output mask is HI `1`, HII `2`, HeI `4`, HeII `8`, HeIII `16`, electron
number density `32`; default `63`. All values are dimensionless fractions
except `electron_density`, in free electrons per cubic metre. Each group has
its own status, optional value and relative diagnostic. Omitted groups carry
no value or diagnostic. A batch `ok` means processing was admitted; per-row
and per-group statuses remain authoritative.

## Resources and evidence

Hard bounds are 65536 rows/row solves, 256 root iterations per valid row,
4000000 cumulative charge evaluations and 1 GiB explicit result payload.
Configured caps may be smaller. Every charge evaluation, including the initial
bracket evaluation and a final refused solve, consumes actual work. Row/batch
counters distinguish solves, root iterations and charge evaluations.

Row/byte bounds precede source scans and allocation. Valid physical rows consume
one solve; invalid physical rows preserve their cause and consume none. Exhausted
work tags affected later valid rows while retaining prior results. The payload
covers the batch header and row array; borrowed inputs, stack, allocator
overhead and RSS are excluded. There is one row allocation per nonempty batch.
Exact source scalar bits and order are owned by result rows. Ordinary
`std::bad_alloc` propagates with RAII cleanup; no server/cache or per-row FFI
is needed.

Owner regressions test separate positivity, nuclei conservation and shared
neutrality; H/He degeneracy identities including HeI's factor four; all 63
nonempty masks; neutral availability below wide electron representation;
subnormal/cast refusal; domain/degeneracy boundaries; source bits/lifetime;
rounding, cumulative quotas, payload thresholds and allocation exceptions.
Independent comparisons allocate `2e-12` relative without an absolute floor
and with a separate positivity gate; references must refine within 5% of that
allocation. Their independent polynomial/high-precision ancestry and actual
execution are separate from these owner identities. Numerical agreement does
not qualify equilibrium as a kinetic recombination history or a cosmic closure.
