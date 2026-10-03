# Homogeneous ground-state hydrogen equilibrium

The native C++20 SDK `irred/hydrogen_equilibrium.hpp` evaluates an ordered batch
of supplied temperatures and physical total hydrogen nuclei densities. It
returns separately represented ionized and neutral fractions under a declared
pure-hydrogen ground-state Saha model. There is no CLI/C ABI route, cosmic
ionization history, helium calculation, recombination kinetics or drag
prediction in this consumer.

## Equation and physical domain

Let `T` be temperature in kelvin and `n_H` the physical total number of hydrogen
nuclei per cubic metre, including both protons and neutral atoms. Charge
neutrality gives `n_e=n_p=x*n_H`; neutral density is `y*n_H`, with `x+y=1`.
For nonrelativistic ideal Maxwell-Boltzmann species in chemical equilibrium with
zero photon chemical potential, the fixed ground-state model is

```
n_Qe = (2*pi*m_e*kB*T/h^2)^(3/2)
K = n_Qe*exp(-chi/(kB*T))/n_H
x^2/y = K,  x+y = 1.
```

The supplied electron/bare-ion/neutral electronic count convention is 2/1/2,
giving a prefactor of one. The neutral electronic count follows `g=2J+1` for
the ASD ground `2S_(1/2)` state; free-electron and bare-ion counts are model
conventions. Common nuclear factors are omitted consistently, and the
translational mass ratio `m_p/m_H` is approximated by one. Ground-state Saha
references are [Hirata, equation 24](http://www.tapir.caltech.edu/~chirata/ph217/lec06.pdf)
and [Sales, Carvalho and Souza, equation 6](https://link.springer.com/article/10.1140/epjc/s10052-023-11671-z).
The scalar-source audit below does not qualify partition truncation, LTE
applicability, nuclear-spin or mass corrections. No cosmological density mapping,
photon temperature history, expansion rate or helium fraction is inferred here.

Every input must be finite, with `1<=T<=100000 K` and `n_H>0`. A conservative
nondegeneracy gate requires `n_H/n_Qe<=0.001` before evaluating fractions. Total
nuclei density bounds the free-electron density from above, so this is sufficient
for the declared dilute electron gas approximation. It excludes denser inputs
even where a smaller actual free-electron density might admit another treatment.
It does not certify ground-state partition truncation, LTE, ideal-gas accuracy or
excited-state effects over the temperature range. These remain separate model
assumptions; a numerical admission is not a physical validation.

## Atomic assets and shared definitions

The atomic assets identity is
`SI2019-exact-h-kB-eV-CODATA2022-me-NIST-ASD5.12-H1-fixed`:

- Electron mass is `9.1093837139e-31 kg`, fixed to the central CODATA 2022
  value. The official electron-mass page explicitly labels its standard
  uncertainty as `2.8e-40 kg`. Its concise `(28)` notation gives the last
  digits of that uncertainty; it is distinct from ASD's theory marker.
- Ground-state hydrogen ionization energy is `13.598434599702 eV`, from
  NIST ASD 5.12, with displayed uncertainty `1.2e-11 eV`. ASD encloses the
  energy in parentheses: its official legend identifies an ab-initio value
  or one otherwise not derived from evaluated experimental data. The selected
  query does not identify a hydrogen isotope.
- The shared `quantities.hpp` owns exact SI definitions
  `h=6.62607015e-34 J s`, `kB=1.380649e-23 J/K` and
  `eV_J=1.602176634e-19 J/eV`. These definitions are not atomic measurements.
  Their exact defining decimal values do not imply exact finite binary storage.

The accepted source binding is the [seven-scalar Prospector serialization at
commit 669093656632a7600fa2a4d883455ccd9a9c3da8](https://github.com/sprajs/prospector/blob/669093656632a7600fa2a4d883455ccd9a9c3da8/register/source-data/hhe-atomic-central-asd512-codata2022-v1.json),
SHA256 `d5c189e4623b5eadc2f44f5035dc959e95ff85df69fc199cf8fc0f57ca9c1839`.
Its [bounded source review](https://github.com/sprajs/prospector/blob/669093656632a7600fa2a4d883455ccd9a9c3da8/register/reviews/2026-10-03-atomic-source-claims.json)
records original snapshot hashes and selected locators.
The seven central decimal values and four nonzero uncertainty decimal values
match the current `long double` literals in `hydrogen_equilibrium.hpp`,
`hydrogen_helium_equilibrium.hpp` and `quantities.hpp`. This comparison checks
decimal transcription; it does not establish platform bit patterns, a binary64
conversion policy or the numerical adequacy of a representation.

The checked ASD help does not establish standard uncertainty, a sigma
multiplier, coverage probability or a probability law for the three ionization
energies. CODATA's electron-mass standard uncertainty does not supply those
missing ASD meanings. ASD says its displayed eV uncertainties already include
its conversion-factor contribution; no second conversion uncertainty is added.
No joint covariance among these seven scalars is established. Repeated query
rows and CODATA page/PDF representations are duplicate source facts.

The header exposes the fixed atomic central values and uncertainties for native
consumers. Atomic uncertainty is excluded from numerical diagnostics. Ground-state
truncation, ideal-gas/equilibrium assumptions and the translational mass
approximation are also excluded. They cannot be added to arithmetic error as if
they were independent measurement noise. The source binding supports selected
scalar facts and conditional ground electronic counts; it does not qualify a
full source-derived Saha closure or recombination history. No external
implementation, partition table or source asset is bundled or copied.

## Stable separate fractions

Production forms `log(K)=log(n_Qe)-chi/(kB*T)-log(n_H)` with wide intermediates.
The quantum-density equation is owned once and also serves the admission gate.
Two algebraic branches avoid overflow and complement cancellation:

```
log(K)<=0:
  q=exp(log(K)/2), s=hypot(q,2)
  x=2*q/(q+s)
  y=2/(2+q*q+q*s)
log(K)>0:
  v=exp(-log(K)), s=sqrt(1+4*v)
  x=2/(1+s)
  y=2*v/(1+2*v+s)
```

Both fractions are computed directly. Production never derives a tiny fraction
by subtracting the rounded larger one from one. Binary64 values may round to
one while the separately stored partner remains strictly positive. These are
independently rounded representations of the mathematical pair; their binary64
sum is not an exact conservation identity or a representation of a zero partner.

The arithmetic profile requires `FE_TONEAREST`, strict floating-point compilation
and long double with at least 64 mantissa bits and exponent range 16384. A baseline
64-wide-epsilon allowance accompanies logs, quantum density and binding ratio.
Sensitivity of `ln(x)` to `ln(K)` is `y/(2-x)`; the magnitude for `ln(y)` is
`x/(2-x)`. Each requested output's diagnostic adds this propagated allowance,
wide arithmetic and binary64 cast loss. Admission requires at most `1e-13`
relative. This is an empirical arithmetic/libm diagnostic, not a universally
certified bound or a propagated atomic/model uncertainty.

Each requested fraction is strictly positive. A value that rounds to zero or
whose cast/sensitivity diagnostic exceeds the admission allowance returns
`conditioning_budget_exceeded` with no value. Positive subnormals are admitted
only if their relative cast loss meets the same budget. A low-temperature tiny
ionized fraction can fail while the separately computed neutral fraction rounds
to one and remains available. If only the neutral group was requested, the tiny
ionized group remains omitted; an omitted or failed partner must not be
interpreted as zero. The same rules apply to a tiny neutral fraction.

The nondegeneracy gate accounts for a conservative wide-arithmetic margin.
An unresolved threshold returns a conditioning failure; a definite violation
returns `outside_domain`. There is no silently clipped boundary.

## Batch, resources and independent evidence

`HydrogenState` carries temperature and density. `evaluate_hydrogen_equilibrium`
copies exact source scalars into result rows in input order. Invalid/nonfinite
rows retain their cause; valid rows have independently requested ionized and
neutral groups. A batch `ok` status means processing was admitted, not that every
row or requested fraction succeeded. Per-group optional values and statuses are
authoritative. Unrequested groups have no value or diagnostic.

Hard bounds are 65536 rows, 65536 cumulative constant-work fraction solves and
1 GiB native result payload. Configured caps can be smaller. Checked row/byte
bounds precede sample scans and allocation. Valid physical rows consume one
solve, including an eventual output refusal; invalid physical rows consume none.
Work exhaustion tags later rows while preserving earlier results. Row/byte
admission failures return an empty batch diagnostic. The payload covers the
batch header and requested row array; borrowed caller spans, stack, allocator
overhead and RSS are excluded. There is one row allocation per nonempty batch;
there is no retained server, per-row FFI or hidden observation copy. Native
allocation failure propagates ordinary `std::bad_alloc` with RAII cleanup.

Named permanent comparisons allocate `2e-12` fractional relative error, with
no absolute floor and a separate positive-output requirement. Four rational
ionization controls at 10000 K use frozen binary64 densities derived from exact
rational targets. Five independently authored 110-digit Decimal logit-equation
bisection controls retain the exact binary64 inputs, including an available
subnormal neutral fraction. 150-digit refinement changes each reference by less
than `1e-100` relative. Ordinary peer tests solve the original log equation by
independent 128/256-step monotone logit bisection, with refinement allocated
`2e-13` relative; they use neither production quadratic branch. Shared physical
assets and CPU libm ancestry are explicit, while Decimal supplies a separate
arithmetic route. Numerical convergence alone does not establish the physical
equilibrium assumptions.

Owner controls challenge temperature/density domains, nondegeneracy, separate
positive fractions, insufficient cast precision, output omission, mixed rows,
source bits/order/lifetime, rounding mode, hard/configured quotas and allocation
exception cleanup. Peer controls also retain monotonic temperature/density and
neutral/ionized limits. A standalone consumer compiles against a fresh installed
header and archive. Ordinary tests require no acquired atomic files or external
engine; the small frozen references and their essential ancestry remain public.
