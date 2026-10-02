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

The electron/proton/neutral degeneracy convention is 2/1/2, giving a prefactor
of one. Factoring nuclear spin into proton/atom counts as 2/2/4 gives the same
ratio. The translational mass ratio `m_p/m_H` is approximated by one. This
equation follows [Hirata's ground-state derivation, equation 24](http://www.tapir.caltech.edu/~chirata/ph217/lec06.pdf)
and [Sales, Carvalho and Souza, equation 6](https://link.springer.com/article/10.1140/epjc/s10052-023-11671-z).
The latter distinguishes equilibrium Saha calculations from the non-equilibrium
processes needed for a recombination history. No cosmological density mapping,
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

- Electron mass is `9.1093837139e-31 kg`, fixed to the central
  [CODATA 2022 value](https://physics.nist.gov/cuu/pdf/all.pdf). Its standard
  uncertainty is `2.8e-40 kg`.
- Ground-state hydrogen ionization energy is `13.598434599702 eV`, from
  [NIST ASD 5.12](https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=on&biblio=on&e_out=0&el_name_out=on&level_out=on&shells_out=on&spectra=H-DS+i&unc_out=on&units=1),
  with supplied uncertainty `1.2e-11 eV`.
- The shared `quantities.hpp` owns exact SI definitions
  `h=6.62607015e-34 J s`, `kB=1.380649e-23 J/K` and
  `eV_J=1.602176634e-19 J/eV`. These definitions are not atomic measurements.

The header exposes the fixed atomic central values and uncertainties for native
consumers. Atomic uncertainty is excluded from numerical diagnostics. Ground-state
truncation, ideal-gas/equilibrium assumptions and the translational mass
approximation are also excluded. They cannot be added to arithmetic error as if
they were independent measurement noise. Original source documents were read;
no external implementation, partition table or source asset is bundled or copied.

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
