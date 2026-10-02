# Finite supplied joint abundance, mass and matter-temperature law

The native C++20 `irred/baryon_abundance_law.hpp` consumer retains an explicitly
supplied finite joint law. Every state contains physical `omega_b`, the He4
effective-rest-mass fraction `Y`, neutral effective H1/He4 masses, and a whole
ordered vector of supplied matter temperatures. Ordered scale factors are fixed.
`prepare_baryon_abundance_law` prepares each existing
[abundance owner](baryon-abundance.md) once. `evaluate_density` returns nuclei
moments; `evaluate_equilibrium(mask)` composes the existing shared-electron
[ground-state LTE consumer](hydrogen-helium-equilibrium.md) and returns requested
LTE moments. There is no CLI or C/Rust ABI route.

## Identity and equations

The unchanged abundance owner is the sole owner of

```
rho_s(a) = rho_critical(H100) * omega_b_s / a^3
n_H_s(a) = rho_s(a) * (1-Y_s) / m_H_s
n_He_s(a) = rho_s(a) * Y_s / m_He_s.
```

`omega_b=Omega_b*h^2` is physical density; no H0 parameter is added. `Y` is the
He4 share of the supplied fixed neutral-effective rest-mass density, rather than
a number fraction. The masses retain the source's neutral electron/binding
convention. Temperature is supplied matter temperature with explicit provenance;
no photon-temperature scaling or thermal evolution is inferred. The same shared
electron density closes all H/He ionization stages. This consumer does not copy
the critical-density, Saha or charge-neutrality equations.

State masses are positive finite normal binary64 relative masses. One compensated
long-double sum normalizes them once, `p_s=w_s/sum_t(w_t)`. Tiny positive
probabilities remain wide numerical coordinates, including probabilities that
would underflow in binary64. Every state belongs to the same complete joint law;
there are no separate marginal laws, Gaussian perturbations, automatic
factorization, rejection sampling, or clipped Y/mass coordinates.

For the entire ordered output vector `f_s`, the supplied-law moments are

```
mu_i = sum_s p_s f_si
C_ij = sum_s p_s (f_si-mu_i) (f_sj-mu_j).
```

These are population moments, without a sample-size correction. Covariance is
computed by a centered second pass, retaining all cross-row/output dependence.
Singular covariance is valid. No Cholesky admission, jitter, PSD projection or
dropped state repairs the law. Zero covariance does not establish independence.

These moments propagate the supplied abundance/mass/temperature law conditional
on the fixed model. They are distinct from numerical diagnostics, a parameter
posterior, an independent measurement or physical qualification. G and atomic
central assets remain fixed: this interface supplies no G or atomic-asset law
coordinates. An extension that varies those assets requires explicit supported
coordinates and corresponding shared owners; a provenance label cannot enable
that propagation. BBN prediction, helium kinetics, evolving binding masses,
expansion feedback and a cosmological thermal history remain separate work.

## Ownership, order and partial results

Inputs are borrowed during preparation. Unique nonempty state and row IDs,
distribution/dependence origins, and each state's source/mass/temperature origins
are required. Each state supplies exactly one temperature per row. Original
scalar bits, IDs and order are retained in one shared immutable source object;
results retain that same object after the original owner is released. Copies
share immutable preparation. Moves invalidate the source owner; self moves
preserve it. Provenance labels are declarations, rather than qualification of
the supplied probability law.

Preparation requires the existing abundance source domain, fixed finite
`0<a<=1`, and finite normal positive temperatures. A structurally admitted source
is retained when later physical admission fails. Mapping admits exact zero
densities at zero omega or Y endpoints. Positive temperatures outside the LTE
operator's narrower `1<=T<=1e5 K` support leave density evaluation useful; LTE
evaluation retains an explicit native refusal. The native dilute-gas gate and
both-positive-nuclei requirements also remain in force. No pure-H or pure-He
operator replaces a refused mixture.

Density axes follow row order, H nuclei before He nuclei. LTE axes follow row
order, then requested native mask bits in ascending order: HI `1`, HII `2`, HeI
`4`, HeII `8`, HeIII `16`, electron density `32`; default `63`. Axes retain a row
index into the immutable row IDs. Density and electron means have units m^-3;
fractions are dimensionless. Covariance units are the product of the associated
axis units, including m^-6 for two number-density axes.

Attempts retain one native batch per state in input order, with all native rows,
mapping results, requested group statuses/optional values and actual solver
counters. Density moments require both mapped nuclei groups. LTE moments require
both mapped nuclei and precisely the requested LTE groups. Any required failure
withholds the entire moments payload, preserving successful attempts and
available partners. Failure causes follow state/row/output order. Unrequested
groups have no value or diagnostic; they neither add an axis nor block requested
groups. Required states are never omitted or renormalized after failure.

## Numerical sensitivity and truthful zero spread

The arithmetic profile requires nearest rounding, strict floating-point flags,
and long double with at least 64 mantissa bits and exponent range 16384. Existing
relative parent diagnostics remain at most `1e-15` for mapping and `1e-13` for
combined LTE. For a parent value `f_si` with relative diagnostic `delta_si`,
the inherited absolute allowance is `e_si=abs(f_si)*delta_si`.

Wide normalization, summation/product arithmetic and actual binary64 cast loss
are charged. Errors propagate deterministically and linearly, rather than as
independent random noise. Mean sensitivity includes weighted inherited errors
and weight-normalization sensitivity. With centered `d_si=f_si-mu_i` and
`E_si=e_si+mean_error_i` plus subtraction arithmetic, covariance sensitivity
includes

```
sum_s p_s (abs(d_si)*E_sj + abs(d_sj)*E_si + E_si*E_sj),
```

and weight, product, summation and output-cast allowances. Returned absolute
diagnostics are rounded upward. The default mean allocation is `1e-12*abs(mu_i)`;
the covariance allocation is `1e-8*sqrt(C_ii*C_jj)`. Policies may tighten these
allocations, but cannot silently relax the named contract. These are empirical
numerical diagnostics, not universally certified bounds or additional physical
uncertainty. Law spread and numerical sensitivity occupy separate result arrays.

Zero variance needs a sufficient mathematical input witness: a single state;
identical relevant raw inputs across states; or an absent density axis in every
state. H density depends on omega/Y/mH, He density on omega/Y/mHe, and LTE on the
full abundance state plus that row's T. Temperature-only variation therefore
witnesses constant density axes. Witnessed constant axes have exact zero
covariance rows/columns and zero covariance sensitivity: shared deterministic
forward error does not manufacture physical spread.

Equal rounded outputs from varying relevant inputs are not a constant witness.
A nonconstant axis with unresolved/nonpositive variance refuses moments. Signed
off-diagonal cancellation is evaluated against the geometric variance scale and
does not establish independence. Positive mean/variance or signed nonzero
covariance underflow is refused, never replaced by zero. Positive subnormals are
accepted only when their measured cast loss and upward representable diagnostic
meet the same allocation. Overflow/nonfinite projections also refuse. No
absolute floor admits a positive reference as zero.

## Combined resources and evidence

Hard bounds are 1024 states, 64 rows, 65536 state-row attempts, 384 output axes,
2^27 upper-triangle state/product terms, 65536 cumulative solves, 256 iterations
per solve, four million cumulative charge evaluations and 1 GiB combined explicit
native payload. Policies may lower these. Dimensions, checked products and
string/payload bounds precede scalar scans/allocation. A conservative product
bound is checked even when constant axes later eliminate products.

The conservative preparation envelope charges copied strings as
`max(32,length+1)` bytes, with checked overflow, including the supported
libstdc++ empty-string capacity growth for lengths 16 through 29. This envelope
is checked before allocation. The source separately records the conservative
`preparation_native_payload_bound` and actual capacity-based
`retained_native_payload_bytes`; a later actual-size check does not replace the
admission envelope. Child LTE temperature-origin copies use the same conservative
envelope in the evaluation bound.

The combined byte bound charges immutable retained source, attempts, IDs/strings,
moments, wide scratch, and simultaneous child input/result temporaries. It
excludes caller input storage, scalar stack frames, allocator bookkeeping,
shared-pointer control-block allocation and RSS. Native child bounds also apply.
Every child receives the remaining global solve/evaluation quota; failed solves
consume actual work and later attempts retain refusals. Results report mapping
rows, solves, root iterations, charge evaluations, actual covariance products and
the checked combined payload bound. Preparation records state preparations and
normalization terms. `std::bad_alloc` propagates with ordinary RAII cleanup.

Owner controls cover exact rational weighted moments, inverse-mass and a^-3/a^-6
scaling, Y-induced anticorrelation, shared cross-row dependence, all 63 LTE masks,
constant input witnesses, equal-rounded varying-state refusals, absent species,
positive subnormal moments and underflow refusals. They also challenge state
permutation, duplicated support, mass normalization, tiny probabilities, global
work/payload quotas, source lifetime/bits/order, moves and every measured
allocation-failure site. Mapping agreement with the same owner remains transport
evidence. Independent exact-binary high-precision abundance/charge-root/moment
controls use density mean `2e-15` relative, LTE mean `5e-13` relative and covariance
`2e-8` of geometric reference variance, with reference refinement consuming at
most 5% of each allocation and separate positivity/zero-witness gates. Shared SI
and atomic facts are declared ancestry, not independent physical evidence.
Matched-quality performance comparisons must include preparation, repeated
evaluation, actual solver/product work, native payload and RSS; retained
preparation alone is not a demonstrated workload speedup.
