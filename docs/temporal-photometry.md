# Finite sampled spectral-time photometry

The native C++20 SDK prepares immutable pooled spectral-time grids and observed
optical passbands, then predicts a coarse batch of declared observer exposures.
The model identity is
`finite_bilinear_rest_spectral_time_zero_outside_full_observer_exposure_mean`.
It composes the existing [sampled radiometry](photometry.md) operator; there is no
new wavelength-density equation, C/Rust ABI, CLI operation, template reader or
training-data qualification in this slice.

## Equations and physical identity

A `TemporalGrid` supplies an ordered rest-time axis τ in seconds, an ordered
rest-wavelength axis λr in metres, and a rectangular time-major array of spectral
luminosities Lλ in W/m. Cell `(i,j)` is stored at
`i*wavelength_count+j`. Within every declared time/wavelength cell the model is
bilinear: at each wavelength knot luminosity is linear in rest time, and each
instantaneous spectrum is linear in rest wavelength. It is zero outside both
finite supports. Endpoint values specify interior limits, with no extrapolation.
This declared zero-outside model does not establish the actual source beyond a
measured, empirical or training-domain interval.

A `TemporalBand` supplies the existing observed-wavelength `SampledPassband`:
metres and dimensionless optical transmission T in [0,1], linear in wavelength
and zero outside finite support. Its entire support B defines incident band flux,
including portions where T=0. Optical transmission enters collected energy and
photons once; it is not quantum efficiency, electronic gain or ADU conversion.

Each `TemporalExposure` declares luminosity distance DL, redshift z, collecting
area A, source epoch e on the observer time axis and an observer interval [a,b]
with b>a. The observer epoch corresponds to rest phase τ=0. Times share one
caller-declared clock origin and second unit; this operator does not convert UTC,
time standards or observer frames. The assumed propagation is

```
r = 1+z
t_observer = e+r*τ
λ_observer = r*λ_rest
Fλ(λ_observer,t_observer) = Lλ(λ_observer/r,(t_observer-e)/r)/(4π DL² r).
```

Distance, redshift, area and passband are fixed during each exposure. This is a
supplied propagation convention, not a time-evolving geometry or fitted
cosmology. Define Δt=b−a and the **full requested exposure** mean spectrum

```
Lbarλ(λr) = (1/Δt) integral from a to b of Lλ(λr,(t_observer-e)/r) dt_observer.
```

The operation integrates this mean at wavelength knots by splitting the
observer interval at mapped rest-time knots. Nonnegative endpoint trapezoids
integrate the declared linear time functions exactly in algebra. Their integral
remains linear in wavelength, so one call to `evaluate_sampled` on `Lbarλ`, the
fixed passband and observer duration Δt supplies each exposure's observables:

- `mean_flux_watt_per_square_metre`: `(1/Δt) ∫dt ∫B Fλ dλ`, in W/m².
- `energy_joule`: `A ∫dt ∫B Fλ T dλ`, in joules.
- `expected_photons`: `A/(hc) ∫dt ∫B λ Fλ T dλ`, a transmitted photon expectation.

The shared sampled operator owns wavelength clipping, propagation, optical
weights, SI constants and photon conversion. There is no additional redshift
time factor multiplying energy. Mean flux divides by the entire observer
interval, including its explicitly zero parts. It is never renormalized by the
covered duration.

## Source ownership and domains

Use `prepare_temporal` from `irred/temporal_photometry.hpp`. Grid and band inputs
carry unique nonempty IDs (within their respective pools, at most 256 bytes) and
nonempty source/calibration origins (at most 1024 bytes). Origins are preserved
caller declarations, not independent source verification. Grids have at least
two finite strictly increasing rest times, which may be negative phases, and at
least two positive finite strictly increasing rest wavelengths. Every luminosity
is finite and nonnegative. Passband axes follow the same wavelength admission;
transmission is finite in [0,1]. All arrays are admitted before a zero shortcut.

Preparation checks shape, hard/configured resource limits and then scientific
arrays before copying each admitted pool once. The immutable owner retains exact
axis/value bits, matrix order, IDs and origins. Copies share that same source
storage. Moves leave their source invalid; self move preserves the owner. Its
`grid`, `band` and `retained_payload_bytes` getters expose owned views and bounds.
A view remains valid while a prepared owner sharing the storage remains alive;
it does not acquire another ownership lifetime on its own. Invalid preparation
returns a precise status with no usable source payload.

`evaluate_temporal` receives ordered exposure rows referring to grid/band indices
and retains their exact scalar/index source values, including the epoch. All
scalars must be finite; DL>0, z≥0, A≥0 and b>a. A zero-length exposure is refused:
its full-interval mean is undefined. Invalid references or scalar domains produce
row failures. Other rows remain in their original order. Independent requested
output masks are the existing flux/energy/photon bits. Unrequested conversions
are omitted; a spectral photon or flux failure does not erase an independently
available collected quantity.

Each row reports `TemporalCoverage::full`, `partial` or `no_overlap`. Coverage
refers to the supplied time support and the entire requested exposure. Wide
`observer_duration_second`, `covered_observer_second` and `covered_fraction`
diagnostics avoid converting a positive tiny coverage to binary64 zero.
`unassessed` is used when support inclusion cannot be established. Exact touching
is zero measure. The declared zero outside time support makes a known absent
intersection a mathematical zero; it does not infer an unobserved source history.

## Numerical allocation and refusals

Named analytic and independent controls use2×10^-12 relative plus 10^-300
absolute for each final observable, with a separate positivity guard. The floor
cannot accept a silently zero positive result. The temporal mean's outward
interval width plus its binary64 conversion has a 1×10^-13 relative admission
allowance; this leaves margin below the shared spectral operator's 1.8×10^-12
interval-admission threshold. The frozen allocation is for the named controls,
not a certified universal all-domain accuracy bound.

Mapped time endpoints use error-free sums and FMA product residuals under
FE_TONEAREST with a wide intermediate profile of at least 64 mantissa bits and
max exponent 16384. Outward endpoints must establish positive segment width,
overlap, interpolation conditioning and full-versus-partial support inclusion.
An unresolved edge returns `conditioning_budget_exceeded`, rather than deleting
a tiny positive segment or reporting zero/full coverage. Adjacent binary-exact
time knots have a positive control; redshift-induced thin overlap and a large
epoch hiding positive phase duration have permanent refusal controls.

The observer duration and averaged spectral knot values must be representable
as binary64 with the declared rounding allowance before calling the shared
operator. Positive underflow is refused, not replaced by zero. A duration larger
than binary64 capacity is an overflow. This composition may refuse an intermediate
mean/duration even where a separately rescaled final integral could be
representable; it does not implement an alternate widened spectral equation.
The independent sampled groups retain their existing typed output failures.

The independent peer uses original bilinear polynomial antiderivatives and direct
observer-time/frequency Gauss-Legendre integration. It reconstructs Lν=λr²Lλ/c,
Fν=rLν(rν)/(4πDL²), optical T(c/ν) and photon weight 1/(hν), splits every time and
source/passband spectral knot, and compares 32/64-point refinements within
2×10^-13 relative. It does not call production interpolation, integration or
radiometric equations in its reference. Exact SI definitions and system pi are
shared ancestry and cannot establish independent physical validation. Inserted
exactly linear time/wavelength knots preserve the model as a separate refinement
control. Original failed controls are preserved in local evidence.

Source-model/interpolation uncertainty is excluded. Supplied distance/redshift
and epoch uncertainty is excluded. Passband/calibration uncertainty is excluded;
there is no joint calibration distribution in this consumer. These are separate
physical/model limitations, not terms hidden within the arithmetic allowance.
No stochastic detector response, count realization, selection, trained SN
spectral template, light-curve fitting or inference is supplied.

## Resource and batch contract

Preparation's hard limits are 4096 grids, 4096 bands, 65536 total pooled time/source-
wavelength/passband knots, 4,000,000 grid luminosities and 1 GiB admitted native
payload. Configured preparation caps can be smaller. Bounds and checked products
precede array scans and allocation. Payload includes the prepared owner/native
storage, pool-owner objects and exact ID/origin/axis/value vector payloads.
Borrowed caller spans, shared-ownership control-block bookkeeping, allocator
overhead and RSS are excluded. Shared storage is charged once by its preparation
lifetime, not once per exposure or copied owner.

Evaluation hard limits are 65536 rows, 4,000,000 cumulative conservative segment
work and 1 GiB result/scratch payload. Per row, work is charged as
`(time_count−1)*wavelength_count + wavelength_count + passband_count−3`.
This is admission work for all declared temporal cells plus the shared spectral
segment bound, not a measured speed claim. Scheduled attempts consume it even
when a later conditioning check fails. Exhaustion gives a tagged row `work_limit`
and preserves prior rows; bad references/scalars receive their own causes.
A configured row or result/scratch-byte failure returns an empty batch diagnostic.

One reused mapped-time interval array, one temporal integral array and one
binary64 mean-spectrum array serve the batch. Result rows and their fixed scalar
storage remain charged even when some outputs are omitted. Existing prepared
source storage is owned/charged separately and is not copied into results.
Native allocation failures propagate ordinary `std::bad_alloc`; scoped failure
sweeps verify preparation/evaluation RAII cleanup. No exception crosses a C ABI
because this slice has no C ABI route.

This is an independently usable native SDK consumer. Passing its ordinary CI
controls establishes the named deterministic numerical and ownership evidence,
without qualifying the source model, passband, propagation or interpretation.
