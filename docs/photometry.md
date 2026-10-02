# Deterministic supplied-spectrum photometry

`photometry.predict` implements this deterministic calculation as one coarse input batch. The first slice predicts incident band-integrated flux, collected energy and expected transmitted photons from a steady, isotropic source. It does not generate noise, select objects, simulate electronics, infer a luminosity or reproduce a survey.

The rest-frame spectral luminosity is a constant Lλ in W per metre on a finite interval [λr0, λr1], and zero elsewhere. The observed rectangular passband is [λo0, λo1] in metres. Inputs supply luminosity distance DL in metres, redshift z, observer exposure t in seconds, collecting area A in m² and constant optical transmission T. Distance and redshift are declared propagation inputs, not a cosmology fitted here. The assumed mapping is λo=(1+z)λr, dt_obs=(1+z)dt_rest and Fλ=Lλ/[4π DL²(1+z)] on the redshifted source support. See [Hogg's wavelength flux relation](https://ned.ipac.caltech.edu/level5/Hogg/Hogg7.html). No alternate geometry automatically inherits this mapping.

Let [a,b] be the intersection of observed passband and redshifted source support. A provably empty intersection has zero measure. For positive width:

- Incident band flux is Fλ(b−a), in W/m².
- Collected energy is A T t Fλ(b−a), in joules.
- Expected transmitted photons are A T t Fλ(b−a)(a+b)/(2hc).

The last expression uses the observed photon energy hc/λo. Area, optical transmission and observer exposure enter once. T is optical transmission, not detector quantum efficiency, electrons per photon or ADU gain. Rest elapsed time t/(1+z) is a convention, not another factor multiplying the prediction. The [BIPM defining constants](https://www.bipm.org/en/measurement-units/si-defining-constants) are exact SI 2019 c=299792458 m/s and h=6.62607015×10^-34 J s, identified as `si_2019_radiometric_definitions`. Numerical storage of their exact decimal definitions is still rounded to the calculation's floating-point type.

All fields must be finite; wavelengths and distance are positive, intervals strictly increasing, z≥0, Lλ/A/t≥0 and 0≤T≤1. Full scientific admission precedes zero shortcuts. Requested outputs are independent: zero area makes collected quantities zero but does not erase incident flux. An unrequested photon conversion cannot invalidate energy. Mathematical zero and positive values that cannot be represented in binary64 are distinct; the latter fail, including underflow.

## Predeclared verification

Named positive controls allocate relative error 2×10^-12 plus an absolute floor 10^-300, with the independent reference allocation 2×10^-13 relative. The floor cannot accept silent zero or catastrophic thin-band loss: positivity/representability and relative edge conditioning are checked separately. Named controls span z=0–10, wavelengths 10^-7–2×10^-6 m, DL=1–10^26 m, Lλ=10^-10–10^40 W/m, A=0–100 m², t=0–10^5 s and T=0–1. These cases do not certify every combination or arbitrary finite inputs. Supplied distance/calibration uncertainty is not numerical rounding; a small distance perturbation contributes approximately twice its relative uncertainty to flux.

Controls use direct wavelength algebra and independent frequency-coordinate integration: Lν=cLλ/ν² and Fν=(1+z)Lν((1+z)ν)/(4πDL²). Composite Simpson refinement tests the nonpolynomial photon integrand Fν/(hν), reverses frequency endpoints and compares energy/photons with the analytic result. Scaling, full-support bolometric conservation, source/observer time, clipping, unit-equivalent inputs and invalid domains are separate gates.

Exact binary input semantics matter at tangency. With s=2^-20, e=2^-52, z=e, source [s,s(1+e)] and observed band starting at s(1+2e), the true positive overlap is 2^-124 m. A second source ending at s with z=2^-65 and an observed band starting at s has positive overlap 2^-85 m. Wide rounding alone can erase either overlap. Outward bounds from error-free sums and FMA product residuals must establish overlap sign and acceptable width conditioning; otherwise return a typed conditioning failure. Exact touching at z=1 and observed lower edge 2s is zero. Adjacent representable endpoints, disjoint bands and thin bands inside the support are permanent adversaries.

Native and boundary tests must preserve source bits/order, independent requested states, output ownership, count/byte limits, allocation-failure cleanup, mixed valid/failed rows and record/source identities. Runtime qualification remains unqualified; passing these named numerical controls is not a stochastic recovery, astrophysical model or inference claim.

Run the synthetic example after building:

```bash
target/release/irred run tests/fixtures/photometry.json /tmp/irred-photometry-example
```

Each requested group carries its own numerical status. Default numerical-contract assurance can pass while interpretation remains unqualified. The operation is available through one-shot `run` and stateless `photometry_predict` stream requests; it does not create a retained observation handle.

The byte quota admits successful combined native/wrapper allocated calculation payload. A fixed bounded empty work-limit diagnostic owner is excluded from that admission budget and may be allocated even when the quota is zero. Borrowed caller storage, scalar call frames, allocator overhead and RSS are also excluded.

## Native sampled-spectrum and passband operator

The independently usable C++20 library also provides `evaluate_sampled` in
`irred/sampled_photometry.hpp`. The same operator is exposed by `photometry.predict` with source model `piecewise_linear_rest_luminosity_observed_optical_passband` through one coarse C/Rust ABI call.
`SampledSpectrum` borrows rest-frame wavelength and Lλ arrays (metres and W/m);
`SampledPassband` borrows observed wavelength and optical-transmission arrays.
Each axis has at least two finite, strictly increasing positive wavelengths;
source luminosities are finite and nonnegative and transmission lies in [0,1].
Scalar distance, redshift, area and observer exposure follow the rectangular
operator's domain and propagation convention.

Both arrays declare exact piecewise-linear functions in wavelength, with zero
outside their finite supports. Endpoint samples specify their interior limits;
there is no extrapolation. Samples do not establish the true continuous source
or instrument response between measurements. The operator does not normalize
transmission, infer a spectrum, accept Fν samples or propagate sample/calibration
covariance. Calibration is fixed. Source, distance, calibration and interpolation
model uncertainty remain outside the numerical allocation.

With r=1+z, use Fλ(λ)=Lλ(λ/r)/(4πDL²r). If B is the entire finite passband
support, the three independent outputs are:

- Incident flux: integral over B of Fλ dλ, in W/m², without transmission weights.
- Collected energy: At integral over B of Fλ T dλ, in joules.
- Expected transmitted photons: At/(hc) integral over B of λ Fλ T dλ.

T is optical transmission and enters collected quantities once. It is neither
quantum efficiency nor electronic gain. Passband regions with zero transmission
still belong to B for incident flux. Source support clipping is applied separately.
Observer exposure receives no additional redshift factor.

The implementation walks the merged source/passband intervals without allocation.
Within each interval it integrates the linear source and quadratic source-times-
transmission product analytically. Multiplying that product by observed wavelength
produces a cubic. Nonnegative Bernstein coefficients and their exact integral
weights avoid subtraction of nearly equal polynomial antiderivatives. There is
no empirical quadrature-error estimate or claim of a certified universal bound.
Outward endpoint/arithmetic checks reject unresolved overlaps and excessive
conditioning; they distinguish exact touching from redshift-induced positive
thin overlap. Positive binary64 underflow and overflow are typed failures.
Each requested result has its own availability and numerical status; unrequested
conversions cannot invalidate another output. Scientific admission precedes zero
shortcuts.

`SampledPolicy` requests the existing flux/energy/photon mask and bounds total
source-plus-passband knots and merged segments. The hard total-knot cap is 65536;
the conservative segment admission count is total knots minus three. Limits are
checked before walking arrays. A policy/work-limit failure has no output payload.
Arrays remain caller-owned and need only live through the synchronous call;
`SampledResult` owns scalar outputs and retains no source views. There is no
cache, retained observation object or per-row FFI.

Owner controls use independent constant/linear polynomial antiderivatives,
clipping, refinement by inserting exactly linear knots, inverse-distance/area/
exposure scaling, bolometric redshift conservation, zero transmission, invalid
arrays, quotas, omitted outputs, unresolved binary endpoint adversaries and
isolated output underflow/overflow. The named allocation is 2×10^-12 relative
plus 10^-300 absolute, with a separate positivity guard; the independent
frequency-coordinate refinement allocation is 2×10^-13 relative. These are
synthetic numerical controls with shared SI definitions, not measurements,
calibration validation, inference or astrophysical qualification.

A native [finite shared passband-calibration law](photometry-calibration.md) composes the sampled operator with an explicitly supplied joint distribution of valid transmission states. Its expected-signal covariance retains cross-band calibration dependence. Source and distance remain fixed; numerical sensitivity, calibration spread and future photon/noise simulation have separate meanings.

## Pooled sampled ingestion and records

The sampled request uses schema version 2 and the same propagation and SI constants
identities as the rectangular request. Supply ordered `spectra`, `passbands` and
`exposures`. Each spectrum supplies `id`, `source_role`, `provenance`,
`rest_wavelength_metre` and `luminosity_watt_per_metre`; each passband supplies
`id`, `source_role`, `provenance`, `observed_wavelength_metre`,
`optical_transmission`, `calibration` and `calibration_provenance`. Array pair
lengths must match. Wavelengths are metres, luminosity is W/m in the rest frame,
and optical transmission is dimensionless in the observed frame. No Fν or
observed flux interpretation is inferred from array shape.

Source roles are caller-declared `measured`, `fitted_summary`,
`calibration_asset` or `synthetic_control`. Provenance declarations are retained,
not independently verified. `calibration` is `fixed` or
`declared_uncertainty_excluded`: both evaluate supplied fixed values, and the
latter explicitly records that declared calibration uncertainty is excluded.
Neither propagates sample covariance or adds noise. IDs must be unique within
each pool, nonempty and at most 256 UTF-8 bytes. Provenance strings are nonempty
and at most 1024 UTF-8 bytes. The JSON schema's string limits count characters;
runtime byte admission is additionally required for non-ASCII declarations.

Exposure rows refer to `spectrum_index` and `passband_index` (zero-based pool
indices) and supply `luminosity_distance_metre`, `redshift`,
`collecting_area_square_metre` and `observer_exposure_second`. They preserve
order and binary64 source values. A bad reference is a row `invalid_input`
failure, not a silently substituted spectrum. Scientifically invalid arrays,
scalars, unresolved conditioning and positive unrepresentable results carry the
native cause. Requested groups remain independent and unrequested groups are
omitted; per-row work-limit failure exposes failed requested groups without a
fabricated value.

`resource_policy` requires `maximum_rows`, `maximum_native_bytes`,
`maximum_total_samples`, `maximum_samples` and `maximum_segments`. Hard limits
are 4096 spectra, 4096 passbands, 65536 exposure rows, 65536 pooled wavelength
knots and 1 GiB admitted native payload. The last two knot/segment policy fields
bound each native source/passband evaluation, with the existing conservative
merged-segment count of pair knots minus three. Aggregate configured row,
pooled-knot or byte exhaustion returns an owned empty `work_limit` diagnostic;
per-row exhaustion preserves other rows. Hard bounds and descriptor lengths are
checked before native sample scans or allocation. CLI requests are bounded to
16 MiB, with bounded individual deserialization arrays; aggregate knot admission
precedes wire-descriptor allocation and the native call. Raw duplicate JSON
members are rejected before intermediate values can erase them.

The additive revision-2 C route is `irred_sampled_photometry_evaluate` followed
by `irred_sampled_photometry_result_view` and
`irred_sampled_photometry_result_destroy`. Caller curve buffers are borrowed
only through the synchronous call. Its result owner copies each pooled curve
once, keeps ID bytes and wavelengths/values immutable, and exports source and
row views whose lifetime ends at destroy. Repeated exposure rows use indices;
there is no per-row FFI, repeated matrix transfer, cache or retained session
source. The independent C++ sampled API remains allocation-free and borrowed.

Byte admission covers the native owner, pooled curve owner objects, exact ID byte storage, two double arrays per curve, exported curve
descriptors and result rows. It excludes borrowed Rust inputs/wire descriptors,
scalar frames, allocator overhead and RSS. The fixed empty diagnostic owner is
excluded and can be allocated at quota zero. Output omission avoids native
unrequested conversions; fixed row wire storage remains charged. The result
owner contains no copied spectrum per exposure.

Run the synthetic pooled fixture after building:

```sh
target/debug/irred run tests/fixtures/sampled_photometry.json /tmp/irred-sampled-example
```

For a stream, wrap the same request in
`{"action":"photometry_predict","request":REQUEST}`. Its transient byte policy
is capped by remaining session retained-byte allowance. The resolved scientific
specification preserves the requested policy; effective runtime policy and
actual resources are execution evidence. One-shot and stream share operation
code and produce the same numerical values with equal effective allowances.

Records retain exact original request bytes by SHA-256, the resolved ordered
pools, caller roles and provenance, row sources, group checks and actual compiled
method `analytic_piecewise_linear_merged_bernstein` and arithmetic
`binary64_storage_longdouble_intermediate`. Native transport controls establish
exact parity with `evaluate_sampled`, immutable source lifetime, malformed
references/descriptors, quotas, allocation/exception containment and omissions.
Independent native constant/linear antiderivative and frequency-coordinate
controls retain the predeclared allocations above. These checks do not qualify
observations, interpolation, calibration, cosmology or inference.
