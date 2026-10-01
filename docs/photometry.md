# Deterministic rectangular-band photometry

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

Each requested group carries its own numerical status. Default numerical-contract assurance can pass while interpretation remains unqualified. The operation is available through the one-shot `run` interface; it does not create a retained observation handle.

The byte quota admits successful combined native/wrapper allocated calculation payload. A fixed bounded empty work-limit diagnostic owner is excluded from that admission budget and may be allocated even when the quota is zero. Borrowed caller storage, scalar call frames, allocator overhead and RSS are also excluded.

## Native sampled-spectrum and passband operator

The independently usable C++20 library also provides `evaluate_sampled` in
`irred/sampled_photometry.hpp`. This extension has no CLI operation or C/Rust ABI.
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
