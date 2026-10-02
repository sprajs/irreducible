# Finite shared optical calibration law

The native C++20 `evaluate_calibration` in `irred/photometry_calibration.hpp`
predicts collected energy and expected transmitted photon signal for a fixed
sampled rest-wavelength spectrum, fixed luminosity distance and redshift, and
per-band fixed observed wavelength grids, collecting areas and observer exposures.
It calls `evaluate_sampled`; its piecewise-linear interpolation, propagation,
optical transmission, SI conventions and scientific domain are unchanged.
There is no CLI or C/Rust ABI for this bounded native consumer.

Each joint state contains every band's transmission array. All state masses are
positive finite normal binary64 relative masses. One long-double sum normalizes
them once: p_s=w_s/sum(w). The expectation is mu_i=sum_s p_s f_si and the
population calibration covariance is C_ij=sum_s p_s(f_si-mu_i)(f_sj-mu_j).
There is no sample-size correction. Shared state membership preserves cross-band
and energy/photon dependence. Rank-deficient positive semidefinite covariance is
valid and receives no jitter. This is variance of expected signal, not photon shot
noise, detection/electronic noise, source, distance or population uncertainty.
Incident flux is fixed under transmission changes and is not a stochastic axis.

`CalibrationResult::normalized_state_mass` retains the owner's wide normalized
weights in input state order, including very small positive states. The default
`CalibrationAggregation::population_moments` retains the existing full moment
calculation and sensitivity gates. An explicit `state_resolved_only` request
returns the same required state-band predictions and weights while omitting
population moments and their aggregate sensitivity gates. This supports the
[joint detector composition](optical-detector.md), whose likelihood uses each
conditional signal directly. It does not substitute unresolved rounded spread
for a zero calibration covariance; every required state/output still must pass.

Axes follow input band order, energy before photons when both are requested.
Energy means use joules, photon means use expected counts; covariance units are
products of the corresponding axis units. Source, calibration, distribution and
dependence origins and unique nonempty band/state IDs are supplied explicitly.
Origins are declarations, not evidence that calibrations are measured or correct.
Every state-band sampled result remains in state then band order. Any required
failed output prevents the whole aggregate; no state is dropped or renormalized.

Inputs are borrowed only during the synchronous call; results own IDs and values.
Hard limits are 1024 states, 64 bands and 1048576 cumulative source/passband knots
across calls (including one extra source length). Each parent call also retains
its 65536 knot cap. Declared dimensions, cumulative knots, covariance size and a
conservative allocation payload estimate are checked before sample scans or
allocation. Policy quotas can lower these limits. The byte estimate includes
owned string characters, statuses, moments, axes and wide scratch; caller storage,
allocator overhead, scalar frames and RSS are excluded. Allocation failure returns
a work-limit status. Empty admission diagnostics have no aggregate payload.

The named empirical forward allocation is 3e-12 abs(f) per parent scalar. Wide
summation arithmetic estimates are added; deterministic sensitivities propagate
linearly rather than as independent random errors. Centered covariance sensitivity
uses |d_i| e_j + |d_j| e_i + e_i e_j, with centered e including mean sensitivity.
Default mean sensitivity is 1e-10 relative and covariance sensitivity is 1e-8
relative to sqrt(C_ii C_jj). These are empirical diagnostics, not rigorous error
bounds or certification of the calibration distribution. Physical calibration
spread and numerical sensitivity are reported separately. Unresolved tiny spread
is refused, including equal rounded outputs from changing transmissions.
Exact input transmission equality, a single state, or zero source/area/exposure
can witness constant axes. Mathematical zero differs from nonzero output or
sensitivity underflow; nonzero subnormal, overflow and nonfinite casts fail.

Owner synthetic controls predeclare mean absolute 1e-280 plus relative 2e-10;
covariance absolute 1e-280 plus relative 2e-8 against geometric variance scale.
Independent analytic polynomial references have no quadrature truncation and share
SI/pi ancestry only; the absolute floor never admits a nonzero reference as zero.
Controls cover correlated energy/photon axes, anticorrelated bands, exact constant
inputs, linear spectra/transmissions, state permutation, mass/distance scaling,
zero source/area, partial failures, hostile masses/IDs, quotas and unresolved spread.
Independent frequency-coordinate refinement is a separate peer gate. Numerical
agreement does not establish astrophysical or inference qualification.
