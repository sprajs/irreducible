# Supplied continuous scalar CMB transfers

The C++ SDK projects a supplied continuous scalar source into finite-support
temperature and E-mode transfer amplitudes. It owns the radial projection and
declared time interpolation. The producer owns photon, species, metric,
scattering and polarization physics, its gauge, signed initial mode and
normalization. This consumer does not construct those source terms.

Include `irred/continuous_cmb_projection.hpp` and link `libirred_core.a`.
`ContinuousCmbSource` contains strictly increasing positive `k_mpc_inverse`,
strictly increasing nonnegative `eta_mpc`, one `observer_eta_mpc` at or after
the last source time, and four equally shaped channel arrays. Their order is
eta-major, k-minor: `channel[eta_index*k.size()+k_index]`. Every channel has
units of inverse conformal Mpc per the same dimensionless signed initial mode.
Nonempty `producer_id`, `signed_mode_id` and `normalization_id` are opaque
provenance supplied by the caller; the SDK does not infer or verify a physical
mode from their spelling. All four channel arrays must be present and finite.

With \(x=k(\eta_0-\eta)\), the stored source channels define

\[
 \Theta_\ell(k)=\int_{\eta_a}^{\eta_b}
 \left[T_0j_\ell(x)+T_1j'_\ell(x)
 +T_2\frac{3j''_\ell(x)+j_\ell(x)}{2}\right]d\eta,
\]

\[
 E_\ell(k)=\int_{\eta_a}^{\eta_b}P_{\rm src}
 \sqrt{\frac{3(\ell-1)\ell(\ell+1)(\ell+2)}{8}}
 \frac{j_\ell(x)}{x^2}\,d\eta.
\]

Primes act on the dimensionless argument. The temperature convention is the
outward angular phase \(e^{+ix\mu}\) with angular split source
\(T_0+i\mu T_1-P_2(\mu)T_2\). The E kernel follows the positive scalar CLASS
convention when `polarization` is its actual
\(P_{\rm src}=\sqrt6\,g\,P_{\rm internal}\). Other producers must supply
equivalently defined channels and record their conversion. No Kelvin factor,
primordial amplitude or spectrum normalization is applied by this operator.

The first supported domain is temperature \(0\leq\ell\leq64\), E-mode
\(2\leq\ell\leq64\), and \(0\leq x\leq512\) throughout the supplied
support. Regular endpoint kernels retain their limits, including the
quadrupole/E \(\ell=2\) value \(1/5\). Every time cell uses the linear
interpolant of its original endpoint source values. There is no interpolation
in k, source differentiation, integration by parts or deletion of finite
boundary terms. Nothing assigns zero outside supplied support, accounts for
omitted early/late tails or renormalizes missing survival mass.

Move a source into `prepare_continuous_cmb_projection` once. Validation and
payload checks occur before acquisition. On success, the move-only prepared
owner retains the original buffers through one shared immutable source;
callers must relinquish mutable aliases. `project_continuous_cmb` accepts an
ordered multipole span and temperature/E output mask. Rows preserve
multipole-major, original k order. Each result retains that same source,
including evaluation refusals, and may outlive the prepared owner. The
requested mask controls the evaluated channels and reported amplitudes.

Vector adaptive Simpson integrates the declared linear source with phase
splitting. Its nested time estimate is empirical. A separate conditional
arithmetic estimate retains positive interpolation/assembly scales, actual
node counts and the measured final cast loss, including signed cancellation.
The default per-output gate is
\(10^{-9}+10^{-5}|X|\); passing it concerns these two estimates only.
`radial_error_estimate` and `source_grid_error_estimate` remain absent
(`std::nullopt`), not zero. Time refinement supplies neither a standard-library
Bessel certificate nor error for a different smooth producer source.

Defaults bound source acquisition to 256 k values, 16,384 times, 1,048,576
array cells and 64 MiB of counted source payload. Evaluation defaults allow
128 requested multipoles, two million radial-node attempts, recursion depth
20 and 128 MiB of counted source/result/numerical scratch payload. Capacity
and size multiplication checks precede numerical work. Payload accounting
excludes allocator metadata, RSS and arbitrary copies of earlier results.
These are bounds for the declared consumer, not process-memory certificates.

Invalid axes, shapes, provenance, arithmetic profile, masks and policy values
refuse. Domain, work, payload, unresolved refinement, nonnormal nonzero
products and narrowing failures also refuse. A failed row withholds its
dependent amplitudes while retaining attempted work, source identity and the
completed-cell prefix diagnostics; it is never dropped or renormalized into a
different source. Inspect result and row status before using optional outputs.

The fresh Release build at `aae8bded2f35af171cafb72ba3a078aa49dbbf07`
passed four native contracts: continuous projection, shared thermal conformal
epoch, old perfect-fluid transfer and its independent peer. Projection controls
cover analytic finite boundaries, independent angular Legendre/spin integration
and refinement, regular small-phase limits, signed sources, smooth-versus-linear
source-grid refinement, ownership and hostile input refusals. This is scoped
engineering and synthetic numerical evidence. The separate installed SDK
consumer checks a constant-monopole \(\ell=1\) finite-boundary integral against
an independently written trigonometric \(j_0\), masks, lifetime and refusals.
Both this installed consumer and the existing broad SDK consumer passed against
the same fresh installed archive with conservative floating-point flags. The
four native tests and two installed consumers remain separate execution records;
all six passed before publication.

Actual same-k CLASS source/transfer matching and support/truncation/precision
controls remain pending. Complete native photon/species/metric/opacity and
polarization source closure also remains open. This operator does not compute
\(C_\ell\), primordial k integration, lensing, a full CMB prediction or an
observational likelihood. See the broader
[primary CMB prerequisites](primary-cmb-projection.md).
