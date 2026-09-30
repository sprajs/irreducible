# Fixture provenance and budgets

These are small, intentional native test assets. Their expected values are independent of the production C++ implementation. Full comparison runs and generator environments stay local; the derivations below preserve the meaning of the committed fixtures.

## Constants and high-precision values

`foundations_oracles.hpp` freezes 43 reference values. Independent Python standard-library Decimal calculations at 70 and 100 digits agreed at 50 displayed significant digits. Pi was constructed using Machin's identity, `pi = 16 atan(1/5) - 4 atan(1/239)`, summing the alternating power series with next-term truncation bounds. This is a precision comparison with truncation bounds, not a general interval-arithmetic proof. No production conversion code or astronomy engine was imported.

The constant set `SI-IAU-definitions-v1` uses:

| Quantity | Definition | Source |
| --- | --- | --- |
| Speed of light | Exactly 299792458 m/s | [BIPM metre definition](https://www.bipm.org/en/si-base-units/metre) |
| Astronomical unit | Exactly 149597870700 m | [IAU 2012 resolution B2, recommendation 1](https://iauarchive.eso.org/static/resolutions/IAU2012_English.pdf) |
| Parsec | Exactly `648000 au / pi` | [IAU 2015 resolution B2, note 4](https://iauarchive.eso.org/static/resolutions/IAU2015_English.pdf) |
| Megaparsec | `10^6 pc` | SI prefix applied to the preceding definition |

These are unit definitions, not measured cosmological parameters. The fixture `H0 = 70 km/s/Mpc` is a chosen synthetic input, giving `H0 = 70000/Mpc` in inverse seconds and `c/H0 = 4282.7494 Mpc`.

The header retains both decimal strings and binary64 rounded values. A decimal input close to a singularity and its rounded binary64 input can denote materially different mathematical problems. Compare the represented input, not a superficially similar decimal label.

## Analytic numerical cases

- `integral_0^1 exp(x) dx = e - 1`; `integral_0^1 1/(1+x) dx = log(2)`.
- `log Gamma(1/2) = log(pi)/2`; `log Gamma(5) = log(24)`.
- The toy constant-q integral is `I(z,q) = integral_0^z (1+t)^(-1-q) dt`. For nonzero q it equals `(1-exp(-q*log(1+z)))/q`; at q=0 it is `log(1+z)`. Frozen points use z in `{0.001, 0.5, 3}` and q in `{-1, -0.5, -1e-10, 0, 1e-10, 0.5, 1}`. These are mathematical fixtures, not an accepted cosmological model.
- For `C = [[4,1],[1,9]]` and residual `(2,-3)`, adjugate algebra gives `det(C)=35`, `C^-1=[[9,-1],[-1,4]]/35`, `C^-1 r=(3/5,-2/5)` and `r^T C^-1 r=12/5`. The normalized Gaussian log density is `-(2 log(2 pi) + log(35) + 12/5)/2`.

The rational cases in [gaussian_exact_cases.hpp](../../../tests/verification/gaussian_exact_cases.hpp) use explicit adjugate/scalar algebra (initially cross-checked with Python Fraction), independent of production factorization. For a flat common offset, its precision is `11/35`, fitted value `7/11`, and projected chi-square `25/11`. A proper common nuisance `t ~ N(1/3,2)` gives marginal covariance `[[6,3],[3,11]]`, determinant `57`, quadratic `1175/513`, posterior mean `77/171` and variance `70/57`. Prepared fixture values do not themselves establish acceptance of a future consumer.

## Predeclared budgets

| Fixture consumer | Budget | Scope/rationale |
| --- | --- | --- |
| Quantity conversions | Relative `1e-13` | Synthetic ten-significant-digit consumer; sampled Mpc values `1e-12` through `1e12`, normal binary64 outputs |
| log1p/expm1 points | Absolute `1e-27` + relative `3e-15` | Named stable-function points; compare identical represented inputs |
| log-gamma half/integer | Absolute `1e-14` | Named analytic identities |
| Toy integrals | Absolute `1e-12` + `1e-10 * abs(reference)` | Sampled q/z cases and refinement, not arbitrary integrands |
| Tiny Gaussian log density | Absolute `1e-10` | Synthetic small-system consumer: quadratic/logdet allocations `4e-11` each, normalization `2e-11` |

The tiny Gaussian contract targets n <= 16, condition number <= 1e6 and whitened residual norm <= 10; named fixtures exercise parts of that intended scope, not a proof over the whole domain. Larger or ill-conditioned consumers need their own propagated budgets.

## Discrepancies and limits

A cancellation regression integrates `x^6 - rounded(1/7)` on `[0,1]`. A local relative-error acceptance rule could accept a wrong near-zero integral; the correction checks the final global estimate against the unchanged absolute/relative policy. Preserve the regression instead of adjusting its expected value or tolerance to the old answer.

A sufficiently narrow unsampled peak can evade an adaptive quadrature estimator entirely. Tests retain that limitation; a finite estimate is not a certified bound. Shared libm calls, shared constants and similar algorithms are disclosed dependencies, not independent confirmations.

See [testing](../../../docs/testing.md) for actual commands and [scientific contracts](../../../docs/scientific-contracts.md) for qualification rules. These notes describe fixture ancestry and intended budgets; they are not run receipts or blanket scientific qualification.

## Large synthetic rank-one covariance

`test_large_spd` uses exactly dyadic `C = I + 11^T/64` and `r_i = ((i mod 11)-5)/4`. Sherman–Morrison gives `C^-1 r = r - 1 sum(r)/(64+n)`; the determinant lemma gives `det(C)=1+n/64`. Integer sums produce the reference quadratic without calling the production factorization or reduction. The exact infinity-norm condition number is `1+(2n-2)/64`.

The ordinary suite runs `test_large_spd --size64`. Its default n=1590 execution is explicit and expensive, not an ordinary test workload. The fixed assembled log-density absolute budget is `1e-8`, intended for a synthetic consumer reporting roughly seven decimal places. Component checks allocate `2e-9` each to quadratic and log determinant, `4e-9` to normalization, and `1e-10` to maximum solution component error. Backward residual and forward sensitivity are diagnostics, not proved bounds. Reference logarithms share standard-library ancestry; differences below `1e-12` are not asserted as proven accuracy digits.

This tests the named rank-one family only. It is not a survey covariance, a likelihood reproduction, or qualification of arbitrary large SPD matrices. The numerical ABI tests separately compare direct native and coarse batch values/statuses/diagnostics, and explicitly disclose their shared numerical-core ancestry.

## Observation preservation fixtures

The synthetic ASCII example retains two distinct measurement IDs for one repeated event, an explicitly ordered `[[4,1],[1,9]]` covariance, and magnitude-labelled synthetic values. Native adversaries separately preserve an asymmetric/singular full matrix at the preparation layer; that layer claims no SPD validation or likelihood. Selection is strictly zHD > 0.01, with exact masks and source indices; no event deduplication occurs.

`handcrafted_fits.hpp` constructs independent big-endian BINTABLE bytes: signed integer LENGTH with TSCAL=0.5, TZERO=1, integer null sentinel, explicit metre unit, repeated event IDs and raw QUALITY bits. AUX is unsupported auxiliary content retained in original bytes, not scientific output. The narrow codec rejects scaled QUALITY to avoid silently changing bits. Disabled codecs return unavailable. These fixtures qualify neither a universal source reader nor a physical dataset.

## Probability and Gaussian contracts

Native statistics tests distinguish finite log densities, genuine outside support, invalid inputs and numerical failures. Named normal/Poisson identities and a selected standard normal on strictly positive support challenge normalization. Gaussian cases use exact 2x2/3x3 adjugate algebra, including full precision converted before marginal selection and an explicitly zero-centered-complement conditional target.

The normalized tiny-system density budget remains `1e-10`; q/logdet/normalization allocations are `4e-11`, `4e-11`, `2e-11`. The proper-prior case with mean `1/3` and variance `2` is independently checked by fixed-panel direct integration with analytic tail accounting, separate from production covariance assembly. Its integration check uses a `1e-9` normalization budget; system-libm ancestry is disclosed. Prior records describe original generative assumptions, not conditional latent posteriors. Ordered-ID mismatches and repeated latent identities reject rather than silently reusing a prior or permuting a response.

An initial profile coefficient reference transcription of `1/11` was incorrect: exact inverse algebra gives `(1^T C^-1 r)/(1^T C^-1 1)=7/11`. The native fixture uses the derived value; its mathematical inputs and tolerances did not change. These cases do not qualify the released 1590-row likelihood.

Selected observation covariance is checked and factored as the declared principal block; the full raw source matrix is retained but its full joint probability validity is not assessed. Durable native cases place asymmetry outside versus inside selection, retain source bytes, and require full precision validation/inversion before marginal selection. This is a scope declaration, not silent symmetrization or a repair of an invalid selected block.

The Gaussian CLI smoke fixture uses C=[[4,1],[1,9]] and residuals [1,2], [0,0]. The exact inverse gives q=21/35=0.6 and q=0; logdet=ln35 and normalization=2ln(2pi). It is a transport/control fixture, not independent math qualification or a W01 comparison. The independent native statistics cases challenge the probability algebra. The bounded Gaussian interface gate covers hostile allocation/layout/ownership and CLI identity/metadata/failure regressions. These transport tests share the native core and do not replace its separate independent mathematical references.

Background fixtures derive de Sitter, Einstein-de Sitter and constant-q distance/clock limits analytically, with exact zero-redshift/scaling/convention checks. Independent long-double fixed-panel quadrature challenges the production adaptive Simpson path and checks refinement before spending the mixed abs2e-15+rel2e-10 integral budget and 1e-8 mag shape budget for named z>=1e-4 cases. These are bounded late-time radiation-free model cases, not age/CMB or actual W01 score evidence. An initial denormal tolerance test expected depth_limit; the executed work-budget witness instead yields work_limit after three callbacks. Its failure is preserved locally and the permanent native case now checks the actual earlier work gate without changing scientific tolerance.

Cached offset profiling retains one Gaussian factor and the solved design/Gram term. `test_profile_operator_hostile.cpp` challenges it against the uncached calculation, ordered IDs, prior exclusion and a tighter conditioning limit. Gaussian ABI/CLI regressions require an OK numerical cause on finite profile scores and preserve `conditioning_budget_exceeded` on failed solves; profile backward/sensitivity fields describe the adjusted-residual solve and are empirical diagnostics, not certified score-error bounds. No performance advantage is claimed by these checks.

The background C ABI tests compare parameter-major/query-minor batches to direct native calls; that parity shares the compiled physics core and is transport evidence only. Hostile tests cover descriptor/resource guards, retained source values, mixed failures, one global callback budget and five allocation-failure points with cleanup. CLI tests check equation/unit/convention identities and absent physical payloads on failed rows. A first CLI test helper expected the observation-specific `result.calculation` wrapper; background returns its calculation directly in `result`. Correcting that structural test assumption changed no physical equation or tolerance.

## Explicit arithmetic and original-input supernova cases

`longdouble_cpu_v1` is a named native policy, not a silent replacement for the binary64 baseline. It retains promoted factor/solve arithmetic, checks intermediate range failures and preserves the caller's floating-point environment. Owner and independent tests cover exact SPD fixtures, invalid policies, range/conditioning failures and cached profile agreement. Empirical sensitivity remains a screen, not an interval proof. Normalized and profiled quadratic payloads reject nonzero values that cannot be represented normally in binary64.

The native supernova consumer retains selected observations and one profile operator, then evaluates model batches. Its prediction is `5 log10((1+zHEL) I(zHD))`; H0 is absorbed by a free additive magnitude offset. It reports a relative profile score, with no normalized-density, evidence, absolute-calibration or H0-identification claim. Synthetic native controls and hostile cases remain explicitly synthetic.

`ldlt_reference.hpp` supplies an independent unit-LDLT algorithm on the same canonical binary64 input promoted to long double. The ordinary exact 3x3/dyadic checks need no external assets. Optional W01 checks use seven frozen LCDM/constant-q points and three quadrature policies, with an assembled absolute score budget of `1e-6`: reference `2e-7`, factor/solve `3e-7`, background `5e-7`, plus a magnitude screen of `1e-8`. Floating-point cast allowances and current selected-matrix inverse-norm diagnostics are included; shared libm/reference authorship limits remain explicit.

`w01_historical.hpp` pins the original table/covariance hashes, read-only historical implementation identity, software versions and seven results. The original GL48 full and Cholesky-stable score routes are compared separately. The historical 40-knot compressed route is approximate and fails its precise-reference allocation at some points; it is retained as a discrepancy rather than substituted for the direct calculation. `tests/w01_reference.rs` verifies original asset hashes with existing SHA-256 infrastructure before invoking the optional native harness. No survey files or external Python environment are required by ordinary tests.

The supernova ABI/CLI fixtures are generated release-profile-shaped transport controls, not released survey assets or independent observations. The CLI control declares three rows with one excluded by the strict zHD cut; its selected covariance is `[[4,1],[1,9]]`, while the retained full matrix deliberately fails full probability validity. In the de Sitter limit the shape is `(1+zHEL)*zHD`; the independent two-row profile formulas give `q=(r0-r1)^2/11` and offset `(8*r0+3*r1)/11`. The fixed `1e-10` CLI comparison tolerance challenges transport and profile composition, not arbitrary survey qualification. Tests retain original bytes and validate altered same-path identities, explicit precision, mixed failures, empty batches, numerical causes and omitted arrays.
