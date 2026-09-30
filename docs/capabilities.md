# Capabilities and limits

This page describes the public command surface at this documentation revision. Development moves quickly: `describe --json` and source at the exact commit take precedence. A native C++ function is not automatically a CLI operation.

| Surface | Available now | Limit |
| --- | --- | --- |
| `describe --json`, `version --json` | Product, build, ABI/schema, commands and capability metadata | No automatic scientific qualification |
| `fixture.checked_i64_add.v1` | Checked integer batches and infrastructure fault tests | Engineering fixture, not science |
| `quantity.convert.v1` | Tagged scale conversions with role/frame/convention metadata | No physical model or frame transformation; finite results currently unaccepted |
| `numerics.scalar_batch.v1` | Compensated sum, log-sum-exp, log1p, expm1, positive log-gamma | Bounded interface; finite results currently unaccepted |
| C++ numerical library | Integration and dense SPD/Gaussian building blocks with native tests | Specific sampled contracts, not unrestricted domain qualification |
| `observations.prepare.v1` | Typed local ASCII observations, full uncertainty retention, immutable prepared selection | Two explicit profiles; no matrix repair, likelihood, independence claim or registered qualification |
| C++ statistics library | Scalar normalized densities, retained Gaussian factors, explicit marginal/conditional operations and offset treatments | Tested scalar/2x2/3x3 fixtures; no qualified large survey likelihood |
| `statistics.gaussian_batch.v1` | Retained normalized densities, proper latent priors or explicit offset profile scores | Bounded native/CLI interface; completed outputs remain unaccepted |
| `background.parameter_query_batch.v1` | Compiled model-array × shared-query background batches | Bounded native/CLI interface; late-time models only, no survey likelihood |
| C++ supernova consumer | Retained observations/profile factor and compiled model batches returning relative offset-profile scores | Named synthetic/original-input comparisons; no normalized density or inference claim |
| `supernova.profile_batch.v1` | Retained released-profile observations and compiled relative offset-profile parameter batches | Named native/transport checks; arbitrary requests remain unqualified, no normalized density or inference claim |
| Local run store | Objects, resolved specifications, attempt receipts, hashes | Replay/cache/resume service not implemented |

The scientific qualification list in CLI discovery is currently empty. Native comparisons and regression tests are useful development evidence, but cannot be promoted into blanket acceptance of a scientific request. The integer fixture may be accepted as an engineering run while its discovery metadata still says `unqualified`.

## Direction of travel

Development follows the calculations agents need for cosmology: shared physics, physical models and the readers/numerical operations they depend on. Each new consumer needs validated prerequisites and an error budget. Fitting and orchestration can stay in external tools; internal statistical methods are added when a particular calculation needs them. This is not a checklist for a complete analysis suite.

GPU/distributed execution, broad cosmological analyses and a demonstrated performance advantage remain future work. Historical research results are not established merely by importing their ideas or documenting them.

## Known numerical limits

Adaptive quadrature can miss a narrow unsampled feature. Small SPD examples do not qualify arbitrary large or ill-conditioned matrices. Agreement through a shared system library is not independent numerical evidence. Domain rejection, including supported quantity roles and underflow policy, is part of the contract. See [testing](testing.md) and [fixture provenance](../cpp/tests/fixtures/README.md).

Observation profiles are `pantheon_plus_released_v1` and `gaussian_fixture_v1`. The former keeps released fitted summaries distinct from raw measurements and retains zHD, zCMB and zHEL separately. Compiled selection uses strictly zHD > 0.01. Repeated event IDs are preserved; unique measurement IDs and declared uncertainty axes determine row order. Unknown calibration/dependence remains unknown. A dedicated synthetic FITS codec test exists, but FITS is unavailable in the product build.

The native Gaussian API requires ordered measurement IDs for residual and design vectors. Full precision is converted before marginal selection; slicing precision alone would instead mean conditioning. Proper named latent priors retain their mean, variance, response, row IDs and independence assumption. Profiling returns a fitted coefficient and projected score, not a normalized density or evidence. Native comparisons are scoped prerequisites, not automatic runtime qualification.

For prepared observation Gaussians, covariance validation is scoped to the selected principal block. It makes no full-source joint probability claim. A precision input requires full inversion before marginal selection; neither path silently symmetrizes the source matrix.

Native `irred::cosmology::Background` implements named radiation-free flat LCDM and kinematic constant-q models. Bounded analytic and independent fixed-panel comparisons cover late-time z=0..5, H0=40..100 km/s/Mpc and named parameter fixtures. Outputs distinguish geometric same-redshift distances from the explicit released zHD/zHEL luminosity convention; the dimensionless shape contains no absolute H0 information. Native CPL is available through the explicit required-parameter `prepare_cpl` factory; existing v1 CLI/ABI rows admit only LCDM and constant-q. Curvature, age and predicted sound horizons remain absent. Quadrature estimates are empirical diagnostics. Native CPL comparison coverage does not automatically qualify arbitrary parameter requests.

`statistics.gaussian_batch.v1` exposes retained normalized Gaussian densities or single-offset profile scores, with ordered residual IDs/units and explicit proper-prior assumptions. Bounded native ABI allocation/ownership and CLI adversarial checks passed. Profile scores are not normalized densities; proper priors and profiling are mutually exclusive. Mathematical evidence remains limited to the separately tested native domains, and no actual W01 likelihood qualification is claimed.

Native dense solves have explicit `binary64_legacy_v1` and `longdouble_cpu_v1` policies. The latter is tested on native fixtures and the named original-input comparison; it does not silently change existing Gaussian/scalar CLI arithmetic. Precision policy, sensitivity screen, ordered row identities and score conventions belong to each calculation. Optional original-input comparisons use pinned hashes and permanent fixtures; no local development receipts are runtime dependencies.

The CPU triangular solve stores a transpose in the previously unused allocated upper triangle. Both precision policies keep the same multiplication/subtraction order. On one GCC 16.2.1 CPU host, three serial interleaved Release trials of the fixed original-data seven-point consumer gave median preparation 10.383→6.149 s, retained batch 0.3480→0.3056 s and total 12.028→7.644 s, with all quality checks retained and no material RSS increase. These numbers describe that fixture/compiler/machine, not a general performance guarantee.
