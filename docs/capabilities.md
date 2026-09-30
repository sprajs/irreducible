# Capabilities

Use the executable's `describe --json` output for its exact build, ABI revision and request schema. Irreducible has one current interface; obsolete model-specific routes are removed together with their callers. Historical records retain their original identities and meanings.

| Operation | Calculation | Scope |
| --- | --- | --- |
| `fixture.checked_i64_add` | Overflow-checked integer addition | Exact integer contract |
| `quantity.convert` | Typed physical conversions | Explicit units, roles, frames and conventions |
| `numerics.scalar_batch` | Compiled scalar numerical methods | Method-specific arithmetic and domain |
| `cosmology.sound_horizon` | Conditional comoving sound horizon | Supplied drag redshift; flat pressureless matter + massless radiation + Lambda |
| `background.evaluate` | Requested flat-FLRW expansion and projections | LCDM, constant q, CPL or fixed five-bin q |
| `observations.prepare` | Immutable typed source preparation | Structural checks; no probability or lineage upgrade |
| `statistics.gaussian` | Normalized Gaussian density, offset profile or proper latent prior | Explicit ordered residuals and nuisance assumptions |
| `supernova.profile` | Conditional single-offset magnitude profile | Same expansion types, with no source effect or an explicit grey magnitude effect |
| `bao.density` | Conditional normalized free-ruler Gaussian density | Same expansion types and explicit H0rd |
| `photometry.predict` | Deterministic flat-spectrum rectangular-passband prediction | Supplied distance/redshift; incident flux, collected energy and transmitted photon expectation |

Background requests select expansion, radial, luminosity_shape, clock, physical or kinematics groups. Dimensionless E and clock integrals remain usable when a requested dimensional H or lookback-time projection fails. Observer and physical-scale failures do not erase independent outputs. E-only and DH-only calculations avoid radial integration. Exact redshift bits share nodes within a model batch; repeated model specifications across model rows currently recompute their geometry under the same global work budget.

Observation preparation distinguishes measured quantities, released fitted summaries and synthetic controls. The generic magnitude-covariance profile does not confer Pantheon release identity, calibration, standardization or independence. Source hashes, axis declarations and selection remain explicit. A retained process owns the same immutable observation object through its dependent consumers, even after the parent handle is released; factors are prepared once and reused across evaluations.

Gaussian covariance validation applies to the selected principal block. Precision input requires the full precision factor and inversion before marginal selection: slicing precision would mean conditioning. Ordered residual and response IDs must match. A proper named latent prior retains mean, variance, response and independence declaration. An offset profile is a relative score, not normalized density or evidence.

Native supernova comparisons cover the named original seven-point, CPL four-point, fixed-q six-point and conditional grey twelve-point cases on the pinned 1590-row sample. The grey cases do not reproduce a historical 1820-row analysis. Stable/full historical fixed-q comparisons pass their allocations; five compressed approximations fail the reference allocation and remain withheld. These facts do not qualify a posterior, an arbitrary supplied sample or every point in the model domain.

Native BAO comparisons include the named eleven-point and fixed-q twenty-four-point cases on pinned thirteen-row inputs. The latter use an optional test-only directed GMP/MPFR interval oracle. The unchanged log-density budget is 1e-8, allocated as reference 2e-9, factor/solve 3e-9 and background 5e-9. Earlier conservative failures remain historical evidence; the directed certificate establishes the named cases, with a narrow margin at the highest-residual endpoint. Production has no GMP/MPFR dependency. SN and BAO remain separate conditional consumers; no joint likelihood or source independence is inferred.

Required numerical checks that pass satisfy default `numerical_contract` assurance and exit 0. Explicit `qualified` assurance exits 6 when applicable named evidence is absent. Completed scientific failures exit 2. Interpretation remains unqualified; no arbitrary request receives automatic scientific qualification.

Measured mirrored dense preparation improved one pinned fixture on one compiler/machine; that historical measurement does not establish performance of the consolidated interface. Optimization follows matched-quality measurements, not the presence of a wider type or more threads. Samplers, priors over cosmological parameters, smoothing campaigns and whole-domain certificates remain outside these capabilities.

The [conditional sound horizon](sound-horizon.md) uses an exact compact scale-factor interval, with explicitly supplied photon/baryon fractions and drag-redshift provenance. It predicts neither thermal history nor drag epoch; no existing BAO likelihood qualification is inherited.

See [photometry](photometry.md) for the bounded deterministic projection.
