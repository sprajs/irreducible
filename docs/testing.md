# Testing

The native tests run from a source checkout using the fixtures included in the repository. External comparison software is not required.

## Optional directed interval reference

The test-only `IRRED_REFERENCE_GMP_MPFR` profile requires installed dynamic GMP/GMPXX and MPFR libraries. The integrated product build forces it OFF; these libraries are not production dependencies. Missing headers or dynamic libraries make this explicit profile unavailable. No libraries or original datasets are bundled.

The following separate reference build and no-assets controls have been executed with GMP/GMPXX 6.3.0 and MPFR 4.2.2:

```sh
.build-tools/bin/cmake -S cpp -B /tmp/irred-reference-profile-check -DCMAKE_BUILD_TYPE=Release -DIRRED_REFERENCE_GMP_MPFR=ON
.build-tools/bin/cmake --build /tmp/irred-reference-profile-check --target test_bao_piecewise_interval --parallel 1
.build-tools/bin/ctest --test-dir /tmp/irred-reference-profile-check -R '^bao_piecewise_interval_controls$' --output-on-failure
```

The executable target is excluded from the default build, so build it explicitly before running its control test. Controls challenge exact rational and directed interval arithmetic; they do not establish an original-data scientific gate. Original-data mode requires an external SHA guard and separate component-budget review. Pin actual library/header versions and hashes, compiler flags, reference source and raw input identities for that comparison. The current host's MPFR headers declare LGPL 3 or later; GMPXX headers declare LGPL 3 or later or GPL 2 or later. Review the installed notices for other enabled environments.

## Ordinary checks

After [dependency preparation](getting-started.md), run:

```sh
python3 tools/build.py
.build-tools/bin/ctest --test-dir build/native --output-on-failure
cargo test --locked --offline -j4
python3 tools/check_docs.py
python3 tools/check_install.py
```

C++ tests exercise scientific/numerical contracts and hostile ABI cases. Rust tests cover safe wrappers, request parsing, records and the CLI. `check_docs.py` checks tracked publication paths and relative Markdown file links; it does not validate equations, external websites or every Markdown feature.

## Exclusive build-identity test

This test changes source comments and CMake flags, then rebuilds and restores them. Run it in an isolated checkout, with no other build using that tree:

```sh
cargo test --locked --offline -j4 --test build_identity --no-run
```

Cargo prints the compiled test executable path. Run **that exact path**, followed by:

```text
--ignored --exact source_receipt_flags_and_cache_identity
```

Do not launch it through a nested Cargo test invocation. Ordinary Cargo testing intentionally ignores it. It checks that source changes affect build identity, mutable local files do not, unsupported overrides fail, and poisoned CMake flags are reset. Preserve the original checkout if the process is forcibly terminated before restoration.

## Fixture ancestry

[Native fixture provenance](../cpp/tests/fixtures/README.md) records the essential constants, derivations and budgets. Keep small explicit expected values and independent algorithms with the tests. Preserve the distinction between exact algebra, high-precision comparisons and shared-library checks.

Full external-tool output and development receipts stay local. Once a comparison becomes a durable regression, routine tests should not require rerunning the external software. Do not delete original historical research or accepted records merely to make the public tree smaller.

A passing test is scoped evidence. It does not qualify arbitrary inputs, every supported toolchain, cosmological inference or a performance claim. Quadrature missing a narrow unsampled feature is a known limitation; a finite estimate is not a proof of accuracy.

## Observation and optional codec tests

Ordinary native tests include immutable observation ownership, exact masks/row order, uncertainty asymmetry retention, hostile descriptors and scoped allocation-failure cleanup. Rust ingestion/CLI tests retain raw source bytes, challenge same-path mutations, semantic failures and nonfinite/missing distinctions. They use small synthetic inputs and require no acquired survey files.

The product build deliberately reports FITS unavailable. A separate optional native codec profile requires installed CFITSIO and can be tested with `-DIRRED_TEST_CFITSIO=ON` in a separate CMake build directory. It covers one handcrafted synthetic BINTABLE profile, not universal FITS support. The wrapper resets this test-only option to OFF for its supported product build. CFITSIO is file-format infrastructure; its installed version/library identity and NASA permissive notice must be reviewed for any enabled deployment.

`check_install.py` installs to a fresh temporary prefix, rejects the old generic include root, and compiles/runs the durable `test_installed_consumer.cpp` against the installed `libirred_core.a`. It checks independent library usability and cleans its temporary installation.

## Explicit Release profile

Build and test the optimized profile with the same native and Rust cases:

```sh
python3 tools/build.py --profile release
.build-tools/bin/ctest --test-dir build/native-release --output-on-failure
cargo test --release --locked --offline -j4
target/release/irred describe --json
```

Debug remains the default and uses `build/native`, `target/debug` and `build/build-manifest.json`. Release uses `build/native-release`, `target/release` and `build/build-manifest-release.json`. Cargo selects the matching native archive and embeds the matching manifest; do not manually cross-link profiles. Both retain strict floating-point options and panic containment. Release adds C++ `-O3 -DNDEBUG` and the pinned standard Cargo release settings; its source/build/profile identity differs from Debug. Native checks remain active under NDEBUG. These build checks establish no workload speedup.

After building both profiles in an isolated checkout, run the durable identity/negative witness:

```sh
cargo test --locked --offline -j4 --test profile_identity -- --ignored
```

This checks each executable's describe/version manifest, identical source digests with distinct profile/build/native-archive identities, and a deliberately wrong native fixture compiled with NDEBUG. That fixture must fail, proving its required checks were executed. It compiles one temporary native test, cleans it after success, and never mutates production source.

## Native and public regression suites

Ordinary CTest includes the current expansion, retained supernova, BAO and Gaussian boundary suites, independent analytic/quadrature/cofactor controls, ownership/allocation failures and payload accounting. Ordinary Cargo tests include current one-shot and retained-stream requests, source identities, quota failures and immutable run records. Request schema and C ABI revision are 2; operation names are canonical, while scientific method and source identifiers retain their declared identities. Scientific qualification remains bounded to named tested cases.

`test_w01_reference` is an optional independent LDLT condition/reconstruction control, not a retired consumer interface. Existing native fixtures retain historical direct/stable and compressed-route discrepancies separately. No expected values are regenerated during builds.

## Optional original-data native regressions

Original datasets are not bundled and are unnecessary for ordinary CI. Explicit ignored Rust guards hash each original file before and after executing a prebuilt native harness, print its executable hash, and verify its named suite and case count. No Rust physical equations are used.

Build the optional native harnesses in the selected profile:

```sh
.build-tools/bin/cmake --build build/native-release --target test_supernova_reference test_bao_reference test_w01_reference --parallel 2
```

For supernova, set `IRRED_W01_TABLE`, `IRRED_W01_COVARIANCE` to the exact assets pinned in `cpp/tests/fixtures/w01_historical.hpp`, and `IRRED_W01_NATIVE_HARNESS` to `build/native-release/test_supernova_reference`. Then run the explicit named29 guard:

```sh
cargo test --release --locked --offline -j2 --test w01_reference released_assets_and_native_twenty_nine_points -- --ignored --nocapture
```

This explicitly enables `IRRED_SN_ALL_NAMED=1` inside the guard: historical seven controls, CPL four, fixed-five-bin six and grey twelve. All use the current 1590 selected source. Grey12 is a conditional application to these inputs, not a replay of the historical 1820-source workflow. The default historical-comparison guard deliberately selects only grey12 and cannot stand in for all29.

For BAO, set `IRRED_BAO_MEAN`, `IRRED_BAO_COVARIANCE` to the exact assets pinned in `cpp/tests/fixtures/bao_reference.hpp`, and `IRRED_BAO_NATIVE_HARNESS` to `build/native-release/test_bao_reference`:

```sh
cargo test --release --locked --offline -j2 --test bao_reference released_assets_and_native_eleven_point_comparison -- --ignored --nocapture
```

The directed fixed-five-bin24 oracle is a separate optional native target `test_bao_piecewise_interval`, enabled by the existing CMake GMP/MPFR test option. GMP/MPFR are test-only dependencies and are never required by production. Point the separate `IRRED_BAO_INTERVAL_HARNESS` variable at that executable and run:

```sh
cargo test --release --locked --offline -j2 --test bao_reference released_assets_and_native_directed_twentyfour_point_comparison -- --ignored --nocapture
```

The guard requires all24 distinct point records and a zero-failure directed suite marker. Independent references and original-input comparisons establish named numerical evidence; they supply no posterior, joint-probe or parameter-campaign qualification.

Use `python3 tools/check_install.py --profile release` for a fresh installed consumer against the chosen Release build. The helper never chooses a profile merely because its files exist.

## Durable comparison coverage

| Accepted external or historical comparison | Durable native coverage |
| --- | --- |
| Smooth/CPL analytic limits and independent fixed-panel quadrature | Background owner and independent hostile tests |
| Fixed-five-bin historical expansion facts, independent split quadrature, boundary conventions | Piecewise owner and hostile tests; requested groups retain separate outcomes |
| Historical selected covariance, direct/stable profile targets, approximate compressed discrepancy | Frozen W01 headers, independent LDLT controls and named29 harness |
| CPL4, fixed-q6 and Pantheon1590 grey12 fixed-point allocations | Named29 harness; fixtures retain their separate source/model lineage |
| Extracted original-class fixed-coefficient basis/profile comparison | Grey native fixture assertions; native geometry ancestry remains explicit |
| Released BAO11 independent quadrature/LDLT and fixed-q24 directed reference | BAO11 harness and optional GMP/MPFR directed24 oracle |
| Correlated tiny Gaussian/profile controls and input permutations | Ordinary statistics/SN/BAO independent native peers |

Optional original-data hash guards test acquired source identity as well as native fixture agreement. Ordinary native tests retain small accepted numerical facts and provenance without requiring Python, historical environments, original datasets, Cobaya or GMP/MPFR. Test-only external algorithms remain explicitly distinct from production kernels. An approximate historical route that failed a precise budget is retained as a discrepancy, rather than promoted to a new expected value.

Ordinary consumer controls include diagonal weighted common-intercept elimination and a correlated 3×3 BAO cofactor density, with source order and failure checks. Historical performance receipts describe their recorded builds; they are not current speed claims.

The separate ignored `original_cli` suite checks representative current CLI ingestion against native transcripts on the original sources: one LCDM supernova point and one BAO point, including raw identities, selection/order and stored records. The SN native reference reader uses CID/survey/row labels while current acquisition uses raw-table-SHA/row IDs; exact numerical parity and independently audited current source ordering are distinct from full metadata-object equality. The BAO reference reader likewise uses its historical labels while current acquisition binds raw-mean-SHA/row/observable IDs; both namespace differences are disclosed rather than normalized. This is transport parity for those cases; it does not establish a full parameter campaign. Work priorities and remaining capabilities are listed in [the roadmap](roadmap.md) and [the gaps](gaps.md).

To repeat the representative public check, first capture stdout from the matching native guards above. The executable and native harnesses must use the same selected Release archive and numerical settings; the test checks the recorded harness digest. Supply all eight path variables explicitly:

```sh
IRRED_W01_TABLE=/path/to/Pantheon+SH0ES.dat \
IRRED_W01_COVARIANCE=/path/to/Pantheon+SH0ES_STAT+SYS.cov \
IRRED_SN_NATIVE_OUTPUT=/path/to/guarded-native-grey12-stdout.log \
IRRED_W01_NATIVE_HARNESS="$PWD/build/native-release/test_supernova_reference" \
IRRED_BAO_MEAN=/path/to/desi_gaussian_bao_ALL_GCcomb_mean.txt \
IRRED_BAO_COVARIANCE=/path/to/desi_gaussian_bao_ALL_GCcomb_cov.txt \
IRRED_BAO_NATIVE_OUTPUT=/path/to/guarded-native-bao11-stdout.log \
IRRED_BAO_NATIVE_HARNESS="$PWD/build/native-release/test_bao_reference" \
cargo test --release --locked --offline -j2 --test original_cli -- --ignored --nocapture
```

## Extending numerical and inference coverage

Preserve a first divergence before changing a model, convention, input or expectation. Required outputs cannot be relabelled optional after failure. Deterministic dependency biases are not independent random errors to combine in quadrature; observable and log-density sensitivity must justify their allocations. Refinement agreement is empirical evidence, not an automatic certified bound or proof of asymptotic convergence.

When adding differential equations or interpolation, test event location, dense output, stiffness/inner solves and derivative accuracy separately. An iterative solver's true residual need not equal its recurrence residual. Tail range and resolution both require checks. Optimization KKT conditions do not calibrate statistical tails; nominal sampling diagnostics do not prove support or identifiability.

Fixed-truth recovery differs from prior-predictive simulation-based calibration. Shared simulator/reference ancestry can conceal the same physics defect. Reweighting weight-tail collapse, absent support and unresolved posterior-predictive tails are failures to report rather than repaired results. These requirements guide future implemented consumers; they do not imply an inference engine currently exists.

Wide long-double reference suites assert at least 64 mantissa bits. The tested Linux/GCC platform meets this requirement; a target where long double equals double is a known unsupported qualification profile, not an observed run failure and not a reason to weaken budgets. Portable CPU refers to the baseline implementation, not proof that every host/compiler arithmetic profile is validated. Ubuntu CI adds engineering evidence only for its actual run.

## CI checks and build artifacts

The `CI` workflow separates native C++ GCC/Clang tests, Rust module unit tests, CLI integration tests, Release application/package checks, and documentation/generated contracts. Native jobs build independently with CMake; Rust jobs build a linked native baseline. Cargo metadata supplies the test target inventory. The current Rust crate has a binary with module unit tests and separate integration targets; a future library also receives unit and documentation tests. Original-asset and exclusive mutation checks remain explicitly ignored in ordinary runs, not claimed as historical reproductions.

Ubuntu 24.04 x86_64 is the CI binary scope. GCC and Clang checks are engineering evidence within the Linux wide-arithmetic contract, not blanket qualification on other platforms. Matrix failures do not cancel the other compiler's evidence. Local checks still share the four-job budget.

Cargo.toml owns the deliberate SemVer source version, currently `0.1.0-dev`; Cargo.lock agrees. Remain at `0.1.0-dev` until an intentional release decision. Version bumps are reviewed source changes: edit Cargo.toml and update the root package entry in Cargo.lock, never make bump commits from CI. For deliberate releases before 1.0, choose a patch for compatible corrections and a minor version for new capabilities or interface/meaning changes; document breaking changes rather than adding speculative compatibility routes. No automatic release or tag is created. The product version is separate from request/ABI revisions, equation/method identifiers and the content-based build_id. Discovery/version and the build manifest retain that source version; none confers scientific qualification. A non-version CI execution label and separate metadata record run/attempt, checked-out SHA and PR head/base: a tested PR merge commit is not its source head. Changing a run attempt does not replace the existing build hash. See [Cargo's version field](https://doc.rust-lang.org/cargo/reference/manifest.html#the-version-field), [SemVer pre-releases](https://semver.org/spec/v2.0.0.html), and [GitHub's contexts](https://docs.github.com/en/actions/reference/workflows-and-actions/contexts).

The Release job tests optimized native assertions and CLI integration, checks a fresh native installation, then uploads a tarball containing the executable, native library/public headers, build manifest, CI metadata, upstream runtime notices and checksums. This is a scoped CI artifact, not a portable installer or published release. Missing payload or inconsistent version/build/checkout metadata fails packaging. Post-download smokes check the actual GitHub artifact separately.

Uploaded CI artifacts expire after 14 days. They do not replace durable Git branch backups or immutable scientific receipts. Dependency and Rust-library notices are copied faithfully with a hashed inventory; their presence is provenance, not a new licence interpretation. The package refuses a host outside the declared Ubuntu scope.

Finite passband-calibration regressions retain analytic constant/linear radiometry, independently integrated frequency-coordinate controls and separately accumulated weighted moments. Refinement occupies at most 5% of the named comparison allocation. Shared-state covariance, permutations, constant witnesses, unresolved spread, partial failures and resource quotas are challenged without original calibration assets or external engines.

Conditional BAO regressions compare independent direct-redshift/sqrt(a) fixed-panel quadrature and cofactor/LDLT normalized Gaussian calculations. They retain separate ratio and 1e-8 density-component allocations, reference refinement, H0 cancellation, projection refusal, positive diagnostic underflow, output omission, lifetime and global callback quotas. Original released assets remain optional local comparisons, with exact source hashes and conditional fitted-summary roles.

Thermal-relic regressions retain exact massless and nonrelativistic limits,
collisionless continuity, independent high-precision momentum integrals and a
matched CLASS background grid. Explicit temperature/statistical-weight and
constant conversions prevent external defaults from masquerading as agreement.
The prepared owner's copies, moves, self moves and failed source use are tested;
the original moved-source failure remains preserved in local research evidence.
These controls cover species states and E/H, without qualifying distances,
recombination, drag, perturbations or a Planck posterior.

Conditional estimator-variance controls compare analytic contrasts, independent
cofactor/KKT calculations and summed deterministic Gaussian-noise responses.
They challenge covariance/parameter permutations, unit scaling, retained
ownership, quotas and positive unrepresentable outputs. Synthetic H0 sampling
quantiles use independent high-precision inversion of the exact binary64
probabilities, including unequal complementary endpoint tails. They test the
explicit fixed-design sampling law, without establishing a posterior,
observational H0 uncertainty or confidence coverage. The original failed-wrapper
move-status witness is preserved separately from successful-owner controls.

## Thermal observable composition

Ordinary native tests exercise the retained thermal physical mapping and flat
distance/supplied-drag ruler consumer, its independent analytic/high-precision
and matched CLASS controls, masks, limits and owner lifetimes. See
[thermal observables](thermal-observables.md) for the equations, physical inputs,
frozen allocations and reference ancestry. External engines, acquired sources
and Python high-precision libraries are comparison tooling; ordinary CI consumes
small declared facts and independent native controls. These tests add no CLI
operation or full-reference observational qualification.

## Coupled hydrogen temperature and finite optical histories

Ordinary native tests retain the prescribed-temperature history controls and
add an independently authored, coupled Radau IIA reference for the pure-H
Compton/adiabatic variant. Named physical models and near-initial stiff-layer
queries have separate Kelvin, depth, rate, visibility and survival allocations.
Constant-coefficient thermal limits and exact finite-opacity controls check
the reference; mesh refinement and residuals qualify its comparisons. Common
atomic/source constants and long-double/libm ancestry remain explicit.
Finite visibility integrates to one minus the boundary survival, without
renormalization. Quotas, omitted outputs, partial refusal, source ownership and
allocation failures are separate engineering gates. Read
[the history guide](recombination-drag.md) before composing this native SDK.
These tests do not qualify helium, a physical drag epoch, present-day CMB
visibility, observational uncertainty or a full LambdaCDM recombination model.

## Explicit relics in the conditional pure-H history

The direct-SI reference independently integrates FD density with composite
GL16/32 quadrature and an explicit momentum tail, then advances the declared
history with resolved coupled Radau stages. Three evolved mass/temperature
profiles and the prescribed-temperature initial boundary test the applicable
outputs and conditional roots. Independent opacity/visibility quadrature retains
survival mass. Split-weight and three-positive-species availability are distinct
owner controls, not extra independent history comparisons. Fixed-density H0
controls test the retained expansion law with explicit cancellation diagnostics.
All momentum work is charged; the default cap can refuse these profiles.

A linker-only forwarding observer of the unchanged native integration entry
point verifies actual momentum callbacks at each measured provider/history
allocation site. Source acquisition precedes momentum work and successful
transfer cannot allocate. Linux GCC/Clang with GNU-compatible ELF linker wrapping is the scope of
this fault instrumentation; the installed scientific library uses ordinary linkage.
Successful values/counters and affected distance/ruler/BAO consumers are checked
separately from failure accounting.

## Retained joint Gaussian predictive controls

Original covariance-route/cofactor and exact rational controls check predictive
mean and full covariance. Direct prior-times-training-times-future integration
checks three future vectors, including the predictive mean; independent native
joint-density quadrature checks normalization. Null/rank-deficient training,
shared/repeated responses, internally correlated future noise, coherent parameter
permutations and unit Jacobians challenge the fixed synthetic contract. Quotas,
lifetime, all measured allocation sites and independent noise/event declarations
remain separate gates. See [Gaussian predictive](gaussian-predictive.md); these
controls do not qualify observational coverage or a released finite-box target.

## Addressed Gaussian recovery campaigns

Ordinary CI runs `gaussian_simulation_contract` and
`gaussian_simulation_peer_contract` for replay, original mathematical facts,
source order, retained-factor colouring and hostile resource/lifetime cases.
The `gaussian_recovery_contract` and `calibration_recovery_contract` tests also
run the frozen two-seed campaigns in ordinary CTest and required compiler CI
under [the simulation contract](gaussian-simulation.md).
Preserve the raw receipt for every attempt and failed cell. Ideal
Gaussian moments/coverage, finite-grid approximation, arithmetic and empirical
discrepancy have separate allocations. Seeded comparisons are not an
independence certificate or a universal coverage guarantee.
Each campaign emits its complete attempt ledger; the named local runs produce
approximately 43 MB and 68 MB of stdout. CTest capture adds this payload to
the executable's memory use. Runtime and peak-memory comparisons must include
logging/capture separately from the retained numerical factors.
