# Eight new native models, 2026-10-09

This campaign adds eight distinct compiled physical models. The first six share
the original scientific build below; Plummer and Bianchi I have separate later
build identities. It advances new model
coverage rather than changing historical cosmology runs or their source/runtime
pins. Each guide declares its own assumptions and domain:

| Model | Native outputs and physical scope |
| --- | --- |
| [Exponential quintessence](quintessence.md) | Canonical scalar, conserved dust/radiation, flat GR E/H, composition and optional past distances |
| [DGP growth](dgp-growth.md) | Self-accelerating flat dust branch, E/H, matter fraction, quasistatic effective coupling and growing mode |
| [NFW halo](nfw-halo.md) | Untruncated spherical enclosed mass, surface density and conditional thin-lens observables |
| [Decaying matter](decaying-matter.md) | Forward finite-anchor GR parent decay into massless daughters, proper time and composition |
| [Curved FLRW](curved-flrw.md) | Conserved dust/radiation/Lambda, nonflat H and null distances, turnaround/antipode refusal |
| [Hernquist sphere](hernquist-sphere.md) | Finite-mass SI density, mass, Newtonian potential, acceleration, circular speed and potential Hessian |
| [Plummer sphere](plummer-sphere.md) | Finite-core SI dynamics, projected mass/density and conditional thin-lens observables |
| [Bianchi I](bianchi-i.md) | Anisotropic perfect-fluid GR expansion, signed axis rates, shear and directional photon redshift |

The numerical diagnostics are empirical, with their precise estimator scope in
each guide. A passing synthetic control does not certify global accuracy or an
observed-data likelihood. Hernquist's original-paper equation access remains
unverified; its declared-density Newtonian derivation is independently checked.

## Production build and installed consumers

The scientific build pins source commit
`775d9ab05fc4b14b128d7c895e31eabb90b282f5`, clean tracked state, GCC 14.2.0,
CMake 3.31.6, Rust 1.98.0 and portable CPU Release arithmetic. The build driver
used one job with `--profile release --native-tests core`; six model consumers
were separately compiled with strict floating-point flags against only the
fresh installed include directory and the complete production archive.

- All 12 new registered model/example contracts passed against the full archive.
- All six external installed consumers compiled and ran successfully.
- Complete installed `libirred_core.a` SHA256:
  `657e4f8f6144f626f32883e705b5640571ad342870024d92d52a3985fddf1fed`.
- Build-manifest SHA256:
  `34bce63651700233847fb1b1e14e81ebca3f2fcd5e3bda8c0e405062750202e2`.
- The registered inventory contains 148 tests; all 12 additions remain full-only.
  The full installed suite contains 13 consumers, with routine membership unchanged.

Worker evidence additionally includes analytic limits, independent formulations,
quadratures/refinement, hostile inputs and resource/ownership controls. Quintessence
passed 177 strict checks and sanitizers, DGP 259 strict checks, and decaying matter
240 strict checks and sanitizers. NFW, curved FLRW and Hernquist include independent
quadrature or high-precision fixtures. Worker component-install receipts are
separate from the complete production SDK gate recorded above.

The first installation dispatch check failed because its historical count was
seven; the corrected count and explicit full-only checks now pass. Curved tiny
positive distances that underflow to zero are refused. The original small-radius
Hernquist finite-difference failure is retained alongside the corrected reference
step; its original tolerance was preserved.

Six scientific workers were requested as GPT-6.1 Sol: four existing workers were
reused and two were newly launched. Launch reasoning effort was inherited rather
than explicitly selected; backend model identity was not independently attested.
Independent content review covered all six models and the complete integrated
43-file diff. Shared local compiler/compute concurrency was bounded by four jobs.
Historical source, scientific result, storage and runtime pins remain unchanged.

## Final eight-model production build

After the first six-model merge, Plummer and Bianchi I were independently
implemented and reviewed. Bianchi I passed 701 strict Release checks and
ASan/UBSan; Plummer passed strict Release, sanitizer, shell/line-of-sight and
high-precision fixture controls. Original source-access failures remain explicit
in their guides. These do not establish optical distances for Bianchi I or a
stellar distribution function for Plummer.

The final clean scientific source is
`980252384b7999333626b39b5ade6d7e670b802f`. Its complete native/Rust Release
production build passed with one job. All 16 new registered model/example
contracts passed against its complete archive. Eight external consumers then
compiled and ran against only the fresh installed headers and archive. The
registered inventory contains 152 tests, the 16 additions remain full-only, and
the full installed suite contains 15 consumers. Routine membership is unchanged.

- Final complete installed archive SHA256:
  `a9d490de79f83e81ccdc71005878c148f6110c954f1f82bc65df5b1cd6373cff`.
- Final build-manifest SHA256:
  `dd8687fbbf18b11279c8ca48f7cd33317b74fe8e45cf42c552cf1553ac4a19f4`.
- Models merged through [PR65](https://github.com/sprajs/irreducible/pull/65),
  [PR66](https://github.com/sprajs/irreducible/pull/66) and
  [PR67](https://github.com/sprajs/irreducible/pull/67), after independent complete
  diff review and seven passing required checks on each exact candidate.

## Actual compiled parameter sweeps

[Reproducible PR28](https://github.com/sprajs/reproducible/pull/28) adds an actual
compiled consumer of the original six-model production SDK. Each run contains
54 synthetic/conditional parameter cases: quick has 280 rows, broader 630 rows.
All four endpoint/reference/refusal controls passed in both corrected attempts.
The numerical qualification is explicitly **rejected or incomplete**: quick has
four failures; broader 24. These consist of scalar constraint diagnostics beyond
the frozen comparison allowance and refined NFW projected/lens-output conditioning refusals.
No allowance was weakened, refused row dropped or incomplete result relabeled.
The original mislabeled Hessians, corrected new attempts, dirty-source refusal
and successful model predictions are retained. Six-panel plots label the failed
qualification prominently. This is a parameter-response experiment, with no
observations, likelihood or inferred cosmological constraints. The historical
CLASS/runtime/scientific pins and this campaign's original SDK source pin775 remain
unchanged.

[Prospector PR39](https://github.com/sprajs/prospector/pull/39) records targeted
equation/source checks, exact PDF hashes where acquired, physical-scope limits
and source-access failures. These are bounded source reviews, with no invented
full-paper stage or theory-viability promotion. Its exact source metadata S3
pins and fresh readback proof are recorded in its dossier; raw papers were
excluded from publication while rights remain unresolved.

## Persistent evidence after transport recovery

Storage startup passed with the prepared verified AWS CLI2.37.11 and current
restricted identity. Subsequent upstream tunnel errors and actual create-only
PUT HTTP503 failures were preserved. After transport recovered, the unchanged
immutable engine selections were uploaded with exact-version object readback
and manifest-last completion. Both were restored into fresh directories and
compared with the originals: 119 files, 1295061 bytes, zero hash differences.
Credentials, profiles, proxies and TLS settings were unchanged. The historical
setup identity failure remains distinct from this runtime's passing startup and
later temporary transport failure.

The first six-model selection has 54 files, 437754 bytes, selection SHA256
`214412ea7d7b6c4bbc55d27b0664dfb5069f5bdd555018999c5a45c569f9dd1b`.
The Plummer/Bianchi followup has 65 files, 857307 bytes, selection SHA256
`db6cfad31d16dc0d1fd973ff246b91d94c156aa89814f2f9f99db7a5c2855265`.
They retain original worker/integration logs and receipts, failed controls,
corrected independent high-precision fixtures, exact build/installed identities,
source-access failures and publication diagnostics. Reconstructable binaries
and rights-unverified raw papers are excluded. The original first-stage version
and failed attempts/partial remote objects remain preserved; no eviction was used.

These are exact scientific evidence pins, distinct from a mutable catalog:

### First six models

- Manifest: `s3://research-data-436908790672-eu-west-2/reproducible/handoffs/irreducible/native-models-20261009t0200z/versions/0204af3da5360005c628ea8025f6353092e0adf93f29f9f32dd3f217cbbbd41e/manifest.json`
- SHA256: `c0d130ffd7f793036b02ecf1aeb22faff3cb73117a64624949d6fb877b32cbb7`
- VersionId: `gK0arALOqaEZWSB9TdYkSjOEKlYjkaKe`
- Format: `research-named-manifest/v1`.

### Plummer and Bianchi I

- Manifest: `s3://research-data-436908790672-eu-west-2/reproducible/handoffs/irreducible/native-model-followups-20261009t0220z/versions/2eb1873ad6abf0baea2520adf35d8da24fc3bc77cf112fa4bad2d8fc32bbda32/manifest.json`
- SHA256: `235d631696d68da99c370570ea57b34fd9438bc8aced5985e3937c88c9678684`
- VersionId: `JFqnoXUZA9.YYcFl2t2CvTxibrqvRiYR`
- Format: `research-named-manifest/v1`.

Restore with the named transport's exact manifest URI, SHA256 and VersionId
into a new directory. Use `list --collection` for completed uncataloged uploads.
The separately published Prospector source metadata and Reproducible compiled
sweep keep their own immutable pins and verified custody. The source-paper
copying limitations remain unchanged. Historical science/runtime/input pins
are not replaced by this new model campaign.

The code merges passed CI on their exact main commits: PR65
`95f9b1f3d933132b7fb82c031b38b37d6ee116a2` (run37872807287), PR66
`2545550781ee97d1eb3cfa6fc23c235a26ea269f` (run37873470252), and PR67
`36c942601a590ede6a4cc2fbd0cc8b01eb2e91c1` (run37874008720). Prospector
PR39 and Reproducible PR28 also passed exact-main CI. Bulk preservation does not
promote the numerical sweep's rejected/incomplete qualification.
