# New native model campaign, 2026-10-09

This campaign adds six distinct compiled physical models. It advances new model
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
