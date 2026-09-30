# Scientific contracts

Irreducible makes changing an assumption possible and its consequences inspectable. Skepticism applies to this implementation as much as to a published package.

## Before implementation

State the equation, model identity, input semantics and numerical domain. Identify the consumer: an observable, log-likelihood or final estimand. Allocate its error budget and explain how dependency errors propagate. Design the independent comparison before optimizing.

Each fixture needs reconstructible ancestry: source or derivation, units/conventions, input representation, reference uncertainty and tolerance rationale. Agreement is incomplete evidence if implementations share an algorithm, constants or system math library.

## Observations and assumptions

Distinguish observations, fitted summaries, empirical training/calibration assets, physical assumptions and synthetic truth. Preserve original values and ordering; record transformations. Track repeated sources, overlap, selection and shared calibrators. Unknown covariance stays unknown. Never silently promote a posterior summary to an independent datum.

Matching dimensions do not justify a physical-role or frame conversion. A parsec definition, for example, does not specify whether a distance is luminosity, angular-diameter or comoving distance.

## Qualification and discrepancies

A test provides evidence for named cases and domains. Qualification additionally requires scope, consumer budget, source/build/model/data identity and relevant comparisons. Keep execution, numerical, inference and interpretation separate.

On a discrepancy, preserve the input, old output and failure locally. Reduce it to a useful regression, identify the cause, and explain the correction. Do not add jitter, discard rows, change conventions or relax tolerances without scientific justification and explicit review. Withhold the affected claims and continue independent work.

Keep the regression and concise derivation in public tests. Large oracle dumps, external-tool environments and development receipts stay local. Native tests should remove routine dependence on external comparison software while preserving the meaning and ancestry of expected values.

## Reproduction

Use **exact** when required original inputs and procedures are reproduced, **approximate** for stated substitutions, **conditional** for unresolved conditions, or **blocked** when a necessary asset/capability is missing. Explain the label; never silently substitute author inputs.

For posterior claims, check support, error, diagnostics and sensitivity for the actual estimand. Weak reweighting or an unresolved numerical screen is not a corrected cosmological result. Historical results from another research tree are not automatically results of Irreducible.

## Data and probability semantics

A quantity's role, frame, observer, calibration and measure matter beyond its dimensions. Signed measured flux is valid even where a logarithmic magnitude is undefined. Preserve missing/nonfinite masks; zero is not a replacement for missing data. Source and covariance axes must stay ordered through selection. A covariance principal block represents selected marginal uncertainty; slicing a precision matrix does not. Unknown cross-survey dependence is not a zero block.

Profile scores, relative targets, proper-prior marginal densities and evidence are different statistical objects. A free-offset SN shape profile carries no absolute calibration or H₀ information. An improper flat nuisance measure does not provide normalized evidence. Include coordinate-transform Jacobians where the measure requires them; MAP changes with coordinates. Diagnose rank, identifiability, missing support, consequential tails and nonexistent moments before reporting an estimand. Mathematical zero density is a support result, never a replacement for numerical failure.

Future spectral and detector operators must declare Fλ/Fν and redshift/time Jacobians, source/rest versus foreground/observer attenuation, and photon or energy weighting. Detector area and photon-energy conversion enter exactly once. Preserve absolute flux/conservation checks before introducing a free amplitude. Censoring, retained truncation and population count/intensity likelihoods need their own normalization; selection is not a residual weight.

Decoded format size needs its own allocation bound; compressed file size is insufficient. A supported FITS/HDF/image profile must state scaling, nulls, units, time/frame conventions and decoding scope. Operator uncertainty is separate from fixed linear propagation C′ = W C Wᵀ. These are contracts for future consumers, not claims that all formats/operators are implemented.
