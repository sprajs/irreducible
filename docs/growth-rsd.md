# Conditional supplied-amplitude growth and synthetic RSD

The native `irred/growth_amplitude.hpp` consumer composes the existing
[GR growth](gr-growth.md) owner with an explicitly supplied fixed sigma8
amplitude. The separate `irred/rsd.hpp` consumer scores ordered synthetic
f sigma8 rows with a retained full Gaussian covariance. Both are standalone
C++20 APIs, with no CLI or C ABI route.

The compatible background is radiation-free, flat, pressureless LCDM with
1e-6 <= Omega_m <= 1 and 1e-8 <= a <= 1. D/a tends to 1 at early times;
D(1) is not silently set to 1. The source supplies sigma8_ref and a_ref,
where a_ref has the same scale-factor domain. Its compiled convention declares
the dimensionless linear total pressureless-matter RMS in a comoving top-hat
of radius 8 h^-1 Mpc. Nonempty amplitude identity/provenance and the explicit
fixed-supplied uncertainty treatment are required.

The shared amplitude equation owner computes

    sigma8(a) = sigma8_ref D(a) / D(a_ref),
    f_sigma8(a) = f(a) sigma8(a).

This propagates a supplied normalization; it does not calculate sigma8 from
a primordial spectrum or transfer function. No H0 enters at fixed Omega_m
and supplied amplitude identity. Changing the physical meaning of the
8 h^-1 Mpc smoothing scale is outside this conditional calculation.

`amplitude_sigma8` and `amplitude_f_sigma8` request outputs independently.
Each row preserves its original scale factor, optional values, cause and
absolute numerical diagnostic. The reference is evaluated once per batch.
Exact a == a_ref rows cancel the identical D numerator and denominator:
sigma8 equals the supplied amplitude exactly, while f_sigma8 retains the
reference f diagnostic. Other rows share the same reference-D diagnostic.
Product/ratio endpoint propagation combines the empirical D/f diagnostics
with arithmetic and storage effects, without treating them as independent
random errors. The reference denominator must remain resolved positive.

The inherited growth relative allowance remains 1e-10; composed values use
1e-12 + 5e-10 |value| by default. Diagnostics are empirical numerical estimates,
not certified universal bounds or posterior variances. Positive outputs and
positive diagnostics require normal binary64 representation; overflow,
underflow or unresolved values are refused. A finite exact zero amplitude
is a declared unperturbed control. It produces exact zeros without quadrature
only after source, query, background and arithmetic validation.

`rsd::prepare_density` admits synthetic-control rows only. Source values,
ordered unique row IDs, event IDs, exact covariance axes, amplitude convention,
full covariance units and table/order/calibration/dependence identities are
explicit. Event IDs may repeat; neither their presence nor full within-source
covariance declares independence from another probe. Observed synthetic Gaussian
values may be negative. All rows remain present, and covariance must be exactly
symmetric SPD. There is no jitter or row removal.

The move-only prepared owner retains source records and one Gaussian factor.
The disposable input matrix is released after preparation; the Gaussian's
required covariance/factor storage remains owned and charged. Model batches
reuse that state without covariance readback, copying or refactorization.
For r = y - mu, the density in the product of dimensionless f_sigma8 coordinates is

    log p = -1/2 [r^T C^-1 r + log det C + n log(2 pi)].

Output bits request normalized density, predictions and residuals separately.
Requested arrays retain ordered optional row values and failure causes.
A failed row withholds density while preserving other useful requested rows.
Successful predictions and residuals survive a density projection refusal.
Batch finite status reports admission, not success of every model.

Projection uses the existing shared density helper. With eps containing
prediction, subtraction and storage diagnostics, it estimates

    ||C^-1 r||inf ||eps||1
      + 1/2 estimated||C^-1||inf ||eps||inf ||eps||1.

The inverse norm estimate uses the retained factor condition estimate and
covariance infinity norm. The default absolute log-density projection allowance
is 1e-8, separately enforced from prediction admission and factor/solve
sensitivity. Numerical diagnostics are never added to observation covariance.
Arithmetic must match the retained Gaussian and satisfy the supported wide
host profile and round-to-nearest environment.

Preparation bounds rows, matrix elements, metadata strings and simultaneous
payload before allocation. Evaluation bounds model/row counts, copied source
strings, requested arrays and serial growth/residual/solve scratch together.
Reference plus query callbacks, including failed callbacks, consume the same
original per-model cap and the consumer whole-batch cap. Existing growth
per-point/depth/provider limits apply. Payload excludes borrowed inputs,
preexisting retained owners, allocator overhead, stack frames and RSS.

Permanent owner controls include exact Einstein-de Sitter limits, original
3x3 covariance cofactors, source/axis failures, amplitude scaling, endpoint
cancellation, omissions, quotas, moves and projection refusal. Separately
authored high-precision quadrature and original GR ODE references qualify named
non-EdS f_sigma8 cases with relative 1e-8 comparison allocation; reference
refinement occupies at most 5%. Density components retain absolute 1e-8
comparison allocation. Shared physical equations remain shared ancestry.

These are bounded numerical synthetic controls. Radiation or massive-relic
perturbations, AP/window/estimator mapping, released fitted-compression validity,
transfer prediction, measured H0, parameter posterior and observational
qualification require separate contracts and evidence.
