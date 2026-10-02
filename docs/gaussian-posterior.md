# Proper Gaussian parameter posterior

The native `irred/gaussian_posterior.hpp` calculation conditions an explicit
proper Gaussian parameter prior on a fixed linear observation model. It is
available through the C++ SDK, with no CLI/C ABI operation. For explicitly
independent beta~N(m,S) and epsilon~N(0,C), r=X beta+epsilon gives

    P = S^-1 + X^T C^-1 X,
    V = P^-1,
    mu = m + V X^T C^-1 (r-Xm).

The returned covariance V and conditional mean mu define a normalized Gaussian
parameter density in the caller's declared product of parameter-coordinate
Lebesgue measures. `log_density` includes the determinant and p log(2 pi).
This is a parameter posterior under the supplied proper prior. It is distinct
from the relative design profile, the conditional estimator sampling law and
the observed-residual density obtained by marginalizing a calibration prior.
It does not fit a nonlinear cosmological model or reproduce a released posterior.

The source Gaussian owns the supplied observation covariance and its retained
factor. `ParameterPrior` specifies ordered unique parameter IDs, units, shared
nuisance IDs (a unique subset), mean, full SPD covariance S, prior/design IDs,
residual unit, parameter measure and an explicit noise/prior independence
declaration. Unit/measure labels are provenance, not physical conversions.
At least two parameters are supported. Observation/design axes must match
exactly. Prior-bearing or mean-shifted source Gaussians are rejected, preventing
an accidental second application of the same prior. A null or rank-deficient
design is permitted because S is proper; no jitter or implicit regularization
is introduced. Zero-width box support and improper flat priors are unsupported.

Preparation whitens design columns with the retained C factor and identity
columns with a temporary S factor. Wide Gram sums form P, with diagnostics for
whitening, accumulation and final binary64 storage. The existing precision-input
Gaussian preparation produces V and its retained normalized-density factor.
The owner retains the design, whitened design, precision and prior. All candidate
allocations and checks finish before it consumes the source; failed preparation
leaves that source usable. Moves invalidate the old owner and self move preserves
it. Evaluation uses retained factors and matrix-vector operations without
refactorization, copying the source covariance or scanning it again.

Finite normal binary64 inputs are required; zero is allowed. Arithmetic requires
round-to-nearest and long double with at least 64 mantissa bits/exponent range
16384. Conditioning, underflow/overflow or unrepresentable positive diagnostics
cause refusal. Element/payload limits apply before allocations. The preparation
work admission uses the conservative dimension count
32*(n*n*p+n*p*p+8*p*p*p), with a default 100-million ceiling; this bounds this
implementation's declared work scope and is not a measured wall-time guarantee.
The payload estimates cover owned/transient matrix, vector and metadata storage,
excluding borrowed input, allocator overhead, stack and RSS.

The posterior covariance's relative numerical error estimate propagates the
precision diagnostic through its inverse norm. Mean outputs carry numerical
absolute-error estimates separately from physical posterior variance. Density
admission propagates the mean/subtraction estimates through P and adds covariance
quadratic/log-determinant sensitivity. Default maximum scaled forward sensitivity
is 1e-10. These are empirical arithmetic/conditioning estimates, not rigorous
continuous or universal error certificates. Numerical noise is not added to S
or C. A tighter unattainable request fails without changing its allocation.

Permanent owner controls use exact rational two-dimensional cofactor references
derived before implementation, plus a separately composed Bayes density identity,
null/rank-deficient designs, parameter-unit Jacobian, exact axis order, proper
prior rejection, quotas and owner moves. The fixed named mean/covariance and
density comparison allocation is 2e-12*(1+absolute reference). Shared factor and
Gaussian kernels in the Bayes comparison are ancestry; they are not independent
evidence. A separately derived covariance-conditioning/cofactor reference and
direct joint prior-times-likelihood quadrature now test normalization, moments,
rank-deficient/null and underdetermined designs, axis permutations and unit
Jacobians. Named mean/covariance controls use 2e-12*(1+absolute reference),
log density uses 2e-11*(1+absolute reference); reference refinement consumes at
most 5% of those allocations. Two initial reference failures remain preserved.
These controls do not qualify arbitrary inputs or physical priors. No observational calibration uncertainty,
confidence coverage or systematic-error claim follows from these controls.
