# Bounded detector counts and censoring

The native C++ `irred/detector.hpp` and `irred/detector_selection.hpp` APIs
compose supplied expected transmitted photons with a declared linear detector
and the same measurement law for simulation and threshold likelihoods. They
have no CLI/C ABI route. The supplied photon arrival process must be explicitly
Poisson; a deterministic realized photon count or arbitrary thermal photon
statistics cannot be substituted for its expectation. Quantum efficiency is a
fixed wavelength-independent electron conversion probability, applied after
optical transmission. Source, distance and calibration uncertainty stay separate.

With N expected transmitted photons, efficiency eta, supplied background B in
electrons, dark rate d in electrons/s, observer exposure t, read RMS sigma in
electrons, gain g in electrons/ADU and bias b in ADU, the model is

    lambda = eta N + B + d t,
    K ~ Poisson(lambda), R ~ Normal(0,sigma^2), independently,
    Y = (K+R)/g + b,
    E(Y) = lambda/g+b, Var(Y) = (lambda+sigma^2)/g^2.

This initial slice supports lambda and sigma from 0 through 64 and finite normal
binary64 inputs (zero allowed). It excludes saturation, digitization, pileup,
cosmic rays, uncertain background subtraction and measured-camera qualification.
The linear Poisson/read model follows the equations in [EMVA 1288 Release 4.0
Linear](https://www.emva.org/wp-content/uploads/EMVA1288Linear_4.0Release.pdf);
this calculation is not an EMVA instrument characterization or compliance claim.
The returned Poisson mean is the stored binary64 model coordinate used by both
simulation and likelihood. Unrepresentable positive means or variances fail.
Bias cancellation is admitted only within the fixed 2e-12 relative arithmetic
allocation or for an FMA-verified exact zero; unresolved positive measurements
cannot become fabricated zeros.

The caller owns a uint64 seed and each (uint64 stream,uint64 sample) address.
The counter is (sample low, sample high, stream low, stream high); the seed is
the key for Philox4x32 with ten rounds. An original transcription of the
[published Random123 algorithm](https://www.thesalmons.org/john/random123/papers/random123sc11.pdf)
uses its mathematical constants; no library implementation is copied. Three
factual known-answer vectors are pinned to Random123 commit
9545ff6413f258be2f04c1d319d99aaef7521150. Uniforms are (word+0.5)/2^32;
word0 drives inverse Poisson CDF, words1/2 drive Box–Muller Gaussian read noise.
All four words and the address are returned. No mutable generator state exists.
Duplicate addresses in a batch fail; intentional replay across calls is allowed.
Sharding and order leave each address unchanged. Distinct counters and empirical
cross-moment tests are not a proof of statistical independence. Uniform
quantization and finite normal tails are explicit discretization assumptions.
Integer words replay exactly; transcendental outputs require the same supported
round-to-nearest wide arithmetic and libm identity, not cross-platform bit equality.

Detection is Y>=threshold. With positive sigma, a detected ADU value has the
normalized continuous mixture density sum_k Poisson(k;lambda) Normal_ADU(Y;k/g+b,
sigma/g). At sigma=0 it uses an integer electron count and a discrete mass.
A nondetection reports only its threshold and log P(Y<threshold). The default
joint record likelihood retains both detections and nondetections. An explicitly
selected-only detected sample divides by P(detected); it cannot accept a
nondetection. Neither measure silently discards undetected rows. Payload and
measure mismatches fail. Structural zeros have an explicit zero_probability
flag and no finite log value; unresolved positive probabilities fail instead.

Likelihood sums both Gaussian/Poisson tails directly rather than subtracting a
rare detection from one. At most 256 Poisson terms are admitted. A geometric
omitted-tail estimate and floating-point diagnostics gate each result before
normalization/logarithms; default probability allocation is 5e-13 absolute plus
2e-11 relative, with a separate scaled log sensitivity gate. These numerical
diagnostics are not rigorous universal certificates or observational noise.
Configured term exhaustion preserves failed rows; no Gaussian approximation,
jitter or tail replacement is used. Batch row count is at most 65536 with explicit
payload admission. Bounds cover owner/vector payload and address-sort scratch,
excluding borrowed inputs, stack, allocator overhead and RSS. Successful other
rows survive invalid detector inputs in simulation; likelihood shares one fixed
detector source across its batch.

Permanent controls include analytic moments, Poisson masses/censoring, the
zero-Poisson Gaussian limit and its ADU Jacobian, structural support, generator
known answers, replay/sharding, resource and arithmetic failures, and fixed seeded
ensemble mean/variance/CDF/cross-moment checks. A small synthetic censored Poisson
recovery retains nondetections and compares the native joint objective with an
independently solved analytic score; detected-only raw averages show the expected
selection bias. This is a conditional synthetic control, not posterior inference,
coverage, astrophysical population recovery or qualification of measured data.
