# Synthetic SIS point-source lens

`irred/sis_thin_lens.hpp` provides one named synthetic thin-lens consumer: an
axial point source, a singular isothermal sphere (SIS) template, an explicit
nonnegative mass sheet, and a circular Gaussian PSF integrated over rectangular
pixels. The lens is centered at (0,0); the source is at (beta,0). Angles are
radians, sigma_v is km/s, relative observer arrival times are seconds, and flux
uses the caller's fixed linear unit. This is a synthetic prediction, not an
arbitrary astrophysical lens model, observed stellar-kinematics interpretation,
lens-mass inference or parameter posterior.

Preparation retains an owned copy of the existing early/late massless-radiation
source and computes only its two comoving distances once. It requires that
provider's compatible flat pressureless-matter, positive massless-radiation and
nonnegative Lambda state. Its supplied-drag provenance travels with the source,
but no ruler is evaluated or used; a zero sound-callback quota succeeds.
Thermal massive relics and the radiation-free growth identity are distinct.

For chi_l and chi_s from that provider, angular distances are
Dd=chi_l/(1+zl), Ds=chi_s/(1+zs), Dds=(chi_s-chi_l)/(1+zs), in physical Mpc.
The base SIS angle is thetaE=4 pi (sigma_v/c)^2 Dds/Ds. Its explicit transformed
potential is psi(theta)=lambda thetaE |theta|+(1-lambda)|theta|^2/2, with beta
already the transformed source coordinate. Lambda=1 is the pure SIS;
when lambda differs from one, sigma_v denotes the base SIS template parameter.

The positive image is theta+=beta/lambda+thetaE. A negative image exists only
for beta<lambda thetaE and is theta-=beta/lambda-thetaE. Signed magnifications
are mu±=(1±lambda thetaE/beta)/lambda^2, parities are their signs, and image
fluxes are source_flux times |mu|. The positive image arrives first with relative
delay zero; the second arrives after

    Delta_t=(1+zl) Dd Ds/(c Dds) 2 thetaE beta.

The authoritative native SI/IAU Mpc and c definitions convert this to seconds.
The SIS, Jacobian and Fermat equations are given in
[Narayan and Bartelmann (1996), sections 3.1–3.3](https://arxiv.org/abs/astro-ph/9606001).
The explicit sheet implements the transformation discussed in
[Schneider and Sluse (2013)](https://arxiv.org/abs/1306.0901).
All code and reference calculations here are original; no external code or
assets are copied.

The bounded slice admits 0<zl<zs<=20, 0<sigma_v<=0.01c,
1e-6<=lambda<=1, 0<beta<=0.1, and 0<=PSF sigma<=0.1 radians. Beta and its
distance from lambda thetaE must exceed a 1e-6 relative topology margin and
propagated distance uncertainty. Unresolved near-coincident redshifts and
caustics are refused. The largest image angle, including its propagated
diagnostic, must be at most 0.1 radians; beta/lambda must remain in this
tangent-plane slice. Rectangles have strictly increasing axes with finite
bounds of absolute value at most 0.2 radians.

For positive PSF width, each image contributes its flux times the product of
normal-CDF differences along the two axes. A tail-aware erf/erfc difference
avoids subtracting two CDFs near one. At width zero, the explicitly supported
point limit assigns flux to half-open rectangles [xlo,xhi) × [ylo,yhi).
An image position whose diagnostic interval intersects an x boundary is refused;
y=0 is exact and follows the half-open convention. Positive Gaussian flux that
underflows is refused rather than replaced by zero.
This includes wide joint-probability, weighted-term and diagnostic products:
positive marginals cannot silently become a zero pixel or a dropped image.

Images can be projected with separate position, magnification/parity, flux and
delay masks. Pixel failures leave those retained predictions available.
Each requested scalar has an optional value, a cause and an empirical absolute
error estimate. Distance subtraction and ratios propagate into thetaE, positions,
magnifications and delays. Pixel diagnostics additionally propagate position
uncertainty using endpoint-density derivative bounds over the retained interval.
Platform erf/erfc arithmetic diagnostics are estimates, not certified universal
bounds. The default projection allowance is 1e-8 relative; narrower budgets can
refuse individual outputs. No jitter, topology repair or dropped image is used.

Preparation reports actual geometry quadrature callbacks. Pixel evaluation
reports attempted image-pixel terms and actual CDF endpoint evaluations
(four per positive-width term, none at width zero). Work limits retain preceding
successful rows and count failed attempted terms. Default limits are 4096 pixels,
8192 terms and 16 MiB concurrent owned payload; hard limits are 65536 pixels,
65536 origin bytes and 1 GiB payload. Admission precedes allocation. Payload
accounting includes conservative temporary provider/source ownership, excluding
borrowed inputs, stack recursion, allocator metadata and RSS. Copying owns the
source; moving invalidates the old owner, and self move preserves it. Supported
arithmetic requires round-to-nearest and the native wide host profile.

Permanent controls use an elementary radiation-only distance, an original
potential/Jacobian derivation, and independent fixed-panel Gaussian-density GL8
integration with 64/128-panel refinement. Named angle/magnification/image-flux
controls allocate 2e-12 relative, delays 2e-11, and pixels 2e-12 absolute in
source-flux units. Reference refinement must use at most 5% of the pixel budget;
distances inherit 1e-9 Mpc plus 2e-11 relative. Finite partitions test total flux,
and synthetic image, source-flux and delay recovery is explicitly conditional.

The degeneracy controls transform beta to lambda beta0 and source_flux to
lambda^2 source_flux0: images and pixel fluxes stay fixed while delays scale by
lambda. With compatible fixed present-day fractions, changing H0 by k keeps
angles fixed and divides delays by k. Combining H0->lambda H0 with the sheet
transformation also preserves the delay. Distances alone therefore do not
identify the lens mass or break this synthetic degeneracy.
