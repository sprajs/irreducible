# Finite-mass Hernquist sphere with supplied Newtonian coupling

`irred/hernquist_sphere.hpp` supplies the native C++ model
`finite-mass-Newtonian-Hernquist-supplied-G-SI/v1`. Its independent positive
parameters are total mass M in kg, scale radius a in metres, and explicitly
supplied Newtonian gravitational coupling G in m³ kg⁻¹ s⁻². Radius is in metres.
The model assigns no implicit solar-mass conversion or fixed G asset; uncertainty
in supplied physical parameters is separate from numerical diagnostics.

The declared spherical density and its Newtonian operators are

    rho(r) = M a / [2π r (r+a)³]
    M(<r) = M r²/(r+a)²
    Phi(r) = -G M/(r+a), with Phi(infinity)=0
    inward radial acceleration = -G M/(r+a)²
    circular speed squared = G M r/(r+a)².

Mass is kg, density kg/m³, potential and squared speed m²/s², and acceleration
m/s². The acceleration is the signed outward radial component; its negative
sign means inward attraction. The eigenvalues of the Cartesian potential
Hessian are -2GM/(r+a)³ in the radial direction and GM/[r(r+a)²] in each of the
two tangential directions, in s⁻². They are potential curvatures; the tidal
acceleration tensor is their negative.

The density is the finite-mass spherical profile associated with
[Hernquist (1990), ApJ 356, 359–364](https://doi.org/10.1086/168845).
The operators here are independently derived from the explicitly declared
profile: integrate 4πr²rho for enclosed mass, use spherical Newtonian gravity
Phi'=GM(<r)/r², fix the zero of potential at infinity, and differentiate that
potential for acceleration and curvature. This gives
Phi''+2Phi'/r=4πG rho. The central cusp is rho∝1/r, the outer density is ∝r⁻⁴,
and the enclosed mass tends to the supplied finite total M. This is a distinct
physical identity from an untruncated NFW halo. The original-paper PDF endpoint timed out and its alternate full-page endpoint
returned an empty response during targeted source review. Original-paper equation
locators remain unverified; the bibliography alone does not establish a
completed original-source audit. All implementation and reference code is
original; no third-party code or data are copied.

At the exact origin the operator retains M(<0)=0, Phi(0)=-GM/a and the continuous
limit of circular speed squared, zero. Density diverges; the Cartesian force
has no unique direction, and a Cartesian potential Hessian does not exist.
Those four scalars are returned with `singular` status. No fictitious central
force vector or finite tidal tensor is supplied. The finite central potential
does not remove the density cusp.

Inputs must be finite with M,a,G and relative tolerance positive and radius
nonnegative. The arithmetic contract requires round-to-nearest, at least 64
long-double mantissa bits and exponent range at least that of the usual
extended format. All finite binary64 input radii are considered, but a requested
nonzero output must be representable as a normal binary64 number. Overflow and
positive underflow are refused. Default relative tolerance is 1e-12. The
implementation uses positive sums, scaled radius ratios and long-double
products; no difference of nearly equal terms or iterative solver is used.
Each scalar carries its own arithmetic diagnostic and admission status.
Aggregate status reports any refused scalar while retaining other admitted
outputs. Diagnostics include wide arithmetic and binary64 projection, and are
not certified bounds or physical-parameter uncertainties.

The model supplies classical Newtonian spherical dynamics only. It has no
relativistic interior, orbit distribution, velocity anisotropy, stellar-population
or kinematic data interpretation, projected lens, cosmological background,
image formation, mass fit or parameter posterior. Its finite-mass potential is
a useful explicit object-dynamics prerequisite; a supplied profile does not
establish that an observed object follows it.

The scientific regression uses direct density integrals for mass and the
interior/exterior Newtonian shell potential. In dimensionless M=a=G=1 controls, the interior radial change of
variable r=expm1(u) resolves the central cusp; the exterior integral ends at
u=50, retaining an omitted dimensionless potential tail bounded by exp(-100).
The tested dimensionless radii range from 1e-6 to 1e6, where that tail is below
1e-37 relative. Mesh doubling must change mass by at most 1e-10 and potential
by at most 1e-9, 5% of their 2e-9 and 2e-8 comparison allocations. Independent 70- and 90-digit mpmath density-shell integrals agreed to 1e-50
relative and provide permanent literals at r/a=0.1,1,10,1e-12,1e12, checked
against the reported arithmetic diagnostics and literal conversion uncertainty.
The high-precision force and curvatures use shell gravity and Poisson, not
production rational formulas. Rational identities at r=a, central and Kepler
limits, independent shell gravity and
Poisson identities, finite-difference derivatives, SI scaling, origin partial
results, invalid inputs, tighter-budget failures, overflow, underflow and
arithmetic-mode refusals challenge the model. Tiny-radius finite-difference
steps use 0.01r below r/a=0.01 to resolve binary64 potential differences while
keeping the original derivative budgets; elsewhere they use 1e-5r. Strict
Release-mode scientific tests, AddressSanitizer/UndefinedBehaviorSanitizer tests
and a fresh installed component-archive SDK consumer passed. Full product
registration and integration checks remain the integration owner's gate.

A standalone SDK consumer can use:

```cpp
#include <irred/hernquist_sphere.hpp>
// Synthetic parameters; G is explicitly supplied, without propagated uncertainty.
auto value = irred::gravity::evaluate_hernquist_sphere(
    {2e41, 3e19, 6.67430e-11}, 3e19);
// Check each required scalar.status before consuming its value.
```
