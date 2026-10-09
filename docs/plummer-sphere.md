# Finite-core Plummer sphere and conditional thin lens

`irred/plummer_sphere.hpp` supplies the native C++ model
`finite-core-Newtonian-Plummer-SI-supplied-G-and-lens-geometry/v1`.
Its independent positive parameters are total mass M in kg, scale radius b in
metres and explicitly supplied Newtonian coupling G in m³ kg⁻¹ s⁻². All radii
and lens distances are metres. This model has a finite central density and
finite total mass, unlike a cusped Hernquist sphere or an untruncated NFW halo.
There is no implicit solar-mass conversion, cosmological state or fixed G asset.

For s=sqrt(r²+b²), the declared spherical profile and Newtonian operators are

    rho = 3 M b²/(4π s⁵)            [kg/m³]
    M(<r) = M r³/s³                [kg]
    Phi = -G M/s, Phi(infinity)=0  [m²/s²]
    a_r = -G M r/s³                [m/s²]
    v_c² = G M r²/s³               [m²/s²]
    Phi'' = G M (b²-2r²)/s⁵        [s⁻²]
    Phi'/r = G M/s³                [s⁻²].

Acceleration is the signed outward radial component, negative for inward
attraction. The last two scalars are the radial and two identical tangential
eigenvalues of the Cartesian potential Hessian. The tidal acceleration tensor
is their negative. At the origin M(<0)=0, a_r=0, v_c²=0, rho=3M/(4πb³),
Phi=-GM/b, and the potential Hessian is the isotropic GM/b³ times the identity.
All central outputs are supported. A radial scalar acceleration of zero
corresponds to the unique zero Cartesian force there.

The profile is associated with [Plummer (1911), MNRAS 71, 460–470](https://doi.org/10.1093/mnras/71.5.460).
Targeted access to that original DOI returned HTTP 403; original-paper equation
locators have not been verified. This is an explicitly specified Plummer-density
model with independently derived Newtonian and projected operators, not a
completed original-paper reproduction audit. Integrate 4πr²rho for enclosed
mass, apply spherical Newtonian Phi'=GM(<r)/r² and the zero at infinity, then
differentiate to obtain force and curvatures. Their trace obeys
Phi''+2Phi'/r=4πG rho. All implementation and comparison code is original;
no external implementation or paper assets are copied.

At cylindrical projected radius R, define s=sqrt(R²+b²). Direct line-of-sight
integration and cylindrical mass integration give

    Σ = M b²/[π(R²+b²)²]       [kg/m²]
    M_2D(<R) = M R²/(R²+b²)    [kg]
    meanΣ = M/[π(R²+b²)]        [kg/m²]
    ΔΣ = M R²/[π(R²+b²)²]       [kg/m²].

The mean is the continuous central limit at R=0, where Σ=meanΣ=M/(πb²),
M_2D=0 and ΔΣ=0. Projected cylindrical mass includes material outside a sphere
of radius R; it is a different quantity from the three-dimensional enclosed
mass. In particular, at R=b the cylindrical mass is M/2 and the spherical
mass is M/(2sqrt(2)).

`plummer_thin_lens` requires positive supplied critical surface density C in
kg/m² and angular-diameter lens distance D_l in metres. The caller owns the
source/lens geometry, consistency of C with supplied G, and validity of the
thin-lens approximation. No cosmological distance or redshift law is inferred.
Its reduced deflection and angular potential are

    κ = Σ/C, γ_t = ΔΣ/C
    α = meanκ R/D_l                     [radians]
    ψ = M/[2π C D_l²] ln(1+R²/b²)        [radians²], ψ(0)=0.

Angles obey R=D_l |θ|, and β=θ-α. The angular potential Hessian eigenvalues
are meanκ (b²-R²)/(b²+R²) radially and meanκ tangentially. These are potential
curvatures; the lens Jacobian eigenvalues are one minus them. The gradient and
Hessian use the same angular potential. At the origin α=ψ=γ_t=0 and the
angular Hessian is isotropic. At R=b its radial eigenvalue is exactly zero.
This conditional operator predicts no image topology, arrival-time distance,
observed stellar kinematics, named lens, shear catalogue, matter spectrum,
mass fit or posterior. Synthetic SI parameters do not establish a physical
object or qualify observational assumptions.

Inputs must be finite, M,b,G,C,D_l and relative tolerance positive, and radius
nonnegative. Default relative tolerance is 1e-12. The arithmetic contract
requires round-to-nearest, at least 64 long-double mantissa bits and extended
exponent range of at least 16384. All finite binary64 input radii are considered;
overflow and nonzero outputs not representable as normal binary64 are refused.
The implementation uses wide hypot and normalized radius/core ratios. The
small angular potential uses log1p without cancelling a central constant, and
ΔΣ uses its positive product rather than meanΣ-Σ. The lens radial curvature
factors (b-R)(b+R), retaining its exact zero and both adjacent floating-point
radii. The Newtonian radial curvature has a cancellation diagnostic near
r=b/sqrt(2), where a tight relative allocation refuses only that scalar while
retaining admitted scalar outputs.

Each scalar has its own status and absolute arithmetic diagnostic. Aggregate
status reports any refused output while preserving other admitted outputs.
Exact known central zeros and the lens radial zero remain exact; they are not
positive underflows. Diagnostics include wide arithmetic, libm and binary64
projection allowances. They are empirical, not certified enclosures or
uncertainty in supplied parameters. Physical parameter uncertainty and the
applicability of classical gravity are caller-owned.

Permanent scientific controls use direct density quadratures for spherical
mass, Newtonian shell potential, line-of-sight Σ and exterior-cylinder meanΣ.
Finite spherical and line-of-sight intervals use tangent substitutions to
resolve both the core and the tail. The cylindrical exterior uses r=R cosh(u)
and stops at u=32. For the tested R/b≥1e-6, its omitted meanΣ is bounded by
3/[8π(R cosh(32))⁴], below 2e-31 relative. Mesh doubling differences must be
≤5% of the mass, potential, Σ and meanΣ comparison budgets (2e-9,2e-9,2e-9,
2e-8 respectively). Independent 90- and 110-digit density-shell, line-of-sight
and cylindrical integrals agreed to 1e-40 relative, supplying permanent literal
fixtures at r/b=0.1,1,10,1e-12,1e12 for all seventeen outputs. Their potential
uses a Green-function integral of the independently line-of-sight-checked
surface density, rather than the production logarithmic expression. Fixtures
are compared to returned diagnostics and binary64 literal conversion uncertainty.
The initial 70-/90-digit oracle failed the same 1e-40 refinement allocation
for tiny shear and the extreme-radius shell integral; increased precision
resolved it without weakening that allocation. Original failed controls are
retained in the comparison evidence, and no oracle dependency is needed by
permanent C++ tests.

Controls also test Poisson against inherited absolute curvature diagnostics,
SI scaling, shell force, derivative identities at numerically resolved radii,
central and Kepler limits, exact central zeros, adjacent lens-curvature zeros,
Newtonian near-zero refusal, invalid parameters, arithmetic mode, overflow,
underflow and tighter-budget failures. Binary64 central potential differences
and far-tail enclosed-mass differences are unresolved at the smallest/largest
radii; their derivative checks use intermediate radii, while independent
integrals, analytic limits and high-precision fixtures cover the extremes.

Strict Release-mode scientific tests, AddressSanitizer/UndefinedBehaviorSanitizer
tests and a fresh installed component-archive SDK consumer passed. Full product
registration and integration checks remain the integration owner's gate.

A standalone SDK consumer can use:

```cpp
#include <irred/plummer_sphere.hpp>
irred::gravity::PlummerSphere sphere{2e41, 3e19, 6.67430e-11};
auto lens = irred::gravity::plummer_thin_lens(sphere, 3e19, 20, 3e24);
// Check each required scalar.status before reading values.
```
