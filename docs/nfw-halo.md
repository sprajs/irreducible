# Spherical NFW halo and conditional thin lens

`irred/nfw_halo.hpp` supplies the native C++ model
`spherical-untruncated-NFW-supplied-critical-density/v1`. Its independent positive
parameters are rho_s in Msun/Mpc³ and r_s in physical Mpc. All radii are physical
Mpc. There is no implicit redshift, concentration, critical versus mean cosmic
density convention, truncation or cosmological geometry.

For x=r/r_s, the density is rho_s/[x(1+x)²] and the enclosed three-dimensional
mass is 4π rho_s r_s³ [ln(1+x)-x/(1+x)], in Msun. The mass at zero is exactly
zero. The untruncated halo has logarithmically divergent total mass; the API
returns only finite enclosed mass, never a finite total halo mass.

The thin-lens projection at x=R/r_s is

    Σ = 2 rho_s r_s f(x)
    meanΣ = 4 rho_s r_s g(x)/x²
    ΔΣ = meanΣ - Σ
    A(x) = acosh(1/x)/sqrt(1-x²)       (x<1)
           1                         (x=1)
           acos(1/x)/sqrt(x²-1)       (x>1)
    f(x) = (1-A(x))/(x²-1), f(1)=1/3
    g(x) = ln(x/2)+A(x), g(1)=1-ln(2).

These surface densities have units Msun/Mpc². The projected cylindrical mass
can be reconstructed as πR² meanΣ. It includes matter outside the sphere of
radius R and must not be confused with the three-dimensional enclosed mass.

`nfw_lens` requires independently supplied positive Σcrit in Msun/Mpc² and
angular-diameter lens distance D_l in physical Mpc. The caller owns their
physical thin-lens geometry, gravity and source/lens distance conventions.
The operator returns κ=Σ/Σcrit, tangential γ=ΔΣ/Σcrit and the reduced deflection
α=meanκ R/D_l in radians. It also returns the radial and tangential eigenvalues
of the angular lens-potential Hessian, 2κ-meanκ and meanκ. These are **potential**
curvatures; the lens Jacobian eigenvalues are one minus them. Deflection points
radially outward in the lens equation β=θ-α. This conditional mass projection
predicts no image topology, time delays, observed stellar kinematics, named
lens, shear catalogue, matter power spectrum or parameter posterior.

The original source equations are [Navarro, Frenk and White (1996/1997),
Eq. 1](https://arxiv.org/abs/astro-ph/9611107), [Bartelmann (1996), section 2](https://arxiv.org/abs/astro-ph/9602053),
and [Wright and Brainerd (2000), Eqs. 10–14](https://arxiv.org/abs/astro-ph/9908213).
Wright–Brainerd Eq. 11 supplies Σ, Eq. 13 meanΣ and Eq. 12 their shear relation.
Their named g_< and g_> shear functions are different from the projected-mass
function g used here. All implementation and comparison code is original;
no third-party implementation or licensed data are copied.

The admitted nonzero radius domain is 1e-8≤x≤1e8. Projection at R=0 is singular
and refused. Outside this bounded numerical slice is refused, rather than
extrapolated. Positive outputs that cannot be represented as normal binary64
values and overflowed outputs are refused. Inputs must be finite, parameters
and relative tolerance positive. The arithmetic contract requires round-to-nearest
and at least 64 mantissa bits in long double. Default relative tolerance is
1e-10. Each scalar carries status and an absolute arithmetic diagnostic;
aggregate status reports any refused output while retaining other admitted
scalar outputs. A near-zero radial curvature can fail the relative allocation.
These diagnostics are empirical arithmetic and truncation allowances, not
certified libm enclosures or physical uncertainty.

For x<0.01, original logarithmic power series through x⁶ in f and x⁸ in g
avoid loss of g from cancellation. With L=ln(2/x), the leading terms are
f=L-1+x²(3L/2-5/4)+… and g=x²(L/2-1/4)+…; hence ΔΣ→rho_s r_s.
The next projected term is allocated as 4x⁸(L+1). For |x-1|<0.001,
f=1/3-2t/5+13t²/35-20t³/63+61t⁴/231, t=x-1;
g is obtained by integrating g'=xf from its exact value at one. The omitted
f term is allocated by |t|⁵. Elsewhere the closed form uses long double with
cancellation scales for g at small x and f near one. The binary64 projection
and inherited projection diagnostics are included in lens admission.

Permanent tests compare independent line-of-sight and cylindrical-shell
Simpson quadratures. For the former z=R sinh(u); for the latter exterior
r=R cosh(u) resolves the shell cusp. Both integrate to u=32. The omitted Σ
fraction is bounded by approximately 4 exp(-64)/x² relative to its shape;
for the tested x≥1e-5 this is below 1e-16. The exterior cylindrical mean tail
is likewise bounded by 8 exp(-64)/x². Mesh doubling differences must be ≤5%
of the respective 2e-9 and 2e-8 comparison budgets. Separate 60-digit mpmath
integrals to u=150 supply portable literal fixtures at both domain endpoints,
at one and on both sides of one, and at ten, checked to 1e-13 relative.
The fixtures integrate the density directly, including the interior sphere,
rather than evaluating the production closed forms. Finite-difference mass
and deflection identities, branch joins, analytic limits, tight-budget failures,
invalid inputs, arithmetic-mode and overflow refusals challenge the model.
The fixture calculation dependency is not needed by permanent C++ tests.

A standalone SDK consumer can use:

```cpp
#include <irred/nfw_halo.hpp>
auto lens = irred::lensing::nfw_lens({1e15, 0.2}, 0.2, 3e15, 1000);
// Check lens.status and each required scalar.status before reading values.
```
