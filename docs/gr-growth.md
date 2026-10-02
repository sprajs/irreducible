# Bounded pressureless GR growing mode

The native `irred/gr_growth.hpp` consumer retains a compatible radiation-free
`cosmology::Expansion` with its `LCDM` identity. It admits 1e-6<=Omega_m<=1,
flat nonnegative Lambda and 1e-8<=a<=1. Other expansion hypotheses are refused.
This explicit pressureless-only approximation contains no radiation or thermal
massive relics. The distinct early/late supplied-ruler provider keeps its
positive-radiation admission; it is not substituted for this background.

For linear pressureless GR perturbations the growing mode is

    D(a)=(5/2) Omega_m E(a) integral_0^a da'/(a'^3 E(a')^3),
    f(a)=d log D / d log a.

Normalization is D/a->1 at early times, not D(1)=1. The source is
[Heath (1977)](https://doi.org/10.1093/mnras/179.3.351); the expression and
matter-era normalization are reviewed in
[Carroll, Press and Turner (1992), section 3.5](https://w0.ned.ipac.caltech.edu/level5/Carroll/Carroll3_5.html).
These equations are implemented originally; no external code/assets are copied.

The existing late scalar E equation and the growth consumer's scaled a^3 E^2
coordinate share a private compiled LCDM helper. Existing public late domains
and scalar arithmetic are preserved. The bounded growth domain reaches earlier
scale factors solely under its declared pressureless approximation.

With Q=(Omega_Lambda/Omega_m)a^3 and a'=a t^2, production evaluates
J=integral_0^1 2t^4/(1+Qt^6)^(3/2)dt, then
D=(5/2)a sqrt(1+Q)J and f=-3/[2(1+Q)]+1/[(1+Q)^(3/2)J].
This removes the singular endpoint from the original integral. Omega_m=1 returns
exact D=a and f=1 without callbacks. The lower matter bound limits the narrow
transition that must be resolved; unsupported lower densities are refused.

An output mask requests D and/or f. Both are dimensionless, and each retains
its own optional scalar, cause and empirical absolute error estimate. Their
numerical diagnostics propagate quadrature, subtraction and storage effects;
they are not posterior variance or certified universal bounds. The default
relative numerical allowance is 1e-10; integration reserves a factor of 32.
Unattainable budgets refuse outputs without changing the allowance.

Point, callback, depth and owned-payload limits precede output allocation.
The hard limits are 65536 points and 1 GiB simultaneous owned payload, with
smaller caller limits by default. Counts include failed quadrature callbacks.
The payload excludes the retained background, borrowed input, allocator metadata,
recursion stack and RSS. Copies own their background; moves invalidate the old
owner and self move preserves it. Arithmetic requires round-to-nearest and the
supported wide host profile.

Permanent controls compare an originally written independent GR ODE RK4 mesh
and early growing-series refinements with an exact Einstein–de Sitter limit.
Named D/f comparisons allocate 1e-8 relative, with reference refinement at most
5% of that allocation. Low matter, early times, cancellation in f, invalid
closures, independent mixed failures, quotas, omissions and ownership are tested.
The shared background change also requires affected background/SN/BAO regressions.

No redshift-space distortion likelihood, scale-dependent transfer function,
thermal/radiation perturbation closure, observational inference or arbitrary
cosmological accuracy is supplied by this bounded calculation.
