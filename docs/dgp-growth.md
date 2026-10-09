# Flat self-accelerating DGP dust background and quasistatic growth

This compiled model has its own identity, distinct from GR LCDM. It evaluates
flat self-accelerating Dvali-Gabadadze-Porrati expansion and linear pressureless,
subhorizon quasistatic growth. Inputs are H0 in km/s/Mpc and the present dust
fraction Omega_m0; flat closure fixes Omega_rc=(1-Omega_m0)^2/4. Radiation,
massive relics, curvature, nonlinear screening, full five-dimensional
perturbations and CMB are absent. The caller must establish k/(aH)>>1
applicability: no supplied wavenumber or dataset is qualified by this calculation.

The source equations are the modified Friedmann closure in
[Deffayet (2001), astro-ph/0010186](https://arxiv.org/abs/astro-ph/0010186)
and quasistatic growth in
[Koyama and Maartens (2006), astro-ph/0511634](https://arxiv.org/abs/astro-ph/0511634v1):

    E=sqrt(Omega_m0 a^-3+Omega_rc)+sqrt(Omega_rc),
    Omega_m(a)=Omega_m0 a^-3/E^2, h_N=d ln H/d ln a=-3 Omega_m(a)/(1+Omega_m(a)),
    beta=1-2 H r_c[1+h_N/3], mu=1+1/(3 beta),
    D_NN+(2+h_N)D_N-(3/2)mu Omega_m(a)D=0.

Flat closure has H0 r_c=1/(1-Omega_m0), with c=1. Algebraically equivalent
production mu=1-(1-Omega_m(a)^2)/[3(1+Omega_m(a)^2)] avoids divergent beta.
Omega_m0=1 is an explicitly named infinite-r_c Einstein-de Sitter endpoint:
E=a^-3/2, mu=1, D=a, f=1; finite beta is not defined there.
The self-accelerating branch has a known ghost caveat:
[Gorbunov, Koyama and Sibiryakov (2006), hep-th/0512097](https://arxiv.org/abs/hep-th/0512097v1).
Numerical success supplies no viability, fitting or CMB claim.

## Frozen contract, declared before implementation

Admit finite 0.05<=Omega_m0<=1, 0<H0<=1000 and arbitrary input points 1e-4<=a<=1.
This is a bounded mathematical dust domain, not a claim of negligible radiation
at its lower endpoint. E, Omega_m, mu, D and f are dimensionless; H is km/s/Mpc.
D approaches the matter normalization D/a->1, not D(1)=1. Optional masks avoid
unrequested growth; duplicate/unordered points retain their original order.

Production evolves g=D/a and v=g_N with adaptive step-doubled RK4 in N=ln a.
At a_init=1e-8 define x=(1-Omega_m0)a_init^(3/2)/(2 sqrt(Omega_m0)). The finite
regular-mode approximation g=1-11x/12, v=-11x/8 follows by inserting
h_N=-3/2+(3/2)x and mu Omega_m=1-(8/3)x into the g equation. O(x^2) initial
truncation is retained as a separate numerical allowance; it is not an exact
finite-endpoint normalization. Default D/f relative allowance is 1e-8. Local
step diagnostics receive a 1/128 allocation; sum, initial truncation and
binary64 storage/roundoff diagnostics remain empirical, not certified bounds.
Failed error/work budgets withhold outputs without changing allowances.

Independent references use explicit midpoint plus Richardson extrapolation
on fixed N meshes of 32768/65536/131072 steps, evolving original D,D_N and
computing E/beta directly from unsimplified source equations. This differs from
production RK4 state, step selection and mu algebra. Named D/f comparisons
allocate 2e-7 relative; independent mesh and a_init/10 refinements occupy <=5%
of that allocation. Production default/refined comparisons allocate 1e-8.
Analytic controls check exact EdS, E(1)=1, modified Friedmann closure, mu limits,
and early g/f coefficients. Adversarial controls cover unsupported parameters,
nonfinite/mixed points, invalid/empty masks, quotas, impossible tolerances,
rounding and copy/move ownership. These are synthetic mathematical controls.

Hard ceilings before allocation: 65536 points, 1 GiB row payload, 200000 step
trials per point and 4000000 total. Defaults: 4096 points, 16 MiB, 20000 point
trials and 200000 total. Rejected trials count; callbacks and trial diagnostics
are retained. Background-only outputs perform no ODE calls. Native SDK only:
no CLI/ABI, transfer function, sigma8, likelihood or observational result.

Source review verified Koyama/Maartens v1 page 3 equations 32 and 34 and page 2
quasistatic assumptions; original PDF SHA256
`c3713bf257b6b545e2787799effa563050572963860f9e10cefe21598d8f3c5a`.
Ghost caveat source v1 PDF SHA256
`f32ad9d3beb5e38ca98a5eb1098f01b7293273f4a980e2e03837cbea594a9bc7`.
These targeted equation/abstract checks do not imply full-paper qualification.
