# Bounded perfect-fluid matter transfer

This workstream freezes a standalone C++ scalar perturbation consumer for a
**declared self-interacting perfect radiation fluid, pressureless CDM and smooth
Lambda** in flat GR. Its homogeneous density uses the retained
[thermal background](thermal-neutrino.md) equation, but its perturbation closure
is distinct. A zero-baryon photon gas would free-stream; this model does not
reinterpret it as a perfect fluid. Baryons, collisionless massless/massive relics,
anisotropic stress and the photon hierarchy are excluded. This is a NEXT-14
prerequisite, not a standard matter transfer, ordinary sigma8, Planck base LCDM
or a [primary CMB prediction](primary-cmb-projection.md).

The proposed native owner is `irred/linear_transfer.hpp`. Until its native
controls have run, the equations and acceptance below are a frozen implementation
contract, not evidence of an implemented calculation. No CLI/ABI route is needed
for the bounded SDK experiment.

## Physical identity and equations

Conformal time has length units (Mpc), comoving k has units Mpc^-1, and
`mathcalH=a H/c` has units Mpc^-1; H is supplied in km/s/Mpc and
c=299792.458 km/s. Adopt Fourier derivatives `nabla^2 -> -k^2` and

\[
 ds^2=a^2[-(1+2\phi)d\eta^2+(1-2\phi)d\mathbf x^2].
\]

There is one scalar potential because the admitted fluids have zero anisotropic
stress. With N=ln a, x=k/mathcalH, V_i=mathcalH theta_i/k^2,
F_c=Omega_c a/(a^4 E^2), F_r=Omega_r/(a^4 E^2),
F_Lambda=Omega_Lambda a^4/(a^4 E^2), and
L=d ln mathcalH/dN=-1+F_c/2+2 F_Lambda, evolve

\[
\begin{aligned}
 \phi_N&=-\phi+\tfrac32(F_c V_c+\tfrac43F_r V_r),\\
 (\delta_c)_N&=-x^2 V_c+3\phi_N,&(V_c)_N&=(L-1)V_c+\phi,\\
 (\delta_r)_N&=-\tfrac43x^2 V_r+4\phi_N,&(V_r)_N&=L V_r+\delta_r/4+\phi.
\end{aligned}
\]

These are the conservation and Einstein momentum equations in regular scaled
velocity coordinates. The implementation evolves the cancellation-safe state
`Delta=delta_c+3V_c`, `S=delta_c-3delta_r/4`, `dV=V_r-V_c`, `V_c`, `phi`.
Writing F=F_c+4F_r/3, its equivalent equations are

\[
\begin{aligned}
 Delta_N&=-x^2V_c+6F_r dV,&S_N&=x^2dV,\\
 (dV)_N&=L dV+(Delta-S)/3,\\
 (V_c)_N&=(L-1)V_c+phi,&phi_N&=-phi+\tfrac32F V_c+2F_r dV.
\end{aligned}
\]

Recover `delta_c=Delta-3V_c`, `delta_r=4(Delta-3V_c-S)/3`.
This retains the small superhorizon comoving density without subtracting two
order-one final values. The independently checked Hamiltonian constraint is

\[
 C=x^2\phi+3(\phi_N+\phi)+\tfrac32(F_c\delta_c+F_r\delta_r)=0.
\]

In the evolved variables it is exactly
`C=x^2 phi+(3/2) F Delta-2 F_r S+6 F_r dV`. Report the maximum residual
at retained RK4 endpoints divided by the sum of absolute magnitudes of these
terms; its default gate is 1e-6. Intermediate RK stages are not accepted states. It is a consistency
gate, not a substitute for independent evolution or a certified forward bound.
The common background supplies `a^4 E^2` at each stage. Admitted input has
Omega_b=0 and no explicit thermal species. Nonphoton and photon homogeneous
radiation fractions are summed only under the explicit supplied perfect-fluid
closure. Zero radiation is a pressureless analytic limit. Zero CDM permits a
tracer limit for equation tests but does not define a matter-variance consumer.

The primary source for the gauge and conservation equations is
[Ma and Bertschinger (1995), equations 23, 43 and 64](https://arxiv.org/abs/astro-ph/9506072).
The [original reduced-system analysis](https://doi.org/10.1093/mnras/stx1662)
provides a second source for radiation/CDM limits; its potential sign differs
and must be translated. No upstream source code or text is copied.

## Primordial normalization and finite start

The explicitly named initial scalar is uniform-density curvature
`zeta=-phi+delta_rho/[3(rho+P)]` on superhorizon scales. A unit growing mode has
phi_rad=-2/3 and phi_matter=-3/5. In the Malik-Wands comoving-curvature convention
R tends to -zeta; the power is identical but signed cross transfers differ. Do
not silently rename this initial scalar R.

For r=Omega_c a_i/Omega_r and x_i=k/mathcalH_i in radiation domination, the
regular adiabatic expansion with phi0=-2/3 is

\[
\begin{aligned}
 \phi/\phi_0&=1-r/16-x_i^2/30,\\
 \delta_c/\phi_0&=-3/2-3r/16-7x_i^2/20,\\
 \delta_r/\phi_0&=-2-r/4-7x_i^2/15,\\
 V_c/\phi_0&=1/2+r/16-x_i^2/120,\\
 V_r/\phi_0&=1/2+r/16-x_i^2/20.
\end{aligned}
\]

Terms r^2, r x_i^2 and x_i^4 are omitted. Use S=0,
`dV=-phi0 x_i^2/24`; solve the Hamiltonian constraint for Delta exactly,
`Delta=(-x_i^2 phi-6 F_r dV)/[(3/2)F]`. This adjusts Delta only within the
omitted orders rather than leaving an inconsistent finite-start constraint.
The leading radiation limit is `Delta=-3 phi0 x_i^2/8`. Admit only r<=1e-4,
x_i<=1e-3 and F_Lambda<=1e-10. Refuse an unsupported start instead of moving it
silently. The radiation-free limit starts with phi=-3/5,
V_c=2phi/3 and delta_c=-2phi-(2/3)x_i^2 phi; its early Lambda fraction must be
<=1e-10. In pure Einstein-de Sitter this is the exact growing solution at every k.

Every requested point evolves separately from a_i, a_i/2 and a_i/4. The last two
changes, multiplied by eight and summed with time refinement diagnostics,
measure residual finite-start sensitivity. This does not turn asymptotic initial
conditions into exact initial data or certify arbitrary starts.

## Transfer and band-limited amplitude

The matter statistic is the CDM comoving density
Delta_c=delta_c+3V_c. Its returned signed transfer is T_c=Delta_c/zeta; CDM is the
entire admitted matter component. The Einstein-de-Sitter control is
T_c=(2/5)k^2 a/(H0/c)^2. Pressureless D/f alone does not supply the radiation-era
shape; this owner evolves radiation and metric variables explicitly.

A separate consumer declares a primordial law supported **only** on its finite
[k_min,k_max] interval:

\[
 \mathcal P_\zeta(k)=A_s(k/k_*)^{n_s-1},\qquad
 P_c(k,a)=\frac{2\pi^2}{k^3}\mathcal P_\zeta(k)T_c(k,a)^2,
\]

\[
 \sigma_{R,\mathrm{band}}^2=\int_{\ln k_{min}}^{\ln k_{max}}
 \mathcal P_\zeta(k)T_c(k,a)^2 W(kR)^2 d\ln k,
 \quad W(y)=3(\sin y-y\cos y)/y^3.
\]

A_s, n_s and k_* are supplied primordial inputs, not a supplied final sigma8.
For the explicitly named `sigma8_band`, R=8/h Mpc and h=H0/100. The law is zero
outside its declared support. The result retains both band endpoints and an
explicit band-limited identity. There is no estimated standard-spectrum tail or
ordinary full-support sigma8 output. W uses its regular Taylor expansion near
zero; negative oscillatory W is squared without clipping.

## Numerical and independent acceptance

Before execution, freeze portable strict floating point and wide state arithmetic.
Production uses fixed RK4 in N with phase-resolved steps bounded by
h_N<=h_max/(1+x). Run h_max, h_max/2 and h_max/4 at the earliest admitted start;
combine eight times the maximum of the last two changes with the separate start
diagnostic and arithmetic allowance. Default output acceptance is abs1e-8 +
rel3e-5. Work counts every RHS call, including every refused and discarded
refinement; no retry receives a fresh quota. Point/order and total work/payload
caps remain explicit. Missing requested outputs are withheld on refusal.
Refinement estimates remain empirical.

The variance consumer uses nested composite Simpson quadrature in ln k. The
4n-node grid evaluates each required transfer once and reuses its ordered values
for n/2n/4n comparisons. Propagate transfer diagnostics additively through
`2 abs(T) e_T+e_T^2`; retain window-quadrature sensitivity separately. Default
variance acceptance is abs1e-10 + rel1e-3. This native path predicts only the
finite declared model and law, with no observational likelihood.

Required durable controls are: exact Einstein-de-Sitter Delta_c/metric and
primordial-amplitude scaling; the pure-radiation analytic potential
`3 phi0(sin s-s cos s)/s^3`, s=k eta/sqrt(3); its separately integrated CDM tracer;
adiabatic/conservation/Einstein residuals; independent conformal-time ODE with
refined RK or a different high-order algorithm; separate initial/time/k/window
refinements; input/order/work/payload/move/mask adversaries. A matched CLASS/CAMB
comparison needs that same perfect-fluid closure and is a separate gate; a
standard photon/neutrino run cannot qualify this toy by resemblance.
