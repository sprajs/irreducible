# Nonflat conserved GR FLRW geometry

`irred/curved_flrw.hpp` supplies the distinct compiled identity
`GR/nonflat-conserved-dust-radiation-lambda-FLRW/v1`. Its expanding homogeneous
GR state conserves pressureless matter and massless radiation separately and
retains a nonnegative cosmological constant. It adds actual spatial curvature
and curved null propagation, not a renamed flat fluid. It does not supply
species transfer functions, perturbations, lens masses, recombination, opacity,
CMB spectra or an observational likelihood.

Inputs are H0 in km/s/Mpc, present matter/radiation/Lambda fractions; curvature
is derived once as Ok=1-Om-Or-OL. H0 must be 1–1000, each supplied fraction
0–10, and queried redshift 0–10000. The bounded normalization has E(0)=1:

    E(z)^2 = Or (1+z)^4 + Om (1+z)^3 + Ok (1+z)^2 + OL
    chi(z) = integral_0^z dz'/E(z')
    D_C = (c/H0) chi
    D_M = (c/H0) sinh(sqrt(Ok) chi)/sqrt(Ok)       Ok>0
        = D_C                                    Ok=0
        = (c/H0) sin(sqrt(-Ok) chi)/sqrt(-Ok)      Ok<0
    D_A = D_M/(1+z), D_L = (1+z) D_M

Distances are Mpc. c=299792.458 km/s is the exact SI speed of light converted
to these units. Photon conservation, metric null propagation and distance
reciprocity are assumptions of these outputs. Peculiar velocities, observer
redshift corrections and lens-source two-epoch distance are outside this slice.
Do not subtract two angular distances to create a lens distance.

[Hogg, astro-ph/9905116v1](https://arxiv.org/abs/astro-ph/9905116v1),
*Distance measures in cosmology*, Eqs.13–15,17,20 supplies the dust/Lambda
expansion and radial/transverse/distance-duality equations. PDF SHA256 is
`cfde8333f9247c94d1b846441f81f6e60338755783999a77cc1631df20da98fe`.
That version calls curvature Omega_R; this API calls it Omega_k and reserves
Omega_r for genuine radiation. The added a^-4 term follows independent radiation
energy conservation; no source radiation or perturbation parity is claimed.
The coordinate-dependent spatial-curvature remark in that early version is
not adopted: intrinsic spatial curvature cannot be erased by a coordinate
change. Equations are implemented originally, with no copied source code/assets.

Admission checks the complete light path, not only its endpoint: the quartic
E² has at most one positive interior minimum, found from its analytic derivative.
Nonpositive or roundoff-unresolved minima are refused. Closed geometry refuses
the first antipode sqrt(-Ok)chi>=pi, including its numerical diagnostic margin.
It does not continue onto another image/winding branch or cross a turnaround.
A loitering model can have positive endpoint H² yet fail this admission.
Einstein-static H0=0 is incompatible with this H0-normalized expanding model and
is explicitly refused; it is not represented by adjusting fractional densities.

Adaptive quadrature returns an empirical radial error estimator. The curvature
Jacobian propagates it into the transverse distance with arithmetic allowances;
it is not a certified bound. Near-flat curvature uses an entire power series
through (Ok chi²)^4, avoiding sqrt(Ok) division. Angular/luminosity errors inherit
the transverse diagnostic divided/multiplied by (1+z). Their factors are fixed
query inputs, not redshift measurement uncertainty. Numerical tolerances are
caller policy, not physical uncertainties. Failed rows retain statuses with
zero unaccepted value fields; consumers must inspect status before use.

The batch admits at most 65,536 points, policy-limited output bytes, per-point
and shared callback quotas and bounded recursion. Owned payload accounts for
batch/rows, not allocator/RSS. Model copies retain independent state; moving
invalidates the source, self move preserves the destination. The supported
arithmetic profile uses round-to-nearest with at least 64 long-double mantissa
bits, no fast math and no floating-point contraction.

Permanent controls use exact Milne D_L=(c/H0)(z+z²/2), Einstein–de Sitter and
pure radiation distances, independent composite Simpson integration in ln(1+z),
flat-limit continuity from both curvature signs, an interior forbidden band,
closed loitering antipode refusal, invalid/static inputs, quotas and moves.
`test_curved_flrw_sdk.cpp` is a standalone C++ consumer printing bounded
nonflat predictions. Native model access does not qualify a released probe or
any cosmological interpretation.
