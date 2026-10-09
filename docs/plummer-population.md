# Isotropic collisionless Plummer population

`irred/plummer_population.hpp` provides the native C++ model
`self-gravitating-isotropic-collisionless-Plummer-mass-traces-tracer-SI/v1`.
It adds a stationary collisionless phase-space distribution to the existing
[finite-core Plummer sphere](plummer-sphere.md). Supplied SI parameters are
mass M in kg, core scale b in metres and Newtonian coupling G in m³ kg⁻¹ s⁻².
This is an isolated, self-gravitating, one-component isotropic model. Mass
traces the tracer population. It supplies neither an arbitrary tracer embedded
in a different potential nor a distribution function for Hernquist or NFW.

The existing compiled sphere owns density, potential, mass and projection.
Preparation consumes those operators, retains rho(r), relative potential
Psi(r)=-Phi(r), enclosed mass and projected surface density, and never
reimplements the gravity laws. The ordered supplied radius grid serves two
separate operators: a three-dimensional local radius r and a projected radius
R having the same supplied numeric coordinate. The local state at r does not
represent a line-of-sight average at R. Repeated coordinates retain their
original order. Enclosed mass is retained with its native status as provenance;
only density, potential and projected surface density gate the population.
An unrelated native tidal-curvature refusal does not invalidate these fields.

For specific binding energy epsilon=Psi-v²/2 and speed v in m/s, the mass
phase-space density is

    f(epsilon) = K epsilon^(7/2) for epsilon>0, otherwise 0
    K = 24 sqrt(2)/(7π³) b²/(G⁵ M⁴),

with units kg s³ m⁻⁶. The local vector-velocity density is f/rho, in s³ m⁻³;
isotropy makes it depend only on the vector's magnitude. The local speed
probability density is 4πv²f/rho, in s m⁻¹. Speed zero has positive f and
vector density but exactly zero speed density. The escape speed is sqrt(2Psi),
in m/s. It is returned with a numerical diagnostic, not used as a rounded
threshold. There is no Gaussian approximation, random generator or synthetic
sample emission.

The exact local one-axis velocity variance is Psi/6, in m²/s². The projected
mass-weighted line-of-sight second moment is 3πPsi(R)/64, also in m²/s². The
latter is a density-weighted integral along the entire line of sight, rather
than the local variance at R. It becomes light-weighted only under an additional
constant mass-to-light assumption. No measured stellar kinematics, luminosity
weighting, anisotropy, binary stars, selection, aperture or PSF, observed mass
constraint, full projected velocity PDF or parameter posterior is supplied.
The local one-axis PDF is also outside this bounded implementation.

The inversion and normalization are reconstructible mathematics. The declared
Plummer profile gives rho(Psi)=C Psi⁵, C=3b²/(4πG⁵M⁴). In the spherical isotropic
Eddington inversion, rho'(0)=0 and

    f(epsilon) = [1/(sqrt(8)π²)] integral_0^epsilon
                 rho''(Psi)/sqrt(epsilon-Psi) dPsi.

Since rho''=20C Psi³, B(4,1/2)=32/35 gives K above. Independently, integrating
4πv²f over 0<v<sqrt(2Psi) recovers rho because B(3/2,9/2)=7π/256. The normalized
variable q=v/sqrt(2Psi) has density
512 q²(1-q²)^(7/2)/(7π), and q² has a Beta(3/2,9/2) law. Thus
E[v²]=Psi/2, E[v⁴]=5Psi²/14, and one-axis variance is E[v²]/3=Psi/6.
Combining rho sigma_r² with the line-of-sight density integral gives
3πPsi(R)/64. The steady isotropic Jeans equation is

    d(rho sigma_r²)/dr = -rho G M(<r)/r².

[Plummer (1911)](https://doi.org/10.1093/mnras/71.5.460) supplies the historical
profile context and [Eddington (1916)](https://doi.org/10.1093/mnras/76.7.572)
the inversion context. Targeted access to both original DOI sources returned
HTTP 403; neither original-paper equation locators nor a complete paper
reproduction are claimed. The mathematics above specifies the implemented
conditional model. All implementation and comparison code is original.

Preparation owns a copy of the supplied source and ordered radius states,
computes the log coefficient once and calls two native operators per radius.
Defaults admit at most 4096 radii, 8192 native evaluations and 4 MiB requested
retained payload. Every required radius must be admitted; failed preparations
retain available diagnostics and withhold velocity batches without dropping
or renormalizing radii. The default relative allocation is 1e-10; native
prerequisites use one sixteenth of it. Copies own their states independently;
moves explicitly invalidate the source object. The C++ SDK stays useful without
a CLI or ABI operation.

Velocity batches borrow requests only during evaluation, preserve all input
rows and own the resulting source parameters, radius coordinates, speeds and
output diagnostics. Each request references a retained radius index. They make
no additional native gravity calls. Defaults permit 65536 rows and velocity
evaluations, with 32 MiB requested retained payload; even invalid rows count
against the original work allowance. Payload arithmetic is checked before
allocation, and allocation failure returns a work-limit refusal. The caps bound
requested retained native payload, not allocator metadata or process RSS.
An empty velocity batch is supported; an empty preparation grid is invalid.

Both stages require finite valid inputs, positive parameters and tolerance,
round-to-nearest, at least 64 wide mantissa bits and extended exponent range.
Nonzero outputs not representable as normal binary64 are refused. Inside
positive support, logarithmic evaluation prevents avoidable intermediate
underflow; a genuinely unrepresentable positive DF or PDF remains a refusal,
even if other scalar outputs are representable. Each scalar has its own status;
aggregate status reports any refusal without discarding admitted scalar fields.

Support uses a wide energy difference and an **empirical diagnostic interval**
from inherited native potential error plus wide arithmetic. If its upper endpoint
is nonpositive, the row declares outside support and returns exact zero f and
PDFs. If its lower endpoint is positive, it evaluates the bound law and propagates
energy and density sensitivity using log1p/expm1. An interval straddling zero
returns ambiguous-support conditioning refusal. Rounded escape and adjacent
floating-point speeds therefore remain ambiguous when the inherited diagnostic
cannot resolve their sign. This is an explicit numerical admission rule, not
a certified enclosure, proof of physical support or uncertainty in supplied
M,b,G. No arbitrary rounded-escape cutoff silently converts positive support
to zero.

Permanent scientific controls integrate every admitted interior velocity node
with a theta-midpoint rule, v=sqrt(2Psi) sin(theta). No unresolved endpoint is
sampled, dropped or replaced by zero. A declared 1e-7 point allocation resolves
near-tail nodes; pooled arithmetic diagnostics must stay below 1e-10. The
256/512-node normalization, recovered density and second/fourth moments must
refine within 1e-10 before comparison at 2e-9. Independently, 128/256-node
native density-weighted LOS integrals must refine within 1e-10 before the same
2e-9 comparison. These gates allocate at most 5% of the downstream budgets to
mesh differences; they are not asymptotic convergence certificates.
Native enclosed mass and a centered density-pressure derivative challenge the
Jeans equation at a separate 2e-7 allocation. Nine 70-/90-digit Beta-inversion,
density-shell and LOS fixtures agree to 1e-50 and check reported diagnostics,
including literal conversion uncertainty. SI scaling, center, zero speed,
ambiguous escape and adjacent speeds, admitted outside support, positive bound
underflow, invalid rows, tight budgets, resource caps and copy/move ownership
provide adversarial controls. No Python oracle is required by permanent tests.

The compiled `example_plummer_population` consumer emits a synthetic SI table
from one owned preparation and one velocity batch. Ambiguous endpoints are
retained as `NA`, rather than printed as zero or omitted. It is a deterministic
population-law readout, not a Monte Carlo sample or observed dataset.
