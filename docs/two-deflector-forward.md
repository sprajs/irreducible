# Fixed two-deflector extended-source pixels

The native source in `irred/two_deflector_forward.hpp` implements a bounded
synthetic question: what ordered electron means follow from two smooth
deflectors on one lens plane, an extended steady source, an affine detector
cutout and a supplied finite shift PSF? Fresh ordinary native controls and an
installed SDK joint Gaussian-stage control have passed for the fixed synthetic
cases. The pixel diagnostics remain empirical; the separate directed interval
and complete scene-to-log-density gates remain open. It does not qualify
B1608+656,
ACS FLT pixels, an elliptical-density SPLE, an observed velocity dispersion,
exactly four images, a noise law or parameter inference.

Use dimensionless tangent coordinates `u=theta/theta_scale`; `theta_scale` is
in radians. Counterclockwise angles rotate a right-handed x/y plane. For
component j, `v=u-center`,

    Q = R(angle) diag(1,q^-2) R(angle)^T
    d = sqrt(core^2 + v^T Q v)
    Psi_j = strength*d
    alpha_j = strength*Q*v/d.

These are elliptical **potentials** with positive cores. Both components share
one lens plane. With the traceless shear matrix
`Gamma=[[shear_1,shear_2],[shear_2,-shear_1]]`,

    Psi0 = Psi_0 + Psi_1 + u^T Gamma u/2
    PsiLambda = Lambda*Psi0 + (1-Lambda)*|u|^2/2
    b(u) = Lambda*(u-alpha_0-alpha_1-Gamma*u).

The physical potential is `theta_scale^2*PsiLambda`. The convergence of each
component is
`strength*(core^2*tr(Q)+det(Q)*|v|^2)/(2*d^3)`, which is nonnegative. The
admitted `0<Lambda<=1` also gives a nonnegative sheet. Supplied strengths have
no inferred stellar-kinematics or cosmological mass convention. Positive cores
can retain central images. The brightness calculation does not enumerate
stationary points and makes no image-completeness claim.

The pulled-back physical source brightness is

    S(b) = peak*exp(-(b-source_center)^T P (b-source_center)/2)
    P = R(source_angle) diag(major_width^-2,minor_width^-2) R^T.

The source center and widths are in u coordinates; peak is electrons per
second per radian squared. This is an explicitly supplied steady observer-band
electron-brightness law. Lensing conserves its surface brightness: no extra
magnification factor belongs in a pixel integral. No spectrum, area, throughput
or QE map is inferred from this input, and no source epoch variation is hidden.

Detector coordinates map to u as `origin+A*x`. Each pixel is a strictly
increasing axis-aligned rectangle in detector coordinates; A may rotate/shear
it. Positive PSF states shift the observer/image field after lensing:
`I(u)=sum_h w_h S(b(u-delta_h))`. States are kept in supplied order. Each mass
is a positive integer multiple of 2^-32 and their integer sum must be 2^32.
This finite law has exact normalization, with no tolerance, dropped state or
renormalization. A singleton PSF has mass 1. Its finite shifts do not claim a
measured or sampled Gaussian PSF.

For the full exposure T and uniform sky B in the same brightness units,

    mean_i = T*theta_scale^2*abs(det(A)) *
             (integral_pixel I(origin+A*x) dx dy + B*pixel_area).

The result is in electrons. The angular Jacobian enters once. The source field
is not truncated or renormalized to the supplied cutout; omitted image flux
remains omitted. This affine geometry is distinct from ACS TAN-SIP and its
lookup distortion. Original SCI/ERR/DQ/WCS products, calibrated source masks,
lens-light/dust/PSF operators and a measured joint noise law need separate
source-qualified consumers.

Preparation consumes one `TwoDeflectorScene` by value, retains its exact
supplied fields and origin, and prepares Q0/Q1/P once. Pass `std::move(scene)`
to transfer an acquired scene; an lvalue argument copies in the caller before
native admission. The owner is move-only. Moving invalidates the old owner;
self move keeps it valid. PSF/cutout metadata are copied into fixed arrays.
`means(span<ForwardPixelRectangle>,policy)` borrows the ordered rectangle batch
and returns a move-only exact-count row array through `rows()`. Every row keeps
its original index, status, optional mean, optional computed diagnostic and
actual work. A malformed batch has a causal status before output allocation;
a row failure keeps the attempted prefix and does not erase other rows.
An empty required batch never produces a successful empty scientific result.

The fixed domain is theta scale [1e-6,1e-4] radians, exposure (0,3600] seconds,
strength [0,4], core [0.02,1], q [0.5,1], center coordinate magnitudes at most4,
source widths [0.02,4] with major>=minor, angles [-pi,pi], shear norm<=0.3,
Lambda [0.25,1], peak/sky [0,1e15]. Decimal endpoint constants are the native
binary64 values. All inputs must be finite. Detector coordinate magnitudes
are at most1e9; A coefficients at most16, origin at most8, PSF shifts at most1.
A is nonsingular with a guarded condition ratio at most100; every shifted
footprint stays within |u_x|,|u_y|<=8. These are synthetic domains, not priors
on a measured lens.

Positive-source pixels must be resolved before quadrature. Each shifted
footprint has center u_c and enclosing radius rho. Let
`rmin=max(0,|u_c-center_j|-rho)` and
`dmin=sqrt(core_j^2+rmin^2)`, since lambda_min(Q)=1. Then

    L = abs(Lambda)*(1 + sum_j strength_j/(q_j^2*dmin_j) + norm(Gamma))
    delta = L*rho + recorded map arithmetic allowance
    V = (|b(u_c)-source_center|*delta + delta^2/2)/minor_width^2.

Require guarded V<=1/4 for every state. This resolved-cell domain is not a
quadrature certificate. The source-zero branch still validates the entire
scene and every footprint, but omits unused maps/quadrature. Sky and a true
zero law are computed directly; positive numerical underflow is a refusal.

The mean uses the existing scalar adaptive Simpson integrator, nested y inside
x. Its original field-sample cap covers preflight maps and all attempted PSF
terms across all rows, including failures. Primitive callback limits never
reset that cap. Returned work separately counts started rows, comparisons,
PSF preflights, field terms, deflector evaluations, inner/outer starts and
callbacks, and causal failed starts. Required callback failure causes survive
secondary NaN detection by the scalar integrator.

For requested `absolute+relative*abs(mean)` electrons, the literal native
allocations are inner quadrature35%, outer quadrature35%, field/PSF arithmetic
20%, projection/area/exposure5%, with5% reserved inside the same total for an
independent reference. The positive absolute floor proposes quadrature work;
all final contributors and their outward stored sum must pass the final
allocation checks before a mean is available. Inner error uses the maximum
observed inner diagnostic times outer length, alongside outer error. Field
arithmetic includes map cancellation, prepared rotations, exponent sensitivity,
PSF accumulation and callback rounding; final projection includes output
rounding. Every sample carries affine-coordinate product magnitudes through a
global map-sensitivity bound, so cancellation before the bounded u coordinate
cannot erase that diagnostic. Inner/outer contributors include the returned
binary64 integral's larger adjacent half spacing (computed wide at subnormals)
and a depth-dependent wide-arithmetic allowance. Their internal quadrature
proposals use30% each, leaving room inside the same35% final allocations.
A positive returned outer integral must stay positive before sky is combined.
Projection separately propagates retained determinant, width/height, area,
sky, prefactor and final casting diagnostics. These are empirical diagnostics, including standard-library
primitives, and do not certify unsampled integrands or universal libm error.
An unavailable diagnostic differs from a computed zero. A final budget refusal
retains its computed contributors but withholds the mean.

Defaults and hard work limits are4096 pixels,9 PSF states,8million whole field
terms,500000 per pixel and depth20. Origin size is at most65536 bytes. Default
requested native payload is4MiB; its hard bound is16MiB. Admission precedes
the single exact-count row allocation. The conservative source charge includes
two owner and two result headers, incoming scene, two actual origin capacities
plus terminators, and all row objects, covering return without optional
elision. Borrowed input, recursion stack, allocator bookkeeping/cookies and
process RSS are excluded and must be recorded separately in process controls.
No covariance/factor or source/pixel arrays are copied to cross an ABI; this
consumer has no CLI or C ABI registration.

The permanent controls retain a no-deflection Gaussian/CDF limit, a genuine
two-offset-component reference, source/exposure linearity, PSF linearity,
component-label exchange, reflection and mass-sheet invariance. The MST scales
source center and both widths by Lambda while keeping peak, geometry, PSF and
sky fixed. Every pixel integrand stays fixed; unlensed source rate scales by
Lambda^2. All actual source widths still satisfy the domain. No H0 constraint
follows. A later requested time-delay route must use compatible retained
distance diagnostics and the same potential; it is absent from this mean-only
owner.

The ordinary native run passed the independent principal-coordinate positive
GL8 reference with 8/16-panel refinement, charging its empirical refinement and
arithmetic allowance inside the same total error budget and 5% reservation.
It also passed the analytic limits, invariances and adversarial domain, work,
payload and numerical-refusal controls. Private source regressions share the
implementation's ancestry. This passage does not supply a high-precision
enclosure or a universal quadrature or standard-library error certificate.

A synthetic observation caller can prepare `statistics::Gaussian` once with
fixed SPD electron covariance, ordered pixel IDs and explicit electron product
measure/provenance, then evaluate the complete residual vector `data-means`.
It must withhold the joint density if any required mean fails. A product of
scalar marginal densities loses supplied covariance. Gaussian noise here is a
separate declared synthetic assumption; ERR, seed/address independence,
selection conditioning and a posterior are not supplied by this forward law.
The installed caller's independent rank-one inverse/determinant control checks
the Gaussian stage at the same supplied native-rounded residual vector. It
does not propagate independent forward-mean or residual-subtraction uncertainty
into a complete scene-to-log-density interval; that gate remains open.
The installed control retained all16 ordered means and passed the joint
Gaussian inverse/determinant and exact-zero correlation witnesses with one
retained covariance factor. Its receipts retain both scene preparations and all
base, refused and exact-zero forward batches, including structural work with
zero field samples.

All source and numerical controls are original; no third-party code or image
assets are copied. Measured B1608 Paper I/SPLE, calibrated pixels, operator
order/noise, delay PDF/sign/order, kinematics and environment remain open.
Numerical acceptance of this synthetic law cannot close those gates.
