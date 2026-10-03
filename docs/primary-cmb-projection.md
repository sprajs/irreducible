# Primary CMB projection prerequisites

NEXT-19 remains open. A primary temperature prediction requires a closed
photon/baryon/species/metric perturbation owner, an initial primordial law and a
source-required thermal/ionization/visibility history. The bounded
[perfect-fluid CDM transfer](linear-transfer.md) excludes the photon hierarchy;
it supplies neither a temperature transfer nor a CMB spectrum.

A narrower [continuous scalar projector](continuous-cmb-projection.md) now
integrates supplied finite-support split sources in the C++ library. Its bounded
analytic/angular controls exercise transfer geometry, including polarization
sign and regular radial limits. The caller still owns the physical source,
signed initial mode, normalization and missing history/support errors. This
prerequisite neither evolves the hierarchy nor integrates a primordial spectrum.

The future first consumer must freeze one gauge and metric convention, scalar
initial mode, Fourier units and primordial normalization with its hierarchy.
Species density definitions and distributions must equal the retained thermal
background. Photon/baryon Thomson momentum exchange must cancel in total
momentum, photon polarization/collision terms need their own closure, and any
collisionless species needs momentum/multipole truncation and refined tails.
Tight coupling is an approximation with its own switch/refinement diagnostics;
an untracked switch cannot replace a hierarchy calculation.

A line-of-sight temperature transfer has the explicit convention

\[
 \Theta_\ell(k)=\int_0^{\eta_0}S_T(k,\eta)
 j_\ell[k(\eta_0-\eta)]\,d\eta,\qquad
 C_\ell^{TT}=4\pi\int\mathcal P_\zeta(k)
 |\Theta_\ell(k)|^2 d\ln k.
\]

The source definition must specify intrinsic photon temperature, gravitational
redshift, Doppler velocity, time-dependent metric and scattering/polarization
terms with signed gauge conventions. The
[original line-of-sight formulation](https://arxiv.org/abs/astro-ph/9603033)
is a source for that contract; this document does not implement its source or
projection. Kelvin spectra require the explicitly supplied monopole temperature
exactly once; dimensionless temperature spectra are a distinct output identity.

Present [pure-H histories](recombination-drag.md) and finite-endpoint Thomson
visibility retain omitted boundary survival mass. They cannot provide the full
source-required present-day last-scattering/reionization history. Helium,
multilevel/temperature/rate closure and endpoint/tail error must qualify before
composition. A physical drag epoch is a different weighted optical-depth
criterion and cannot substitute for last-scattering visibility.

Independent gates must separately refine thermal/visibility tails, initial
conditions, time, wavenumber, momentum/multipole hierarchy and projection. Test
Einstein/conservation/tight-coupling limits, a separately implemented hierarchy
and line-of-sight route, then pinned CLASS/CAMB temperature spectra under exactly
matched species/recombination/primordial inputs. Shared atomic/rate/source
ancestry must be declared. Numerical agreement is not source or observational
qualification.

A measured primary likelihood additionally requires exact released assets,
versions, beam/calibration/instrument/foreground/estimator definitions, ordered
covariance and primary/lensing/shared-source overlap. Those prerequisites are
separate from one closed synthetic projection. No Planck base LCDM, CMB lensing,
posterior or measured-likelihood capability is claimed here.
