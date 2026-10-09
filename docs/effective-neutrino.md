# Effective thermal neutrino radiation partition and GR acceleration

The standalone C++20 `irred/effective_neutrino.hpp` API adds a compiled physical
closure: total early-time Neff is partitioned between explicitly declared
thermal neutrino/antineutrino pairs and the remaining massless nonphoton fluid.
It retains that mapped source, present relic normalization and flat Lambda
closure. The genuinely new `ThermalBackground::evaluate_stress_energy` coarse
batch supplies scaled total and per-species density/pressure, equation of state
and flat-GR deceleration from the same existing FD momentum owner. These APIs
have no CLI/C ABI route and supply no BBN, recombination, physical drag,
perturbations, growth or CMB spectra.

## Physical source and conventions

`EffectiveNeutrinoModel` requires explicit finite H0>0 in km/s/Mpc, physical
omega_b>=0, omega_cdm>=0, Tcmb>0 Kelvin and total Neff>=0. Up to sixteen species
supply nonnegative mass in eV, positive temperature ratio Tnu/Tgamma and
positive degeneracy. Degeneracy counts neutrino/antineutrino pairs: one pair is
CLASS `deg_ncdm=1`, native `statistical_weight=2`. No mass hierarchy or mass-sum
conversion is inferred. Baryon and CDM densities exclude these relics.

For r0=(4/11)^(1/3) and Fnu=(7/8)r0^4, the exact massless FD energy-density
moment in [thermal neutrino](thermal-neutrino.md) gives

    Ni = degeneracy_i * (temperature_ratio_i/r0)^4
    Nur = Neff - sum Ni
    omega_ur = omega_gamma * Fnu * Nur.

`omega_gamma` here is the physical photon density from the existing authoritative
`map_thermal_physical_model`, using its fixed SI/G convention; the first map's
wide witness supplies this conversion without copying a blackbody equation.
Kelvin species temperatures are ratio*Tcmb, and native weights are
2*degeneracy. A negative Nur is refused without clipping, tolerance or species
removal. Flat Lambda is the existing present-energy-density closure; a
normalization whose numerical diagnostic reaches the closure boundary remains
a numerical refusal. It is not a new physical exclusion.

The source-linked first external reference realization is total Neff=3.046,
Tcmb=2.7255 K, one mass=0.06 eV, ratio=0.71611, degeneracy=1. Its Nur is
2.0327983725164924, rather than 2.046. The source contract is Prospector's
`references/lcdm/class-planck-primary-six-parameter-v1.json`; the existing
Planck VI selected-source review and pinned CLASS explanatory/input conventions
motivate this named effective thermal realization. It is not proof of exact
CAMB equality or a nonthermal neutrino-decoupling distribution. The historical
ratio r0 with Nur=2.046 remains a different source realization.

Diagnostics retain the ideal partition, its arithmetic estimate, actual early
Neff measured from the retained background's emitted scaled radiation coordinate,
today's physical relic density omega_nu=Omega_species(today)*h^2 and Lambda.
The present relic density includes kinetic energy. No independently conserved
dust/radiation split or BBN prediction is supplied. Scalar physical-map wide,
operation and cast witnesses remain available. Copies own independent storage;
moves transfer all source/diagnostics/background and invalidate the old owner;
self move preserves the owner. Admitted sources and failed normalization
callbacks remain inspectable after a closure or work refusal.

## Same-state stress energy and acceleration

With D=a^4*rho/rhocrit0 and Q=a^4*p/rhocrit0, both dimensionless,

    D = Omega_gamma + Omega_ur + (Omega_b+Omega_cdm)*a
        + Omega_Lambda*a^4 + sum C_i*I_rho(m_i*a/T_i0)
    Q = (Omega_gamma+Omega_ur)/3 - Omega_Lambda*a^4
        + sum C_i*I_pressure(m_i*a/T_i0)
    C_i = g_i*T_i0^4/(2*pi^2*rhocrit0)
    w_total = Q/D
    q = -a*addot/adot^2 = (1+3*w_total)/2.

These are the zero-chemical-potential collisionless FD moments and flat-GR
acceleration equation, with pressureless baryons/CDM and pLambda=-rhoLambda.
The compiled shared momentum kernel evaluates density and pressure together;
no finite-difference pressure, new relic interpolation table or second source
normalization is introduced. At a=1 pressure still requires momentum work even
though the existing E/H owner has its exact present-day normalization shortcut.
Per-species density/pressure are scaled by a^4/rhocrit0, not eV^4.

Finite 0<=a<=1 is admitted. At a=0, positive D0 gives the declared analytic
radiation-limit w=1/3 and q=1; these are not finite physical rho or H. A
radiation-free D0=Q0=0 refuses these ratios. Ordered invalid rows survive as
explicit failures. Batch status reports processing admission; row status reports
scaled stress admission, while requested w/q retain their own statuses.
Successful scaled stress can remain inspectable when an output budget refuses
w or q. A failed species never yields a partial aggregate or a fabricated q.

## Numerical estimates and resources

D's estimate includes inherited present-normalization error times a^4, scaled
FD density diagnostics and arithmetic. Q uses the same closure estimate, FD
pressure diagnostics and arithmetic scaled by the sum of positive radiation,
Lambda and species-pressure centers. This preserves the diagnostic at signed
pressure cancellation. Require D-eD>0. The ratio estimate is

    ew = (eQ + abs(Q/D)*eD)/(D-eD) + arithmetic
    eq = 1.5*ew + acceleration assembly arithmetic.

Reporting adds the measured signed-to-binary64 cast loss. The default w/q
admission allocation is **2e-10 absolute + 2e-10 relative**. Absolute allocation
is essential at zero crossings and for cold pressure. Estimates are empirical
numerical diagnostics, not certified continuous bounds, random uncertainties or
physical-input errors. Existing E/H and distance/ruler budgets are unchanged.

The effective source consumer adds a conservative 128*epsilon(binary64)
component-conversion estimate. Scalar map arithmetic/casts and the two species
temperature conversions are covered: FD density's logarithmic temperature
response lies in [3,4], and pressure's in [4,5], so their temperature-cast losses
remain below this allowance. Flat closure propagates every mapped non-Lambda
component direction at today's normalization, separately from the runtime FD
normalization diagnostic. The Neff subtraction's absolute estimate is
128*epsilon(wide)*(Neff+sum Ni); the photon-normalized remainder estimate and
physical-density reporting loss propagate separately, including into Lambda.
This conservatively preserves cancellation rather than assigning a relative
error to an arbitrarily small residual.

Stress batches use the retained momentum method and existing per-moment and
whole-batch ceilings. Both density and pressure callbacks, failed nodes and
fallback work count. Preparation callbacks remain separately exposed; no
implicit larger work limit is installed. Conservative payload admission covers
owned model/map/background temporaries, row/species capacities and preparation
storage; borrowed input owners, allocator metadata, recursion stack and RSS are
excluded. No cache/server or observation readback is needed.

`observable_request` produces the same explicit physical source for the existing
retained thermal distance/ruler and BAO consumers. Supplied finite z_drag and
nonempty origins remain required; this composition earns no physical drag or
released-compression qualification. Its distance/ruler preparation retains its
own original numerical diagnostics and does not inherit the stress-only source
conversion estimates automatically.

## Evidence and limits

Owner regressions cover pure Lambda, matter, radiation and mixed polynomial
states, pressure cancellation, exact radiation limits, inherited closure,
ordered failures, source partition/double-count refusal, casts/rounding, work and
payload limits, retained method and copy/move lifetimes. Independent direct
infinite-momentum mpmath1.3.0 60/90-digit controls test transition/cold species;
refinement must occupy <=5% of the frozen 2e-10+2e-10*abs(reference) allocation.
A five-point log-H continuity diagnostic uses its own **8e-7** discretization
allocation and factor-two refinement <=4e-8; it is not an independent2e-10
pressure reference. The independently installed archive/header consumer
composes actual source preparation, acceleration and an existing distance.

The matched CLASSv3.3.0 background driver at source commit
`0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18` independently checks D,Q,w,q and
scaled per-species density/pressure at a=0.001,0.01,0.1,0.5,1. It explicitly
maps the native SI/G photon density, dimensionless mass/temperature and FD
phase-space amplitude into the existing runtime background owner without
changing its source or object files. qmax=80,bins4096/8192 and background
10000/20000 with tolerance1e-12 give maximum reference-refinement fraction
0.010225 of the frozen2e-10abs+2e-10relative allocation, below5percent. The
native direct-momentum consumer matches every coordinate within that unchanged
allocation. Immutable custody pins belong to the package evidence record.

Matched CLASS comparison requires explicit physical constant/temperature/state
mapping. Different fixed constants or source defaults are model-convention
changes and cannot be hidden by widening numerical budgets. These finite
controls establish their named numerical scope, not arbitrary-domain accuracy,
observational qualification, a posterior or native full standard cosmology.
