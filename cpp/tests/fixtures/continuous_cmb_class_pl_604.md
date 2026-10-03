# Original CLASS source table and independent PL projection references

This fixture contains generated numerical facts and independently calculated
geometry references. It freezes one illustrative supplied CLASS point, rather
than a likelihood best fit, observation, spectrum or physical validation.

The producer is [CLASS 3.3.0](https://github.com/lesgourg/class_public/tree/0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18),
revision `0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18`, by Julien Lesgourgues,
Thomas Tram, Nils Schoeneberg and the CLASS contributors. Its pinned README
allows free use with a publication citation to
[CLASS II: Approximation schemes](https://arxiv.org/abs/1104.2933).
No top-level standard LICENSE/COPYING file was found at this revision; no
upstream SPDX or public-domain status is assigned. The independently authored
test, wrapper and optional reference code use Irreducible's BSD-3-Clause license.
The fixture includes generated source values, not copied CLASS/HyRec code,
atomic assets, MPFR/GMP binaries, official likelihood products or observational
data. The producer used HyRec2020 by Yacine Ali-Haimoud and Christopher Hirata;
its source and assets remain separate upstream ingredients.

## Supplied physical identity

The input fixes H0=67.32117 km/s/Mpc, physical omega_b=.02238280,
omega_cdm=.1201075, Tcmb=2.7255 K, N_ur=2.046 and one .06 eV ncdm species,
with Tncdm/Tcmb=.7137658555036082, degeneracy1 and ksi0. Curvature, scalar-field
and fluid densities are zero. YHe=.2454006 is supplied. Recombination is HyRec;
reionization is reio_camb with tau_reio=.05430842, width.5 and exponent1.5,
He-II redshift3.5 and width.5. The scalar initial mode is synchronous adiabatic,
with actual curvature_ini=1. The producer's analytic primordial input is
A_s=2.100549e-9, n_s=.9660499, alpha_s=0 and k_pivot=.05/Mpc.

The exporter runs input through transfer. It does not run harmonic, lensing
or output modules, despite broader product switches in the input. The fixture
channels and dimensionless transfer references are per unit signed curvature.
No primordial spectrum, Tcmb or microkelvin factor is applied.

## Exact source order and ancestry

The header contains all604 original increasing conformal times and four
original source channels, without interpolation, resampling, dropped endpoints
or a zero-outside rule. Row order is `[eta_Mpc,T0,T1,T2,Psrc]`; all source channels
are per conformal Mpc. This is source-column position0, actual source-k index93,
with exact k=`0x1.4709f9265abbap-8` /Mpc and observer eta=
`0x1.ba323d47c7a8dp+13` Mpc. The actual CLASS source type indices are
`[2,3,0,1]`, not ascending type order. Psrc uses CLASS's positive
sqrt6 times visibility times internal scattering quadrupole convention.
The two actual multipoles are2 and19; requested20 selected stored19.

The permanent law is the native linear interpolation within each of603 original
time cells. The original CLASS producer's continuum/source-grid, support and
transfer discretizations are separate. See the
[projection guide](../../../docs/continuous-cmb-projection.md) for the common
Fourier, angular, units, kernels and endpoint conventions.

| Retained member | SHA256 |
| --- | --- |
| Pinned CLASS README | `24738c5700f7f0754ad29f14874d35b7fd5a314b8d6cb7e6729b93b7158c1bf9` |
| Exporter source v2 | `da85c020b90bfdc6ba254fb2f3ae972cd454436f1dd9e527b6fcfa25c6b8d3ab` |
| Supplied input INI | `8d23c8483882a950c99a55ce4bd6caeb772d749a5388f6bc54f01ba2421d6879` |
| Exporter actual object | `80f9fa23815d9cca4769d1f36cc5d8f60c0896b3bf5afbca6be594776dad525f` |
| Exporter actual binary | `41f5ed66a9c027d36476ca628f658c4cd3eb6112bd47aa49eebd3edcb010a7ca` |
| Exporter build receipt | `a291cdfbeb7ceef4a94f3e9e0a27e0ef9d0f13f0be7d41a64c4af3ac7f815881` |
| Original raw export,315719 bytes | `2c871e216f2b599eb6f31543ce1678b731a66a77a539ec65b41595fdcf25332b` |
| Structural transport source | `3871b6dca23d064ab3c4830d6d89e62eeb218fa1189f128e3e685472e8b63d1e` |
| Original coarse transport,65985 bytes | `6ca4a504758e01b7009f31b557e31cef174cc3b2ac8dd454823b600d64e19d0f` |
| Independently authored original angular caller | `f4ab4d994c3e6a07c4c0fc4741115437745adb1e6a016be1da07f53316e9811d` |
| Original angular reference/native raw log | `c6ed18f881dd79673647f3264d7e12663948dc3a823095c3a91a9c99b0a5ed7f` |
| Original run receipt | `5a1d05f2e3743b9ad03ec1bdfeebea845525320fc3b9dccc073991e04fa5b7c0` |
| Independent source review | `bec6f42ab64aefdce61284c347e8b3c4d09a28ac4292cf05ba69b22805b63ae5` |
| Independent actual-evidence review | `01da272e51254fa301d0ea83998f934e2a492cecf7fe91a398e3642acd208858` |

The original coarse transport also retains eight historical CLASS final-q-table
interpolants. They are unused by this fixture's independent expectations.
All original interpolated-q and direct-q discrepancies and execution records
remain preserved. Passing this fixture supplies no attribution of those
remaining discrepancies.

## Earned reference and portable conversion

The four original expectations are the aggregate MPFR256/N256 real intervals
in the retained raw log. Their source is the independent analytic PL time-cell
integral followed by angular Legendre/spin projection, with actual emitted-rule
moment defects, a real-phase Taylor remainder and propagated MPFR rounding.
MPFR192/N128,N256 refinement and MPFR256/N256 agreed with common overlapping
enclosures; imaginary intervals contained zero. No reference Bessel, native
Simpson callback or CLASS transfer interpolation was used. The standalone
[reference source](../reference/continuous_cmb_angular_mpfr.cpp) retains this
method, including the signed Doppler and regular observer-endpoint controls.
Its input adapter now uses these exact compiled source literals; the numerical
body is preserved from the original caller.

Every original256-bit center and earned radius is retained as a hexadecimal
string beside a portable lower/upper binary64 bit pair. Source transcription
checked the exact serialized dyadic inequalities `lower <= center-radius` and
`upper >= center+radius`. The pair uses the two outward neighbors of the rounded
center. This enlarges the enclosure to include reference conversion loss; the
tiny MPFR radius is never assigned to a rounded binary64 center. The four
synthetic references have the same treatment. Default CI decodes exact bit
patterns and does not need decimal parsing, MPFR or external files.

Acceptance retains the original per-output allocation
`1e-9 + 1e-5*abs(native)`. A test-local nearest-binary64 guard evaluates an upper
bound of `max(abs(native-lower),abs(native-upper))+native_time+native_arithmetic`
against a downward allocation. Each basic assembly operation is padded outward
by one adjacent binary64 value under explicit IEC559/53-bit/nearest admission.
Native quadrature/arithmetic estimates stay separate from the independently
bounded reference. Radial/source-grid optional fields remain absent. This local
test guard is not a production-wide interval or libm certificate.

The original independent run passed all four combined geometry gates and the
signed finite-endpoint/regular1/5 controls in31.978s including compilation.
It charged459344269 explicit MPFR calls,771936 vector-cell visits and6453 native
nodes. The new permanent regression subsequently passed its own fresh strict
GCC 16.2.1 C++20 Release build and native CTest contract at
`93a14f559affeac0a9d0a68ace7f29b662e2de45`. Its 0.010663 s test checked the retained
four-output brackets and sign/endpoint/ownership/refusal controls without an
MPFR dependency; the existing continuous-projection contract also passed.
The adapted optional reproducer remains separately uncompiled and unexecuted.
The original MPFR reference run, this native regression and integrated CI
retain distinct execution identities and qualification scopes.
