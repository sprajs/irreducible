# Retained supplied-source CMB geometry regression

The permanent native regression compares the continuous projector with four
independently bounded temperature/E-mode transfer facts for one original
604-row supplied source. Its input and output roles follow the
[projection guide](continuous-cmb-projection.md). The fixture preserves the
original source order, one exact k, multipoles2/19, original observer and every
finite endpoint. It supplies numerical geometry evidence for that piecewise
linear source; it does not construct photon/species/metric or opacity physics.

[Fixture provenance](../cpp/tests/fixtures/continuous_cmb_class_pl_604.md)
records the illustrative CLASS3.3.0 physical point, upstream attribution,
exporter/input/raw-output identities and the independent reference ancestry.
The expected values come from analytic PL time-cell integration and angular
Legendre/spin integration at MPFR192/256 bits, with an earned finite-rule moment,
Taylor and rounding bound. They are independent of the native Bessel/Simpson
algorithm and the CLASS final-q interpolator. They are not native golden values.

The default regression uses portable outward binary64 brackets of the original
MPFR intervals, including reference conversion loss. It preserves the original
`1e-9+1e-5*abs(native)` combined allocation. It checks all four outputs,603
completed cells, source/order/masks/lifetime, retained failed work, signed
constant Doppler endpoints and the positive quadrupole/E observer limit.
It needs the compiled C++ SDK and standard library only, with conservative
floating-point flags and the native admitted arithmetic profile. No CLASS,
MPFR/GMP or external fixture files are required by ordinary CI.

The separate optional [high-precision source](../cpp/tests/reference/continuous_cmb_angular_mpfr.cpp)
reproduces the independently authored analytic PL/angular calculation from the
same compiled inputs. It requires explicit MPFR/GMP linking and execution. Its
original work limits are600M explicit MPFR calls,1M attempted vector cells,
three1024-node synthetic calls and one default2M-node native call; its historical
preparation ceiling20000 is an explicit override of default16384, with604 rows
inside both. It is not a default CI workload. The adapted optional reference
has not been built or run; the original independently executed method supplies
the retained reference ancestry.

The registered `continuous_cmb_class_pl_contract` and existing
`continuous_cmb_projection_contract` both passed a fresh GCC 16.2.1 C++20 Release
build at `93a14f559affeac0a9d0a68ace7f29b662e2de45`, with `-fno-fast-math` and
`-ffp-contract=off`. The new retained-source contract took 0.010663 s; both CTest
cases took 0.35 s, and the serial configuration/build/test transaction took 63.592 s.
This is scoped local validation of the unchanged numerical policy and supplied
PL geometry. Current integrated CI and publication remain separate gates.

The original independent same-PL control passed its complete four-output budget
and sign/endpoint controls. Original CLASS comparisons remain separately
preserved: final-q interpolation differences were large, while a distinct direct-q
control retained0.1271%/0.1933% temperature and2.9358%/0.4389% E differences.
Their remaining source/time/radial/support cause has no accepted bound.
This regression supplies no smooth-source-grid or omitted-support certificate,
full CLASS match, primordial k integration, C_l, lensing prediction, official
likelihood or observational/joint inference qualification. Production
`radial_error_estimate` and `source_grid_error_estimate` remain absent.
