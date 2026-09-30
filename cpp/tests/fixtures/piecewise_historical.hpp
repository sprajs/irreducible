#pragma once
#include <array>
namespace irred::test_reference {
// Frozen factual outputs: original lib/cosmology.py qbins calculator, Python
// 3.12.14/NumPy2.2.6, read-only/bytecode-disabled, BLAS/OMP threads1.
// SourceSHA256
// 6c1d9c8cf13e7126aa31b6dc7976b39ee6c34b030dc4735cc60f3ffbf2209e1b. Fixed
// edges0,.1,.3,.6,1,2.5; right-bin q at internal edges, final last bin.
// Calculator endpoints qi=-3/2 are admitted here; historical strict prior
// excludes them. No sampler/prior/campaign acceptance. Shared libm/equations.
// Budgets E abs2e-14+rel2e-12, I abs2e-15+rel2e-10, q exact assignment.
struct PiecewiseHistorical {
  std::array<double, 5> q;
  double z, E, I, assigned_q;
};
inline constexpr std::array<PiecewiseHistorical, 18> piecewise_historical{{
    {{0, 0, 0, 0, 0}, 0, 1, 0, 0},
    {{0, 0, 0, 0, 0}, 0.1, 1.1, 0.09531017980432493, 0},
    {{0, 0, 0, 0, 0}, 0.3, 1.3, 0.262364264467491, 0},
    {{0, 0, 0, 0, 0}, 0.6, 1.6, 0.4700036292457356, 0},
    {{0, 0, 0, 0, 0}, 1, 2, 0.6931471805599454, 0},
    {{0, 0, 0, 0, 0}, 2.5, 3.5, 1.252762968495368, 0},
    {{-0.4, -0.4, -0.2, 0.1, 0.3}, 0, 1, 0, -0.4},
    {{-0.4, -0.4, -0.2, 0.1, 0.3},
     0.1,
     1.0588528529217844,
     0.09715029563521174,
     -0.4},
    {{-0.4, -0.4, -0.2, 0.1, 0.3},
     0.3,
     1.1704854282221195,
     0.2766257670858054,
     -0.2},
    {{-0.4, -0.4, -0.2, 0.1, 0.3},
     0.6,
     1.3819976853697582,
     0.512095939679986,
     0.1},
    {{-0.4, -0.4, -0.2, 0.1, 0.3},
     1,
     1.7664783943032025,
     0.7675780631675353,
     0.3},
    {{-0.4, -0.4, -0.2, 0.1, 0.3},
     2.5,
     3.656434459784413,
     1.3508429426223894,
     0.3},
    {{-3, 2, -3, 2, -3}, 0, 1, 0, -3},
    {{-3, 2, -3, 2, -3}, 0.1, 0.8264462809917356, 0.11033333333333344, 2},
    {{-3, 2, -3, 2, -3}, 0.3, 1.3641641467609638, 0.2993510848126233, -3},
    {{-3, 2, -3, 2, -3}, 0.6, 0.9005614875101675, 0.5739193718527777, 2},
    {{-3, 2, -3, 2, -3}, 1, 1.7589091552932956, 0.8937198563219169, -3},
    {{-3, 2, -3, 2, -3}, 2.5, 0.574337683361076, 2.546022359412469, -3},
}};
} // namespace irred::test_reference
