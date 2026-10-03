// Standalone installed-header/archive consumer; no private/in-tree header.
#include <irred/supplied_shell_projection.hpp>
#include <iostream>
#include <utility>
int main() {
  using namespace irred::projection;
  using S = irred::numerics::Status;
  const double k[]{0, 1}, chi[]{0, 2}, amplitudes[]{.25, .5, 1, -.25};
  const std::string_view ids[]{"first", "second"};
  const AtomRole roles[]{AtomRole::radial_atom, AtomRole::boundary_response};
  SuppliedShellView input;
  input.k_mpc_inverse = k; input.chi_mpc = chi; input.amplitudes = amplitudes;
  input.k_ids = ids; input.shell_ids = ids; input.roles = roles;
  input.source_identity = "synthetic installed finite atoms";
  input.mode_origin = "declared scalar unit-zeta; physical closure absent";
  input.ordering_provenance = "literal two-by-two k-major table";
  const unsigned ell[]{0, 1, 8, 0};
  auto owner = acquire_supplied_shell_projection(input);
  const auto *original = owner.source();
  auto result = owner.evaluate({ell});
  auto retained = std::move(result);
  owner = PreparedProjection{};
  if (!original || retained.status() != S::ok || retained.source() != original ||
      result.source() || retained.rows().size() != 8 ||
      retained.rows()[0].signed_value != .75 ||
      retained.rows()[1].signed_value != 0 ||
      !retained.rows()[4].absolute_numerical_radius ||
      retained.work_after().total != retained.work_before().total + retained.work_delta().total)
    return 1;
  std::cout << "PASS standalone supplied-shell public consumer\n";
}
