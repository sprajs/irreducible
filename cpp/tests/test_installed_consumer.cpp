// Standalone installed-library linkage contract; no survey/science qualification.
#include <irred/numerics.hpp>
#include <cmath>
int main() {
 const auto result=irred::numerics::log1p_checked(0.5);
 return result.status==irred::numerics::Status::ok && std::abs(result.value-0.4054651081081643819780131154643491)<1e-14 ? 0 : 1;
}
