# Test-only prerequisite; no reference/physical-drag qualification is implied.
add_executable(test_conditional_drag_reference_certificate
  tests/test_conditional_drag_reference_certificate.cpp)
target_link_libraries(test_conditional_drag_reference_certificate PRIVATE irred_core)
target_compile_options(test_conditional_drag_reference_certificate PRIVATE
  -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off -fno-elide-constructors)
add_test(NAME conditional_drag_reference_certificate_contract
  COMMAND test_conditional_drag_reference_certificate)

add_executable(test_conditional_drag_reference_constants
  tests/test_conditional_drag_reference_constants.cpp)
target_link_libraries(test_conditional_drag_reference_constants PRIVATE irred_core)
target_compile_options(test_conditional_drag_reference_constants PRIVATE
  -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off -fno-elide-constructors)
add_test(NAME conditional_drag_reference_constants_contract
  COMMAND test_conditional_drag_reference_constants)
