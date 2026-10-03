# Conditional supplied late-depth interval under exact emitted working inputs.
target_sources(irred_core PRIVATE src/conditional_hydrogen_helium_drag.cpp)
foreach(conditional_hhe_test IN ITEMS conditional_hydrogen_helium_drag conditional_hydrogen_helium_drag_peer)
  add_executable(test_${conditional_hhe_test} tests/test_${conditional_hhe_test}.cpp)
  target_link_libraries(test_${conditional_hhe_test} PRIVATE irred_core)
  target_compile_options(test_${conditional_hhe_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${conditional_hhe_test}_contract COMMAND test_${conditional_hhe_test})
endforeach()
# Observe owning returns in the SDK harness without optional copy elision.
target_compile_options(test_conditional_hydrogen_helium_drag PRIVATE -fno-elide-constructors)
