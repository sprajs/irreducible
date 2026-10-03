# Root owns the single include in the common CMakeLists.txt.
target_sources(irred_core PRIVATE src/supplied_shell_projection.cpp)
foreach(shell_test IN ITEMS supplied_shell_projection supplied_shell_projection_installed)
  add_executable(test_${shell_test} tests/test_${shell_test}.cpp)
  target_link_libraries(test_${shell_test} PRIVATE irred_core)
  target_compile_options(test_${shell_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${shell_test}_contract COMMAND test_${shell_test})
  set_tests_properties(${shell_test}_contract PROPERTIES TIMEOUT 240)
endforeach()
# Independent reference/checker execution is an explicit root-leased qualification
# step against a saved native --reference-cases corpus, not an unbounded CTest.
