# Distinct supplied linear Kaiser/geometric AP mean, no released likelihood.
target_sources(irred_core PRIVATE src/windowed_linear_power.cpp)
foreach(power_test IN ITEMS windowed_linear_power windowed_linear_power_peer)
  add_executable(test_${power_test} tests/test_${power_test}.cpp)
  target_link_libraries(test_${power_test} PRIVATE irred_core)
  target_compile_options(test_${power_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${power_test}_contract COMMAND test_${power_test})
endforeach()
