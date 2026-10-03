target_sources(irred_core PRIVATE src/ideal_acoustic.cpp)
add_executable(test_ideal_acoustic tests/test_ideal_acoustic.cpp)
target_link_libraries(test_ideal_acoustic PRIVATE irred_core)
target_compile_options(test_ideal_acoustic PRIVATE
  -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME ideal_acoustic_contract COMMAND test_ideal_acoustic)
