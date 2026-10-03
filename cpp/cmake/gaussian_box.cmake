# Included by the integration owner after shared CMake coordination.
target_sources(irred_core PRIVATE src/gaussian_box.cpp)
add_executable(test_gaussian_box tests/test_gaussian_box.cpp)
target_link_libraries(test_gaussian_box PRIVATE irred_core)
target_compile_options(test_gaussian_box PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME gaussian_box_contract COMMAND test_gaussian_box)
