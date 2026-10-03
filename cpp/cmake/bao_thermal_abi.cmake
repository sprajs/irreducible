# One coarse thermal consumer; shared scientific owners remain unchanged.
target_sources(irred_core PRIVATE src/bao_thermal_abi.cpp)
add_executable(test_bao_thermal_abi tests/test_bao_thermal_abi.cpp)
target_link_libraries(test_bao_thermal_abi PRIVATE irred_core)
target_compile_options(test_bao_thermal_abi PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME bao_thermal_abi_contract COMMAND test_bao_thermal_abi)
# The source-only installed consumer is compiled against an actual installed
# exported target by the separately allocated package-validation controller.
