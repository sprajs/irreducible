# Coarse synthetic predictive consumer; native equation owners are unchanged.
target_sources(irred_core PRIVATE src/gaussian_predictive_abi.cpp)
add_executable(test_gaussian_predictive_abi tests/test_gaussian_predictive_abi.cpp)
target_link_libraries(test_gaussian_predictive_abi PRIVATE irred_core)
target_compile_options(test_gaussian_predictive_abi PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME gaussian_predictive_abi_contract COMMAND test_gaussian_predictive_abi)
# Same forwarding hook as native owner controls; original arithmetic unchanged.
if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
  target_compile_definitions(test_gaussian_predictive_abi PRIVATE IRRED_PREDICTIVE_WORK_OBSERVATION)
  target_link_options(test_gaussian_predictive_abi PRIVATE "-Wl,--wrap=_ZN5irred8numerics8choleskyESt4spanIKdLm18446744073709551615EEmmNS0_10ArithmeticE")
endif()
