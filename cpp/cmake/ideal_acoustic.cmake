target_sources(irred_core PRIVATE src/ideal_acoustic.cpp)
add_executable(test_ideal_acoustic tests/test_ideal_acoustic.cpp)
target_link_libraries(test_ideal_acoustic PRIVATE irred_core)
target_compile_options(test_ideal_acoustic PRIVATE
  -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME ideal_acoustic_contract COMMAND test_ideal_acoustic)

# Independent portable TRACE controls; external arithmetic remains opt-in.
add_executable(test_ideal_acoustic_peer tests/test_ideal_acoustic_peer.cpp)
target_link_libraries(test_ideal_acoustic_peer PRIVATE irred_core)
target_compile_options(test_ideal_acoustic_peer PRIVATE
  -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME ideal_acoustic_peer_contract COMMAND test_ideal_acoustic_peer)
set_tests_properties(ideal_acoustic_contract ideal_acoustic_peer_contract
  PROPERTIES TIMEOUT 240)

option(IRRED_IDEAL_ACOUSTIC_MPFR_REFERENCE
  "Enable optional high-precision ideal-acoustic reference controls" OFF)
if(IRRED_IDEAL_ACOUSTIC_MPFR_REFERENCE)
  find_path(IRRED_ACOUSTIC_MPFR_INCLUDE_DIR mpfr.h REQUIRED)
  find_path(IRRED_ACOUSTIC_GMP_INCLUDE_DIR gmp.h REQUIRED)
  find_library(IRRED_ACOUSTIC_MPFR_LIBRARY NAMES mpfr REQUIRED)
  find_library(IRRED_ACOUSTIC_GMP_LIBRARY NAMES gmp REQUIRED)
  add_executable(test_ideal_acoustic_mpfr_reference EXCLUDE_FROM_ALL
    tests/test_ideal_acoustic_peer.cpp)
  target_compile_definitions(test_ideal_acoustic_mpfr_reference PRIVATE
    IRRED_IDEAL_ACOUSTIC_MPFR_REFERENCE)
  target_include_directories(test_ideal_acoustic_mpfr_reference PRIVATE
    ${IRRED_ACOUSTIC_MPFR_INCLUDE_DIR} ${IRRED_ACOUSTIC_GMP_INCLUDE_DIR})
  target_compile_options(test_ideal_acoustic_mpfr_reference PRIVATE
    -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  target_link_libraries(test_ideal_acoustic_mpfr_reference PRIVATE irred_core
    ${IRRED_ACOUSTIC_MPFR_LIBRARY} ${IRRED_ACOUSTIC_GMP_LIBRARY})
  add_test(NAME ideal_acoustic_mpfr_reference_contract
    COMMAND test_ideal_acoustic_mpfr_reference)
  set_tests_properties(ideal_acoustic_mpfr_reference_contract PROPERTIES TIMEOUT 720)
endif()
