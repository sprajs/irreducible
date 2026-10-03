# Explicit perfect-radiation-fluid/CDM slice, not a photon/neutrino hierarchy.
target_sources(irred_core PRIVATE src/linear_transfer.cpp)
foreach(transfer_test IN ITEMS linear_transfer linear_transfer_peer thermal_conformal_epoch)
  add_executable(test_${transfer_test} tests/test_${transfer_test}.cpp)
  target_link_libraries(test_${transfer_test} PRIVATE irred_core)
  target_compile_options(test_${transfer_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${transfer_test}_contract COMMAND test_${transfer_test})
endforeach()
