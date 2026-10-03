# Pure explicitly supplied zero-mass FD radiation/CDM/Lambda; no photons/baryons.
target_sources(irred_core PRIVATE src/massless_fd_transfer.cpp)
foreach(fd_test IN ITEMS massless_fd_transfer massless_fd_trial massless_fd_peer installed_massless_fd)
  add_executable(test_${fd_test} tests/test_${fd_test}.cpp)
  target_link_libraries(test_${fd_test} PRIVATE irred_core)
  target_compile_options(test_${fd_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${fd_test}_contract COMMAND test_${fd_test})
  set_tests_properties(${fd_test}_contract PROPERTIES TIMEOUT 300)
endforeach()
