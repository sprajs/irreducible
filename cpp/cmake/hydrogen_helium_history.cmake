# Root includes this fragment after the existing equilibrium/history owners.
target_sources(irred_core PRIVATE src/hydrogen_helium_history.cpp)
foreach(hhe_history_test IN ITEMS hydrogen_helium_history hydrogen_helium_history_peer)
  add_executable(test_${hhe_history_test} tests/test_${hhe_history_test}.cpp)
  target_link_libraries(test_${hhe_history_test} PRIVATE irred_core)
  target_compile_options(test_${hhe_history_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${hhe_history_test}_contract COMMAND test_${hhe_history_test})
endforeach()
