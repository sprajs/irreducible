# Conditional HeII/He history with explicit supplied H-neutral/thermal driver.
target_sources(irred_core PRIVATE src/helium_escape_history.cpp)
foreach(helium_escape_test IN ITEMS helium_escape_history helium_escape_history_peer)
  add_executable(test_${helium_escape_test} tests/test_${helium_escape_test}.cpp)
  target_link_libraries(test_${helium_escape_test} PRIVATE irred_core)
  target_compile_options(test_${helium_escape_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${helium_escape_test}_contract COMMAND test_${helium_escape_test})
  set_tests_properties(${helium_escape_test}_contract PROPERTIES TIMEOUT 300)
endforeach()
