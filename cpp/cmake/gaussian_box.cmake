# Included by the integration owner after shared CMake coordination.
target_sources(irred_core PRIVATE src/gaussian_box.cpp src/box_enclosure.cpp
  src/retained_gaussian.cpp src/gaussian_box_heldout.cpp)
add_executable(test_gaussian_box tests/test_gaussian_box.cpp)
target_link_libraries(test_gaussian_box PRIVATE irred_core)
target_compile_options(test_gaussian_box PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME gaussian_box_contract COMMAND test_gaussian_box)

add_executable(test_gaussian_box_heldout tests/test_gaussian_box_heldout.cpp)
target_link_libraries(test_gaussian_box_heldout PRIVATE irred_core)
target_compile_options(test_gaussian_box_heldout PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME gaussian_box_heldout_contract COMMAND test_gaussian_box_heldout)
add_executable(test_gaussian_box_heldout_allocation tests/test_gaussian_box_heldout_allocation.cpp)
target_link_libraries(test_gaussian_box_heldout_allocation PRIVATE irred_core)
target_compile_options(test_gaussian_box_heldout_allocation PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME gaussian_box_heldout_allocation_contract COMMAND test_gaussian_box_heldout_allocation)

add_executable(test_gaussian_box_heldout_peer tests/test_gaussian_box_heldout_peer.cpp)
target_link_libraries(test_gaussian_box_heldout_peer PRIVATE irred_core)
target_compile_options(test_gaussian_box_heldout_peer PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
add_test(NAME gaussian_box_heldout_peer_contract COMMAND test_gaussian_box_heldout_peer)
add_executable(test_gaussian_box_heldout_installed tests/test_gaussian_box_heldout_installed.cpp)
target_link_libraries(test_gaussian_box_heldout_installed PRIVATE irred_core)
target_compile_options(test_gaussian_box_heldout_installed PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
# Genuine runtime parse/schema checks; declaration tokens are deliberately
# synthetic and never establish installed SDK/compiler/binary provenance.
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME gaussian_box_heldout_wire_contract
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_gaussian_box_heldout_wire.py
    $<TARGET_FILE:test_gaussian_box_heldout_installed>)
