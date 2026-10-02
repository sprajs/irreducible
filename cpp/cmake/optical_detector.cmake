# Root includes this fragment after the existing photometry/detector owners.
target_sources(irred_core PRIVATE src/optical_detector.cpp)
foreach(optical_detector_test IN ITEMS optical_detector optical_detector_peer)
  add_executable(test_${optical_detector_test} tests/test_${optical_detector_test}.cpp)
  target_link_libraries(test_${optical_detector_test} PRIVATE irred_core)
  target_compile_options(test_${optical_detector_test} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
  add_test(NAME ${optical_detector_test}_contract COMMAND test_${optical_detector_test})
endforeach()
