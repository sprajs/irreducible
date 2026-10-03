# One private owner; no public header/ABI or existing physical consumer changes.
target_sources(irred_core PRIVATE src/thermal_elementary_postcheck.cpp)

# First qualified arithmetic implementation request is deliberately narrower
# than the conditional mathematics. Other toolchains receive causal refusal.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux"
   AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_SIZEOF_VOID_P EQUAL 8
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
  file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/src/thermal_elementary_postcheck.cpp" _elementary_cpp)
  file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/src/thermal_elementary_postcheck.hpp" _elementary_hpp)
  file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_thermal_elementary.cpp" _elementary_fixture)
  file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/tests/thermal_elementary_wire.hpp" _elementary_wire)
  string(SHA256 _elementary_source "${_elementary_cpp}${_elementary_hpp}${_elementary_fixture}${_elementary_wire}")
  set_source_files_properties(src/thermal_elementary_postcheck.cpp PROPERTIES COMPILE_OPTIONS "-fstack-usage")
  foreach(_target test_thermal_elementary test_thermal_elementary_fixture test_thermal_elementary_no_elide)
    add_executable(${_target} tests/test_thermal_elementary.cpp)
    target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off -fstack-usage)
    target_compile_definitions(${_target} PRIVATE
      IRRED_ELEMENTARY_SOURCE_SHA256="${_elementary_source}" IRRED_ELEMENTARY_WRAP_ALLOCATIONS)
    target_link_options(${_target} PRIVATE -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc -Wl,--wrap=free)
    if(_target STREQUAL "test_thermal_elementary_no_elide")
      target_sources(${_target} PRIVATE src/thermal_elementary_postcheck.cpp)
      target_compile_options(${_target} PRIVATE -fno-elide-constructors)
    else()
      target_link_libraries(${_target} PRIVATE irred_core)
    endif()
    if(NOT _target STREQUAL "test_thermal_elementary")
      target_compile_definitions(${_target} PRIVATE IRRED_ELEMENTARY_FIXTURE_ONLY)
    endif()
  endforeach()
  add_test(NAME thermal_elementary_native_contract COMMAND test_thermal_elementary)
  add_test(NAME thermal_elementary_fixture_contract COMMAND test_thermal_elementary_fixture --facts)
  add_test(NAME thermal_elementary_no_elide_contract COMMAND test_thermal_elementary_no_elide --facts)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME thermal_elementary_exact_reference_contract
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/thermal_elementary_reference.py"
      "$<TARGET_FILE:test_thermal_elementary_fixture>" "$<TARGET_FILE:test_thermal_elementary>")
    add_test(NAME thermal_elementary_lifetime_frame_contract
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/thermal_elementary_frames.py"
      "${CMAKE_CURRENT_BINARY_DIR}" "$<TARGET_FILE:test_thermal_elementary_fixture>"
      "$<TARGET_FILE:test_thermal_elementary_no_elide>")
  endif()
endif()
