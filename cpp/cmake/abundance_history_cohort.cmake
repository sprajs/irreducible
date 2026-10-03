# Public-only local SDK proof caller. Root owns adding this fragment's include.
# No source is added to irred_core and no public API/ABI is introduced.
add_executable(test_abundance_history_cohort tests/test_abundance_history_cohort.cpp)
target_link_libraries(test_abundance_history_cohort PRIVATE irred_core)
target_compile_options(test_abundance_history_cohort PRIVATE
  -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off -fno-elide-constructors)
add_test(NAME abundance_history_cohort_source_contract COMMAND test_abundance_history_cohort)
add_test(NAME abundance_history_cohort_history_contract COMMAND test_abundance_history_cohort --history)
add_test(NAME abundance_history_cohort_allocation_contract COMMAND test_abundance_history_cohort --fault-history)
set_tests_properties(abundance_history_cohort_history_contract abundance_history_cohort_allocation_contract
  PROPERTIES LABELS "scientific;abundance-history-cohort")

add_executable(test_installed_abundance_history_cohort tests/test_installed_abundance_history_cohort.cpp)
target_link_libraries(test_installed_abundance_history_cohort PRIVATE irred_core)
target_compile_options(test_installed_abundance_history_cohort PRIVATE
  -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off -fno-elide-constructors)
add_test(NAME abundance_history_cohort_sdk_composition_contract COMMAND test_installed_abundance_history_cohort)
add_test(NAME abundance_history_cohort_wire_argument_refusal COMMAND test_installed_abundance_history_cohort --unexpected)
set_tests_properties(abundance_history_cohort_wire_argument_refusal PROPERTIES WILL_FAIL TRUE)
set_tests_properties(abundance_history_cohort_sdk_composition_contract
  PROPERTIES LABELS "scientific;abundance-history-cohort")

# Deliberate isolated caller flag mutation. The guard must refuse before any
# native work. This target does not mutate the library or earn its no-elision gate.
if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
  add_executable(test_abundance_history_cohort_unsafe_fp tests/test_abundance_history_cohort.cpp)
  target_link_libraries(test_abundance_history_cohort_unsafe_fp PRIVATE irred_core)
  target_compile_definitions(test_abundance_history_cohort_unsafe_fp PRIVATE IRRED_COHORT_EXPECT_REJECTED_FP=1)
  target_compile_options(test_abundance_history_cohort_unsafe_fp PRIVATE
    -Wall -Wextra -Wpedantic -ffast-math -ffp-contract=off -fno-elide-constructors)
  add_test(NAME abundance_history_cohort_unsafe_fp_contract COMMAND test_abundance_history_cohort_unsafe_fp)
endif()
