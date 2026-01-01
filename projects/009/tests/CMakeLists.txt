cmake_minimum_required(VERSION 3.20)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(FAST_DIR "${CMAKE_CURRENT_LIST_DIR}/cases_fast")
set(EXP_DIR  "${CMAKE_CURRENT_LIST_DIR}/cases_experimental")

file(GLOB FAST_CASES CONFIGURE_DEPENDS "${FAST_DIR}/*.json")
file(GLOB EXP_CASES  CONFIGURE_DEPENDS "${EXP_DIR}/*.json")

function(add_009_case case_file label)
  get_filename_component(case_name "${case_file}" NAME_WE)

  add_test(
    NAME "009.${label}.${case_name}"
    COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_LIST_DIR}/run_case.py"
            --exe "$<TARGET_FILE:009>"
            --case "${case_file}"
  )

  set_tests_properties("009.${label}.${case_name}" PROPERTIES
    LABELS "${label}"
    ENVIRONMENT "LC_ALL=C"
  )
endfunction()

foreach(f IN LISTS FAST_CASES)
  add_009_case("${f}" "fast")
endforeach()

foreach(f IN LISTS EXP_CASES)
  add_009_case("${f}" "experimental")
endforeach()
