cmake_minimum_required(VERSION 3.20)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(CASE_DIR "${CMAKE_CURRENT_LIST_DIR}/cases")
file(GLOB CASE_FILES CONFIGURE_DEPENDS "${CASE_DIR}/*.json")

foreach(case_file IN LISTS CASE_FILES)
  get_filename_component(case_name "${case_file}" NAME_WE)

  add_test(
    NAME "008.${case_name}"
    COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_LIST_DIR}/run_case.py"
            --exe "$<TARGET_FILE:008>"
            --case "${case_file}"
            --emit-wkt
  )

  # 失敗時にstdout/stderrを見やすくするため（任意）
  set_tests_properties("008.${case_name}" PROPERTIES
    ENVIRONMENT "LC_ALL=C"
  )
endforeach()
