cmake_minimum_required(VERSION 3.20)
project(009_s2_vs_cgal_nef_s2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(nlohmann_json REQUIRED CONFIG)

find_library(S2_LIB NAMES s2 REQUIRED)
add_library(s2::s2 UNKNOWN IMPORTED)
set_target_properties(s2::s2 PROPERTIES IMPORTED_LOCATION "${S2_LIB}")

# CGAL (Nef_S2を使う)
find_package(CGAL REQUIRED COMPONENTS Core)

add_executable(009 main.cpp)
target_link_libraries(009 PRIVATE
  nlohmann_json::nlohmann_json
  s2::s2
  CGAL::CGAL
  CGAL::CGAL_Core
)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
