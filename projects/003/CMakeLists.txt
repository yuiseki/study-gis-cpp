cmake_minimum_required(VERSION 3.20)
project(geojson_geodesic_metrics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(nlohmann_json REQUIRED CONFIG)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GEOGRAPHICLIB REQUIRED IMPORTED_TARGET geographiclib)

add_executable(003 main.cpp)
target_link_libraries(003 PRIVATE
  PkgConfig::GEOGRAPHICLIB
  nlohmann_json::nlohmann_json
)
