cmake_minimum_required(VERSION 3.20)
project(mini03_geodesic_buffer_geojson LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(nlohmann_json REQUIRED CONFIG)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GEOGRAPHICLIB REQUIRED IMPORTED_TARGET geographiclib)

add_executable(002 main.cpp)
target_link_libraries(002 PRIVATE PkgConfig::GEOGRAPHICLIB nlohmann_json::nlohmann_json)
