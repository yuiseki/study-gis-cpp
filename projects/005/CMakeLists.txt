cmake_minimum_required(VERSION 3.20)
project(gdal_geodesic_metrics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GDAL CONFIG REQUIRED)  # GDAL::GDAL :contentReference[oaicite:5]{index=5}
find_package(nlohmann_json REQUIRED CONFIG)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GEOGRAPHICLIB REQUIRED IMPORTED_TARGET geographiclib)

add_executable(005 main.cpp)
target_link_libraries(005 PRIVATE
  GDAL::GDAL
  PkgConfig::GEOGRAPHICLIB
  nlohmann_json::nlohmann_json
)
