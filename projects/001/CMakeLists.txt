cmake_minimum_required(VERSION 3.20)
project(mini02_proj_geos_geojson LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PROJ REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)

# まずはGEOSのCMake Config（Ubuntuのlibgeos-devが入れてくれることが多い）を優先
find_package(GEOS CONFIG QUIET)

if (GEOS_FOUND)
  set(GEOS_TARGET GEOS::geos_c)  # C API（安定ABI）
else()
  # 無ければpkg-configにフォールバック
  find_package(PkgConfig REQUIRED)
  # 注意: pkg-configのモジュール名は geos（geos_c ではない）:contentReference[oaicite:2]{index=2}
  pkg_check_modules(GEOS REQUIRED IMPORTED_TARGET geos)
  set(GEOS_TARGET PkgConfig::GEOS)
endif()

add_executable(001 main.cpp)
target_link_libraries(001 PRIVATE
  PROJ::proj
  ${GEOS_TARGET}
  nlohmann_json::nlohmann_json
)
