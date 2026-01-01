cmake_minimum_required(VERSION 3.20)
project(008_geos_vs_boost_overlay LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(nlohmann_json REQUIRED CONFIG)

# Boost.Geometry はヘッダオンリー（リンク不要）だが、ヘッダ探索のため find_package だけしておく
find_package(Boost REQUIRED)

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

add_executable(008 main.cpp)
target_link_libraries(008 PRIVATE
  ${GEOS_TARGET}
  nlohmann_json::nlohmann_json
)
target_include_directories(008 PRIVATE ${Boost_INCLUDE_DIRS})
