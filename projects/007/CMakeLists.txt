cmake_minimum_required(VERSION 3.20)
project(007_distance_strategies LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PROJ REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)

# Boost.Geometry はヘッダオンリーなのでリンク不要（ただしヘッダが見つかる保証のため find_package だけ）
find_package(Boost REQUIRED)

# GeographicLib: まずCMake packageを試し、ダメなら include/library を直接探す
find_package(GeographicLib QUIET CONFIG)

if(TARGET GeographicLib::GeographicLib)
  set(GEOGRAPHICLIB_TARGET GeographicLib::GeographicLib)
else()
  find_path(GEOGRAPHICLIB_INCLUDE_DIR
    NAMES GeographicLib/Geodesic.hpp GeographicLib/Constants.hpp
  )
  find_library(GEOGRAPHICLIB_LIBRARY
    NAMES GeographicLib Geographic
  )
  if(NOT GEOGRAPHICLIB_INCLUDE_DIR OR NOT GEOGRAPHICLIB_LIBRARY)
    message(FATAL_ERROR "GeographicLib not found. Install GeographicLib dev package or set GeographicLib_DIR / CMAKE_PREFIX_PATH.")
  endif()

  add_library(geographiclib_unknown UNKNOWN IMPORTED)
  set_target_properties(geographiclib_unknown PROPERTIES
    IMPORTED_LOCATION "${GEOGRAPHICLIB_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GEOGRAPHICLIB_INCLUDE_DIR}"
  )
  set(GEOGRAPHICLIB_TARGET geographiclib_unknown)
endif()

add_executable(007 main.cpp)
target_link_libraries(007 PRIVATE
  PROJ::proj
  nlohmann_json::nlohmann_json
  ${GEOGRAPHICLIB_TARGET}
)

target_include_directories(007 PRIVATE ${Boost_INCLUDE_DIRS})
