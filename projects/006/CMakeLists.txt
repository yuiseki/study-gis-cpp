cmake_minimum_required(VERSION 3.20)
project(gdal_reproject_rewrite LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GDAL CONFIG REQUIRED)
find_package(PROJ REQUIRED CONFIG)

add_executable(006 main.cpp)
target_link_libraries(006 PRIVATE GDAL::GDAL PROJ::proj)
