cmake_minimum_required(VERSION 3.20)
project(gdal_inspect_min LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GDAL CONFIG)
if(NOT GDAL_FOUND)
    find_package(GDAL REQUIRED)
endif()

add_executable(004 main.cpp)
target_link_libraries(004 PRIVATE GDAL::GDAL)
