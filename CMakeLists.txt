cmake_minimum_required(VERSION 3.15)
project(Matrix)

set(CMAKE_CXX_STANDARD 14)
add_subdirectory(googletest-master)

set(MAIN_SOURCES main.cpp)
include_directories(lib/include)
include_directories(googletest-master/googletest/include)
add_executable(Matrix ${MAIN_SOURCES})
target_link_libraries(Matrix gtest gtest_main)