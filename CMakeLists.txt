cmake_minimum_required(VERSION 3.15)

# Set the project name to your project name, my project isn't very descriptive
project(acm CXX)
include(cmake/StandardProjectSettings.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC_LIBRARY "${PROJECT_NAME}Library")
set(MAIN_BINARY "${PROJECT_NAME}")
# set(TEST_BINARY "test${PROJECT_NAME}")
# set(BENCH_BINARY "bench${PROJECT_NAME}")

# Project warning section, this library will be empty if USE_WARNINGS is OFF
add_library(project_warnings INTERFACE)
option(USE_WARNINGS "To turn on standard warnings" OFF)
if(USE_WARNINGS)
    # standard compiler warnings
    include(cmake/CompilerWarnings.cmake)
    set_project_warnings(project_warnings)
endif()

# creating the main executable and include the CMakeList file from src to include and link all the files
add_executable(${MAIN_BINARY} src/main.cpp)
add_library(${SRC_LIBRARY})
target_link_libraries(${SRC_LIBRARY} PUBLIC project_warnings)
include(src/CMakeLists.txt)

# Testing section
option(ENABLE_TESTING "Enable Test Builds" OFF)
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# Benchmark section
option(ENABLE_BENCHMARK "Enable Benchmark Builds" OFF)
