cmake_minimum_required(VERSION 3.10)
project(tphase)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Detecting C++ default compiler
if (NOT UNIX)
    message(FATAL_ERROR "Non Unix system not supported yet")
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Default C++ compiler: GCC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message(STATUS "Default C++ compiler: LLVM Clang")
else ()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} currently not supported by CMake build")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

# Set source and include directories
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Use GLOB to collect source and header files
file(GLOB SOURCE_FILES "${SOURCE_DIR}/*.cpp")
file(GLOB INCLUDE_FILES "${INCLUDE_DIR}/*.h")

# Include directories
# 更明确的包含目录设置
include_directories(
    ${INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/htslib
)

# HTSLIB configuration
set(HTSLIB ${CMAKE_CURRENT_SOURCE_DIR}/lib)
link_directories(${HTSLIB})

# Add executable and link libraries
add_executable(tphase ${SOURCE_FILES} ${INCLUDE_FILES} main.cpp)
target_link_libraries(tphase PUBLIC z hts)
message(STATUS "Source files: ${SOURCE_FILES}")
message(STATUS "Include files: ${INCLUDE_FILES}")
message(STATUS "Include directories: ${INCLUDE_DIR}")