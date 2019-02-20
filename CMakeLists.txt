# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(paf2sxs)
# We build using c++11
set(CMAKE_CXX_STANDARD 11)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3")

# We use OpenMP for parallelism
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# gzip file reader, header only
ExternalProject_Add(gzipreader
  GIT_REPOSITORY "https://github.com/gatoravi/gzip_reader.git"
  GIT_TAG "0ef26c0399e926087f9d6c4a56067a7bf1fc4f5e"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(gzipreader SOURCE_DIR)
set(gzipreader_INCLUDE "${SOURCE_DIR}")

set(CMAKE_BUILD_TYPE Release) 

# set up our target executable and specify its dependencies and includes
add_executable(paf2sxs
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  ${CMAKE_SOURCE_DIR}/src/paf.cpp
  ${CMAKE_SOURCE_DIR}/src/cigar.cpp
  )
add_dependencies(paf2sxs tayweeargs)
add_dependencies(paf2sxs gzipreader)
target_include_directories(paf2sxs PUBLIC
  "${tayweeargs_INCLUDE}"
  "${gzipreader_INCLUDE}")
target_link_libraries(paf2sxs
  z)
set(CMAKE_EXE_LINKER_FLAGS "-static")
