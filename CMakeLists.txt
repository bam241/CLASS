# ----------------------------------------------------------
# Project name
# ----------------------------------------------------------

project(Class CXX)

# ----------------------------------------------------------
# Minimum Cmake version
# ----------------------------------------------------------

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# ----------------------------------------------------------
# Verbose mode for DEBUG only
# ----------------------------------------------------------

set(CMAKE_VERBOSE_MAKEFILE off)

# ----------------------------------------------------------
# Test CXX compiler and C++11 comptability
# ----------------------------------------------------------

INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
IF(COMPILER_SUPPORTS_CXX11)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
  MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
ENDIF()

# ----------------------------------------------------------
# ROOT integration
# ----------------------------------------------------------

# defining ROOT location
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# Defining needed classes
find_package(ROOT REQUIRED COMPONENTS TMVA)

message("------- ROOT PATH DETECTION -------")

if (ROOT_FOUND)
  message(STATUS "ROOT has been found")
else ()
  message(FATAL_ERROR "ROOT Not found")   
endif()

message("-- Include path   : ${ROOT_INCLUDE_DIRS}")
message("-- Libraries list : ${ROOT_LIBRARIES}")
message("-- Libraries DIR  : ${ROOT_LIBRARY_DIR}")
message("-- Definitions    : ${ROOT_DEFINITIONS}")
message("-- CXX Flags      : ${ROOT_CXX_FLAGS}")
message("-- CC Flags       : ${ROOT_CC_FLAGS}")
message("------------------------------------")

# Use ROOT file

INCLUDE(${ROOT_USE_FILE})

# Add directories
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

# Build Decay.idx File
set(decay_dir "${CMAKE_CURRENT_SOURCE_DIR}/DATA_BASES/DECAY/ALL")
execute_process(COMMAND bash "-c" "sed -e 's:PATHTOBASE:${decay_dir}:g' ${decay_dir}/Decay.in > ${decay_dir}/Decay.idx")

# Set CLASS variable
# SET(CYCLASS_INCLUDE_DIRS $ROOT_INCLUDE_DIRS )
SET(LIBS ${LIBS} ${ROOT_LIBRARIES})

# Add compilation flags
# add_definitions(${ROOT_CXX_FLAGS})

# ----------------------------------------------------------
# Variables definitions
# ----------------------------------------------------------

set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)

# print all
MESSAGE( STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY :" ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# ----------------------------------------------------------
# In source
# ----------------------------------------------------------

message("------------------------------------")
message( STATUS "in source ")
add_subdirectory(source)
message("------------------------------------")
