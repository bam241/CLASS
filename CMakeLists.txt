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

# => Which define Following Variables : 
# 
# ROOT_INCLUDE_DIRS       PATH    include directories for the ROOT installation
# ROOT_LIBRARY_DIR        PATH    libraries directory of ROOT
# ROOT_LIBRARIES          LIST    Libraries to link against. The actual list of libraries is composed using the COMPONENTS listed in the find_package(...)
# ROOT_DEFINITIONS        STRING  Compile definitions needed to use ROOT
# ROOT_CXX_FLAGS          STRING  C++ compiler flags used to build ROOT
# ROOT_CC_FLAGS           STRING  C compiler flags used to build ROOT
# ROOT_<library>_LIBRARY  PATH    Full path for each of the ROOT libraries listed in COMPONENTS
# ROOT_<command>_CMD      PATH    Full path for each ROOT executable (rootcling, root, hadd, etc.)
# ROOT_<option>_FOUND     BOOL    True for each enabled build option (e.g. cocoa, python, xrootd, etc.)
# ROOT_FOUND              BOOL    True if the ROOT package has been found
# ROOT_USE_FILE           PATH    Path to a CMake module which makes use of the previous variables and loads modules with useful macros or functions such as ROOT_GENERATE_DICTIONARY
# 

message("")
message("------- ROOT PATH DETECTION -------")
message("")

if (ROOT_FOUND)
  message(STATUS "ROOT has been found")
  message("")
else ()
  message(FATAL_ERROR "ROOT Not found")   
  message("")
endif()

message("-- Include path   : ${ROOT_INCLUDE_DIRS}")
message("")
message("-- Libraries list : ${ROOT_LIBRARIES}")
message("")
message("-- Libraries DIR  : ${ROOT_LIBRARY_DIR}")
message("")
message("-- Definitions    : ${ROOT_DEFINITIONS}")
message("")
message("-- CXX Flags      : ${ROOT_CXX_FLAGS}")
message("")
message("-- CC Flags       : ${ROOT_CC_FLAGS}")
message("")

message("------------------------------------")
message("")

# Use ROOT file
INCLUDE(${ROOT_USE_FILE})

# Add directories
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

# Add compilation flags
# add_definitions(${ROOT_CXX_FLAGS})

# ----------------------------------------------------------
# Variables definitions
# ----------------------------------------------------------

SET(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)

# print all
MESSAGE( STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY :" ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# ----------------------------------------------------------
# In source
# ----------------------------------------------------------

message("")
message("------------------------------------")
message( STATUS "in source ")
add_subdirectory(source)
message("------------------------------------")
message("")
