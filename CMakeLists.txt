# ----------------------------------------------------------
# Project name
# ----------------------------------------------------------

project(Class CXX)

# ----------------------------------------------------------
# Minimum Cmake version
# ----------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

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
FIND_PACKAGE(ROOT REQUIRED COMPONENTS TMVA Gui)

MESSAGE("------- ROOT PATH DETECTION -------")

IF (ROOT_FOUND)
  MESSAGE(STATUS "ROOT has been found")
ELSE ()
  MESSAGE(FATAL_ERROR "ROOT Not found")   
ENDIF()

MESSAGE("-- Include path   : ${ROOT_INCLUDE_DIRS}")
MESSAGE("-- Libraries list : ${ROOT_LIBRARIES}")
MESSAGE("-- Libraries DIR  : ${ROOT_LIBRARY_DIR}")
MESSAGE("-- Definitions    : ${ROOT_DEFINITIONS}")
MESSAGE("-- CXX Flags      : ${ROOT_CXX_FLAGS}")
MESSAGE("-- CC Flags       : ${ROOT_CC_FLAGS}")
MESSAGE("------------------------------------")

# Use ROOT file

INCLUDE(${ROOT_USE_FILE})

# Add directories
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

# ----------------------------------------------------------
# Variables definitions
# ----------------------------------------------------------

set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)

# print all
MESSAGE( STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY :" ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# ----------------------------------------------------------
# Add Directories
# ----------------------------------------------------------

# Build Decay.idx File
add_subdirectory(DATA_BASES/DECAY/ALL)

# Add source
ADD_SUBDIRECTORY(source)

# Add gui
ADD_SUBDIRECTORY(gui)

################################
# GTest
################################

enable_testing()
add_subdirectory(GTest)
