cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

PROJECT(Class)
#set the project name

SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)
#set the cmake module find path

# check for and enable c++11 support
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
IF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
    MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
ENDIF()
# Find ROOT
# root is required
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -pthread -stdlib=libc++ -std=c++1z -m64  -stdlib=libc++ -std=gnu++11 -I/opt/local/libexec/root6/include/root")


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
FIND_PACKAGE(ROOT REQUIRED COMPONENTS Gui Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread pthread m dl TMVA)
MESSAGE("--    ROOT Root: ${ROOTSYS}")
MESSAGE("--    ROOT Include directory: ${ROOT_INCLUDE_DIR}")
MESSAGE("--    ROOT Library directories: ${ROOT_LIBRARY_DIRS}")
MESSAGE("--    ROOT Libraries: ${ROOT_LIBRARIES}")
MESSAGE("--    ROOT Root Use File : ${ROOT_USE_FILE}")
INCLUDE(${ROOT_USE_FILE})


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
ADD_DEFINITIONS(${ROOT_CXX_FLAGS})
SET(CYCLASS_INCLUDE_DIRS $ROOT_INCLUDE_DIRS )
set(LIBS ${LIBS} ${ROOT_LIBRARIES})

SET(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
export(PACKAGE CLASS)

MESSAGE( STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY :" ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
ADD_SUBDIRECTORY(source)


