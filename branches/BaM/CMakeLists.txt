PROJECT(Class)
#set the project name

SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)
#set the cmake module find path


# Find ROOT
# root is required

FIND_PACKAGE(ROOT REQUIRED COMPONENTS Gui Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread pthread Wl rpath m dl)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
ADD_DEFINITIONS(${ROOT_CXX_FLAGS})
SET(CYCLASS_INCLUDE_DIRS $ROOT_INCLUDE_DIRS )
set(LIBS ${LIBS} ${ROOT_LIBRARIES})
MESSAGE("--    ROOT Root: ${ROOTSYS}")
MESSAGE("--    ROOT Include directory: ${ROOT_INCLUDE_DIR}")
MESSAGE("--    ROOT Library directories: ${ROOT_LIBRARY_DIRS}")
MESSAGE("--    ROOT Libraries: ${ROOT_LIBRARIES}")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
export(PACKAGE CLASS)

MESSAGE( STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY :" ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
ADD_SUBDIRECTORY(source)

#i Please enter the commit message for your changes. Lines starting

