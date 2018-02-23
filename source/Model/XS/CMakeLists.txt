# ----------------------------------------------------------
# Define header dir and includes 
# ----------------------------------------------------------

# include directories 
SET( CLASS_HEADERS_DIR_XS ${CMAKE_CURRENT_SOURCE_DIR}  CACHE INTERNAL "")

LIST( APPEND CLASS_SRC_FILE_XS "${CMAKE_CURRENT_SOURCE_DIR}/XSM_CLOSEST.cxx")
LIST( APPEND CLASS_SRC_FILE_XS "${CMAKE_CURRENT_SOURCE_DIR}/XSM_MLP.cxx")

SET( CLASS_SRC_FILE_XS ${CLASS_SRC_FILE_XS} PARENT_SCOPE)

message("------- CLASS_SRC_FILE VARIABLE -------")
MESSAGE( STATUS "CLASS_SRC_FILE_XS=" ${CLASS_SRC_FILE_XS})
