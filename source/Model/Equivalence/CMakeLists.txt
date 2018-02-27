# ----------------------------------------------------------
# Define header dir and includes 
# ----------------------------------------------------------

# include directories 
SET( CLASS_HEADERS_DIR_EQUIVALENCE ${CMAKE_CURRENT_SOURCE_DIR}  CACHE INTERNAL "")

LIST( APPEND CLASS_SRC_FILE_EQUIVALENCE "${CMAKE_CURRENT_SOURCE_DIR}/EQ_OneParameter.cxx")

SET( CLASS_SRC_FILE_IRRADIATION ${CLASS_SRC_FILE_EQUIVALENCE} PARENT_SCOPE)

MESSAGE( STATUS "CLASS_SRC_FILE_EQUIVALENCE=" ${CLASS_SRC_FILE_EQUIVALENCE})
