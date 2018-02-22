# ----------------------------------------------------------
# Define header dir and includes 
# ----------------------------------------------------------

# include directories 
SET( CLASS_HEADERS_DIR_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")



INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Equivalence:${CMAKE_CURRENT_SOURCE_DIR}/Irradiation:${CMAKE_CURRENT_SOURCE_DIR}/XS)
