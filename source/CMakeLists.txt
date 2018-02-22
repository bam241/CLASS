# ----------------------------------------------------------
# In include
# ----------------------------------------------------------

message( STATUS "in source/include ")
add_subdirectory(include)

# ----------------------------------------------------------
# In model
# ----------------------------------------------------------

message( STATUS "in source/Model ")
add_subdirectory(Model/Equivalence)
add_subdirectory(Model/Irradiation)
add_subdirectory(Model/XS)

# ----------------------------------------------------------
# In External
# ----------------------------------------------------------

message( STATUS "in source/external ")
add_subdirectory(external)

# ----------------------------------------------------------
# Includes
# ----------------------------------------------------------

INCLUDE_DIRECTORIES( ${CLASS_HEADERS_DIR_EQUIVALENCE} ${CLASS_HEADERS_DIR_IRRADIATION} ${CLASS_HEADERS_DIR_XS} ${CLASS_HEADERS_DIR_INCLUDE} ${CLASS_HEADERS_DIR_EXTERNAL})

# ----------------------------------------------------------
# In Src
# ----------------------------------------------------------

message( STATUS "in source/src ")
add_subdirectory(src)

# ----------------------------------------------------------
# Compilation
# ----------------------------------------------------------

message( STATUS "Build Libraries CLASSpkg")
ADD_LIBRARY(CLASSpkg SHARED ${CLASS_SRC_FILE} ${CLASS_HEADERS_DIR_EQUIVALENCE} ${CLASS_HEADERS_DIR_IRRADIATION} ${CLASS_HEADERS_DIR_XS})
target_link_libraries(CLASSpkg CLASSpkg_root)
