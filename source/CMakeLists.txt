# ----------------------------------------------------------
# In include
# ----------------------------------------------------------

message("")
message( STATUS "in source/include ")
message("")
add_subdirectory(include)

# ----------------------------------------------------------
# In model
# ----------------------------------------------------------

message("")
message( STATUS "in source/Model ")
message("")
add_subdirectory(Model/Equivalence)
add_subdirectory(Model/Irradiation)
add_subdirectory(Model/XS)

# ----------------------------------------------------------
# In External
# ----------------------------------------------------------

message("")
message( STATUS "in source/external ")
message("")
add_subdirectory(external)

# ----------------------------------------------------------
# Includes
# ----------------------------------------------------------

INCLUDE_DIRECTORIES( ${CLASS_HEADERS_DIR_EQUIVALENCE} ${CLASS_HEADERS_DIR_IRRADIATION} ${CLASS_HEADERS_DIR_XS} ${CLASS_HEADERS_DIR_INCLUDE} ${CLASS_HEADERS_DIR_EXTERNAL})

# ----------------------------------------------------------
# In Src
# ----------------------------------------------------------

message("")
message( STATUS "in source/src ")
message("")
add_subdirectory(src)

# ----------------------------------------------------------
# Compilation
# ----------------------------------------------------------

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("inc_dirs = ${inc_dirs}")

message("")
message( STATUS "Build Libraries CLASSpkg")
message("")
ADD_LIBRARY(CLASSpkg SHARED ${CLASS_SRC_FILE} ${CLASS_HEADERS_DIR_EQUIVALENCE} ${CLASS_HEADERS_DIR_IRRADIATION} ${CLASS_HEADERS_DIR_XS})
target_link_libraries(CLASSpkg CLASSpkg_root)
