# ----------------------------------------------------------
# Modify template and build Decay.idx
# ----------------------------------------------------------
message(STATUS "COMMAND : sed -e 's:PATHTOBASE:${CMAKE_CURRENT_SOURCE_DIR}:g' Decay.in > Decay.idx")
execute_process(COMMAND bash "-c" "sed -e 's:PATHTOBASE:${CMAKE_CURRENT_SOURCE_DIR}:g' ${CMAKE_CURRENT_SOURCE_DIR}/Decay.in > ${CMAKE_CURRENT_SOURCE_DIR}/Decay.idx")
