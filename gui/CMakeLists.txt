INCLUDE_DIRECTORIES( ${CLASS_HEADERS_DIR} )
add_executable (CLASSgui CLASSGui.cxx CLASSPlotElement.cxx CLASSWin.cxx CLASSRead.cxx)
target_link_libraries(CLASSgui CLASSpkg)

