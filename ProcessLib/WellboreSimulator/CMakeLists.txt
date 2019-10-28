append_source_files(SOURCES)

add_library(WellboreSimulator ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS WellboreSimulator LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(WellboreSimulator PUBLIC ProcessLib PRIVATE ParameterLib)

include(Tests.cmake)
