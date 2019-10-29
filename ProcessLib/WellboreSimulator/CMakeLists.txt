append_source_files(SOURCES)

ogs_add_library(WellboreSimulator ${SOURCES})
target_link_libraries(WellboreSimulator PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
