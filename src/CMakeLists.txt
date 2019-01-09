
# create binary
set(SOURCES main.f90 numeric_kinds.f90)
message (STATUS "SOURCES: " ${SOURCES})
add_executable(main ${SOURCES})

# linking
set( CMAKE_VERBOSE_MAKEFILE on )
target_link_libraries(main)
