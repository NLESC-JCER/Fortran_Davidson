
# create binary
set(SOURCES main.f90 numeric_kinds.f90 davidson.f90 sort.f90)
message (STATUS "SOURCES: " ${SOURCES})
add_executable(main ${SOURCES})

# linking
set( CMAKE_VERBOSE_MAKEFILE on )

target_compile_options(main
  PRIVATE
    "-std=f2008"
  )

target_link_libraries(main ${LINEAR_ALGEBRA})
