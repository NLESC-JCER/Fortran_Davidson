
# create binary
set(SOURCES main.f90 numeric_kinds.f90 davidson.f90 benchmark.f90)
message (STATUS "SOURCES: " ${SOURCES})
add_executable(main ${SOURCES})

# linking
set( CMAKE_VERBOSE_MAKEFILE on )


## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-std=f2008")
  list(APPEND Fortran_FLAGS_DEBUG -Wall -Wextra -Wconversion
    -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -fbounds-check
    -ffpe-trap=zero,overflow,underflow -finit-real=nan)

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS_DEBUG -check all -fpe0 -warn -traceback -debug extended)
endif()
  
target_compile_options(main
  PRIVATE
  ${OpenMP_Fortran_FLAGS}
   ${Fortran_FLAGS}
   "$<$<CONFIG:Debug>:${Fortran_FLAGS_DEBUG}>"
   )

target_link_libraries(main ${LINEAR_ALGEBRA})

if(ENABLE_TEST)
  MESSAGE(STATUS "enable test")
  enable_testing()
  add_subdirectory(tests)
endif(ENABLE_TEST)
