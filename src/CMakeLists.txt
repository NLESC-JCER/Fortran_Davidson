
# create binary
set(UTILS array_utils.f90 numeric_kinds.f90 lapack_wrapper.f90)
set(SOURCES main.f90 davidson.f90 benchmark.f90 ${UTILS})
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

target_link_libraries(main PRIVATE ${LINEAR_ALGEBRA} ${OpenMP_Fortran_FLAGS})

# Benchmark free matrix
add_executable(benchmark_free benchmark_free.f90 davidson.f90 ${UTILS})

target_compile_options(main
  PRIVATE
  ${OpenMP_Fortran_FLAGS}
   ${Fortran_FLAGS}
   )

target_link_libraries(benchmark_free PRIVATE ${LINEAR_ALGEBRA} ${OpenMP_Fortran_FLAGS}) 


if(ENABLE_TEST)
  MESSAGE(STATUS "enable test")
  enable_testing()
  add_subdirectory(tests)
endif(ENABLE_TEST)
