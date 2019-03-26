# test

find_package(PythonInterp REQUIRED)

execute_process(
  COMMAND
  ${PYTHON_EXECUTABLE} "-c" "import numpy; import scipy"
  RESULT_VARIABLE _numpy_status
  OUTPUT_VARIABLE _numpy_location
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE  
  )

set(SOURCE_DAVIDSON array_utils davidson lapack_wrapper numeric_kinds)
foreach(PROG ${SOURCE_DAVIDSON})
  list(APPEND NUMPY_TEST_SRC ${CMAKE_SOURCE_DIR}/src/${PROG}.f90)
endforeach(PROG)
# set(NUMPY_TEST_SRC ${CMAKE_SOURCE_DIR}/src/davidson.f90 ${CMAKE_SOURCE_DIR}/src/lapack_wrapper.f90 ${CMAKE_SOURCE_DIR}/src/numeric_kinds.f90)

set(test_cases test_call_lapack test_dense_numpy test_free_numpy test_dense_properties
  test_free_properties test_reorder)

foreach(PROG ${test_cases})
  add_executable(${PROG}
    ${PROG}.f90
    test_utils.f90
    ${NUMPY_TEST_SRC}
    )

  target_compile_options(${PROG}
    PRIVATE
    ${OpenMP_Fortran_FLAGS}
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortran_FLAGS_DEBUG}>"
    )
  
  target_link_libraries(${PROG} PRIVATE ${LINEAR_ALGEBRA} ${OpenMP_Fortran_FLAGS})
endforeach(PROG)

MESSAGE(STATUS "Working directory: " ${CMAKE_SOURCE_DIR}/bin)

add_test(
  NAME test_dense_properties
  COMMAND $<TARGET_FILE:test_dense_properties>
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
  )

add_test(
  NAME test_free_properties
  COMMAND $<TARGET_FILE:test_free_properties>
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
  )

add_test(
  NAME test_reorder
  COMMAND $<TARGET_FILE:test_reorder>
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin  
  )

# Tests requiring python, numpy and scipy
if(${_numpy_status} EQUAL 0)
  add_test(
    NAME numpy_test
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_davidson.py -d $<TARGET_FILE:test_dense_numpy>
    -f $<TARGET_FILE:test_free_numpy>
    )
  add_test(
    NAME lapack_test
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_lapack.py -i $<TARGET_FILE:test_call_lapack>
    )
  
endif()  
