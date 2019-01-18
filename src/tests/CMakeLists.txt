# test

find_package(PythonInterp REQUIRED)

execute_process(
  COMMAND
  ${PYTHON_EXECUTABLE} "-c" "import numpy"
  RESULT_VARIABLE _numpy_status
  OUTPUT_VARIABLE _numpy_location
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE  
  )

set(NUMPY_TEST_SRC ${CMAKE_SOURCE_DIR}/src/davidson.f90 ${CMAKE_SOURCE_DIR}/src/numeric_kinds.f90)

set(test_cases test_numpy test_properties)

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
  
  target_link_libraries(${PROG} ${LINEAR_ALGEBRA})
endforeach(PROG)

MESSAGE(STATUS "Working directory: " ${CMAKE_SOURCE_DIR}/bin)

add_test(
  NAME test_properties
  COMMAND $<TARGET_FILE:test_properties>
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
  )

if(${_numpy_status} EQUAL 0)
  add_test(
    NAME numpy_test
    COMMAND $<TARGET_FILE:test_numpy>
    # COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_davidson.py --executable $<TARGET_FILE:test_numpy>
    )
endif()  
