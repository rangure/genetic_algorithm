project(tests)
list(APPEND Tests test_validity test_simulator test_GA) 

foreach(TEST IN LISTS Tests)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} PUBLIC geneticAlgorithm circuitSimulator)
    if(MPI_FOUND)
    target_link_libraries(${TEST} PRIVATE MPI::MPI_C)
    endif()
    target_include_directories(${TEST} PRIVATE ../includes)
    set_target_properties(${TEST} PROPERTIES
        CXX_STANDARD 17
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/bin")
    add_test(NAME ${TEST} COMMAND ${TEST})
endforeach()

add_test(NAME executable COMMAND "${CMAKE_BINARY_DIR}/bin/Circuit_Optimizer")