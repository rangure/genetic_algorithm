cmake_minimum_required(VERSION 3.10)

# enable_language(CXX)
# set the project name
project(Genetic_Algorithm)
MESSAGE(${CMAKE_C_COMPILER})
MESSAGE(${CMAKE_CXX_COMPILER})
find_package(MPI)
if(MPI_FOUND)
    add_compile_definitions(PARALLEL_ENABLE=1)
    MESSAGE(${MPI_C_LIBRARIES})
    MESSAGE(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
add_definitions(-DOMPI_SKIP_MPICXX)
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(SYSTEM ${MPI_INCLUDE_PATH})


# set the include path
include_directories(include)

# add the main code
add_subdirectory(src)

# add the tests
include(CTest)
enable_testing()
add_subdirectory(tests)


