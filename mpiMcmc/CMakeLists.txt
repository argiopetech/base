cmake_minimum_required(VERSION 2.8)

enable_testing()

file(GLOB MPIMCMC_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(mpiMcmc ${MPIMCMC_FILES})

target_link_libraries(mpiMcmc base9 m -ltcmalloc ${MPI_LIBRARIES})

install(TARGETS mpiMcmc
        DESTINATION bin)

add_executable(mpiMcmc_tests
    tests/mpiFuncs_unittest.cpp
    mpiFuncs.cpp
)

target_link_libraries(mpiMcmc_tests gtest gtest_main base9 ${MPI_LIBRARIES})

add_test(
    NAME mpiMcmc_tests
    COMMAND mpiMcmc_tests --gtest_color=yes
)