cmake_minimum_required(VERSION 2.8)

file(GLOB MPIMCMC_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(mpiMcmc ${MPIMCMC_FILES})

target_link_libraries(mpiMcmc base9 m ${MPI_LIBRARIES})

install(TARGETS mpiMcmc
        DESTINATION bin)

#add_executable(mpiMcmc_tests
#    tests/mpiFuncs_unittest.cpp
#    mpiFuncs.cpp
#)

#target_link_libraries(mpiMcmc_tests gtest gtest_main base9)

#add_test(
#    NAME mpiMcmc_tests
#    COMMAND mpiMcmc_tests --gtest_color=yes
#)