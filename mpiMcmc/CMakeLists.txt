cmake_minimum_required(VERSION 2.8)

file(GLOB MPIMCMC_FILES *.c)

add_executable(mpiMcmc ${MPIMCMC_FILES})

target_link_libraries(mpiMcmc base9 m ${MPI_LIBRARIES})

install(TARGETS mpiMcmc
        DESTINATION bin)
