cmake_minimum_required(VERSION 2.8)

file(GLOB SIMCLUSTER_FILES *.cpp *.hpp)

add_executable(simCluster ${SIMCLUSTER_FILES})

target_link_libraries(simCluster base9 m ${MPI_LIBRARIES})

install(TARGETS simCluster
        DESTINATION bin)
