cmake_minimum_required(VERSION 2.8)

file(GLOB SIMCLUSTER_FILES *.c)

add_executable(simCluster ${SIMCLUSTER_FILES})

target_link_libraries(simCluster base9 m)

install(TARGETS simCluster
        DESTINATION bin)
