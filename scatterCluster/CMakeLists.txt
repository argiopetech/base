cmake_minimum_required(VERSION 2.8)

file(GLOB SCATTERCLUSTER_FILES *.cpp *.hpp)

add_executable(scatterCluster ${SCATTERCLUSTER_FILES})

target_link_libraries(scatterCluster base9 m ${MPI_LIBRARIES})

install(TARGETS scatterCluster
        DESTINATION bin)
