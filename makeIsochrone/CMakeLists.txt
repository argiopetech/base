cmake_minimum_required(VERSION 2.8)

file(GLOB MAKEISOCHRONE_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(makeIsochrone ${MAKEISOCHRONE_FILES})

target_link_libraries(makeIsochrone base9 m ${MPI_LIBRARIES})

install(TARGETS makeIsochrone
        DESTINATION bin)
