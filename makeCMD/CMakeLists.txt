cmake_minimum_required(VERSION 2.8)

file(GLOB MAKECMD_FILES *.cpp *.hpp)

add_executable(makeCMD ${MAKECMD_FILES})

target_link_libraries(makeCMD base9 m ${MPI_LIBRARIES})

install(TARGETS makeCMD
        DESTINATION bin)
