cmake_minimum_required(VERSION 2.8)

file(GLOB MAKECMD_FILES *.c)

add_executable(makeCMD ${MAKECMD_FILES})

target_link_libraries(makeCMD base9 m)

install(TARGETS makeCMD
        DESTINATION bin)
