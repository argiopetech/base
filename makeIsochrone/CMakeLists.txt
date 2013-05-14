cmake_minimum_required(VERSION 2.8)

file(GLOB MAKEISOCHRONE_FILES *.c)

add_executable(makeIsochrone ${MAKEISOCHRONE_FILES})

target_link_libraries(makeIsochrone base9 m)

install(TARGETS makeIsochrone
        DESTINATION bin)
