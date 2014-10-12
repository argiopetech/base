cmake_minimum_required(VERSION 2.8)

file(GLOB MAKECMD_FILES makeCMD.cpp) # *.cpp *.hpp)
file(GLOB MAKECMDS_FILES makeCMDs.cpp) # *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(makeCMD ${MAKECMD_FILES})
add_executable(makeCMDs ${MAKECMDS_FILES})

target_link_libraries(makeCMD base9 m yaml-cpp)
target_link_libraries(makeCMDs base9 m yaml-cpp)

install(TARGETS makeCMD makeCMDs
        DESTINATION bin)
