cmake_minimum_required(VERSION 2.8)

file(GLOB MAKECMD_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(makeCMD ${MAKECMD_FILES})

target_link_libraries(makeCMD base9 m yaml-cpp)

install(TARGETS makeCMD
        DESTINATION bin)
