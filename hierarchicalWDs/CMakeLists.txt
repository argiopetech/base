cmake_minimum_required(VERSION 2.8)

file(GLOB HIERARCHICALWDS_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(hierarchicalWDs ${HIERARCHICALWDS_FILES})

target_link_libraries(hierarchicalWDs base9 m yaml-cpp ${GSL_LOC} ${GSLCBLAS_LOC})

install(TARGETS hierarchicalWDs
        DESTINATION bin)
