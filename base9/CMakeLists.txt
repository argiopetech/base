cmake_minimum_required(VERSION 2.8)

file(GLOB BASE9_FILES *.cpp *.hpp) 

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_library(base9 STATIC ${BASE9_FILES})

target_link_libraries(base9 ${GSL_LOC} ${GSLCBLAS_LOC} yaml-cpp)

add_executable(base9_tests
    tests/base9_unittest.cpp
)
 
target_link_libraries(base9_tests gtest gtest_main)

add_test(
    NAME base9_tests
    COMMAND base9_tests --gtest_color=yes
)

