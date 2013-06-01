cmake_minimum_required(VERSION 2.8)

file(GLOB BASE9_FILES *.c *.cpp) 

add_library(base9 STATIC ${BASE9_FILES})

target_link_libraries(base9 ${GSL_LOC} ${GSLCBLAS_LOC} yaml-cpp)
