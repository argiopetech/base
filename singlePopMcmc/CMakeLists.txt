cmake_minimum_required(VERSION 2.8)

file(GLOB MPIMCMC_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} "/usr/local/include" )

add_executable(singlePopMcmc ${MPIMCMC_FILES})

target_link_libraries(singlePopMcmc base9 m ${GSL_LOC} ${GSLCBLAS_LOC} yaml-cpp)

install(TARGETS singlePopMcmc
        DESTINATION bin)
