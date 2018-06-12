cmake_minimum_required(VERSION 2.8)

file(GLOB MULTIPOPMCMC_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(multiPopMcmc ${MULTIPOPMCMC_FILES})

target_link_libraries(multiPopMcmc base9 m ${GSL_LOC} ${GSLCBLAS_LOC} yaml-cpp)

install(TARGETS multiPopMcmc
        DESTINATION bin)
