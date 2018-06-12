cmake_minimum_required(VERSION 2.8)

file(GLOB SOURCE_FILES *.cpp *.hpp)

# include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(base-runSql ${SOURCE_FILES})

target_link_libraries(base-runSql base9)

install(TARGETS base-runSql
        DESTINATION bin)
