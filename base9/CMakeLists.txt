cmake_minimum_required(VERSION 2.8)

file(GLOB BASE9_FILES *.cpp *.hpp MsRgbModels/*.cpp MsRgbModels/*.hpp WdCoolingModels/*.cpp WdCoolingModels/*.hpp WdAtmosphereModels/*.cpp WdAtmosphereModels/*.hpp IO/*.cpp IO/*.hpp sleef/*.c sleef/*.h sqlite/sqlite3.c sqlite/sqlite3.h)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_library(base9 STATIC ${BASE9_FILES})

set_target_properties(base9 PROPERTIES POSITION_INDEPENDENT_CODE true)

install(TARGETS base9
        DESTINATION lib)
