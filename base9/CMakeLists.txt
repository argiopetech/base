cmake_minimum_required(VERSION 2.8)

file(GLOB BASE9_FILES *.cpp *.hpp MsRgbModels/*.cpp MsRgbModels/*.hpp WdCoolingModels/*.cpp WdCoolingModels/*.hpp) 

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_library(base9 STATIC ${BASE9_FILES})

install(TARGETS base9
        DESTINATION lib)

if (USE_TESTS)
  add_executable(base9_tests
    tests/base9_unittest.cpp
  )
 
  target_link_libraries(base9_tests gtest gtest_main base9)

  add_test(
    NAME base9_tests
    COMMAND base9_tests --gtest_color=yes
  )
endif()
