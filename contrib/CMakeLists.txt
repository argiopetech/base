cmake_minimum_required(VERSION 2.8)

include(ExternalProject)
ExternalProject_Add(
        yamlcpp
        URL https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz
        CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=ON
        )
