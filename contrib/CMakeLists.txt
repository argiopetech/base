cmake_minimum_required(VERSION 2.8)

include(ExternalProject)
FIND_PACKAGE(Boost)

IF (NOT Boost_FOUND)
    ExternalProject_Add(
            Boost
            URL http://downloads.sourceforge.net/project/boost/boost/1.53.0/boost_1_53_0.tar.gz
            UPDATE_COMMAND "./bootstrap.sh"
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ./b2
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ./b2 install
            )

    set_property( TARGET Boost PROPERTY FOLDER "Externals")

    ExternalProject_Get_Property( Boost source_dir )
    ExternalProject_Get_Property( Boost binary_dir )

    set( Boost_INCLUDE_DIRS ${source_dir} )

    set( Boost_LIBRARIES ${binary_dir}/stage/lib )
ENDIF()

ExternalProject_Add(
        yamlcpp
        URL https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz
        CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=ON -DYAML_CPP_BUILD_TOOLS:BOOL=OFF
        )
