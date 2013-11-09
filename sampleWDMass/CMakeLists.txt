cmake_minimum_required(VERSION 2.8)

file(GLOB SAMPLEWDMASS_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(sampleWDMass ${SAMPLEWDMASS_FILES})

target_link_libraries(sampleWDMass base9 m)

install(TARGETS sampleWDMass
        DESTINATION bin)
