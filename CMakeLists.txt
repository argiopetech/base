cmake_minimum_required(VERSION 2.8)

# Initialize CXXFLAGS.
#set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
#set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
#set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

## Compiler-specific C++11 activation.
#if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
#    execute_process(
#        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
#    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
#        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
#    endif ()
#elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
#else ()
#    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
#endif ()

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_C_FLAGS "-O3 -ffast-math -fstrict-aliasing" CACHE STRING 
    "Flags used by the compiler during all build types."
    FORCE)

set(CMAKE_CXX_COMPILER "mpiCC")
set(CMAKE_CXX_FLAGS "-O3 -ffast-math -fstrict-aliasing" CACHE STRING 
    "Flags used by the compiler during all build types."
    FORCE)

if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES Profile)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
    "Reset the configurations"
    FORCE)
endif() 

set(CMAKE_C_FLAGS_DEBUG "-O3 -ffast-math -fstrict-aliasing -g -Wall" CACHE STRING
#set(CMAKE_C_FLAGS_DEBUG "-g -Wall" CACHE STRING
  "Flags used by the C compiler during debug builds."
  FORCE)

set(CMAKE_CXX_FLAGS_DEBUG "-O3 -ffast-math -fstrict-aliasing -g -Wall" CACHE STRING
  "Flags used by the C compiler during debug builds."
  FORCE)

set(CMAKE_C_FLAGS_PROFILE "-O3 -ffast-math -fstrict-aliasing -g" CACHE STRING
  "Flags used by the C compiler during profile builds."
  FORCE)
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-L/usr/local/lib -lprofiler" CACHE STRING
  "Flags used when linking executables during profile builds."
  FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-L/usr/local/lib -lprofiler" CACHE STRING
  "Flags used when linking libraries during profile builds."
  FORCE)

MARK_AS_ADVANCED(
  CMAKE_CXX_FLAGS_PROFILE
  CMAKE_C_FLAGS_PROFILE
  CMAKE_EXE_LINKER_FLAGS_PROFILE
  CMAKE_SHARED_LINKER_FLAGS_PROFILE)

#####################################
# Project specific code starts here #
#####################################

project(Base9 C CXX)

find_library(GSL_LOC
             gsl)

find_library(GSLCBLAS_LOC
	     gslcblas)

MARK_AS_ADVANCED(
  GSL_LOC
  GSLCBLAS_LOC)

# Additional include and link directories
include_directories(base9)

# Additional compiler flags
add_definitions(-Wall -Wno-unused)

# Libraries
add_subdirectory(base9)

# Applications
add_subdirectory(mpiMcmc)
add_subdirectory(makeCMD)
add_subdirectory(makeIsochrone)
add_subdirectory(sampleWDMass)
add_subdirectory(scatterCluster)
add_subdirectory(simCluster)
