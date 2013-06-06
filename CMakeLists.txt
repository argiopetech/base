cmake_minimum_required(VERSION 2.8)

project(Base9 C CXX)

if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES PROFILE)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
    "Reset the configurations"
    FORCE)
endif()

find_package(MPI)

if (NOT (MPI_C_FOUND AND MPI_CXX_FOUND))
    message( STATUS ${MPI_CXX_COMPILER})
    message( FATAL_ERROR "It looks like MPI isn't installed..." )
endif()

set(CMAKE_C_FLAGS                  "-march=native -ffast-math -fstrict-aliasing -Wall ${MPI_C_COMPILE_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG            "-O0 -g")
set(CMAKE_C_FLAGS_MINSIZEREL       "-Os -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE          "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
set(CMAKE_C_FLAGS_PROFILE          "-O3 -g -DNDEBUG")

set(CMAKE_CXX_FLAGS                "-ffast-math -march=native -fstrict-aliasing -std=c++11 -Wall ${MPI_CXX_COMPILE_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS_PROFILE        "-O3 -g -DNDEBUG")

set(CMAKE_EXE_LINKER_FLAGS "${MPI_LINK_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_PROFILE    "-L/usr/local/lib -lprofiler")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-L/usr/local/lib -lprofiler")

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
#    if (NOT (GCC_VERSION VERSION_GREATER 4.6 OR GCC_VERSION VERSION_EQUAL 4.6))
    if (NOT (GCC_VERSION VERION_EQUAL 4.6))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.6 or Clang 3.2 (recommended). If Clang is installed, try 'CC=`which clang` CXX=`which clang++` cmake ..")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__extern_always_inline=inline")# -stdlib=libc++ -D__extern_always_inline=inline")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

MARK_AS_ADVANCED(
  CMAKE_CXX_FLAGS_PROFILE
  CMAKE_C_FLAGS_PROFILE
  CMAKE_EXE_LINKER_FLAGS_PROFILE
  CMAKE_SHARED_LINKER_FLAGS_PROFILE)

#####################################
# Project specific code starts here #
#####################################

find_library(GSL_LOC
             gsl)

find_library(GSLCBLAS_LOC
	     gslcblas)

MARK_AS_ADVANCED(
  GSL_LOC
  GSLCBLAS_LOC)

# Additional include and link directories
include_directories(base9)

# Libraries
add_subdirectory(base9)

# Applications
add_subdirectory(mpiMcmc)
add_subdirectory(makeCMD)
add_subdirectory(makeIsochrone)
add_subdirectory(sampleWDMass)
add_subdirectory(scatterCluster)
add_subdirectory(simCluster)
