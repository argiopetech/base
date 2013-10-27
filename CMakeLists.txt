cmake_minimum_required(VERSION 2.8.10)

project(Base9 CXX C)

list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" )

include(CheckIncludeFiles)

set(CMAKE_REQUIRED_LIBRARIES "m")
check_include_files (math.h HAS_MATH)

include(CheckFunctionExists)
check_function_exists (exp HAS_EXP)
check_function_exists (exp10 HAS_EXP10)
check_function_exists (log HAS_LOG)
check_function_exists (log10 HAS_LOG10)
check_function_exists (pow HAS_POW)

######################
# Set project verion #
######################
set (Base9_VERSION_MAJOR 9)
set (Base9_VERSION_MINOR 3)
set (Base9_VERSION_PATCH 1)

option (USE_TESTS "Enable testing" OFF) 

if (USE_TESTS)
  enable_testing()
endif()

##########################
# Set up the config file #
##########################
configure_file (
  "${PROJECT_SOURCE_DIR}/Base9Config.h.in"
  "${PROJECT_BINARY_DIR}/Base9Config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

include(GetProject)

if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES PROFILE)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
    "Reset the configurations"
    FORCE)
endif()

find_package(MPI)

if (NOT (MPI_CXX_FOUND))
  message( STATUS ${MPI_CXX_COMPILER})
  message( FATAL_ERROR "It looks like MPI isn't installed..." )
endif()

include_directories( ${MPI_CXX_INCLUDE_PATH} )

set(CMAKE_CXX_FLAGS                "-march=native -fstrict-aliasing -std=c++11 -W -Wall -Wno-missing-braces -pedantic ${MPI_CXX_COMPILE_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_PROFILE        "-O3 -g -DNDEBUG" CACHE STRING "" FORCE)

set(CMAKE_EXE_LINKER_FLAGS "${MPI_LINK_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_PROFILE    "-L/usr/local/lib -lprofiler")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-L/usr/local/lib -lprofiler")

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  execute_process( COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
  if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
    message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or Clang 3.2 (recommended). If Clang is installed, try `CC=`which clang` CXX=`which clang++` cmake ..")
  endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__extern_always_inline=inline")

  if (APPLE AND (${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  endif ()

else ()
  message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

MARK_AS_ADVANCED(
  CMAKE_CXX_FLAGS_PROFILE
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

# External projects
get_project( "YamlCpp" )

# Libraries
add_subdirectory(base9)

# Applications
add_subdirectory(mpiMcmc)
#add_subdirectory(makeCMD)
#add_subdirectory(makeIsochrone)
add_subdirectory(sampleWDMass)
#add_subdirectory(scatterCluster)
#add_subdirectory(simCluster)