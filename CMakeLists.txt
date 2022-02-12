cmake_minimum_required(VERSION 2.8.10)

project(Base9 CXX C)

list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" )

set(CMAKE_REQUIRED_INCLUDES "/usr/local/include")

include(CheckIncludeFiles)

check_include_files (math.h HAS_MATH)
check_include_files ("gsl/gsl_blas.h;gsl/gsl_math.h;gsl/gsl_eigen.h;gsl/gsl_statistics.h;gsl/gsl_linalg.h;gsl/gsl_spline.h;gsl/gsl_errno.h" HAS_GSL_BLAS)
check_include_files (getopt.h HAS_GETOPT)
check_include_files (unistd.h HAS_UNISTD)

include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES "m")
check_function_exists (exp exp HAS_EXP)
check_function_exists (exp10 HAS_EXP10)
check_function_exists (log HAS_LOG)
check_function_exists (log10 HAS_LOG10)
check_function_exists (pow HAS_POW)

######################
# Set project verion #
######################
include(GetGitRevisionDescription)
#set (Base9_VERSION )
git_describe(Base9_VERSION)

option (USE_TESTS "Enable testing" OFF)

if (USE_TESTS)
  enable_testing()
endif()

option (USE_TCMALLOC "Enable linking to tcmalloc" OFF)
option (BUILD_NONINTERACTIVE "Builds non-interactive applications for use on e.g., supercomputing clusters" OFF)

##########################
# Set up the config file #
##########################
configure_file (
  "${PROJECT_SOURCE_DIR}/Base9Config.h.in"
  "${PROJECT_BINARY_DIR}/Base9Config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES PROFILE)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
    "Reset the configurations"
    FORCE)
endif()

set(CMAKE_C_FLAGS                "${CMAKE_C_FLAGS} -march=native -fstrict-aliasing -W -Wall -Wno-missing-braces -pedantic -pthread -DENABLE_SSE2 -msse3" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG          "-O0 -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_PROFILE        "-O0 -g -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_OPROFILE       "-O3 -g -DNDEBUG" CACHE STRING "" FORCE)

set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -march=native -fstrict-aliasing -std=c++11 -W -Wall -Wno-missing-braces -pedantic -pthread -fno-rtti -DENABLE_SSE2 -msse3" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_PROFILE        "-O0 -g -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_OPROFILE       "-O3 -g -DNDEBUG" CACHE STRING "" FORCE)

set(CMAKE_EXE_LINKER_FLAGS "-pthread -ldl" CACHE STRING "" FORCE)

set(CMAKE_EXE_LINKER_FLAGS_PROFILE    "-L/usr/local/lib -lprofiler")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-L/usr/local/lib -lprofiler")

if (USE_TCMALLOC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltcmalloc")
endif()

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
  CMAKE_C_FLAGS_PROFILE
  CMAKE_C_FLAGS_OPROFILE
  CMAKE_CXX_FLAGS_PROFILE
  CMAKE_CXX_FLAGS_OPROFILE
  CMAKE_EXE_LINKER_FLAGS_PROFILE
  CMAKE_SHARED_LINKER_FLAGS_PROFILE)

#####################################
# Project specific code starts here #
#####################################
find_library(GSL_LOC
  gsl
  ENV LD_LIBRARY_PATH)

find_library(GSLCBLAS_LOC
  gslcblas
  ENV LD_LIBRARY_PATH)

MARK_AS_ADVANCED(
  GSL_LOC
  GSLCBLAS_LOC)

# Additional include and link directories
include_directories(base9)
include_directories(yaml-cpp/include)

# Libraries

option(YAML_CPP_BUILD_TESTS "Enable testing" OFF)
option(YAML_CPP_BUILD_TOOLS "Enable parse tools" OFF)
option(YAML_CPP_BUILD_CONTRIB "Enable contrib stuff in library" OFF)
add_subdirectory(yaml-cpp)

add_subdirectory(base9)

# Applications
add_subdirectory(singlePopMcmc)
add_subdirectory(makeCMD)
#add_subdirectory(makeIsochrone)
add_subdirectory(sampleMass)
add_subdirectory(sampleWDMass)
add_subdirectory(scatterCluster)
add_subdirectory(simCluster)
add_subdirectory(multiPopMcmc)
add_subdirectory(hierarchicalWDs)
