#!/bin/sh

# Get the number of CPUs on the system
# This is not portable
if [ ! "$NCPUS" ]; then
    unamestr=`uname`
    if [ "$unamestr" = 'Linux' ]; then
        NCPUS=`grep -c ^processor /proc/cpuinfo`
    elif [ "$unamestr" = 'Darwin' ] || [ "$unamestr" = 'FreeBSD' ]; then
        NCPUS=`sysctl -n hw.ncpu`
    else
        NCPUS=1
    fi
fi

OWD=$PWD
BASE=`dirname "${0}"`

cd $BASE

if [ -e ".git" ]; then
    # Setup submodules
    git submodule init
    git submodule update
else
    if [ ! -e "yaml-cpp/CMakeLists.txt" ]; then
        # Manually clone yaml-cpp, removing the directory first
        rm -r yaml-cpp
        git clone https://github.com/argiopetech/yaml-cpp.git --depth 1
    fi
fi

cp cmake/yaml-cpp-CMakeLists.txt yaml-cpp/CMakeLists.txt

cd ./BUILD
if [ "$PREFIX" ]; then
    cmake -DCMAKE_BUILD_TYPE="RELEASE" -DCMAKE_INSTALL_PREFIX=$PREFIX ..
else
    cmake -DCMAKE_BUILD_TYPE="RELEASE" -DCMAKE_INSTALL_PREFIX='/usr/local' ..
fi
make -j"$NCPUS"

if [ "$?" = 0 ] && [ ! "$NO_INSTALL" ]; then
    make install
fi;

cd "$OWD"
