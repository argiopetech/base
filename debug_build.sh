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

if [ ! -e "yaml-cpp/CMakeLists.txt" ]; then
    # Manually clone yaml-cpp, removing the directory first
    rm -r yaml-cpp
    git clone --depth=1 --branch=yaml-cpp-0.6.2 https://github.com/jbeder/yaml-cpp.git
fi

cd ./BUILD
cmake -DCMAKE_BUILD_TYPE="RELWITHDEBINFO" -DCMAKE_INSTALL_PREFIX='.' ..
make -j"$NCPUS"

if [ "$?" = 0 ] && [ ! "$NO_INSTALL" ]; then
    make install
fi;

cd "$OWD"
