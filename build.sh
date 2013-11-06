#!/bin/sh

# Get the number of CPUs on the system
# This is not portable
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
  MULTICPUS="-j`grep -c ^processor /proc/cpuinfo`"
elif [[ "$unamestr" == 'Darwin' ]]; then
  MULTICPUS="-j`sysctl -n hw.ncpu`"
fi

BASE=`dirname ${0}`

pushd $BASE

# Setup submodules
git submodule init
cp cmake/yaml-cpp-CMakeLists.txt yaml-cpp/CMakeLists.txt

pushd ./BUILD
cmake -DCMAKE_INSTALL_PREFIX='.' -DCMAKE_BUILD_TYPE="RELWITHDEBINFO" ..
make $MULTICPUS

make install

popd
popd
