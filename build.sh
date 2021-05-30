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

OWD="${PWD}"
BASE=`dirname "${0}"`

cd "${BASE}"

if [ ! -e "yaml-cpp/CMakeLists.txt" ]; then
    # Manually clone yaml-cpp, removing the directory first
    rm -r yaml-cpp
    git clone --depth=1 --branch=yaml-cpp-0.6.2 https://github.com/jbeder/yaml-cpp.git
fi

cd yaml-cpp
git checkout yaml-cpp-0.6.2

TMPFILE=`mktemp`
cp src/regex_yaml.h "${TMPFILE}"
echo '#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"' > src/regex_yaml.h
cat "${TMPFILE}" >> src/regex_yaml.h
rm "${TMPFILE}"
cd ..

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
