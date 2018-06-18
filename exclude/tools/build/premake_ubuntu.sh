#!/bin/bash

PREMAKE_VERSION=trunk

ROOT=$(pwd)
cd ${ROOT}/3rdparty/bv_files
if [[ ! -f build_visit ]]; then
    echo "please provide build_visit file"
    exit 1
fi
DIR_BUILD=${ROOT}/3rdparty/builds/$PREMAKE_VERSION
DIR_INSTALL=${ROOT}/3rdparty/installs/$PREMAKE_VERSION
mkdir -p $DIR_BUILD
mkdir -p $DIR_INSTALL

ARGS=""
ARGS=${ARGS}" --debug --no-visit "
ARGS=${ARGS}" --fortran --cc gcc --cxx g++ --cxxflag -std=c++98 "
ARGS=${ARGS}" --makeflags -j2 "
ARGS=${ARGS}" --cmake --python --qt --parallel "
ARGS=${ARGS}" --hdf5 --szip --zlib --silo "
ARGS=${ARGS}" --ospray "
ARGS=${ARGS}" --alt-tbb-dir /home/qwu/software/tbb2017_20170604oss "
ARGS=${ARGS}" --alt-ispc-dir /home/qwu/software/ispc-v1.9.1-linux "
ARGS=${ARGS}" --alt-embree-dir /home/qwu/software/embree-3.2.0.x86_64.linux "
ARGS=${ARGS}" --thirdparty-path $DIR_INSTALL "
ARGS=${ARGS}" --installation-build-dir $DIR_BUILD "

if   [[ "$PREMAKE_VERSION" == "trunk" ]]; then
    ARGS=${ARGS}" "
fi

# ------------------------------ 
# VisIt MPI
# ------------------------------ 
export PAR_COMPILER=/usr/bin/mpicc
export PAR_COMPILER_CXX=/usr/bin/mpicxx
export PAR_INCLUDE="-I/usr/lib/x86_64-linux-gnu/openmpi/include"

# ------------------------------ 
# VisIt Environment
# ------------------------------ 

echo "yes" | bash ./build_visit $ARGS
