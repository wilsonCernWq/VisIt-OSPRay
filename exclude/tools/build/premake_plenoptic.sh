#!/bin/bash

if [[ -z "$1" ]]; then
    echo "Error: the visit version is required as the 1st argument"
    echo "       candidate versions: trunk 2.13"
    exit -1
elif [[ "$1" == "trunk" ]]; then
    PREMAKE_VERSION=trunk
elif [[ "$1" == "2.13" ]]; then
    PREMAKE_VERSION=rc2.13
else
    echo "Error: unknown visit version: $1"
    exit -2
fi

ROOT=$(pwd)
cd ${ROOT}/3rdparty/bv_files/${PREMAKE_VERSION}
if [[ ! -f build_visit ]]; then
    echo "please provide build_visit file"
    exit 1
fi
DIR_BUILD=${ROOT}/3rdparty/builds/$PREMAKE_VERSION
DIR_INSTALL=${ROOT}/3rdparty/installs/$PREMAKE_VERSION
mkdir -p $DIR_BUILD
mkdir -p $DIR_INSTALL

ARGS=""
ARGS=${ARGS}" --no-visit --stdout "
ARGS=${ARGS}" --fortran --cc gcc --cxx g++ --cxxflag -std=c++11 "
ARGS=${ARGS}" --makeflags -j4 "
ARGS=${ARGS}" --parallel "
ARGS=${ARGS}" --cmake --qt --python "
ARGS=${ARGS}" --hdf5 --netcdf --szip --silo "
ARGS=${ARGS}" --thirdparty-path $DIR_INSTALL "
ARGS=${ARGS}" --installation-build-dir $DIR_BUILD "

if   [[ "$PREMAKE_VERSION" == "trunk" ]]; then
    #ARGS=${ARGS}" --stdout "
    ARGS=${ARGS}" --alt-tbb-dir /home/qadwu/Work/softwares/tbb/tbb-2018_20180618-install "
    ARGS=${ARGS}" --alt-ispc-dir /home/qadwu/Work/softwares/ispc/ispc-1.9.2-install "
    ARGS=${ARGS}" --alt-embree-dir /home/qadwu/Work/softwares/embree/embree-3.2.0-install "
    ARGS=${ARGS}" --alt-ospray-dir /home/qadwu/Work/projects/ospray/builds/ospray-visit-install/lib/cmake/ospray-1.8.0 "
    #ARGS=${ARGS}" --pidx "
elif [[ "$PREMAKE_VERSION" == "rc2.13" ]]; then
    #ARGS=${ARGS}" --stdout "
    ARGS=${ARGS}" --ispc --embree --tbb "
    ARGS=${ARGS}" --ospray "
    #ARGS=${ARGS}" --pidx "
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
