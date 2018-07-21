#!/bin/bash

if [[ -z "$1" ]]; then
    echo "Error: the visit version is required as the 1st argument"
    echo "       candidate versions: trunk 2.13"
    exit -1
elif [[ "$1" == "trunk" ]]; then
    PREMAKE_VERSION=trunk
elif [[ "$1" == "2.13" || "$1" == "rc2.13" ]]; then
    PREMAKE_VERSION=rc2.13
else
    echo "Error: unknown visit version: $1"
    exit -2
fi

# ------------------------------ 
# configure
# ------------------------------ 
ROOT=$(pwd)
cd ${ROOT}/bv_files/${PREMAKE_VERSION}
if [[ ! -f build_visit ]]; then
    echo "please provide build_visit file"
    exit 1
fi
DIR_BUILD=/data/qwu/VisIt/3rdparty/$PREMAKE_VERSION
DIR_INSTALL=/ssd/qwu/VisIt/3rdparty/$PREMAKE_VERSION
mkdir -p $DIR_BUILD
mkdir -p $DIR_INSTALL

ARGS=""
ARGS=${ARGS}" --debug --no-visit "
ARGS=${ARGS}" --fortran --cc gcc --cxx g++ --cxxflag -std=c++11 "
ARGS=${ARGS}" --makeflags -j8 --parallel "
ARGS=${ARGS}" --cmake --qt --python"
ARGS=${ARGS}" --boost --hdf5 --netcdf --szip --silo "
ARGS=${ARGS}" --alt-tbb-dir /home/sci/qwu/software/tbb2018_20180618oss "
ARGS=${ARGS}" --alt-ispc-dir /home/sci/qwu/software/ispc-v1.9.2-linux "
ARGS=${ARGS}" --alt-embree-dir /home/sci/qwu/software/embree-3.2.0.x86_64.linux "
ARGS=${ARGS}" --alt-ospray-dir /home/sci/qwu/OSPRay/Hastur/install-visit/lib64/cmake/ospray-1.6.0 "
ARGS=${ARGS}" --thirdparty-path $DIR_INSTALL "
ARGS=${ARGS}" --installation-build-dir $DIR_BUILD "

if   [[ "$PREMAKE_VERSION" == "trunk" ]]; then
    ARGS=${ARGS}" --skip-opengl-context-check --llvm --mesagl "
elif [[ "$PREMAKE_VERSION" == "rc2.13" ]]; then
    ARGS=${ARGS}" --slivr --uintah "
fi

# ------------------------------ 
# MPI
# ------------------------------ 
source /opt/intel/parallel_studio_xe_2017.2.050/psxevars.sh
export PAR_COMPILER=/opt/intel/impi/2017.2.174/bin64/mpicc
export PAR_COMPILER_CXX=/opt/intel/impi/2017.2.174/bin64/mpicxx
export PAR_INCLUDE="-I/opt/intel/impi/2017.2.174/include64/"

# ------------------------------ 
# 3rdparties
# ------------------------------ 

# export ICET_FILE=IceT-devel.tar.gz
# export ICET_VERSION=devel
# export ICET_COMPATIBILITY_VERSION=devel
# export ICET_BUILD_DIR=IceT-devel
    
echo "yes" | bash ./build_visit $ARGS
