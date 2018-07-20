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
cd ${ROOT}/bv_files/${PREMAKE_VERSION}
if [[ ! -f build_visit ]]; then
    echo "please provide build_visit file"
    exit 1
fi
DIR_BUILD=/scratch/users/qwu/VisIt/3rdparty/$PREMAKE_VERSION
DIR_INSTALL=/ssd/users/qwu/VisIt/3rdparty/$PREMAKE_VERSION
mkdir -p $DIR_BUILD
mkdir -p $DIR_INSTALL

ARGS=""
ARGS=${ARGS}" --debug --no-visit "
ARGS=${ARGS}" --fortran --cc gcc --cxx g++ --cxxflag -std=c++11 "
ARGS=${ARGS}" --makeflags -j24 "
ARGS=${ARGS}" --parallel "
ARGS=${ARGS}" --cmake --qt --python "
ARGS=${ARGS}" --hdf5 --netcdf --szip --silo "
ARGS=${ARGS}" --alt-boost-dir /usr "
ARGS=${ARGS}" --thirdparty-path $DIR_INSTALL "
ARGS=${ARGS}" --installation-build-dir $DIR_BUILD "

if   [[ "$PREMAKE_VERSION" == "trunk" ]]; then
ARGS=${ARGS}" --stdout "
    ARGS=${ARGS}" --skip-opengl-context-check --mesagl --llvm "
    ARGS=${ARGS}" --ispc --embree --tbb "
    ARGS=${ARGS}" --alt-ospray-dir /home/sci/qwu/OSPRay/Lagodelcielo/install-visit/lib64/cmake/ospray-1.6.0 "
    #ARGS=${ARGS}" --pidx "
    #ARGS=${ARGS}" --uintah "
elif [[ "$PREMAKE_VERSION" == "rc2.13" ]]; then
    ARGS=${ARGS}" --ispc --embree "
    ARGS=${ARGS}" --ospray --alt-tbb-dir /home/sci/qwu/software/tbb2018_20180618oss "
    ARGS=${ARGS}" --alt-pidx-dir /home/sci/qwu/software/Lagodelcielo/PIDX/install "
    ARGS=${ARGS}" --slivr --uintah "
fi

# ------------------------------ 
# VisIt MPI
# ------------------------------ 
export PAR_COMPILER=/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpicc
export PAR_COMPILER_CXX=/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpicxx
export PAR_INCLUDE="-I/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/include"

# ------------------------------ 
# VisIt Environment
# ------------------------------ 

# export UINTAH_FILE=Uintah-2.1.0.tar.gz
# export UINTAH_VERSION=2.1.0
# export UINTAH_COMPATIBILITY_VERSION=2.1
# export UINTAH_BUILD_DIR=Uintah-2.1.0/optimized

# export ICET_FILE=IceT-devel.tar.gz
# export ICET_VERSION=devel
# export ICET_COMPATIBILITY_VERSION=devel
# export ICET_BUILD_DIR=IceT-devel

# export ICET_FILE=IceT-1-0-0.tar.gz
# export ICET_VERSION=1.0.0
# export ICET_COMPATIBILITY_VERSION=1.0.0
# export ICET_BUILD_DIR=IceT-1-0-0

# sudo pip2 install numpy pillow Cython requests pyparsing seedme pil-compat setuptools

echo "yes" | bash ./build_visit $ARGS
