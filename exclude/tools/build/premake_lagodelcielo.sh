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
cd ${ROOT}/Downloads/bv_files/${PREMAKE_VERSION}
if [[ ! -f build_visit ]]; then
    echo "please provide build_visit file"
    exit 1
fi

ARGS=""
ARGS=${ARGS}" --debug --no-visit "
ARGS=${ARGS}" --fortran --cc gcc --cxx g++ --cxxflag -std=c++98 "
ARGS=${ARGS}" --makeflags -j24 "
ARGS=${ARGS}" --parallel "
ARGS=${ARGS}" --cmake --qt --python "
ARGS=${ARGS}" --alt-boost-dir /usr "
ARGS=${ARGS}" --hdf5 --netcdf --szip --silo --uintah"
ARGS=${ARGS}" --ispc --embree --tbb --ospray "
ARGS=${ARGS}" --alt-pidx-dir /home/sci/qwu/software/Lagodelcielo/PIDX/install "
ARGS=${ARGS}" --thirdparty-path /ssd/users/qwu/VisIt "
ARGS=${ARGS}" --installation-build-dir /scratch/users/qwu/VisIt "
    
if   [[ "$PREMAKE_VERSION" == "trunk" ]]; then
    ARGS=${ARGS}" --skip-opengl-context-check "
elif [[ "$PREMAKE_VERSION" == "rc2.13" ]]; then
    ARGS=${ARGS}" --slivr "
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