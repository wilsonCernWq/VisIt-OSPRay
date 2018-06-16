#!/bin/bash

PREMAKE_VERSION=trunk

ROOT=$(pwd)
cd ${ROOT}/3rdparty/bv_files
if [[ ! -f build_visit ]]; then
    echo "please provide build_visit file"
    exit 1
fi
DIR_BUILD=/Users/qwu/Work/visit/3rdparty/builds/$PREMAKE_VERSION
DIR_INSTALL=/Users/qwu/Work/visit/3rdparty/installs/$PREMAKE_VERSION
mkdir -p $DIR_BUILD
mkdir -p $DIR_INSTALL

ARGS=""
ARGS=${ARGS}" --debug --no-visit --static "
ARGS=${ARGS}" --fortran --cc gcc --cxx g++ --cxxflag -std=c++98 "
ARGS=${ARGS}" --makeflags -j2 "
ARGS=${ARGS}" --cmake --python --qt --parallel "
ARGS=${ARGS}" --hdf5 --szip --zlib --silo "
ARGS=${ARGS}" --ospray "
ARGS=${ARGS}" --alt-tbb-dir /Users/qwu/Work/softwares/tbb2018_20180312oss "
ARGS=${ARGS}" --alt-ispc-dir /Users/qwu/Work/softwares/ispc-v1.9.2-osx "
ARGS=${ARGS}" --alt-embree-dir /Users/qwu/Work/softwares/embree-2.17.4.x86_64.macosx "
ARGS=${ARGS}" --thirdparty-path $DIR_INSTALL "
ARGS=${ARGS}" --installation-build-dir $DIR_BUILD "

if   [[ "$PREMAKE_VERSION" == "trunk" ]]; then
    ARGS=${ARGS}" "
fi

# ------------------------------ 
# VisIt MPI
# ------------------------------ 
#export PAR_COMPILER=/usr/local/Cellar/mpich/3.2.1_2/mpicc
#export PAR_COMPILER_CXX=/usr/local/Cellar/mpich/3.2.1_2/mpicxx
#export PAR_INCLUDE="-I/usr/local/Cellar/mpich/3.2.1_2/include/"

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
