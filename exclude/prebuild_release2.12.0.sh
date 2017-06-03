#!/bin/bash

# load module
# module load cmake/3.7.2
# module load gcc mvapich2

ROOT=$(pwd)

# # ------------------------------ 
# # Build OSPRay
# # ------------------------------ 

# cd ${ROOT}

# # build ospray
# mkdir -p ospray
# cd ospray

# wget https://downloads.sourceforge.net/project/ispcmirror/v1.9.1/ispc-v1.9.1-linux.tar.gz
# tar -xzf ./ispc-v1.9.1-linux.tar.gz

# wget https://github.com/01org/tbb/releases/download/2017_U5/tbb2017_20170226oss_lin.tgz
# tar -xzf ./tbb2017_20170226oss_lin.tgz

# wget https://github.com/embree/embree/releases/download/v2.15.0/embree-2.15.0.x86_64.linux.tar.gz
# tar -xzf ./embree-2.15.0.x86_64.linux.tar.gz

# wget https://github.com/wilsonCernWq/ospray/archive/qwu-devel.zip
# unzip ./qwu-devel.zip

# # compile ospray
# OSPROOT=$(pwd)

# # -- TBB component
# export TBB_ROOT=${OSPROOT}/tbb2017_20170226oss
# source ${TBB_ROOT}/bin/tbbvars.sh intel64
# export LD_LIBRARY_PATH=${TBB_ROOT}/lib/intel64/gcc4.7:${LD_LIBRARY_PATH}

# # -- ispc component
# export ISPC_ROOT=${OSPROOT}/ispc-v1.9.1-linux
# export PATH=${ISPC_ROOT}:${PATH}

# # -- embree component
# export EMBREE_ROOT=${OSPROOT}/embree-2.15.0.x86_64.linux
# export embree_DIR=${EMBREE_ROOT}

# # hide embree tbb library
# mkdir ${EMBREE_ROOT}/tbb
# mv ${EMBREE_ROOT}/lib/libtbb* ${EMBREE_ROOT}/tbb 2> /dev/null

# # cmake ospray
# export OSPRAY_BUILD=${OSPROOT}/ospray-qwu-devel/build
# export OSPRAY_INSTALL=${OSPROOT}/ospray-qwu-devel/install

# # build ospray
# cd ospray-qwu-devel
# mkdir -p build
# mkdir -p install
# cd build
# rm ./CMakeCache.txt
# cmake -DOSPRAY_APPS_EXAMPLEVIEWER=OFF -DOSPRAY_APPS_PARAVIEW_TFN_CVT=OFF -DOSPRAY_APPS_VOLUMEVIEWER=OFF -DCMAKE_INSTALL_PREFIX=${OSPRAY_INSTALL} ..
# make -j8 install
# cd ..

# # construct script to source ospray
# cat > thisospray.sh << EOF
# #!/bin/bash
# export PATH=$PATH:${OSPRAY_INSTALL}/bin
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${OSPRAY_INSTALL}/bin/lib64
# source ${EMBREE_ROOT}/embree-vars.sh 
# EOF

# cd ../..

# ------------------------------ 
# Build VisIt
# ------------------------------ 

cd ${ROOT}

# # visit mpi
# export PAR_INCLUDE=-I/uufs/chpc.utah.edu/sys/installdir/mvapich2/2.1/include/
# export PAR_COMPILER=/uufs/chpc.utah.edu/sys/installdir/mvapich2/2.1/bin/mpic++
# export PAR_COMPILER_CXX=/uufs/chpc.utah.edu/sys/installdir/mvapich2/2.1/bin/mpic++

# source /opt/intel/parallel_studio_xe_2016.1.056/psxevars.sh
# export PAR_COMPILER=/opt/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/bin/mpicc
# export PAR_COMPILER_CXX=/opt/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/bin/mpicxx
# export PAR_INCLUDE="-I/opt/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/include"

# wget -O build_visit2_12_0 http://portal.nersc.gov/project/visit/releases/2.12.0/build_visit2_12_0

# echo "yes" | bash ./build_visit2_12_0 --makeflags -j8 --parallel --debug --flags-debug --cxxflags -g --cflags -g --fortran --no-visit --boost --hdf5 --netcdf --silo --szip --slivr --uintah --mesa

# wget http://portal.nersc.gov/project/visit/releases/2.12.0/visit2.12.0.tar.gz

# tar -xzf ./visit2.12.0.tar.gz

cp $(hostname).cmake $(hostname).cmake.bkp
cat >> $(hostname).cmake.bkp << EOF

##
## OSPRay
## -- recommend to remove tbb libraries inside embree and ospray binary folder
##
SET(OSPRAY_USE_EXTERNAL_EMBREE ON)
SET(ospray_DIR ${ROOT}/../ospray-qwu-kepler-23.05.2017/lib64/cmake/ospray-1.3.0)
SET(embree_DIR ${ROOT}/../embree-2.13.0.x86_64.linux/lib/cmake/embree-2.13.0)
SET(TBB_ROOT ${ROOT}/../tbb2017_20160916oss)
SET(ISPC_EXECUTABLE ${ROOT}/../ispc-v1.9.1-linux)
VISIT_OPTION_DEFAULT(VISIT_OSPRAY ON TYPE BOOL)

##
## PIDX
##
SET(PUDX_DIR /home/sci/qwu/software/visitospray-kepler/PIDX/install)

EOF

cp $(hostname).cmake.bkp visit2.12.0/src/config-site/$(hostname).cmake

# cd ${ROOT}
# git clone https://bitbucket.org/WilsonOverCloud/visit-pascal-volumerendering.git changes

# cd ${ROOT}/changes
# git pull origin master

cd ${ROOT}
# bash ${ROOT}/changes/deploy.sh ${ROOT}/visit2.12.0/src ${ROOT}/changes source
# bash ${ROOT}/changes/deploy.sh ${ROOT}/visit2.12.0/src ${ROOT}/changes cmake

CHANGEDIR=/home/sci/qwu/VisIt/visitOSPRayCPU/working
bash ${CHANGEDIR}/deploy.sh ${ROOT}/visit2.12.0/src ${CHANGEDIR} "source"
bash ${CHANGEDIR}/deploy.sh ${ROOT}/visit2.12.0/src ${CHANGEDIR} "cmake"

mkdir -p ${ROOT}/visit2.12.0/build
cd ${ROOT}/visit2.12.0/build

# make visit
rm -f CMakeCache.txt
LIBRARY_PATH=${LIBRARY_PATH}:${ROOT}/ospray/embree-2.15.0.x86_64.linux/lib/cmake/embree-2.15.0/../..
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOT}/ospray/embree-2.15.0.x86_64.linux/lib/cmake/embree-2.15.0/../..
${ROOT}/cmake-3.0.2/bin/cmake -DCMAKE_BUILD_TYPE:STRING="Debug" -DVISIT_CXX_FLAGS:STRING="-g -Wno-deprecated" -DVISIT_C_FLAGS:STRING="-g" -DPIDX_DIR="/home/sci/qwu/software/visitospray-kepler/PIDX/install" ${ROOT}/visit2.12.0/src

make -j16
