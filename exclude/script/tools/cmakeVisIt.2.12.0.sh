#!/bin/sh
rm -f CMakeCache.txt
/home/sci/qwu/VisIt/visit/cmake/3.0.2/linux-x86_64_gcc-4.8/bin/cmake  -DCMAKE_BUILD_TYPE:STRING="Debug" -DVISIT_CXX_FLAGS:STRING="-g -Wno-deprecated" -DVISIT_C_FLAGS:STRING="-g" -DVTK_DIR="/home/sci/qwu/VisIt/visit/vtk/6.1.0/linux-x86_64_gcc-4.8" -DPIDX_DIR="/home/sci/qwu/VisIt/visitOSPRayCPU/PIDX/install" $1
