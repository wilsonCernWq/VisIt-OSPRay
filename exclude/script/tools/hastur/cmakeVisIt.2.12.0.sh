#!/bin/sh
rm -rf CMakeCache.txt CMakeFiles/
/home/sci/qwu/VisIt/visit/cmake/3.0.2/linux-x86_64_gcc-4.8/bin/cmake  -DCMAKE_BUILD_TYPE:STRING="Release" -DVTK_DIR="/home/sci/qwu/VisIt/visit/vtk/6.1.0/linux-x86_64_gcc-4.8/lib/cmake/vtk-6.1" /home/sci/qwu/VisIt/visitOSPRayCPU/visit-trunk/src
