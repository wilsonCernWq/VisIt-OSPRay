#!/bin/sh
rm -f CMakeCache.txt
/home/sci/qwu/VisIt/visit/cmake/3.8.1/linux-x86_64_gcc-5.4/bin/cmake  -DCMAKE_BUILD_TYPE:STRING="Release" -DVTK_DIR="/home/sci/qwu/VisIt/visit/vtk/6.1.0/linux-x86_64_gcc-5.4/lib/cmake/vtk-6.1" /home/sci/qwu/VisIt/visitOSPRayCPU/visit-trunk/src
