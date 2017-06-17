#!/bin/sh
rm -f CMakeCache.txt
/home/sci/qwu/VisIt/visit/cmake/3.8.1/linux-x86_64_gcc-5.4/bin/cmake  -DCMAKE_BUILD_TYPE:STRING="Release" /home/sci/qwu/VisIt/visitOSPRayCPU/visit-trunk/src
