#!/bin/bash

cd /home/sci/qwu/visitOSPRay/visitOSPRayCPU/pascal

if [ "$1" == "source" ]; then 

    cp -r plots/* ../src-visit2.12.0/plots

    cp -r avt/Filters/*.C   ../src-visit2.12.0/avt/Filters
    cp -r avt/Filters/*.h   ../src-visit2.12.0/avt/Filters

    cp -r avt/Plotter/vtk/*.C     ../src-visit2.12.0/avt/Plotter/vtk

    cp -r avt/Pipeline/Data/*.C   ../src-visit2.12.0/avt/Pipeline/Data
    cp -r avt/Pipeline/Data/*.h   ../src-visit2.12.0/avt/Pipeline/Data

elif [ "$1" == "cmake" ]; then 
    
    cp ./hastur.sci.utah.edu.cmake ../src-visit2.12.0/config-site
    cp ./CMakeLists.txt            ../src-visit2.12.0
    cp avt/Filters/CMakeLists.txt  ../src-visit2.12.0/avt/Filters
    cp avt/Plotter/CMakeLists.txt  ../src-visit2.12.0/avt/Plotter
    cp avt/Plotter/vtk/*.in        ../src-visit2.12.0/avt/Plotter/vtk
fi

cd -
