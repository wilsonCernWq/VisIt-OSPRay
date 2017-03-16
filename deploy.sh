#!/bin/bash

cd /home/sci/qwu/visitOSPRay/visitOSPRayCPU/pascal

if [ "$1" == "plot" ]; then 

    cp -v plots/Volume/* ../src-visit2.12.0/plots/Volume    

elif [ "$1" == "filter" ]; then 

    cp -v avt/Filters/*.C   ../src-visit2.12.0/avt/Filters
    cp -v avt/Filters/*.h   ../src-visit2.12.0/avt/Filters

elif [ "$1" == "pipeline" ]; then 

    cp -v avt/Pipeline/Data/*.C   ../src-visit2.12.0/avt/Pipeline/Data
    cp -v avt/Pipeline/Data/*.h   ../src-visit2.12.0/avt/Pipeline/Data

elif [ "$1" == "source" ]; then

    cp -v plots/Volume/* ../src-visit2.12.0/plots/Volume
    cp -v avt/Filters/*.C   ../src-visit2.12.0/avt/Filters
    cp -v avt/Filters/*.h   ../src-visit2.12.0/avt/Filters
    cp -v avt/Pipeline/Data/*.C   ../src-visit2.12.0/avt/Pipeline/Data
    cp -v avt/Pipeline/Data/*.h   ../src-visit2.12.0/avt/Pipeline/Data

elif [ "$1" == "vtkospray" ]; then

    gvfs-trash ../src-visit2.12.0/avt/Plotter/OSPRay/vtkOSPRay/*
    cp -vr avt/Plotter/OSPRay/vtkOSPRay/* ../src-visit2.12.0/avt/Plotter/OSPRay/vtkOSPRay

elif [ "$1" == "cmake" ]; then 
    
    cp -v avt/Plotter/OSPRay/CMake/* ../src-visit2.12.0/avt/Plotter/OSPRay/CMake
    cp -v ./*.txt           ../src-visit2.12.0/
    cp -v avt/Filters/*.txt ../src-visit2.12.0/avt/Filters/

fi

cd -
