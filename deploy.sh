#!/bin/bash

cd /home/sci/qwu/visitOSPRay/visitOSPRayCPU/pascal

if [ "$1" == "source" ]; then 

    cp -r plots ../src-visit2.12.0/

    cp -r avt/Filters/*.cpp ../src-visit2.12.0/avt/Filters
    cp -r avt/Filters/*.C   ../src-visit2.12.0/avt/Filters
    cp -r avt/Filters/*.h   ../src-visit2.12.0/avt/Filters

    cp -r avt/Pipeline/Data/*.cpp ../src-visit2.12.0/avt/Pipeline/Data
    cp -r avt/Pipeline/Data/*.C   ../src-visit2.12.0/avt/Pipeline/Data
    cp -r avt/Pipeline/Data/*.h   ../src-visit2.12.0/avt/Pipeline/Data

elif [ "$1" == "cmake" ]; then 
    
    cp ./*.txt           ../src-visit2.12.0/
    cp avt/Filters/*.txt ../src-visit2.12.0/avt/Filters/

fi

cd -
