#!/bin/bash

# define source location
# VISITSRC=/home/sci/qwu/visitOSPRay/visitOSPRayCPU/src-visit2.12.0
# WORKSRC=/home/sci/qwu/visitOSPRay/visitOSPRayCPU/working

VISITSRC=${1}
WORKSRC=${2}

# shopt -s extglob

if [ "$3" == "source" ]; then 

    cp -vr ${WORKSRC}/plots/*                 ${VISITSRC}/plots
    cp -vr ${WORKSRC}/viewer/main/*           ${VISITSRC}/viewer/main
    cp -vr ${WORKSRC}/engine/main/*           ${VISITSRC}/engine/main
    cp -vr ${WORKSRC}/avt/Filters/*.C         ${VISITSRC}/avt/Filters
    cp -vr ${WORKSRC}/avt/Filters/*.h         ${VISITSRC}/avt/Filters
    cp -vr ${WORKSRC}/avt/Plotter/vtk/*.C     ${VISITSRC}/avt/Plotter/vtk
    cp -vr ${WORKSRC}/avt/Pipeline/Data/*.C   ${VISITSRC}/avt/Pipeline/Data
    cp -vr ${WORKSRC}/avt/Pipeline/Data/*.h   ${VISITSRC}/avt/Pipeline/Data

elif [ "$3" == "cmake" ]; then 
    
    cp -v ${WORKSRC}/config-site/hastur.sci.utah.edu.cmake   ${VISITSRC}/config-site
    cp -v ${WORKSRC}/CMakeLists.txt                          ${VISITSRC}
    cp -v ${WORKSRC}/avt/Filters/CMakeLists.txt              ${VISITSRC}/avt/Filters
    cp -v ${WORKSRC}/avt/Plotter/CMakeLists.txt              ${VISITSRC}/avt/Plotter
    cp -v ${WORKSRC}/avt/Plotter/vtk/*.in                    ${VISITSRC}/avt/Plotter/vtk

fi

# shopt -u extglob
