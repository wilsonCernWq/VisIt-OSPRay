#!/bin/bash

THREAD=${1}

ROOT=$(pwd)
cp ~/VisIt/script/kepler_benchmark.py ${ROOT}
ssh -n qwu@wopr.sci.utah.edu rm /tmp/visit.qiwu.*

for NODE in 1 2 4 8 16
do

    # ----------------------------
    mkdir ${ROOT}/n${NODE}p${THREAD}
    cd ${ROOT}/n${NODE}p${THREAD}
    
    # set variables
    export VISIT_NODES=${NODE}
    export VISIT_THREADS=${THREAD}
    export VISIT_WORKDIR=$(pwd)

    # run
    ~/VisIt/visitOSPRayCPU/build-release2.12.0/bin/visit -nowin -cli -s ${ROOT}/kepler_benchmark.py -withhold-timing-output -timing -debug 1

    # clean up
    mv ~/*.vlog .
    mv ~/*.timings .
    ssh -n qwu@wopr.sci.utah.edu mv /tmp/visit.qwu.* ${ROOT}/n${NODE}p${THREAD}
    unset VISIT_NODES
    unset VISIT_THREADS
    unset VISIT_WORKDIR

    # go back
    cd ${ROOT}
    # ----------------------------

done
