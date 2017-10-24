#!/bin/bash
#COBALT -t 60
#COBALT -n 2
#COBALT --attrs mcdram=cache:numa=quad
#COBALT -A Viz_Support

echo "Starting Cobalt job script"
source /home/qiwu/software/embree-2.15.0.x86_64.linux/embree-vars.sh 
EXE=/home/qiwu/software/VisIt/cooley/visit-trunk/build/bin/visit
SCP=/home/qiwu/software/VisIt/changes/exclude/script/cooley/cooley_local_pidx_benchmark.py

# set data and working directories
# WORKDIR=/home/qiwu/timings/visit/theta/submit
# cd $WORKDIR

# run job
$EXE -np 2 -nn 2 -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
