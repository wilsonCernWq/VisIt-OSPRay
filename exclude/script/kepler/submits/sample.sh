#!/bin/bash
#SBATCH --time=4:00:00 # walltime, abbreviated by -t
#
#
EXE=$HOME/software/Kepler/VisIt/visit-trunk/build/bin/visit
SCP=$HOME/VisIt/working/exclude/script/kepler/kepler_local_pidx_benchmark.py
#
# set data and working directories
NUM_THREADS=$((SLURM_NTASKS / SLURM_NNODES))
NUM_NODES=$SLURM_NNODES
WORKDIR=$HOME/Desktop/timing/08-11/Kepler/debug
cd $WORKDIR
#
# run job
$EXE -np 16 -nn 16 -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP &> \
    job-$SLURM_JOBID-$NUM_NODES-$NUM_THREADS.out

