#!/bin/bash
#SBATCH -n 8
#SBATCH -N 8
#SBATCH -c 68
#SBATCH --time=0:30:00 # walltime, abbreviated by -t
#
#
export LIBRARY_PATH=${LIBRARY_PATH}:/work/04915/qiwu/stampede2/software/VisIt/stampede2/visit/vtk/6.1.0/linux-x86_64_gcc-5.4/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/work/04915/qiwu/stampede2/software/VisIt/stampede2/visit/vtk/6.1.0/linux-x86_64_gcc-5.4/lib
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=68
export KMP_AFFINITY=verbose,granularity=core,compact,1,0
export OSPRAY_THREADS=68
#
#
EXE=/work/04915/qiwu/stampede2/software/VisIt/stampede2/visit-trunk/build/bin/visit
SCP=/work/04915/qiwu/stampede2/software/VisIt/stampede2/changes/exclude/script/stampede2/stampede2_local_pidx_benchmark.py
#
# set data and working directories
NUM_THREADS=$((SLURM_NTASKS / SLURM_NNODES))
NUM_NODES=$SLURM_NNODES
#
# run job
echo "job: <$SLURM_JOBID> node = $NUM_NODES threads-per-node = $NUM_THREADS"
$EXE -np 8 -nn 8 -l ibrun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
