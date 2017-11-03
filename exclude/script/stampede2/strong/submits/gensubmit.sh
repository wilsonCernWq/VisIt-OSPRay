#!/bin/bash

NT=$2 # number of tasks
NN=$3 # number of node
NC=$4 # number of cores per task
NAME=n${NN}p$((NT / NN))
EXE=$WORK/software/VisIt/stampede2/visit-trunk/build/bin/visit
SCP=$1

mkdir -p $NAME
cd $NAME

cat >> sample-$NAME.sh <<EOF
#!/bin/bash
#SBATCH -n ${NT}
#SBATCH -N ${NN}
#SBATCH --time=00:30:00 # walltime, abbreviated by -t
#
# load modules
module load intel
module load impi
module load hdf5
export LIBRARY_PATH=\${LIBRARY_PATH}:/work/04915/qiwu/stampede2/software/VisIt/stampede2/visit/vtk/6.1.0/linux-x86_64_gcc-5.4/lib
export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/work/04915/qiwu/stampede2/software/VisIt/stampede2/visit/vtk/6.1.0/linux-x86_64_gcc-5.4/lib
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=$NC
export KMP_AFFINITY=verbose,granularity=core,compact,1,0
export OSPRAY_THREADS=$NC
#
# run job
echo "job: <\$SLURM_JOBID> node = $NN processes-per-node = $((NT / NN))"
$EXE -np $NT -nn $NN -l ibrun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP 
EOF
sbatch -p normal -A TG-ASC170049 sample-$NAME.sh
cd -
