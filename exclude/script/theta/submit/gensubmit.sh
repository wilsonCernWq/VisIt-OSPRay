#!/bin/bash
NT=$2 # number of tasks
NN=$3 # number of node
NC=$4 # number of cores per task
NAME=n${NN}p$((NT / NN))
EXE=/home/qiwu/software/VisIt/theta/visit-trunk/build/bin/visit
SCP=$1

mkdir -p $NAME
cd $NAME

cat >> sample-$NAME.sh <<EOF
#!/bin/bash
#COBALT -t 60
#COBALT -n ${NN}
#COBALT -A Viz_Support

#
#
export LIBRARY_PATH=/home/qiwu/software/VisIt/theta/visit/python/2.7.11/linux-x86_64_cc/lib/python2.7/lib-dynload:/home/qiwu/software/VisIt/theta/visit/python/2.7.11/linux-x86_64_cc/lib:\${LIBRARY_PATH}:/home/qiwu/software/VisIt/theta/visit/vtk/6.1.0/linux-x86_64_cc/lib:/home/qiwu/software/PIDX/install/theta/lib

export LD_LIBRARY_PATH=/home/qiwu/software/VisIt/theta/visit/python/2.7.11/linux-x86_64_cc/lib/python2.7/lib-dynload:/home/qiwu/software/VisIt/theta/visit/python/2.7.11/linux-x86_64_cc/lib:\${LD_LIBRARY_PATH}:/home/qiwu/software/VisIt/theta/visit/vtk/6.1.0/linux-x86_64_cc/lib:/home/qiwu/software/PIDX/install/theta/lib

export PATH=/home/qiwu/software/VisIt/theta/visit/python/2.7.11/linux-x86_64_cc/bin:\${PATH}

#
#
module unload PrgEnv-intel
module load PrgEnv-gnu
source /home/qiwu/software/embree-2.15.0.x86_64.linux/embree-vars.sh
export CRAYPE_LINK_TYPE=dynamic
export MPICH_MAX_THREAD_SAFETY=multiple
export OSPRAY_THREADS=64
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=256
export KMP_AFFINITY=verbose,granularity=core,compact,1,0
export CRAY_OMP_CHECK_AFFINITY=TRUE
#
# run job
echo "job: <\$SLURM_JOBID> node = $NN processes-per-node = $((NT / NN))"
$EXE -np $NT -nn $NN -l aprun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
EOF
chmod +x sample-$NAME.sh
qsub sample-$NAME.sh
cd -
