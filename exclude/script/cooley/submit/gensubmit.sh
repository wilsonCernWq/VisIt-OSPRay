#!/bin/bash
NT=$2 # number of tasks
NN=$3 # number of node
NC=$4 # number of cores per task
NAME=n${NN}p$((NT / NN))
EXE=/home/qiwu/software/VisIt/cooley/visit-trunk/build/bin/visit
SCP=$1

mkdir -p $NAME
cd $NAME

cat >> sample-$NAME.sh <<EOF
#!/bin/bash
#COBALT -t 480
#COBALT -n ${NN}
#COBALT -A Viz_Support
#
# load modules
source /home/qiwu/software/embree-2.15.0.x86_64.linux/embree-vars.sh 
export LIBRARY_PATH=${LIBRARY_PATH}:/home/qiwu/software/VisIt/cooley/visit/vtk/6.1.0/linux-rhel_6-x86_64_gcc-4.9/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/qiwu/software/VisIt/cooley/visit/vtk/6.1.0/linux-rhel_6-x86_64_gcc-4.9/lib
#
# run job
echo "job: <\$SLURM_JOBID> node = $NN processes-per-node = $((NT / NN))"
$EXE -np $NT -nn $NN -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
EOF
chmod +x sample-$NAME.sh
qsub sample-$NAME.sh
cd -
