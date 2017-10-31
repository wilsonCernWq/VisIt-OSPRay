#!/bin/bash

NT=$2 # number of tasks
NN=$3 # number of node
NC=$4 # number of cores per task
NAME=n${NN}p$((NT / NN))
EXE=$HOME/software/Kepler/VisIt/visit-trunk/build/bin/visit
SCP=$1

mkdir -p $NAME
cd $NAME

cat >> sample-$NAME.sh <<EOF
#!/bin/bash
#SBATCH -n ${NT}
#SBATCH -N ${NN}
#SBATCH --time=4:00:00 # walltime, abbreviated by -t
#
# run job
$EXE -np ${NT} -nn ${NN} -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
EOF
sbatch sample-$NAME.sh
cd -
