#!/bin/bash

NT=$1 # number of tasks
NN=$2 # number of node
NC=$3 # number of cores per task
NAME=n${NN}p$((NT / NN))
EXE=$HOME/software/Kepler/VisIt/visit-trunk/build/bin/visit
SCP=$HOME/VisIt/working/exclude/script/kepler/kepler_local_pidx_benchmark.py

mkdir -p $NAME
cd $NAME

cat >> sample-$NAME.sh <<EOF
#!/bin/bash
#SBATCH -n ${1}
#SBATCH -N ${2}
#SBATCH -c ${3}
#SBATCH --time=4:00:00 # walltime, abbreviated by -t
#
# run job
$EXE -np ${NT} -nn ${NN} -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
EOF
sbatch sample-$NAME.sh
cd -
