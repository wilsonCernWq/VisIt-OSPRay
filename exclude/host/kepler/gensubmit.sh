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
mkdir -p no-icet
cd no-icet
export SLIVR_USE_ICET=0
$EXE -np ${NT} -nn ${NN} -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
cd ..

mkdir -p icet-reduce
cd icet-reduce
export SLIVR_USE_ICET=1
export SLIVR_ICET_STRATEGY=0
$EXE -np ${NT} -nn ${NN} -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
cd ..

mkdir -p icet-tree
cd icet-tree
export SLIVR_USE_ICET=1
export SLIVR_ICET_STRATEGY=1
$EXE -np ${NT} -nn ${NN} -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
cd ..

mkdir -p icet-radixk
cd icet-radixk
export SLIVR_USE_ICET=1
export SLIVR_ICET_STRATEGY=2
$EXE -np ${NT} -nn ${NN} -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
cd ..

mkdir -p icet-bswap
cd icet-bswap
export SLIVR_USE_ICET=1
export SLIVR_ICET_STRATEGY=3
$EXE -np ${NT} -nn ${NN} -l mpirun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP
cd ..

EOF
sbatch sample-$NAME.sh
cd -
