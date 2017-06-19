#!/bin/bash
DEST=$1
mv /home/sci/qwu/slurm-*.out $DEST
mv /home/sci/qwu/*.engine_par.*.vlog $DEST
mv /home/sci/qwu/*.engine_ser.*.vlog $DEST
mv /home/sci/qwu/engine_par.*.timings $DEST
cd /home/sci/qwu/
source /home/sci/qwu/software/Kepler/VisIt/visit-depend.sh deactivate
