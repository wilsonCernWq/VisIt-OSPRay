#!/bin/bash
DEST=$1
cd /home/sci/qwu/
mv $(ls -1 /home/sci/qwu/slurm-*.out | head -n 1) $DEST
mv /home/sci/qwu/*.engine_par.*.vlog $DEST
mv /home/sci/qwu/*.engine_ser.*.vlog $DEST
mv /home/sci/qwu/engine_par.*.timings $DEST
cd /home/sci/qwu/
source /home/sci/qwu/software/Kepler/VisIt/visit-depend.sh deactivate
