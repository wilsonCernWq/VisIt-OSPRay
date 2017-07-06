#!/bin/bash
DEST=$1
mv /gpfs/mira-home/qiwu/*.error $DEST
mv /gpfs/mira-home/qiwu/*.output $DEST
mv /gpfs/mira-home/qiwu/*.cobaltlog $DEST
mv /gpfs/mira-home/qiwu/*.engine_par.*.vlog $DEST
mv /gpfs/mira-home/qiwu/*.engine_ser.*.vlog $DEST
mv /gpfs/mira-home/qiwu/engine_par.*.timings $DEST
cd /gpfs/mira-home/qiwu
