#!/bin/bash
source /home/sci/qwu/software/Kepler/VisIt/visit-depend.sh activate
DEST=$1
mkdir -p $DEST
ssh qwu@wopr.sci.utah.edu mv /tmp/visit.qiwu.* $DEST
cd ${DEST}

