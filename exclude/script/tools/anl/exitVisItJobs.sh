#!/bin/bash
DEST=$1
mv ~/*.error $dest
mv ~/*.output $dest
mv ~/*.cobaltlog $dest
mv ~/*.engine_par.*.vlog $dest
mv ~/*.engine_ser.*.vlog $dest
mv ~/engine_par.*.timings $dest
cd ~
