#!/bin/bash

cd /home/sci/qwu/visitOSPRay/visitOSPRayCPU/pascal

cp -r plots/ ../src-visit2.12.0/
cp -r avt/Filters/avt* ../src-visit2.12.0/avt/Filters
cp -r avt/Filters/img* ../src-visit2.12.0/avt/Filters
cp -r avt/Pipeline/* ../src-visit2.12.0/avt/Pipeline
