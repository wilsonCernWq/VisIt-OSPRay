#!/bin/bash

NT=$1 # number of tasks
NN=$2 # number of node
NC=$3 # number of cores per task
NV=$((NN/1))
NAME=n${NN}p$((NT / NN))

mkdir -p $NAME
cd $NAME

EXE=$WORK/software/VisIt/stampede2/visit-trunk/build/bin/visit
SCP=$(pwd)/input-$NAME.py

cat > input-$NAME.py <<EOF
import os
import math
from subprocess import call

server_path = "./"
client_path = "./"
datainfo = {
    'HOSTNAME': "localhost",
    'FULLPATH': "/scratch/04915/qiwu/runs/weak/data/volume-${NV}/GradientVar.idx",
    'TIMESTEP': 0,
    'VARIABLE': "var_0"
}

def makePlot(atts, useOSPRay = True, usePascal = True, useDefault = True):
    OpenDatabase(datainfo['HOSTNAME'] + ":" + datainfo['FULLPATH'])
    AddPlot("Volume", datainfo['VARIABLE'])
    def drawPlots(VolumeAtts, VolumeType):
        # Splatting, Texture3D, RayCasting, RayCastingIntegration
        # SLIVR, RayCastingSLIVR, OSPRaySLIVR, Tuvok
        print "drawing volume type: " + str(VolumeType)
        VolumeAtts.rendererType = VolumeType
        SetPlotOptions(VolumeAtts)
        DrawPlots()
        SaveWindow()
        # camera positions
        c = GetView3D()
        # front/back views
        N = 10
        for i in range(N):
            angle = float(i) / float(N) * 2 * math.pi
            cc = c
            cc.viewNormal = (0, math.sin(angle),  math.cos(angle))
            cc.viewUp     = (0, math.cos(angle), -math.sin(angle))
            SetView3D(cc)
            DrawPlots()
            SaveWindow()
        for i in range(N):
            angle = float(i) / float(N) * 2 * math.pi
            cc = c
            cc.viewNormal = (math.cos(angle), 0, math.sin(angle))
            cc.viewUp     = (0, 1, 0)
            SetView3D(cc)
            DrawPlots()
            SaveWindow()
    if (useOSPRay):
        drawPlots(atts, atts.OSPRaySLIVR)
    if (usePascal):
        drawPlots(atts, atts.RayCastingSLIVR)
    if (useDefault):
        drawPlots(atts, atts.RayCasting)
    # close all
    DeleteActivePlots()
    CloseDatabase(datainfo['HOSTNAME'] + ":" + datainfo['FULLPATH'])

#-----------------------------------------------------------------------------
# Set Default Option
opt = GetDefaultFileOpenOptions("IDX")
opt['Big Endian'] = 0
opt['Use RAW format'] = 1
opt['Use extra cells'] = 1
SetDefaultFileOpenOptions("IDX", opt)

#-----------------------------------------------------------------------------
# setup VolumeAttribute
# set TF
VolumeAtts = VolumeAttributes()
VolumeAtts.lightingFlag = 0
VolumeAtts.rendererType = VolumeAtts.OSPRaySLIVR  
VolumeAtts.sampling = VolumeAtts.Trilinear  # KernelBased, Rasterization, Trilinear
VolumeAtts.samplesPerRay = 2
VolumeAtts.rendererSamples = 3

#-----------------------------------------------------------------------------
# open remote
makePlot(VolumeAtts, False, True, False)
exit()
EOF

cat > sample-$NAME.sh <<EOF
#!/bin/bash
#SBATCH -n ${NT}
#SBATCH -N ${NN}
#SBATCH --time=00:30:00 # walltime, abbreviated by -t
#
# load modules
module load intel
module load impi
module load hdf5
export LIBRARY_PATH=${LIBRARY_PATH}:/work/04915/qiwu/stampede2/software/VisIt/stampede2/visit/vtk/6.1.0/linux-x86_64_gcc-5.4/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/work/04915/qiwu/stampede2/software/VisIt/stampede2/visit/vtk/6.1.0/linux-x86_64_gcc-5.4/lib
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=$NC
export KMP_AFFINITY=verbose,granularity=core,compact,1,0
export OSPRAY_THREADS=$NC
#
# run job
echo "job: <\$SLURM_JOBID> node = $NN processes-per-node = $((NT / NN))"
$EXE -np $NT -nn $NN -l ibrun \
    -withhold-timing-output -timing \
    -nowin -cli -s $SCP 
EOF
sbatch -p normal -A TG-ASC170049 sample-$NAME.sh
cd -
