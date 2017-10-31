import os
import math
from subprocess import call

#-----------------------------------------------------------------------------
server_path = "./"
client_path = "./"
datainfo = {
    'HOSTNAME': "localhost",
    'FULLPATH': "/scratch/04915/qiwu/runs/weak/data-2/GradientVar_0.idx",
    'TIMESTEP': 0,
    'VARIABLE': "var_0"
}
#-----------------------------------------------------------------------------
# functions
def makePlot(atts, useOSPRay = True, usePascal = True, useDefault = True):
    OpenDatabase(datainfo['HOSTNAME'] + ":" + datainfo['FULLPATH'])
    SetTimeSliderState(datainfo['TIMESTEP'])
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
        N = 100
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
makePlot(VolumeAtts, True, False, False)
exit()
