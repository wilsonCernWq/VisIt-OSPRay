import os
import math
from subprocess import call

#-----------------------------------------------------------------------------
server_path = "/home/sci/qwu/Desktop/timing/tmp/Kepler"
client_path = "./"
datainfo = {
    'HOSTNAME': "wopr.sci.utah.edu",
    'FULLPATH': "/usr/sci/cedmav/data/pidx_uintah/CCVars.idx",
    'TIMESTEP': 229829,
    'VARIABLE': "O2"
}
cmd_prefix = "/home/sci/qwu/VisIt/working/exclude/tools/kepler"
cmd_enter = "source " + cmd_prefix + "/enterVisItJobs.sh " + server_path + "/"
cmd_exit  = "source " + cmd_prefix + "/exitVisItJobs.sh "  + server_path + "/"
#------------------------------------------------------------------------------
# functions
def makeColorControlPoint(color, position):
    cPoint = ColorControlPoint()
    cPoint.colors = color
    cPoint.position = position
    return cPoint

def makeOpacityControlPoint(x, height, width, xBias, yBias):
    oPoint = GaussianControlPoint()
    oPoint.x = x
    oPoint.height = height
    oPoint.width = width
    oPoint.xBias = xBias
    oPoint.yBias = yBias
    return oPoint

def makePlot(machine, atts, numThreads, numNodes, \
             useOSPRay = True, usePascal = True, useDefault = True):
    dirpath = client_path + "n" + str(numNodes) + "p" + str(numThreads)
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)

    machine.GetLaunchProfiles(0).numProcessors = numThreads * numNodes
    machine.GetLaunchProfiles(0).numNodes = numNodes
    machine.GetLaunchProfiles(0).sublaunchPreCmdSet = True
    machine.GetLaunchProfiles(0).sublaunchPreCmd  = cmd_enter + dirpath
    machine.GetLaunchProfiles(0).sublaunchPostCmdSet = True
    machine.GetLaunchProfiles(0).sublaunchPostCmd = cmd_exit  + dirpath

    OpenComputeEngine(machine)
    OpenDatabase(datainfo['HOSTNAME'] + ":" + datainfo['FULLPATH'])
    SetTimeSliderState(datainfo['TIMESTEP'])
    AddPlot("Volume", datainfo['VARIABLE'])
    def drawPlots(VolumeAtts, VolumeType):
        # Splatting, Texture3D, RayCasting, RayCastingIntegration
        # SLIVR, RayCastingSLIVR, RayCastingOSPRay, Tuvok
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
        drawPlots(atts, atts.RayCastingOSPRay)
    if (usePascal):
        drawPlots(atts, atts.RayCastingSLIVR)
    if (useDefault):
        drawPlots(atts, atts.RayCasting)
    # close all
    DeleteActivePlots()
    CloseDatabase(hostname + ":" + database)
    CloseComputeEngine(hostname)
    # clean up data
    call("mv visit*.png " + dirpath, shell=True)

#------------------------------------------------------------------------------
# Set Default Option
opt = GetDefaultFileOpenOptions("IDX")
opt['Big Endian'] = 0
opt['Use RAW format'] = 1
opt['Use extra cells'] = 1
SetDefaultFileOpenOptions("IDX", opt)

#------------------------------------------------------------------------------
# setup VolumeAttribute
# set TF
VolumeAtts = VolumeAttributes()
VolumeAtts.legendFlag = 1
VolumeAtts.lightingFlag = 0
# remove all color control points
for _ in range(5):
    VolumeAtts.colorControlPoints.RemoveControlPoints(0)
VolumeAtts.colorControlPoints.AddControlPoints(makeColorControlPoint((0, 0, 128, 255), 0))
VolumeAtts.colorControlPoints.AddControlPoints(makeColorControlPoint((0, 0, 255, 255), 0.158661))
VolumeAtts.colorControlPoints.AddControlPoints(makeColorControlPoint((0, 255, 255, 255), 0.339156))
VolumeAtts.colorControlPoints.AddControlPoints(makeColorControlPoint((0, 255, 0, 255), 0.509461))
VolumeAtts.colorControlPoints.AddControlPoints(makeColorControlPoint((255, 255, 0, 255), 0.681223))
VolumeAtts.colorControlPoints.AddControlPoints(makeColorControlPoint((255, 0, 0, 255), 1))
VolumeAtts.colorControlPoints.smoothing = VolumeAtts.colorControlPoints.Linear  # None, Linear, CubicSpline
VolumeAtts.colorControlPoints.equalSpacingFlag = 0
VolumeAtts.colorControlPoints.discreteFlag = 0
VolumeAtts.colorControlPoints.categoryName = ""
# opacity
VolumeAtts.opacityAttenuation = 0.215686
VolumeAtts.opacityMode = VolumeAtts.GaussianMode  # FreeformMode, GaussianMode, ColorTableMode
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.05, 0.25, 0.05, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.15, 0.35, 0.0622093, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.25, 0.366667, 0.0581395, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.35, 0.516667, 0.0694768, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.45, 0.6, 0.0677326, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.55, 0.616667, 0.0674419, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.652616, 0.516667, 0.0726745, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.803779, 0.483333, 0.117733, 0, 0))
VolumeAtts.opacityControlPoints.AddControlPoints(makeOpacityControlPoint(0.95, 0.4, 0.0604652, 0, 0))
VolumeAtts.resampleFlag = 1
VolumeAtts.resampleTarget = 50000
VolumeAtts.opacityVariable = "default"
VolumeAtts.compactVariable = "default"
VolumeAtts.useColorVarMin = 1
VolumeAtts.colorVarMin = 0
VolumeAtts.useColorVarMax = 1
VolumeAtts.colorVarMax = 0.15
VolumeAtts.useOpacityVarMin = 0
VolumeAtts.opacityVarMin = 0
VolumeAtts.useOpacityVarMax = 0
VolumeAtts.opacityVarMax = 0
VolumeAtts.smoothData = 0
VolumeAtts.samplesPerRay = 2
# Splatting, Texture3D, RayCasting, RayCastingIntegration, SLIVR, RayCastingSLIVR, RayCastingOSPRay, Tuvok
VolumeAtts.rendererType = VolumeAtts.RayCastingOSPRay  
VolumeAtts.gradientType = VolumeAtts.CenteredDifferences  # CenteredDifferences, SobelOperator
VolumeAtts.num3DSlices = 200
VolumeAtts.scaling = VolumeAtts.Linear  # Linear, Log, Skew
VolumeAtts.skewFactor = 1
VolumeAtts.limitsMode = VolumeAtts.OriginalData  # OriginalData, CurrentPlot
VolumeAtts.sampling = VolumeAtts.Trilinear  # KernelBased, Rasterization, Trilinear
VolumeAtts.rendererSamples = 3
#transferFunction2DWidgets does not contain any TransferFunctionWidget objects.
VolumeAtts.transferFunctionDim = 1
VolumeAtts.lowGradientLightingReduction = VolumeAtts.Lower  # Off, Lowest, Lower, Low, Medium, High, Higher, Highest
VolumeAtts.lowGradientLightingClampFlag = 0
VolumeAtts.lowGradientLightingClampValue = 1
VolumeAtts.materialProperties = (0.4, 0.75, 0, 15)

#------------------------------------------------------------------------------
# open remote
m = GetMachineProfile(datainfo['HOSTNAME'])
makePlot(m, VolumeAtts, 1, 16, True, False, False)
call("mv *.vlog *.timings " + server_path, shell=True)
exit()
