import os
from subprocess import call

#-----------------------------------------------------------------------------
server_path = "./"
client_path = "./"
datainfo = {
    'HOSTNAME': "localhost",
    'FULLPATH': "/work/04915/qiwu/stampede2/data/Brown_MPM_MasterUda.uda/index.xml",
    'TIMESTEP': 0,
    'VARIABLE': "w_qn0/1"
}
#-----------------------------------------------------------------------------
def makePlot(atts, useOSPRay = True, usePascal = True, useDefault = True):
    # open dataset
    OpenDatabase(datainfo['HOSTNAME'] + ":" + datainfo['FULLPATH'])
    SetTimeSliderState(datainfo['TIMESTEP'])
    AddPlot("Volume", datainfo['VARIABLE'])
    def drawPlots(VolumeAtts, VolumeType):
        # Splatting, Texture3D, RayCasting, RayCastingIntegration, SLIVR
        # RayCastingSLIVR, OSPRaySLIVR, Tuvok
        print "drawing volume type: " + str(VolumeType)
        VolumeAtts.rendererType = VolumeType
        SetPlotOptions(VolumeAtts)
        DrawPlots()
        SaveWindow()
        # camera positions
        c = [GetView3D(), GetView3D(), GetView3D(), 
             GetView3D(), GetView3D(), GetView3D()]
        # side views
        c[0].viewNormal = (0, 1, 0)
        c[0].viewUp = (0, 0, -1)
        c[1].viewNormal = (0, 0,-1)
        c[1].viewUp = (0,-1, 0)
        c[2].viewNormal = (0,-1, 0)
        c[2].viewUp = (0, 0, 1)
        c[3].viewNormal = (0, 0, 1)
        c[3].viewUp = (0, 1, 0)
        # front/back views
        c[4].viewNormal = ( 1, 0, 0)
        c[4].viewUp = (0, 1, 0)
        c[5].viewNormal = (-1, 0, 0)
        c[5].viewUp = (0, 1, 0)
        # N front/back views
        for i in range(2):
            SetView3D(c[i % 2 + 4])
            DrawPlots()
            SaveWindow()
        # N side views
        for i in range(4):
            SetView3D(c[i % 4])
            DrawPlots()
            SaveWindow()
    # do plots
    if (useOSPRay):
        drawPlots(atts, atts.OSPRaySLIVR)
    if (usePascal):
        drawPlots(atts, atts.RayCastingSLIVR)
    if (useDefault):
        drawPlots(atts, atts.RayCasting)
    # close all
    DeleteActivePlots()
    CloseDatabase(datainfo['HOSTNAME'] + ":" + datainfo['FULLPATH'])

#-----------------------------------------------------------------------------------
VolumeAtts = VolumeAttributes()
VolumeAtts.freeformOpacity = (0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 19, 36, 63, 71, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173)
VolumeAtts.opacityAttenuation = 0.3529412
VolumeAtts.sampling = VolumeAtts.Trilinear  # KernelBased, Rasterization, Trilinear
VolumeAtts.rendererSamples = 3
#------------------------------------------------------------------------------------
# open remote
makePlot(VolumeAtts, False, True, False)
exit()

