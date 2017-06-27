def makePlot(atts):
    # open dataset
    OpenDatabase('~/data/Brown_MPM_MasterUda.uda/index.xml')
    # setup plot
    AddPlot('Volume','w_qn0/1')
    def drawPlots(VolumeAtts, VolumeType):
        # Splatting, Texture3D, RayCasting, RayCastingIntegration, SLIVR
        # RayCastingSLIVR, OSPRaySLIVR, Tuvok
        print "drawing volume type: " + str(VolumeType)
        VolumeAtts.rendererType = VolumeType
        SetPlotOptions(VolumeAtts)
        DrawPlots()
        # camera positions
        c = [GetView3D(), GetView3D(), GetView3D(), GetView3D(), GetView3D(), GetView3D()]
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
        # N side views
        for i in range(4):
            SetView3D(c[i % 4])
            DrawPlots()
    # do plots
    drawPlots(atts, atts.OSPRaySLIVR)
    drawPlots(atts, atts.RayCastingSLIVR)
    drawPlots(atts, atts.RayCasting)
    #drawPlots(atts, atts.RayCastingIntegration)
    # close all
    DeleteActivePlots()
    CloseDatabase('~/data/Brown_MPM_MasterUda.uda/index.xml')
    CloseComputeEngine()

#--------------------------------------------------------------------------------------------------------
VolumeAtts = VolumeAttributes()
VolumeAtts.freeformOpacity = (0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 19, 36, 63, 71, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173)
VolumeAtts.opacityAttenuation = 0.3529412
VolumeAtts.sampling = VolumeAtts.Trilinear  # KernelBased, Rasterization, Trilinear
VolumeAtts.rendererSamples = 3
#--------------------------------------------------------------------------------------------------------
# open remote
makePlot(VolumeAtts)
exit()

