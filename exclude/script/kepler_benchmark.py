import os

# open remote
m = GetMachineProfile("wopr.sci.utah.edu")
m.GetLaunchProfiles(0).numProcessors = int(os.environ['VISIT_NODES']) * int(os.environ['VISIT_THREADS'])
m.GetLaunchProfiles(0).numNodes = int(os.environ['VISIT_NODES'])
m.GetLaunchProfiles(0).sublaunchPreCmdSet = True
m.GetLaunchProfiles(0).sublaunchPreCmd = "cd " + os.environ['VISIT_WORKDIR']
OpenComputeEngine(m)
OpenDatabase('wopr.sci.utah.edu:~/data/Brown_MPM_MasterUda.uda/index.xml')

# add plot
AddPlot('Volume','w_qn0/1')

# add plot attribute
p = VolumeAttributes()
p.freeformOpacity = (0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 19, 36, 63, 71, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173)
p.rendererType = 6
p.opacityAttenuation = 0.3529412
SetPlotOptions(p)

# draw plot
DrawPlots()

# cameras
c = [GetView3D(), GetView3D(), GetView3D(), GetView3D(), GetView3D(), GetView3D()]
# side views
c[0].viewNormal = (0, 0, 1)
c[0].viewUp = (0, 1, 0)
c[1].viewNormal = (0, 1, 0)
c[1].viewUp = (0, 0, -1)
c[2].viewNormal = (0, 0,-1)
c[2].viewUp = (0,-1, 0)
c[3].viewNormal = (0,-1, 0)
c[3].viewUp = (0, 0, 1)
# front/back views
c[4].viewNormal = ( 1, 0, 0)
c[4].viewUp = (0, 1, 0)
c[5].viewNormal = (-1, 0, 0)
c[5].viewUp = (0, 1, 0)

# 10 side views
for i in range(10):
    SetView3D(c[i % 4])
    SaveWindow()

# 10 front/back views
for i in range(10):
    SetView3D(c[i % 2 + 4])
    SaveWindow()

# exit program nicely
DeleteAllPlots()
CloseDatabase('wopr.sci.utah.edu:~/data/Brown_MPM_MasterUda.uda/index.xml')
CloseComputeEngine('wopr.sci.utah.edu')
exit()
