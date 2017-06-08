import os
# open remote
m = GetMachineProfile("wopr.sci.utah.edu")

m.GetLaunchProfiles(0).numProcessors = 64
m.GetLaunchProfiles(0).numNodes = 4

# m.sshCommandSpecified = True
# m.sshCommand = ("ssh", "wopr.sci.utah.edu", "cd ~/Desktop/timing/05-28/test && bash")

m.GetLaunchProfiles(0).sublaunchPreCmdSet = True
m.GetLaunchProfiles(0).sublaunchPreCmd = "cd " + os.environ['VISIT_WORKDIR']

print m

OpenComputeEngine(m)
OpenDatabase("wopr.sci.utah.edu:~/data/Brown_MPM_MasterUda.uda/index.xml")

# add plot
AddPlot("Volume","w_qn0/1")

# add plot attribute
p = VolumeAttributes()
p.freeformOpacity = (0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 19, 36, 63, 71, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173)
p.rendererType = 6
p.opacityAttenuation = 0.3529412
SetPlotOptions(p)

# draw plot
DrawPlots()

# save screen shot
SaveWindow()

# end session nicely
# DeleteAllPlots()
CloseDatabase("wopr.sci.utah.edu:~/data/Brown_MPM_MasterUda.uda/index.xml")
CloseComputeEngine("wopr.sci.utah.edu")
exit()
