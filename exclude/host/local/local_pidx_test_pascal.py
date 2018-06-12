# Set Default Option
opt =  GetDefaultFileOpenOptions("IDX")
opt['Big Endian'] = 0
opt['Use RAW format'] = 1
opt['Use extra cells'] = 0
SetDefaultFileOpenOptions("IDX", opt)

# open dataset
OpenDatabase('~/data/RawDataExample_IDX/test_raw.idx')

# setup plot
AddPlot('Volume','variable_0')
p = VolumeAttributes()
p.freeformOpacity = (0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173, 195, 173, 151, 130, 108, 86, 65, 43, 21, 0, 0, 21, 43, 65, 86, 108, 130, 151, 173)
# p.opacityAttenuation = 0.2156863
p.rendererType = 5
SetPlotOptions(p)
DrawPlots()
sw = SaveWindowAttributes()
sw.height = 4096
sw.width  = 4096
SetSaveWindowAttributes(sw)
SaveWindow()