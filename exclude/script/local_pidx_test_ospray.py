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
p.opacityAttenuation = 0.03
p.rendererType = 6
SetPlotOptions(p)
DrawPlots()
