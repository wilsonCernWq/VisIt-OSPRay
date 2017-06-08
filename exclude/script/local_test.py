OpenDatabase('~/data/Brown_MPM_MasterUda.uda/index.xml')
AddPlot('Volume','w_qn0/1')
p = VolumeAttributes()
SaveAttribute('~/VisIt/rcslivr_attribute_save.xml', p)
# print p
LoadAttribute('~/VisIt/rcslivr_attribute.xml', p)
# print p
SetPlotOptions(p)
DrawPlots()
