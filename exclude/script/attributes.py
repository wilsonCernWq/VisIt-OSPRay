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

def getOSPRayAtts():
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
    # Splatting, Texture3D, RayCasting, RayCastingIntegration, SLIVR, RayCastingSLIVR, OSPRaySLIVR, Tuvok
    VolumeAtts.rendererType = VolumeAtts.OSPRaySLIVR  
    VolumeAtts.gradientType = VolumeAtts.CenteredDifferences  # CenteredDifferences, SobelOperator
    VolumeAtts.num3DSlices = 200
    VolumeAtts.scaling = VolumeAtts.Linear  # Linear, Log, Skew
    VolumeAtts.skewFactor = 1
    VolumeAtts.limitsMode = VolumeAtts.OriginalData  # OriginalData, CurrentPlot
    VolumeAtts.sampling = VolumeAtts.Rasterization  # KernelBased, Rasterization, Trilinear
    VolumeAtts.rendererSamples = 3
    #transferFunction2DWidgets does not contain any TransferFunctionWidget objects.
    VolumeAtts.transferFunctionDim = 1
    VolumeAtts.lowGradientLightingReduction = VolumeAtts.Lower  # Off, Lowest, Lower, Low, Medium, High, Higher, Highest
    VolumeAtts.lowGradientLightingClampFlag = 0
    VolumeAtts.lowGradientLightingClampValue = 1
    VolumeAtts.materialProperties = (0.4, 0.75, 0, 15)

def getOSPRayAtts():
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
    # Splatting, Texture3D, RayCasting, RayCastingIntegration, SLIVR, RayCastingSLIVR, OSPRaySLIVR, Tuvok
    VolumeAtts.rendererType = VolumeAtts.OSPRaySLIVR  
    VolumeAtts.gradientType = VolumeAtts.CenteredDifferences  # CenteredDifferences, SobelOperator
    VolumeAtts.num3DSlices = 200
    VolumeAtts.scaling = VolumeAtts.Linear  # Linear, Log, Skew
    VolumeAtts.skewFactor = 1
    VolumeAtts.limitsMode = VolumeAtts.OriginalData  # OriginalData, CurrentPlot
    VolumeAtts.sampling = VolumeAtts.Rasterization  # KernelBased, Rasterization, Trilinear
    VolumeAtts.rendererSamples = 3
    #transferFunction2DWidgets does not contain any TransferFunctionWidget objects.
    VolumeAtts.transferFunctionDim = 1
    VolumeAtts.lowGradientLightingReduction = VolumeAtts.Lower  # Off, Lowest, Lower, Low, Medium, High, Higher, Highest
    VolumeAtts.lowGradientLightingClampFlag = 0
    VolumeAtts.lowGradientLightingClampValue = 1
    VolumeAtts.materialProperties = (0.4, 0.75, 0, 15)
