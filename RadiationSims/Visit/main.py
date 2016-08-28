import Image
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/Utilities")
from sigmaprojectionplot import sigmaprojectionplot

## Use full vtk dump to plot an x-y slice of the density, with contours of radiation energy density overplotted as well as arrows for the flux
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def singleplot(pid, savpos, pheight, filedir):
    
    rhocodeast = 29.06
    conv = 28.9 / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    DefineScalarExpression("nhden", "density*28.9/29.06")
    # Create a plot of the scalar field 'density'
    AddPlot("Pseudocolor","nhden")
    
    # Slice the volume to show only the central slice
    AddOperator("Slice")
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topright":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 100000
    psa.colorTableName = "BuGn"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    sa = SliceAttributes()
    sa.originType = sa.Intercept  # Point, Intercept, Percent, Zone, Node
    sa.originIntercept = 0.0
    sa.axisType = sa.ZAxis
    SetOperatorOptions(sa, 1)
    
    AddPlot("Contour", "rad_energy_density", 1, 1)
    ca = ContourAttributes()
    ca.colorType = ca.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    ca.singleColor = (255, 0, 255, 255)
    ca.contourNLevels = 5
    ca.legendFlag = 0
    ca.minFlag = 1
    ca.maxFlag = 1
    ca.min = 10
    ca.max = 100000
    ca.scaling = ca.Log  # Linear, Log
    SetPlotOptions(ca)
    
    AddPlot("Vector", "rad_flux", 1, 1)
    va = VectorAttributes()
    va.useStride = 0
    va.nVectors = 400
    va.scale = 0.1
    va.scaleByMagnitude = 0
    va.autoScale = 1
    va.headSize = 0.4
    va.useLegend = 0
    va.lineWidth = 2
    va.colorByMag = 0
    va.vectorColor = (255, 255, 0, 255)
    va.vectorOrigin = va.Tail  # Head, Middle, Tail
    SetPlotOptions(va)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 1
    denleg.numTicks = 8
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.02
    
    DrawPlots()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationSlice" + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use projected Sigma vtk to plot the surface density, with star particles
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., axproj: direction in which projection is done e.g. xy, xz or yz, binsize: resolution of simulation e.g. 128, 256, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def sigmaprojectionplot(pid, savpos, axproj, binsize, pheight, filedir):
    
    rhocodeast = 29.06
    dxone = 60
    conv = dxone / rhocodeast
    # Clear any previous plots
    # DeleteAllPlots()
    
    # Open the vtk database
    if axproj == "xy":
        sigtxt = "Sigma3"
    if axproj == "xz":
        sigtxt = "Sigma2"
    if axproj == "yz":
        sigtxt = "Sigma1"
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + sigtxt + "." + pid + ".vtk", 0)
    
    # Create a plot of the scalar field 'density'
    DefineScalarExpression("column", sigtxt + "/29.06")
    AddPlot("Pseudocolor","column")
    
    # Bin the data to get column density
    AddOperator("DataBinning", 0)
    dbas = DataBinningAttributes()
    dbas.numDimensions = dbas.Three  # One, Two, Three
    dbas.dim1BinBasedOn = dbas.X  # X, Y, Z, Variable
    dbas.dim1NumBins = binsize
    dbas.dim2BinBasedOn = dbas.Y  # X, Y, Z, Variable
    dbas.dim2NumBins = binsize
    dbas.dim3BinBasedOn = dbas.Z  # X, Y, Z, Variable
    dbas.dim3NumBins = binsize
    dbas.reductionOperator = dbas.Average  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
    dbas.varForReduction = "column"
    if axproj == "xy":
        dbas.dim3NumBins = 1
    if axproj == "xz":
        dbas.dim2NumBins = 1
    if axproj == "yz":
        dbas.dim1NumBins = 1
    dbas.outputType = dbas.OutputOnBins
    SetOperatorOptions(dbas, 1)
    
    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    if axproj == "xy":
        pas.projectionType = pas.XYCartesian
    if axproj == "xz":
        pas.projectionType = pas.XZCartesian
    if axproj == "yz":
        pas.projectionType = pas.ZYCartesian
    SetOperatorOptions(pas, 1)
    
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 1000
    psa.colorTableName = "Spectral"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    # Open the star particle database
    OpenDatabase("localhost:" + filedir + "/RadParGrav.starpar." + pid + ".vtk", 0)
    
    # Plot the stellar masses
    DefineScalarExpression("starmass", "star_particle_mass / 29.06")
    AddPlot("Pseudocolor", "starmass", 1, 0)
    starpsa = PseudocolorAttributes()
    if savpos == "bottomright":
        starpsa.legendFlag = 1
    else:
        starpsa.legendFlag = 0
    starpsa.scaling = starpsa.Log  # Linear, Log, Skew
    starpsa.pointSize = 0.05
    starpsa.pointType = starpsa.Sphere  # Box, Axis, Icosahedron, Point, Sphere
    starpsa.colorTableName = "RdPu"
    starpsa.pointSizePixels = 30
    starpsa.minFlag = 1
    starpsa.maxFlag = 1
    starpsa.min = 50
    starpsa.max = 5000
    SetPlotOptions(starpsa)
    
    AddOperator("Project", 1)
    starpa = ProjectAttributes()
    if axproj == "xy":
        starpa.projectionType = starpa.XYCartesian
    if axproj == "xz":
        starpa.projectionType = starpa.XZCartesian
    if axproj == "yz":
        starpa.projectionType = starpa.ZYCartesian
    SetOperatorOptions(starpa, 1)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 1
    denleg.numTicks = 6
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.02
    
    starleg = GetAnnotationObject(alllegs[1])
    starleg.drawTitle = 0
    starleg.drawMinMax = 0
    starleg.managePosition = 0
    starleg.position = (0.01, 0.92)
    starleg.numberFormat = "%# -9.1g"
    starleg.controlTicks = 1
    starleg.numTicks = 3
    starleg.xScale = 1
    starleg.yScale = 3
    starleg.fontFamily = starleg.Times
    starleg.fontBold = 1
    starleg.fontHeight = 0.02
    
    DrawPlots()
    ResetView()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationColumn" + axproj + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to show scatter plot of velocity as a function of cell density
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, pheight: height of output plot in pixels, default 1024, vtype: type of velocity to plot e.g. vturb or va (for alfven speed) filedir: directory in which input vtk dump is stored and in which output png is plotted
def magscatterplot(pid, pheight, vtype, filedir):
    
    rhocodeast = 29.06
    conv = 28.9 / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    DefineScalarExpression("nhden", "density*28.9/29.06")
    DefineScalarExpression("vturb", "momentum_magnitude/density")
    DefineScalarExpression("va", "cell_centered_B_magnitude/sqrt(4.0*3.14159*density)")
    
    AddPlot("Scatter","nhden")
    sa = ScatterAttributes()
    
    sa.var1 = "nhden"
    sa.var1Role = sa.Coordinate0  # Coordinate0, Coordinate1, Coordinate2, Color, None
    sa.var1MinFlag = 1
    sa.var1MaxFlag = 1
    sa.var1Min = 0.01
    sa.var1Max = 100000
    sa.var1Scaling = sa.Log  # Linear, Log, Skew
    sa.var1SkewFactor = 1

    sa.var2Role = sa.Coordinate1  # Coordinate0, Coordinate1, Coordinate2, Color, None
    sa.var2 = vtype
    sa.var2MinFlag = 1
    sa.var2MaxFlag = 1
    sa.var2Min = 0.01
    sa.var2Max = 100
    sa.var2Scaling = sa.Log  # Linear, Log, Skew
    sa.var2SkewFactor = 1
    
    sa.scaleCube = 1
    sa.pointSize = 0.05
    sa.pointSizePixels = 1
    sa.legendFlag = 1
    sa.pointType = sa.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    
    SetPlotOptions(sa)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 1.5
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.visible = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 1.5
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    aa.legendInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    DrawPlots()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = vtype + "_Scatter" + "_T" + pid
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a slice of the density, with magnetic field overplotted as yellow arrows
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., axproj: direction which is normal to plane e.g. xy, xz or yz, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def singlemagplot(pid, savpos, axproj, pheight, filedir):
    
    rhocodeast = 29.06
    conv = 28.9 / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    DefineScalarExpression("nhden", "density*28.9/29.06")
    # Create a plot of the scalar field 'density'
    AddPlot("Pseudocolor","nhden")
    
    # Slice the volume to show only the central slice
    AddOperator("Slice")
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topright":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 100000
    psa.colorTableName = "BuGn"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    sa = SliceAttributes()
    sa.originType = sa.Intercept  # Point, Intercept, Percent, Zone, Node
    sa.originIntercept = 0.0
    if axproj == "xy":
        sa.originIntercept = 0.0
        sa.axisType = sa.ZAxis
    if axproj == "xz":
        sa.originIntercept = 0.0
        sa.axisType = sa.YAxis
    if axproj == "yz":
        sa.originIntercept = 0.0
        sa.axisType = sa.XAxis
    SetOperatorOptions(sa, 1)
    
    AddPlot("Vector", "cell_centered_B", 1, 1)
    va = VectorAttributes()
    va.useStride = 0
    va.nVectors = 400
    va.scale = 0.1
    va.scaleByMagnitude = 0
    va.autoScale = 1
    va.headSize = 0.4
    va.useLegend = 0
    va.lineWidth = 2
    va.colorByMag = 0
    va.vectorColor = (255, 255, 0, 255)
    va.vectorOrigin = va.Tail  # Head, Middle, Tail
    SetPlotOptions(va)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 1.5
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 1.5
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 4
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.03
    denleg.suppliedValues = (0.01,1,100,10000)
    denleg.suppliedLabels = ("0.01","1","10^2","10^4")
    denleg.drawLabels = denleg.Labels

    DrawPlots()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationSlice" + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a slice of the density, with magnetic field overplotted as yellow arrows and star particles
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., axproj: direction which is normal to plane e.g. xy, xz or yz, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def singlestarmagplot(pid, savpos, pheight, axproj, filedir):
    
    rhocodeast = 29.06
    conv = 28.9 / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the star particle database
    OpenDatabase("localhost:" + filedir + "/RadParGrav.starpar." + pid + ".vtk", 0)
    
    # Plot the stellar masses
    DefineScalarExpression("starmass", "star_particle_mass / 29.06")
    AddPlot("Pseudocolor", "starmass", 1, 0)
    DrawPlots()
    
    # Get the position of the maximum star particle mass
    Query("Max")
    import string
    s = string.split(GetQueryOutputString(), " ")
    maxid = 0
    nextGood = 0
    for token in s:
        if token == "(node":
            nextGood = 1
            continue
        if nextGood == 1:
            maxid = int(token)
            break

    Query("Node Coords", domain=0, element=maxid, use_global_id=0)
    xpar = GetQueryOutputValue()[0]
    ypar = GetQueryOutputValue()[1]
    zpar = GetQueryOutputValue()[2]
    DeleteAllPlots()

    print xpar, ypar, zpar
    
    # Plot the stellar masses
    AddPlot("Pseudocolor", "starmass", 1, 0)
    starpsa = PseudocolorAttributes()
    if savpos == "bottomright":
        starpsa.legendFlag = 1
    else:
        starpsa.legendFlag = 0
    starpsa.scaling = starpsa.Log  # Linear, Log, Skew
    starpsa.pointSize = 0.05
    starpsa.pointType = starpsa.Sphere  # Box, Axis, Icosahedron, Point, Sphere
    starpsa.colorTableName = "RdPu"
    starpsa.pointSizePixels = 30
    starpsa.minFlag = 1
    starpsa.maxFlag = 1
    starpsa.min = 50
    starpsa.max = 5000
    SetPlotOptions(starpsa)
    
    AddOperator("Box", 1)
    dz = 15
    bas = BoxAttributes()
    bas.minx = -30
    bas.maxx = 30
    bas.miny = -30
    bas.maxy = 30
    bas.minz = -30
    bas.maxz = 30
    if axproj == "xy":
        bas.minz = zpar - dz
        bas.maxz = zpar + dz
    if axproj == "xz":
        bas.minz = ypar - dz
        bas.maxz = ypar + dz
    if axproj == "yz":
        bas.minz = xpar - dz
        bas.maxz = xpar + dz
    SetOperatorOptions(bas, 1)
    
    AddOperator("Project", 0)
    starpa = ProjectAttributes()
    if axproj == "xy":
        starpa.projectionType = starpa.XYCartesian
    if axproj == "xz":
        starpa.projectionType = starpa.XZCartesian
    if axproj == "yz":
        starpa.projectionType = starpa.ZYCartesian
    SetOperatorOptions(starpa, 0)
    DrawPlots()
    
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    DefineScalarExpression("nhden", "density*28.9/29.06")
    # Create a plot of the scalar field 'density'
    AddPlot("Pseudocolor","nhden")
    
    # Slice the volume to show only the central slice
    AddOperator("Slice")
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topright":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 100000
    psa.colorTableName = "BuGn"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    sa = SliceAttributes()
    sa.originType = sa.Intercept  # Point, Intercept, Percent, Zone, Node
    if axproj == "xy":
        sa.originIntercept = zpar
        sa.axisType = sa.ZAxis
    if axproj == "xz":
        sa.originIntercept = ypar
        sa.axisType = sa.YAxis
    if axproj == "yz":
        sa.originIntercept = xpar
        sa.axisType = sa.XAxis
    SetOperatorOptions(sa, 1)

    AddPlot("Vector", "cell_centered_B", 1, 1)
    va = VectorAttributes()
    va.useStride = 0
    va.nVectors = 400
    va.scale = 0.1
    va.scaleByMagnitude = 0
    va.autoScale = 1
    va.headSize = 0.4
    va.useLegend = 0
    va.lineWidth = 2
    va.colorByMag = 0
    va.vectorColor = (255, 255, 0, 255)
    va.vectorOrigin = va.Tail  # Head, Middle, Tail
    SetPlotOptions(va)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 1.5
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 1.5
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)

    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[1])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 3
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.03
    denleg.suppliedValues = (0.01,1,100,10000)
    denleg.suppliedLabels = ("0.01","1","10^2","10^4")
    denleg.drawLabels = denleg.Labels

    starleg = GetAnnotationObject(alllegs[0])
    starleg.drawTitle = 0
    starleg.drawMinMax = 0
    starleg.managePosition = 0
    starleg.position = (0.01, 0.92)
    starleg.numberFormat = "%# -9.1e"
    starleg.controlTicks = 0
    starleg.numTicks = 2
    starleg.xScale = 1
    starleg.yScale = 3
    starleg.fontFamily = starleg.Times
    starleg.fontBold = 1
    starleg.fontHeight = 0.03
    starleg.suppliedValues = (100,1000)
    starleg.suppliedLabels = ("10^2","10^3")
    starleg.drawLabels = starleg.Labels

    DrawPlots()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationSlice" + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a slice of the density, with star particles
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., axproj: direction which is normal to plane e.g. xy, xz or yz, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def singlestarnoradplot(pid, savpos, pheight, axproj, filedir):
    
    rhocodeast = 29.06
    conv = 28.9 / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the star particle database
    OpenDatabase("localhost:" + filedir + "/RadParGrav.starpar." + pid + ".vtk", 0)
    
    # Plot the stellar masses
    DefineScalarExpression("starmass", "star_particle_mass / 29.06")
    AddPlot("Pseudocolor", "starmass", 1, 0)
    DrawPlots()
    
    # Get the position of the maximum star particle mass
    Query("Max")
    import string
    s = string.split(GetQueryOutputString(), " ")
    maxid = 0
    nextGood = 0
    for token in s:
        if token == "(node":
            nextGood = 1
            continue
        if nextGood == 1:
            maxid = int(token)
            break

    Query("Node Coords", domain=0, element=maxid, use_global_id=0)
    xpar = GetQueryOutputValue()[0]
    ypar = GetQueryOutputValue()[1]
    zpar = GetQueryOutputValue()[2]
    DeleteAllPlots()

    print xpar, ypar, zpar
    
    # Plot the stellar masses
    AddPlot("Pseudocolor", "starmass", 1, 0)
    starpsa = PseudocolorAttributes()
    if savpos == "bottomright":
        starpsa.legendFlag = 1
    else:
        starpsa.legendFlag = 0
    starpsa.scaling = starpsa.Log  # Linear, Log, Skew
    starpsa.pointSize = 0.05
    starpsa.pointType = starpsa.Sphere  # Box, Axis, Icosahedron, Point, Sphere
    starpsa.colorTableName = "RdPu"
    starpsa.pointSizePixels = 30
    starpsa.minFlag = 1
    starpsa.maxFlag = 1
    starpsa.min = 50
    starpsa.max = 5000
    SetPlotOptions(starpsa)
    
    AddOperator("Box", 1)
    dz = 15
    bas = BoxAttributes()
    bas.minx = -30
    bas.maxx = 30
    bas.miny = -30
    bas.maxy = 30
    bas.minz = -30
    bas.maxz = 30
    if axproj == "xy":
        bas.minz = zpar - dz
        bas.maxz = zpar + dz
    if axproj == "xz":
        bas.minz = ypar - dz
        bas.maxz = ypar + dz
    if axproj == "yz":
        bas.minz = xpar - dz
        bas.maxz = xpar + dz
    SetOperatorOptions(bas, 1)
    
    AddOperator("Project", 0)
    starpa = ProjectAttributes()
    if axproj == "xy":
        starpa.projectionType = starpa.XYCartesian
    if axproj == "xz":
        starpa.projectionType = starpa.XZCartesian
    if axproj == "yz":
        starpa.projectionType = starpa.ZYCartesian
    SetOperatorOptions(starpa, 0)
    DrawPlots()
    
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    DefineScalarExpression("nhden", "density*28.9/29.06")
    # Create a plot of the scalar field 'density'
    AddPlot("Pseudocolor","nhden")
    
    # Slice the volume to show only the central slice
    AddOperator("Slice")
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topright":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 100000
    psa.colorTableName = "BuGn"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    sa = SliceAttributes()
    sa.originType = sa.Intercept  # Point, Intercept, Percent, Zone, Node
    if axproj == "xy":
        sa.originIntercept = zpar
        sa.axisType = sa.ZAxis
    if axproj == "xz":
        sa.originIntercept = ypar
        sa.axisType = sa.YAxis
    if axproj == "yz":
        sa.originIntercept = xpar
        sa.axisType = sa.XAxis
    SetOperatorOptions(sa, 1)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 1.5
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 1.5
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)

    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[1])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 3
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.03
    denleg.suppliedValues = (0.01,1,100,10000)
    denleg.suppliedLabels = ("0.01","1","10^2","10^4")
    denleg.drawLabels = denleg.Labels

    starleg = GetAnnotationObject(alllegs[0])
    starleg.drawTitle = 0
    starleg.drawMinMax = 0
    starleg.managePosition = 0
    starleg.position = (0.01, 0.92)
    starleg.numberFormat = "%# -9.1e"
    starleg.controlTicks = 0
    starleg.numTicks = 2
    starleg.xScale = 1
    starleg.yScale = 3
    starleg.fontFamily = starleg.Times
    starleg.fontBold = 1
    starleg.fontHeight = 0.03
    starleg.suppliedValues = (100,1000)
    starleg.suppliedLabels = ("10^2","10^3")
    starleg.drawLabels = starleg.Labels

    DrawPlots()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationSlice" + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a slice of the density, with star particles, contours of radiation energy density overplotted as well as arrows for the flux
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., axproj: direction which is normal to plane e.g. xy, xz or yz, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def singlestarplot(pid, savpos, pheight, axproj, filedir):
    
    rhocodeast = 29.06
    conv = 28.9 / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the star particle database
    OpenDatabase("localhost:" + filedir + "/RadParGrav.starpar." + pid + ".vtk", 0)
    
    # Plot the stellar masses
    DefineScalarExpression("starmass", "star_particle_mass / 29.06")
    AddPlot("Pseudocolor", "starmass", 1, 0)
    DrawPlots()
    
    # Get the position of the maximum star particle mass
    Query("Max")
    import string
    s = string.split(GetQueryOutputString(), " ")
    maxid = 0
    nextGood = 0
    for token in s:
        if token == "(node":
            nextGood = 1
            continue
        if nextGood == 1:
            maxid = int(token)
            break

    Query("Node Coords", domain=0, element=maxid, use_global_id=0)
    xpar = GetQueryOutputValue()[0]
    ypar = GetQueryOutputValue()[1]
    zpar = GetQueryOutputValue()[2]
    DeleteAllPlots()

    print xpar, ypar, zpar
    
    # Plot the stellar masses
    AddPlot("Pseudocolor", "starmass", 1, 0)
    starpsa = PseudocolorAttributes()
    if savpos == "bottomright":
        starpsa.legendFlag = 1
    else:
        starpsa.legendFlag = 0
    starpsa.scaling = starpsa.Log  # Linear, Log, Skew
    starpsa.pointSize = 0.05
    starpsa.pointType = starpsa.Sphere  # Box, Axis, Icosahedron, Point, Sphere
    starpsa.colorTableName = "RdPu"
    starpsa.pointSizePixels = 30
    starpsa.minFlag = 1
    starpsa.maxFlag = 1
    starpsa.min = 50
    starpsa.max = 5000
    SetPlotOptions(starpsa)

    dz=10
    AddOperator("Box", 1)
    bas = BoxAttributes()
    bas.minx = -30
    bas.maxx = 30
    bas.miny = -30
    bas.maxy = 30
    bas.minz = -30
    bas.maxz = 30
    if axproj == "xy":
        bas.minz = zpar - dz
        bas.maxz = zpar + dz
    if axproj == "xz":
        bas.minz = ypar - dz
        bas.maxz = ypar + dz
    if axproj == "yz":
        bas.minz = xpar - dz
        bas.maxz = xpar + dz
    SetOperatorOptions(bas, 1)
    
    AddOperator("Project", 0)
    starpa = ProjectAttributes()
    if axproj == "xy":
        starpa.projectionType = starpa.XYCartesian
    if axproj == "xz":
        starpa.projectionType = starpa.XZCartesian
    if axproj == "yz":
        starpa.projectionType = starpa.ZYCartesian
    SetOperatorOptions(starpa, 0)
    DrawPlots()
    
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    DefineScalarExpression("nhden", "density*28.9/29.06")
    # Create a plot of the scalar field 'density'
    AddPlot("Pseudocolor","nhden")
    
    # Slice the volume to show only the central slice
    AddOperator("Slice")
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topright":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 100000
    psa.colorTableName = "BuGn"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    sa = SliceAttributes()
    sa.originType = sa.Intercept  # Point, Intercept, Percent, Zone, Node
    if axproj == "xy":
        sa.originIntercept = zpar
        sa.axisType = sa.ZAxis
    if axproj == "xz":
        sa.originIntercept = ypar
        sa.axisType = sa.YAxis
    if axproj == "yz":
        sa.originIntercept = xpar
        sa.axisType = sa.XAxis
    SetOperatorOptions(sa, 1)
    
    AddPlot("Contour", "rad_energy_density", 1, 1)
    ca = ContourAttributes()
    ca.colorType = ca.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    ca.singleColor = (255, 0, 255, 255)
    ca.contourNLevels = 5
    ca.legendFlag = 0
    ca.minFlag = 1
    ca.maxFlag = 1
    ca.min = 10
    ca.max = 100000
    ca.lineWidth = 3
    ca.scaling = ca.Log  # Linear, Log
    SetPlotOptions(ca)
    
    AddPlot("Vector", "rad_flux", 1, 1)
    va = VectorAttributes()
    va.useStride = 0
    va.nVectors = 400
    va.scale = 0.1
    va.scaleByMagnitude = 0
    va.autoScale = 1
    va.headSize = 0.4
    va.useLegend = 0
    va.lineWidth = 2
    va.colorByMag = 0
    va.vectorColor = (255, 255, 0, 255)
    va.vectorOrigin = va.Tail  # Head, Middle, Tail
    SetPlotOptions(va)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 1.5
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 1.5
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)

    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[1])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 3
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.03
    denleg.suppliedValues = (0.01,1,100,10000)
    denleg.suppliedLabels = ("0.01","1","10^2","10^4")
    denleg.drawLabels = denleg.Labels

    starleg = GetAnnotationObject(alllegs[0])
    starleg.drawTitle = 0
    starleg.drawMinMax = 0
    starleg.managePosition = 0
    starleg.position = (0.01, 0.92)
    starleg.numberFormat = "%# -9.1e"
    starleg.controlTicks = 0
    starleg.numTicks = 2
    starleg.xScale = 1
    starleg.yScale = 3
    starleg.fontFamily = starleg.Times
    starleg.fontBold = 1
    starleg.fontHeight = 0.03
    starleg.suppliedValues = (100,1000)
    starleg.suppliedLabels = ("10^2","10^3")
    starleg.drawLabels = starleg.Labels

    DrawPlots()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationSlice" + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot the surface density
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., axproj: direction in which projection is done e.g. xy, xz or yz, binsize: resolution of simulation e.g. 128, 256, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def projectionplot(pid, savpos, axproj, binsize, pheight, filedir):
    
    rhocodeast = 29.06
    dxone = 60
    conv = dxone / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the vtk database
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    
    # Create a plot of the scalar field 'density'
    DefineScalarExpression("column", "density*60/29.06")
    AddPlot("Pseudocolor","column")
    
    # Bin the data to get column density
    AddOperator("DataBinning", 0)
    dbas = DataBinningAttributes()
    dbas.numDimensions = dbas.Three  # One, Two, Three
    dbas.dim1BinBasedOn = dbas.X  # X, Y, Z, Variable
    dbas.dim1NumBins = binsize
    dbas.dim2BinBasedOn = dbas.Y  # X, Y, Z, Variable
    dbas.dim2NumBins = binsize
    dbas.dim3BinBasedOn = dbas.Z  # X, Y, Z, Variable
    dbas.dim3NumBins = binsize
    dbas.reductionOperator = dbas.Average  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
    dbas.varForReduction = "column"
    if axproj == "xy":
        dbas.dim3NumBins = 1
    if axproj == "xz":
        dbas.dim2NumBins = 1
    if axproj == "yz":
        dbas.dim1NumBins = 1
    dbas.outputType = dbas.OutputOnBins
    SetOperatorOptions(dbas, 1)
    
    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    if axproj == "xy":
        pas.projectionType = pas.XYCartesian
    if axproj == "xz":
        pas.projectionType = pas.XZCartesian
    if axproj == "yz":
        pas.projectionType = pas.ZYCartesian
    SetOperatorOptions(pas, 1)
    
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 1000
    psa.colorTableName = "Spectral"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 2
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 2
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 3
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.04
    denleg.suppliedValues = (0.01,1,100)
    denleg.suppliedLabels = ("0.01","1","100")
    denleg.drawLabels = denleg.Labels
    
    DrawPlots()
    ResetView()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationColumn" + axproj + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot the surface density with star particles
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., axproj: direction in which projection is done e.g. xy, xz or yz, binsize: resolution of simulation e.g. 128, 256, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def starprojectionplot(pid, savpos, axproj, binsize, pheight, filedir):
    
    rhocodeast = 29.06
    dxone = 60
    conv = dxone / rhocodeast
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the vtk database
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    
    # Create a plot of the scalar field 'density'
    DefineScalarExpression("column", "density*60/29.06")
    AddPlot("Pseudocolor","column")
    
    # Bin the data to get column density
    AddOperator("DataBinning", 0)
    dbas = DataBinningAttributes()
    dbas.numDimensions = dbas.Three  # One, Two, Three
    dbas.dim1BinBasedOn = dbas.X  # X, Y, Z, Variable
    dbas.dim1NumBins = binsize
    dbas.dim2BinBasedOn = dbas.Y  # X, Y, Z, Variable
    dbas.dim2NumBins = binsize
    dbas.dim3BinBasedOn = dbas.Z  # X, Y, Z, Variable
    dbas.dim3NumBins = binsize
    dbas.reductionOperator = dbas.Average  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
    dbas.varForReduction = "column"
    if axproj == "xy":
        dbas.dim3NumBins = 1
    if axproj == "xz":
        dbas.dim2NumBins = 1
    if axproj == "yz":
        dbas.dim1NumBins = 1
    dbas.outputType = dbas.OutputOnBins
    SetOperatorOptions(dbas, 1)
    
    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    if axproj == "xy":
        pas.projectionType = pas.XYCartesian
    if axproj == "xz":
        pas.projectionType = pas.XZCartesian
    if axproj == "yz":
        pas.projectionType = pas.ZYCartesian
    SetOperatorOptions(pas, 1)
    
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 1000
    psa.colorTableName = "Spectral"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    # Open the star particle database
    OpenDatabase("localhost:" + filedir + "/RadParGrav.starpar." + pid + ".vtk", 0)
    
    # Plot the stellar masses
    DefineScalarExpression("starmass", "star_particle_mass / 29.06")
    
    # Plot the stellar masses
    AddPlot("Pseudocolor", "starmass", 1, 0)
    starpsa = PseudocolorAttributes()
    if savpos == "bottomright":
        starpsa.legendFlag = 1
    else:
        starpsa.legendFlag = 0
    starpsa.scaling = starpsa.Log  # Linear, Log, Skew
    starpsa.pointSize = 0.05
    starpsa.pointType = starpsa.Sphere  # Box, Axis, Icosahedron, Point, Sphere
    starpsa.colorTableName = "RdPu"
    starpsa.pointSizePixels = 30
    starpsa.minFlag = 1
    starpsa.maxFlag = 1
    starpsa.min = 50
    starpsa.max = 5000
    SetPlotOptions(starpsa)

    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    if axproj == "xy":
        pas.projectionType = pas.XYCartesian
    if axproj == "xz":
        pas.projectionType = pas.XZCartesian
    if axproj == "yz":
        pas.projectionType = pas.ZYCartesian
    SetOperatorOptions(pas, 1)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 2
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 2
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 3
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.04
    denleg.suppliedValues = (0.01,1,100)
    denleg.suppliedLabels = ("0.01","1","100")
    denleg.drawLabels = denleg.Labels
    
    DrawPlots()
    ResetView()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SimulationColumn" + axproj + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a spherical slice of the density, with radiation energy density overplotted
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., rad: string indicating the radius at which slice is taken in pc e.g. "15.0", radflt: float indicating the radius at which slice is taken in pc e.g. 15.0, cont: flag indicating whether to plot energy density as contours e.g. "contour" or not, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def singlesphericalplot(pid, savpos, rad, radflt, cont, pheight, filedir):
    
    # Clear any previous plots
    DeleteAllPlots()

    # Open the vtk database
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)

    # Define Displacement operatios
    DefineVectorExpression("c", "coords(mesh)")
    DefineScalarExpression("radius", "sqrt(c[0]^2+c[1]^2+c[2]^2)")
    DefineScalarExpression("phi", "if(gt(c[1],0.),acos(c[0]/sqrt(c[0]^2+c[1]^2)),-acos(c[0]/sqrt(c[0]^2+c[1]^2)))")
    DefineScalarExpression("theta", "acos(c[2]/radius)-0.5*3.14159")
    DefineVectorExpression("sph", "{radius, theta, phi}")
    DefineVectorExpression("sphoffset", "sph-c")
    DefineScalarExpression("zint", "sqrt(1.+cos(theta)*cos(0.5*phi))")
    DefineScalarExpression("xint", "180.*cos(theta)*sin(0.5*phi)/zint")
    DefineScalarExpression("yint", "90.*sin(theta)/zint")
    DefineVectorExpression("cart", "{radius, yint, xint}")
    DefineVectorExpression("cartoffset", "cart-c")
    DefineScalarExpression("nhden", "density*28.9/29.06")
    
    # Create a plot of the scalar field 'density'
    AddPlot("Pseudocolor","nhden")
    
    # Apply a box to exclude the line y = 0
    AddOperator("Box", 1)
    bas = BoxAttributes()
    bas.minx = -30
    bas.maxx = 30
    bas.miny = -0.001
    bas.maxy = 0.001
    bas.minz = -30
    bas.maxz = 30
    bas.inverse = 1
    SetOperatorOptions(bas, 1)
    
    # Apply a displacement operator to warp to hammer projection
    AddOperator("Displace", 1)
    das = DisplaceAttributes()
    das.factor = 1
    das.variable = "cartoffset"
    SetOperatorOptions(das, 1)

# Slice the volume to show only the central slice
    AddOperator("Slice")
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 100000
    psa.colorTableName = "bluehot"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    sa = SliceAttributes()
    sa.originType = sa.Intercept  # Point, Intercept, Percent, Zone, Node
    sa.originIntercept = radflt
    sa.axisType = sa.XAxis
    SetOperatorOptions(sa, 1)
    
    if cont == "contour":
        AddPlot("Contour", "rad_energy_density", 1, 1)
        ca = ContourAttributes()
        ca.colorType = ca.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
        ca.singleColor = (255, 0, 255, 255)
        ca.contourNLevels = 5
        ca.legendFlag = 0
        ca.minFlag = 1
        ca.maxFlag = 1
        ca.min = 10
        ca.max = 100000
        ca.scaling = ca.Log  # Linear, Log
        SetPlotOptions(ca)

# Annotation
# Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 1
    denleg.numTicks = 8
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.02
    
    DrawPlots()
    ResetView()
    
    if cont == "pseudo":
        AddPlot("Pseudocolor", "rad_energy_density", 1, 1)
        ca = PseudocolorAttributes()
        ca.scaling = ca.Log  # Linear, Log, Skew
        if savpos == "bottomleft":
            ca.legendFlag = 1
        else:
            ca.legendFlag = 0
        ca.legendFlag = 0
        ca.minFlag = 1
        ca.maxFlag = 1
        ca.min = 10
        ca.max = 100000
        ca.colorTableName = "RdPu"
        ca.invertColorTable = 0
        ca.opacityType = ca.Constant  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
        ca.opacity = 0.50
        SetPlotOptions(ca)

    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 1
    denleg.numTicks = 6
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.02
    
    DrawPlots()
    ResetView()

    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SphericalSimulationSlice_R" + rad + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight / 2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a spherical slice of the density, with radiation energy density overplotted and star particles
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., rad: string indicating the radius at which slice is taken in pc e.g. "15.0", radflt: float indicating the radius at which slice is taken in pc e.g. 15.0, cont: flag indicating whether to plot energy density as contours e.g. "contour" or not, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def singlesphericalstarplot(pid, savpos, rad, radflt, cont, pheight, filedir):
    
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the vtk database
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    
    # Define Displacement operatios
    DefineVectorExpression("c", "coords(mesh)")
    DefineVectorExpression("starpos", "{0,0,0}")
    DefineVectorExpression("cstar", "c - starpos")
    DefineScalarExpression("radius", "sqrt(cstar[0]^2+cstar[1]^2+cstar[2]^2)")
    DefineScalarExpression("phi", "if(gt(cstar[1],0.),acos(cstar[0]/sqrt(cstar[0]^2+cstar[1]^2)),-acos(cstar[0]/sqrt(cstar[0]^2+cstar[1]^2)))")
    DefineScalarExpression("theta", "acos(cstar[2]/radius)-0.5*3.14159")
    DefineVectorExpression("sph", "{radius, theta, phi}")
    DefineVectorExpression("sphoffset", "sph-cstar")
    DefineScalarExpression("zint", "sqrt(1.+cos(theta)*cos(0.5*phi))")
    DefineScalarExpression("xint", "180.*cos(theta)*sin(0.5*phi)/zint")
    DefineScalarExpression("yint", "90.*sin(theta)/zint")
    DefineVectorExpression("cart", "{radius, yint, xint}")
    DefineVectorExpression("cartoffset", "cart-c")
    DefineScalarExpression("nhden", "density*28.9/29.06")
    
    # Create a plot of the scalar field 'density'
    AddPlot("Pseudocolor","nhden")
    
    # Apply a box to exclude the line y = 0
    AddOperator("Box", 1)
    bas = BoxAttributes()
    bas.minx = -30
    bas.maxx = 30
    bas.miny = -0.001
    bas.maxy = 0.001
    bas.minz = -30
    bas.maxz = 30
    bas.inverse = 1
    SetOperatorOptions(bas, 1)
    
    # Apply a displacement operator to warp to hammer projection
    AddOperator("Displace", 1)
    das = DisplaceAttributes()
    das.factor = 1
    das.variable = "cartoffset"
    SetOperatorOptions(das, 1)
    
    # Slice the volume to show only the central slice
    AddOperator("Slice")
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 100000
    psa.colorTableName = "bluehot"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    sa = SliceAttributes()
    sa.originType = sa.Intercept  # Point, Intercept, Percent, Zone, Node
    sa.originIntercept = radflt
    sa.axisType = sa.XAxis
    SetOperatorOptions(sa, 1)
    
    if cont == "contour":
        AddPlot("Contour", "rad_energy_density", 1, 1)
        ca = ContourAttributes()
        ca.colorType = ca.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
        ca.singleColor = (255, 0, 255, 255)
        ca.contourNLevels = 5
        ca.legendFlag = 0
        ca.minFlag = 1
        ca.maxFlag = 1
        ca.min = 10
        ca.max = 100000
        ca.scaling = ca.Log  # Linear, Log
        SetPlotOptions(ca)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 2
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 2
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    if cont == "pseudo":
        AddPlot("Pseudocolor", "rad_energy_density", 1, 1)
        ca = PseudocolorAttributes()
        ca.scaling = ca.Log  # Linear, Log, Skew
        if savpos == "midleft":
            ca.legendFlag = 1
        else:
            ca.legendFlag = 0
        ca.minFlag = 1
        ca.maxFlag = 1
        ca.min = 10
        ca.max = 100000
        ca.colorTableName = "RdPu"
        ca.invertColorTable = 0
        ca.opacityType = ca.Constant  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
        ca.opacity = 0.50
        SetPlotOptions(ca)
        
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    if savpos == "topleft":
        denleg = GetAnnotationObject(alllegs[0])
    elif savpos == "midleft":
        denleg = GetAnnotationObject(alllegs[1])
    else:
        denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 3
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.04
    if savpos == "topleft":
        denleg.suppliedValues = (0.01,1,100,10000)
        denleg.suppliedLabels = ("0.01","1","10^2","10^4")
    elif savpos == "midleft":
        denleg.suppliedValues = (10,1000,100000)
        denleg.suppliedLabels = ("10","10^3","10^5")
    denleg.drawLabels = denleg.Labels
    
    DrawPlots()
    ResetView()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SphericalSimulationSlice_R" + rad + "_T" + pid + savpos
    swa.width = pheight
    swa.height = pheight / 2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a spherical projection of the density
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., rad: string indicating the interior radius from which projection begins in pc e.g. "5.0", radflt: float indicating the interior radius from which projection begins in pc e.g. 5.0, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def sphericalprojplot(pid, savpos, rad, radflt, pheight, filedir):
    
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the vtk database
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    
    # Define Displacement operatios
    DefineVectorExpression("c", "coords(mesh)")
    DefineScalarExpression("radius", "sqrt(c[0]^2+c[1]^2+c[2]^2)")
    DefineScalarExpression("phi", "if(gt(c[1],0.),acos(c[0]/sqrt(c[0]^2+c[1]^2)),-acos(c[0]/sqrt(c[0]^2+c[1]^2)))")
    DefineScalarExpression("theta", "acos(c[2]/radius)-0.5*3.14159")
    DefineVectorExpression("sph", "{radius, theta, phi}")
    DefineVectorExpression("sphoffset", "sph-c")
    DefineScalarExpression("zint", "sqrt(1.+cos(theta)*cos(0.5*phi))")
    DefineScalarExpression("xint", "180.*cos(theta)*sin(0.5*phi)/zint")
    DefineScalarExpression("yint", "90.*sin(theta)/zint")
    DefineVectorExpression("cart", "{radius, yint, xint}")
    DefineVectorExpression("cartoffset", "cart-c")
    DefineScalarExpression("nhden", "density*28.9/29.06")
    
    # Create a plot of the scalar field 'density'
    DefineScalarExpression("column", "density*(30-"+rad+")/29.06")
    AddPlot("Pseudocolor","column")
    
    # Apply a spherical clipping to remove the inner pcs
    AddOperator("Clip", 1)
    cas = ClipAttributes()
    cas.funcType = cas.Sphere  # Plane, Sphere
    cas.radius = radflt
    cas.sphereInverse = 0
    SetOperatorOptions(cas, 1)
    
    # Apply a spherical clipping to remove the outer pcs
    AddOperator("Clip", 1)
    coas = ClipAttributes()
    coas.funcType = coas.Sphere  # Plane, Sphere
    coas.radius = 30
    coas.sphereInverse = 1
    SetOperatorOptions(coas, 1)
    
    # Apply a displacement operator to warp to hammer projection
    AddOperator("Displace", 1)
    das = DisplaceAttributes()
    das.factor = 1
    das.variable = "cartoffset"
    SetOperatorOptions(das, 1)
    
    AddOperator("DataBinning", 0)
    dbas = DataBinningAttributes()
    dbas.numDimensions = dbas.Three  # One, Two, Three
    dbas.dim1BinBasedOn = dbas.X  # X, Y, Z, Variable
    dbas.dim1NumBins = 256
    dbas.dim2BinBasedOn = dbas.Y  # X, Y, Z, Variable
    dbas.dim2NumBins = 256
    dbas.dim3BinBasedOn = dbas.Z  # X, Y, Z, Variable
    dbas.dim3NumBins = 256
    dbas.reductionOperator = dbas.Average  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
    dbas.varForReduction = "column"
    dbas.dim1NumBins = 1
    dbas.outputType = dbas.OutputOnBins
    SetOperatorOptions(dbas, 1)
    
    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    pas.projectionType = pas.ZYCartesian
    SetOperatorOptions(pas, 1)
    
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 1000
    psa.colorTableName = "Spectral"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 2
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 2
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 2
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.04
    denleg.suppliedValues = (0.01,1,100)
    denleg.suppliedLabels = ("0.01","1","100")
    denleg.drawLabels = denleg.Labels
    
    DrawPlots()
    ResetView()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SphericalSimulationColumn_R" + rad + savpos
    swa.width = pheight
    swa.height = pheight / 2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use full vtk dump to plot a spherical projection of the density with star particles
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, savpos: string indicating position of legends e.g. topright, topleft etc., rad: string indicating the interior radius from which projection begins in pc e.g. "5.0", radflt: float indicating the interior radius from which projection begins in pc e.g. 5.0, pheight: height of output plot in pixels, default 1024, filedir: directory in which input vtk dump is stored and in which output png is plotted
def sphericalprojstarplot(pid, savpos, rad, radflt, pheight, filedir):
    
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the vtk database
    OpenDatabase("localhost:" + filedir + "/RadParGrav_joined." + pid + ".vtk", 0)
    
    # Create a plot of the scalar field 'density'
    DefineScalarExpression("column", "density*(30-"+rad+")/29.06")
    AddPlot("Pseudocolor","column")
    
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 1000
    psa.colorTableName = "Spectral"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    # Define Displacement operatios
    DefineVectorExpression("c", "coords(mesh)")
    DefineScalarExpression("radius", "sqrt(c[0]^2+c[1]^2+c[2]^2)")
    DefineScalarExpression("phi", "if(gt(c[1],0.),acos(c[0]/sqrt(c[0]^2+c[1]^2)),-acos(c[0]/sqrt(c[0]^2+c[1]^2)))")
    DefineScalarExpression("theta", "acos(c[2]/radius)-0.5*3.14159")
    DefineVectorExpression("sph", "{radius, theta, phi}")
    DefineVectorExpression("sphoffset", "sph-c")
    DefineScalarExpression("zint", "sqrt(1.+cos(theta)*cos(0.5*phi))")
    DefineScalarExpression("xint", "180.*cos(theta)*sin(0.5*phi)/zint")
    DefineScalarExpression("yint", "90.*sin(theta)/zint")
    DefineVectorExpression("cart", "{radius, yint, xint}")
    DefineVectorExpression("cartoffset", "cart-c")
    DefineScalarExpression("nhden", "density*28.9/29.06")
    
    # Apply a spherical clipping to remove the inner pcs
    AddOperator("Clip", 1)
    cas = ClipAttributes()
    cas.funcType = cas.Sphere  # Plane, Sphere
    cas.radius = radflt
    cas.sphereInverse = 0
    SetOperatorOptions(cas, 1)
    
    # Apply a spherical clipping to remove the outer pcs
    AddOperator("Clip", 1)
    coas = ClipAttributes()
    coas.funcType = coas.Sphere  # Plane, Sphere
    coas.radius = 30
    coas.sphereInverse = 1
    SetOperatorOptions(coas, 1)
    
    # Apply a displacement operator to warp to hammer projection
    AddOperator("Displace", 1)
    das = DisplaceAttributes()
    das.factor = 1
    das.variable = "cartoffset"
    SetOperatorOptions(das, 1)
    
    AddOperator("DataBinning", 0)
    dbas = DataBinningAttributes()
    dbas.numDimensions = dbas.Three  # One, Two, Three
    dbas.dim1BinBasedOn = dbas.X  # X, Y, Z, Variable
    dbas.dim1NumBins = 256
    dbas.dim2BinBasedOn = dbas.Y  # X, Y, Z, Variable
    dbas.dim2NumBins = 256
    dbas.dim3BinBasedOn = dbas.Z  # X, Y, Z, Variable
    dbas.dim3NumBins = 256
    dbas.reductionOperator = dbas.Average  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
    dbas.varForReduction = "column"
    dbas.dim1NumBins = 1
    dbas.outputType = dbas.OutputOnBins
    SetOperatorOptions(dbas, 1)
    
    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    pas.projectionType = pas.ZYCartesian
    SetOperatorOptions(pas, 1)
    DrawPlots()
    
    # Open the star particle database
    OpenDatabase("localhost:" + filedir + "/RadParGrav.starpar." + pid + ".vtk", 0)
    
    # Plot the stellar masses
    DefineScalarExpression("starmass", "star_particle_mass / 29.06")
    # Plot the stellar masses
    AddPlot("Pseudocolor", "starmass", 1, 0)
    starpsa = PseudocolorAttributes()
    if savpos == "topright":
        starpsa.legendFlag = 1
    else:
        starpsa.legendFlag = 0
    starpsa.scaling = starpsa.Log  # Linear, Log, Skew
    starpsa.pointSize = 0.15
    starpsa.pointType = starpsa.Sphere  # Box, Axis, Icosahedron, Point, Sphere
    starpsa.colorTableName = "RdPu"
    starpsa.pointSizePixels = 30
    starpsa.minFlag = 1
    starpsa.maxFlag = 1
    starpsa.min = 50
    starpsa.max = 5000
    SetPlotOptions(starpsa)
    
    # Apply a spherical clipping to remove the outer pcs
    AddOperator("Clip", 0)
    coas = ClipAttributes()
    coas.funcType = coas.Sphere  # Plane, Sphere
    coas.radius = 30
    coas.sphereInverse = 1
    SetOperatorOptions(coas, 0)
    
    # Apply a displacement operator to warp to hammer projection
    AddOperator("Displace", 0)
    das = DisplaceAttributes()
    das.factor = 1
    das.variable = "cartoffset"
    SetOperatorOptions(das, 0)
    
    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    pas.projectionType = pas.ZYCartesian
    SetOperatorOptions(pas, 0)
    
    # Annotation
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    aa = AnnotationAttributes()
    aa.axes2D.visible = 1
    aa.axes2D.xAxis.title.visible = 0
    aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Times
    aa.axes2D.xAxis.label.font.scale = 2
    aa.axes2D.xAxis.label.font.bold = 1
    aa.axes2D.xAxis.label.font.italic = 0
    aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Times
    aa.axes2D.yAxis.label.font.scale = 2
    aa.axes2D.yAxis.label.font.bold = 1
    aa.axes2D.yAxis.label.font.italic = 0
    aa.axes2D.yAxis.title.visible = 0
    aa.userInfoFlag = 0
    aa.databaseInfoFlag = 0
    SetAnnotationAttributes(aa)
    
    # Set the Legend attributes
    alllegs = GetAnnotationObjectNames()
    denleg = GetAnnotationObject(alllegs[0])
    denleg.drawTitle = 0
    denleg.drawMinMax = 0
    denleg.managePosition = 0
    denleg.position = (0.01, 0.92)
    denleg.numberFormat = "%# -9.1e"
    denleg.controlTicks = 0
    denleg.numTicks = 2
    denleg.xScale = 1
    denleg.yScale = 3
    denleg.fontFamily = denleg.Times
    denleg.fontBold = 1
    denleg.fontHeight = 0.04
    denleg.suppliedValues = (0.01,1,100)
    denleg.suppliedLabels = ("0.01","1","100")
    denleg.drawLabels = denleg.Labels
    
    starleg = GetAnnotationObject(alllegs[1])
    starleg.drawTitle = 0
    starleg.drawMinMax = 0
    starleg.managePosition = 0
    starleg.position = (0.01, 0.92)
    starleg.numberFormat = "%# -9.1e"
    starleg.controlTicks = 0
    starleg.numTicks = 2
    starleg.xScale = 1
    starleg.yScale = 3
    starleg.fontFamily = starleg.Times
    starleg.fontBold = 1
    starleg.fontHeight = 0.03
    starleg.suppliedValues = (100,1000)
    starleg.suppliedLabels = ("10^2","10^3")
    starleg.drawLabels = starleg.Labels
    
    DrawPlots()
    ResetView()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = "SphericalSimulationColumnStar_R" + rad + savpos
    swa.width = pheight
    swa.height = pheight / 2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

## Use fits file and output of Disperse to plot projections of surface density as well as locations of filaments
## Inputs are pid: 4-digit id indicating time of dump e.g. "0003" for RadParGrav_Joined.0003.vtk, filedir: directory in which input vtk dump is stored and in which output png is plotted, nres: resolution of fits file e.g. 256 or 128, savpos: string indicating position of legends e.g. topright, topleft etc., pheight: height of output plot in pixels, default 1024
def filamentplot(pid, filedir, nres, radius, savpos, pheight):
    
    fname = "den{0:02d}".format(pid)
    
    # Clear any previous plots
    DeleteAllPlots()
    
    # Open the vtk database
    OpenDatabase("localhost:" + filedir + "/" + fname + "_proj.fits", 0)
    
    #DefineVectorExpression("c", "coords(mesh)")
    #DefineVectorExpression("cart", "c * 4.0 * 15.0 /"+nres+" - 2.0 * 15.0")
    #DefineVectorExpression("cartoffset", "cart-c")
    DefineScalarExpression("column", "hdu1 / 29.06")
    AddPlot("Pseudocolor","column")
    
    #AddOperator("Displace", 0)
    #das = DisplaceAttributes()
    #das.factor = 1
    #das.variable = "cartoffset"
    #SetOperatorOptions(das, 0)
    
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 1000
    psa.colorTableName = "Spectral"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    OpenDatabase("localhost:" + filedir + "/" + fname + ".vtk", 0)
    AddPlot("Pseudocolor","field_value")
    
    psa = PseudocolorAttributes()
    psa.minFlag = 1
    psa.maxFlag = 1
    psa.scaling = psa.Log  # Linear, Log, Skew
    if savpos == "topleft":
        psa.legendFlag = 1
    else:
        psa.legendFlag = 0
    psa.min = 0.01
    psa.max = 0.1
    psa.colorTableName = "PuRd"
    psa.invertColorTable = 1
    SetPlotOptions(psa)
    
    starname = "/RadParGrav.starpar.{0:04d}.vtk".format(pid)
    
    # Open the star particle database
    OpenDatabase("localhost:" + filedir + starname, 0)
    
    # Plot the stellar masses
    DefineScalarExpression("starmass", "star_particle_mass / 29.06")
    
    # Plot the stellar masses
    AddPlot("Pseudocolor", "starmass", 1, 0)
    starpsa = PseudocolorAttributes()
    if savpos == "bottomright":
        starpsa.legendFlag = 1
    else:
        starpsa.legendFlag = 0
    starpsa.scaling = starpsa.Log  # Linear, Log, Skew
    starpsa.pointSize = 0.02
    starpsa.pointType = starpsa.Sphere  # Box, Axis, Icosahedron, Point, Sphere
    starpsa.colorTableName = "RdPu"
    starpsa.pointSizePixels = 10
    starpsa.minFlag = 1
    starpsa.maxFlag = 1
    starpsa.min = 50
    starpsa.max = 5000
    SetPlotOptions(starpsa)
    
    # Then Project to 2D
    AddOperator("Project", 0)
    pas = ProjectAttributes()
    pas.projectionType = pas.XYCartesian
    SetOperatorOptions(pas, 1)
    
    DefineVectorExpression("c", "coords(mesh)")
    DefineVectorExpression("cart", "(c + 2.0 * "+radius+") * "+nres+" / (4.0 * "+radius+")")
    DefineVectorExpression("cartoffset", "cart-c")
    
    AddOperator("Displace", 0)
    das = DisplaceAttributes()
    das.factor = 1
    das.variable = "cartoffset"
    SetOperatorOptions(das, 0)
    
    DrawPlots()
    ResetView()
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = filedir
    swa.fileName = fname + "_filamentplot" + savpos
    swa.width = pheight
    swa.height = pheight / 2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

plotheight = 1024
plotflag = "fin"

## plotflag is used to set the type of plots we want to output
## Options are:
## fil: filament plots
## sdmag: Surface Densities with magnetic fields for no feedback models with varying magnetic field
## mscat: Scatter plots for turbulent and alfven velocity as a function of density for no feedback models with varying magnetic field
## mag: Density Slices with magnetic fields for no feedback models with varying magnetic field
## him: Density slices and spherical density slices for feedback model with high magnetic field
## lom: Density slices and spherical density slices for feedback model with low magnetic field
## nom: Density slices and spherical density slices for feedback model with no magnetic field
## sde: Surface Density projections and circumcluster projections for models with varying surface density at
## early times
## sdl: Surface Density projections and circumcluster projections for models with varying surface density at
## late times
## fid: Density slices and spherical density slices for fiducial model as a function of time
## hid: Density slices and spherical density slices for high surface density model as a function of time
## lod: Density slices and spherical density slices for low surface density model as a function of time
## fin: Surface Density Projections for fiducial model at very late times as gas is being driven out

if plotflag == "fil":

    #filedir = "/Users/sudhirraskutti/Desktop/Thesis/Code/Disperse/bin/"
    #fname = "den02.fits_c100.up.NDskl.TRIM.ASMB.S006.BRK"
    filedir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/RB0.05"
    for i in xrange(0,11):
        filamentplot(i, filedir, "256.", "15.", "topright", plotheight)
    
if plotflag == "sdmag":

    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/NoFeedback"
    starprojectionplot("0030", "topright", "xy", 256, plotheight, outdir)
    exit()
    AddWindow()
    starprojectionplot("0030", "topright", "yz", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0030", "topleft", "2", 2., 2 * plotheight, outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B5.0"
    starprojectionplot("0003", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    starprojectionplot("0003", "topright", "yz", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0003", "topright", "2", 2., 2 * plotheight, outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B0.5"
    starprojectionplot("0003", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    starprojectionplot("0003", "topright", "yz", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0003", "bottomright", "2", 2., 2 * plotheight, outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B0.1"
    starprojectionplot("0003", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    starprojectionplot("0003", "topright", "yz", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0003", "bottomright", "2", 2., 2 * plotheight, outdir)
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/Visit"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "fsixd"
    swa.width = 4*plotheight
    swa.height = 4*plotheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    
    swa.subWindowAtts.win1.position = (0, 3*plotheight)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (plotheight, 3*plotheight)
    swa.subWindowAtts.win2.size = (plotheight, plotheight)
    swa.subWindowAtts.win3.position = (2*plotheight, 3*plotheight)
    swa.subWindowAtts.win3.size = (2*plotheight, plotheight)
    
    swa.subWindowAtts.win4.position = (0, 2*plotheight)
    swa.subWindowAtts.win4.size = (plotheight, plotheight)
    swa.subWindowAtts.win5.position = (plotheight, 2*plotheight)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (2*plotheight, 2*plotheight)
    swa.subWindowAtts.win6.size = (2*plotheight, plotheight)
    
    swa.subWindowAtts.win7.position = (0, plotheight)
    swa.subWindowAtts.win7.size = (plotheight, plotheight)
    swa.subWindowAtts.win8.position = (plotheight, plotheight)
    swa.subWindowAtts.win8.size = (plotheight, plotheight)
    swa.subWindowAtts.win9.position = (2*plotheight, plotheight)
    swa.subWindowAtts.win9.size = (2*plotheight, plotheight)
    
    swa.subWindowAtts.win10.position = (0, 0)
    swa.subWindowAtts.win10.size = (plotheight, plotheight)
    swa.subWindowAtts.win11.position = (plotheight, 0)
    swa.subWindowAtts.win11.size = (plotheight, plotheight)
    swa.subWindowAtts.win12.position = (2*plotheight, 0)
    swa.subWindowAtts.win12.size = (2*plotheight, plotheight)
    
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "mscat":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B5.0"
    magscatterplot("0003", plotheight, "va", outdir)
    AddWindow()
    magscatterplot("0003", plotheight, "vturb", outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B0.5"
    magscatterplot("0003", plotheight, "va", outdir)
    AddWindow()
    magscatterplot("0003", plotheight, "vturb", outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B0.1"
    magscatterplot("0003", plotheight, "va", outdir)
    AddWindow()
    magscatterplot("0003", plotheight, "vturb", outdir)

    outdir = "/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/Visit"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "fsixe"
    swa.width = 3*plotheight
    swa.height = 2*plotheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 0)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (0, plotheight)
    swa.subWindowAtts.win2.size = (plotheight, plotheight)
    swa.subWindowAtts.win3.position = (plotheight, 0)
    swa.subWindowAtts.win3.size = (plotheight, plotheight)
    swa.subWindowAtts.win4.position = (plotheight, plotheight)
    swa.subWindowAtts.win4.size = (plotheight, plotheight)
    swa.subWindowAtts.win5.position = (2*plotheight, 0)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (2*plotheight, plotheight)
    swa.subWindowAtts.win6.size = (plotheight, plotheight)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "mag":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/NoFeedback"
    singlestarnoradplot("0030", "topleft", plotheight, "xy", outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B5.0"
    singlestarmagplot("0003", "topright", plotheight, "xy", outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B0.5"
    singlestarmagplot("0003", "topleft", plotheight, "xy", outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B0.1"
    singlestarmagplot("0003", "bottomright", plotheight, "xy", outdir)
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/Visit"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "fsixa"
    swa.width = 2*plotheight
    swa.height = 2*plotheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, plotheight)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (plotheight, plotheight)
    swa.subWindowAtts.win2.size = (plotheight, plotheight)
    swa.subWindowAtts.win3.position = (0, 0)
    swa.subWindowAtts.win3.size = (plotheight, plotheight)
    swa.subWindowAtts.win4.position = (plotheight, 0)
    swa.subWindowAtts.win4.size = (plotheight, plotheight)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()


    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/B0.01"
    #projectionplot("0000", "topright", "xy", 128, plotheight, outdir)
    #projectionplot("0001", "topright", "xy", 128, plotheight, outdir)
    #projectionplot("0002", "topright", "xy", 128, plotheight, outdir)
    #projectionplot("0003", "topright", "xy", 128, plotheight, outdir)
    #projectionplot("0004", "topright", "xy", 128, plotheight, outdir)
    #projectionplot("0005", "topright", "xy", 128, plotheight, outdir)

if plotflag == "him":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/RB0.05"
    singlestarplot("0001", "bottomright", plotheight, "xy", outdir)
    AddWindow()
    #starprojectionplot("0002", "topright", "xy", 128, plotheight, outdir)
    singlestarplot("0001", "bottomleft", plotheight, "yz", outdir)
    AddWindow()
    singlesphericalstarplot("0001", "topleft", "8", 8., "pseudo", plotheight, outdir)
    AddWindow()
    singlestarplot("0004", "topright", plotheight, "xy", outdir)
    AddWindow()
    #starprojectionplot("0003", "topleft", "xy", 128, plotheight, outdir)
    singlestarplot("0004", "topleft", plotheight, "yz", outdir)
    AddWindow()
    singlesphericalstarplot("0004", "midleft", "8", 8., "pseudo", plotheight, outdir)
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/Visit"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "fsixf"
    swa.width = 2*plotheight
    swa.height = 5*plotheight/2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 3*plotheight/2)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (0, plotheight/2)
    swa.subWindowAtts.win2.size = (plotheight, plotheight)
    swa.subWindowAtts.win3.position = (0, 0)
    swa.subWindowAtts.win3.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win4.position = (plotheight, 3*plotheight/2)
    swa.subWindowAtts.win4.size = (plotheight, plotheight)
    swa.subWindowAtts.win5.position = (plotheight, plotheight/2)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (plotheight, 0)
    swa.subWindowAtts.win6.size = (plotheight, plotheight / 2)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "lom":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/RB5.0"
    singlestarplot("0003", "bottomright", plotheight, "xy", outdir)
    AddWindow()
    #starprojectionplot("0005", "topright", "xy", 128, plotheight, outdir)
    singlestarplot("0003", "bottomleft", plotheight, "yz", outdir)
    AddWindow()
    singlesphericalstarplot("0003", "topleft", "8", 8., "pseudo", plotheight, outdir)
    AddWindow()
    singlestarplot("0010", "topright", plotheight, "xy", outdir)
    AddWindow()
    #starprojectionplot("0010", "topleft", "xy", 128, plotheight, outdir)
    singlestarplot("0010", "topleft", plotheight, "yz", outdir)
    AddWindow()
    singlesphericalstarplot("0010", "midleft", "8", 8., "pseudo", plotheight, outdir)
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/Visit"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "fsixg-new"
    swa.width = 2*plotheight
    swa.height = 5*plotheight/2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 3*plotheight/2)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (0, plotheight/2)
    swa.subWindowAtts.win2.size = (plotheight, plotheight)
    swa.subWindowAtts.win3.position = (0, 0)
    swa.subWindowAtts.win3.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win4.position = (plotheight, 3*plotheight/2)
    swa.subWindowAtts.win4.size = (plotheight, plotheight)
    swa.subWindowAtts.win5.position = (plotheight, plotheight/2)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (plotheight, 0)
    swa.subWindowAtts.win6.size = (plotheight, plotheight / 2)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "nom":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/RB5000.0"
    singlestarplot("0004", "bottomright", plotheight, "xy", outdir)
    AddWindow()
    #starprojectionplot("0004", "topright", "xy", 128, plotheight, outdir)
    singlestarplot("0004", "bottomleft", plotheight, "yz", outdir)
    AddWindow()
    singlesphericalstarplot("0004", "topleft", "15", 15., "pseudo", plotheight, outdir)
    AddWindow()
    singlestarplot("0008", "topright", plotheight, "xy", outdir)
    AddWindow()
    #starprojectionplot("0008", "topleft", "xy", 128, plotheight, outdir)
    singlestarplot("0008", "bottomleft", plotheight, "yz", outdir)
    AddWindow()
    singlesphericalstarplot("0008", "midleft", "15", 15., "pseudo", plotheight, outdir)
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/Visit"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "fsixh"
    swa.width = 2*plotheight
    swa.height = 5*plotheight/2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 3*plotheight/2)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (0, plotheight/2)
    swa.subWindowAtts.win2.size = (plotheight, plotheight)
    swa.subWindowAtts.win3.position = (0, 0)
    swa.subWindowAtts.win3.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win4.position = (plotheight, 3*plotheight/2)
    swa.subWindowAtts.win4.size = (plotheight, plotheight)
    swa.subWindowAtts.win5.position = (plotheight, plotheight/2)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (plotheight, 0)
    swa.subWindowAtts.win6.size = (plotheight, plotheight / 2)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()


if plotflag == "sdl":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/LowDensity"
    starprojectionplot("0010", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0010", "topleft", "2", 2., 2 * plotheight, outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial"
    starprojectionplot("0037", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0037", "topright", "2", 2., 2 * plotheight, outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/HighDensity"
    starprojectionplot("0010", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0010", "bottomright", "2", 2., 2 * plotheight, outdir)
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "ffive"
    swa.width = 3*plotheight
    swa.height = 3*plotheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 2*plotheight)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (plotheight, 2*plotheight)
    swa.subWindowAtts.win2.size = (2*plotheight, plotheight)
    swa.subWindowAtts.win3.position = (0, plotheight)
    swa.subWindowAtts.win3.size = (plotheight, plotheight)
    swa.subWindowAtts.win4.position = (plotheight, plotheight)
    swa.subWindowAtts.win4.size = (2*plotheight, plotheight)
    swa.subWindowAtts.win5.position = (0, 0)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (plotheight, 0)
    swa.subWindowAtts.win6.size = (2*plotheight, plotheight)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "sde":

    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/LowDensity"
    starprojectionplot("0007", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0007", "topleft", "2", 2., 2 * plotheight, outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial"
    starprojectionplot("0025", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0025", "topright", "2", 2., 2 * plotheight, outdir)
    AddWindow()
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/HighDensity"
    starprojectionplot("0007", "topright", "xy", 256, plotheight, outdir)
    AddWindow()
    sphericalprojstarplot("0007", "bottomright", "2", 2., 2 * plotheight, outdir)

    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial"
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "ffour"
    swa.width = 3*plotheight
    swa.height = 3*plotheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 2*plotheight)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (plotheight, 2*plotheight)
    swa.subWindowAtts.win2.size = (2*plotheight, plotheight)
    swa.subWindowAtts.win3.position = (0, plotheight)
    swa.subWindowAtts.win3.size = (plotheight, plotheight)
    swa.subWindowAtts.win4.position = (plotheight, plotheight)
    swa.subWindowAtts.win4.size = (2*plotheight, plotheight)
    swa.subWindowAtts.win5.position = (0, 0)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (plotheight, 0)
    swa.subWindowAtts.win6.size = (2*plotheight, plotheight)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "fin":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperI/Visit/Fiducial"
    sigmaprojectionplot("0200", "topleft", "yz", 256, plotheight, outdir)
    AddWindow()
    sigmaprojectionplot("0200", "bottomright", "xy", 256, plotheight, outdir)

    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "SimulationColumnAll" + "_T0200"
    swa.width = plotheight
    swa.height = 2*plotheight
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, plotheight)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win1.layer = 0
    swa.subWindowAtts.win1.transparency = 0
    swa.subWindowAtts.win1.omitWindow = 0
    swa.subWindowAtts.win2.position = (0, 0)
    swa.subWindowAtts.win2.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.layer = 0
    swa.subWindowAtts.win2.transparency = 0
    swa.subWindowAtts.win2.omitWindow = 0
    
    SetSaveWindowAttributes(swa)
    SaveWindow()

    DeleteAllPlots()

if plotflag == "ng":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperI/Visit/NG"
    singleplot("0003", "topleft", plotheight, outdir)
    projectionplot("0003", "topleft", "xy", 128, plotheight, outdir)
    projectionplot("0003", "topleft", "yz", 128, plotheight, outdir)
    projectionplot("0003", "topleft", "xz", 128, plotheight, outdir)


if plotflag == "fid":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial"
    singlestarplot("0011", "bottomright", plotheight, "xy", r"$\displaystyle t_{\rm 2}$", outdir)
#    AddWindow()
#    singlesphericalstarplot("0011", "bottomleft", "10", 10., "pseudo", plotheight, outdir)
#    AddWindow()
#    singlesphericalstarplot("0011", "bottomleft", "15", 15., "pseudo", plotheight, outdir)
#    AddWindow()
#    singlesphericalstarplot("0011", "bottomleft", "25", 25., "pseudo", plotheight, outdir)
#    AddWindow()
#    singlestarplot("0025", "topright", plotheight, "xy", outdir)
#    AddWindow()
#    singlesphericalstarplot("0025", "topleft", "10", 10., "pseudo", plotheight, outdir)
#    AddWindow()
#    singlesphericalstarplot("0025", "midleft", "15", 15., "pseudo", plotheight, outdir)
#    AddWindow()
#    singlesphericalstarplot("0025", "bottomleft", "25", 25., "pseudo", plotheight, outdir)
#    
#    swa = SaveWindowAttributes()
#    swa.outputToCurrentDirectory = 0
#    swa.family = 0
#    swa.outputDirectory = outdir
#    swa.fileName = "fone"
#    swa.width = 2*plotheight
#    swa.height = 5*plotheight/2
#    swa.screenCapture = 0
#    swa.quality = 80
#    swa.progressive = 0
#    swa.binary = 0
#    swa.stereo = 0
#    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
#    swa.forceMerge = 0
#    swa.resConstraint = swa.NoConstraint
#    swa.advancedMultiWindowSave = 1
#    swa.subWindowAtts.win1.position = (0, 3*plotheight/2)
#    swa.subWindowAtts.win1.size = (plotheight, plotheight)
#    swa.subWindowAtts.win2.position = (0, plotheight)
#    swa.subWindowAtts.win2.size = (plotheight, plotheight / 2)
#    swa.subWindowAtts.win3.position = (0, plotheight/2)
#    swa.subWindowAtts.win3.size = (plotheight, plotheight / 2)
#    swa.subWindowAtts.win4.position = (0, 0)
#    swa.subWindowAtts.win4.size = (plotheight, plotheight / 2)
#    swa.subWindowAtts.win5.position = (plotheight, 3*plotheight/2)
#    swa.subWindowAtts.win5.size = (plotheight, plotheight)
#    swa.subWindowAtts.win6.position = (plotheight, plotheight)
#    swa.subWindowAtts.win6.size = (plotheight, plotheight / 2)
#    swa.subWindowAtts.win7.position = (plotheight, plotheight/2)
#    swa.subWindowAtts.win7.size = (plotheight, plotheight / 2)
#    swa.subWindowAtts.win8.position = (plotheight, 0)
#    swa.subWindowAtts.win8.size = (plotheight, plotheight / 2)
#    SetSaveWindowAttributes(swa)
#    SaveWindow()
#    DeleteAllPlots()

if plotflag == "lod":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/LowDensity"
    singlestarplot("0003", "bottomright", plotheight, "xy", outdir)
    AddWindow()
    singlesphericalstarplot("0003", "bottomleft", "10", 10., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0003", "bottomleft", "15", 15., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0003", "bottomleft", "25", 25., "pseudo", plotheight, outdir)
    AddWindow()
    singlestarplot("0007", "topright", plotheight, "xy", outdir)
    AddWindow()
    singlesphericalstarplot("0007", "topleft", "10", 10., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0007", "midleft", "15", 15., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0007", "bottomleft", "25", 25., "pseudo", plotheight, outdir)
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "ftwo"
    swa.width = 2*plotheight
    swa.height = 5*plotheight/2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 3*plotheight/2)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (0, plotheight)
    swa.subWindowAtts.win2.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win3.position = (0, plotheight/2)
    swa.subWindowAtts.win3.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win4.position = (0, 0)
    swa.subWindowAtts.win4.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win5.position = (plotheight, 3*plotheight/2)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (plotheight, plotheight)
    swa.subWindowAtts.win6.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win7.position = (plotheight, plotheight/2)
    swa.subWindowAtts.win7.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win8.position = (plotheight, 0)
    swa.subWindowAtts.win8.size = (plotheight, plotheight / 2)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "hid":
    
    outdir = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/HighDensity"
    singlestarplot("0003", "bottomright", plotheight, "xy", outdir)
    AddWindow()
    singlesphericalstarplot("0003", "bottomleft", "10", 10., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0003", "bottomleft", "15", 15., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0003", "bottomleft", "25", 25., "pseudo", plotheight, outdir)
    AddWindow()
    singlestarplot("0007", "topright", plotheight, "xy", outdir)
    AddWindow()
    singlesphericalstarplot("0007", "topleft", "10", 10., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0007", "midleft", "15", 15., "pseudo", plotheight, outdir)
    AddWindow()
    singlesphericalstarplot("0007", "bottomleft", "25", 25., "pseudo", plotheight, outdir)
    
    swa = SaveWindowAttributes()
    swa.outputToCurrentDirectory = 0
    swa.family = 0
    swa.outputDirectory = outdir
    swa.fileName = "fthree"
    swa.width = 2*plotheight
    swa.height = 5*plotheight/2
    swa.screenCapture = 0
    swa.quality = 80
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits  # None, PackBits, Jpeg, Deflate
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    swa.advancedMultiWindowSave = 1
    swa.subWindowAtts.win1.position = (0, 3*plotheight/2)
    swa.subWindowAtts.win1.size = (plotheight, plotheight)
    swa.subWindowAtts.win2.position = (0, plotheight)
    swa.subWindowAtts.win2.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win3.position = (0, plotheight/2)
    swa.subWindowAtts.win3.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win4.position = (0, 0)
    swa.subWindowAtts.win4.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win5.position = (plotheight, 3*plotheight/2)
    swa.subWindowAtts.win5.size = (plotheight, plotheight)
    swa.subWindowAtts.win6.position = (plotheight, plotheight)
    swa.subWindowAtts.win6.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win7.position = (plotheight, plotheight/2)
    swa.subWindowAtts.win7.size = (plotheight, plotheight / 2)
    swa.subWindowAtts.win8.position = (plotheight, 0)
    swa.subWindowAtts.win8.size = (plotheight, plotheight / 2)
    SetSaveWindowAttributes(swa)
    SaveWindow()
    DeleteAllPlots()

if plotflag == "cmf":

    totalim = Image.new('RGB', (2*plotheight, plotheight), "white")

    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperI/Visit/Fiducial/SimulationColumnyz_T0200topleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,0))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperI/Visit/Fiducial/SimulationColumnxy_T0200bottomright.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,plotheight))
    totalim.save("/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/f3a.png")

if plotflag == "cmb":

    totalim = Image.new('RGB', (2*plotheight, 2*plotheight + plotheight/2), "white")
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SimulationSlice_T0011topleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,0))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SphericalSimulationSlice_R10_T0011bottomleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,plotheight))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SphericalSimulationSlice_R15_T0011bottomleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,plotheight+plotheight/2))

    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SphericalSimulationSlice_R25_T0011bottomleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,2*plotheight))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SimulationSlice_T0025topright.png"
    imone = Image.open(imname)
    totalim.paste(imone, (plotheight,0))

    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SphericalSimulationSlice_R10_T0025topleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (plotheight,plotheight))

    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SphericalSimulationSlice_R15_T0025bottomleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (plotheight,plotheight+plotheight/2))

    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SphericalSimulationSlice_R25_T0025bottomleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (plotheight,2*plotheight))
    
    totalim.save("/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/fone.png")


    totalim = Image.new('RGB', (plotheight + plotheight/2, 3*plotheight), "white")

    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/HighDensity/SimulationColumnxytopleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,0))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/HighDensity/SphericalSimulationColumn_R10topleft.png"
    imone = Image.open(imname)
    totalim.paste(imone, (0,plotheight))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SimulationColumnxytopright.png"
    imone = Image.open(imname)
    totalim.paste(imone, (plotheight,0))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/SphericalSimulationColumn_R10topright.png"
    imone = Image.open(imname)
    totalim.paste(imone, (plotheight,plotheight))

    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/LowDensity/SimulationColumnxytopright.png"
    imone = Image.open(imname)
    totalim.paste(imone, (2*plotheight,0))
    
    imname = "/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/LowDensity/SphericalSimulationColumn_R10topright.png"
    imone = Image.open(imname)
    totalim.paste(imone, (2*plotheight,plotheight))
    
    totalim.save("/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/Fiducial/ftwo.png")

