import Image

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

#DeleteAllPlots()
#CloseDatabase("localhost:" + filedir + "/RadParGrav_joined." + sigtxt + "." + pid + ".vtk")
#CloseDatabase("localhost:" + filedir + "/RadParGrav.starpar." + pid + ".vtk")