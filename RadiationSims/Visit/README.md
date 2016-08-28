Package containing all the Python code used to plot density slices and surface density projections from the vtk dumps output from Hyperion simulations. To plot the vtk slices take the following steps

    1)  If the simulation has run on a number of clusters, the vtk files must first be joined. 
    ** Stuff about join here **

    2)  The vtk files must then be stored locally wherever Visit is being run

    3)  main.py contains all the scripts for 
        plotflag controls which type of plot you want to make. Options are:
            fil: filament plots
            sdmag: Surface Densities with magnetic fields for no feedback models with varying magnetic field
            mscat: Scatter plots for turbulent and alfven velocity as a function of density for no feedback models with varying magnetic field
            mag: Density Slices with magnetic fields for no feedback models with varying magnetic field
            him: Density slices and spherical density slices for feedback model with high magnetic field
            lom: Density slices and spherical density slices for feedback model with low magnetic field
            nom: Density slices and spherical density slices for feedback model with no magnetic field
            sde: Surface Density projections and circumcluster projections for models with varying surface density at early times
            sdl: Surface Density projections and circumcluster projections for models with varying surface density at late times
            fid: Density slices and spherical density slices for fiducial model as a function of time
            hid: Density slices and spherical density slices for high surface density model as a function of time
            lod: Density slices and spherical density slices for low surface density model as a function of time
            fin: Surface Density Projections for fiducial model at very late times as gas is being driven out

    4)  main.py can then be run from within Visit by Launching the CLI and then at the terminal entering 
        >> Source("path to main.py")