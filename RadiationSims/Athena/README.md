Folder stores all the files necessary to run simulations of self-gravitating star-forming clouds 
with radiative feedback in Hyperion

Take the following steps to run simulations, assuming that Hyperion has been correctly configured 
on tiger

1)  raduvcloud.c contains all the problem-specific information as well as some analysis and outputs. 
    Move to src/prob/raduvcloud.c

2)  Configure Hyperion (in the main directory) with ./config and the following parameters
    a) Full Hydro:
    --with-problem=raduvcloud --enable-radiation --with-nscalars=5 --with-flux=hllc --with-integrator=vl --with-radflux=hlle_m1 --with-eos=isothermal --enable-fft --with-gravity=fft_obc --enable-starparticle --enable-mpi
        nscalars command optional if using passive scalars, and change radflux depending on which closure is being used e.g. M1, P1, F1
    b) No Radiation:
    --with-problem=raduvcloud --with-gas=hydro --with-flux=hllc --with-integrator=vl --with-eos=isothermal --enable-fft --with-gravity=fft_obc --enable-starparticle --enable-mpi
    c) No Gravity:
    Same as Hydro or No Radiation configures, except in raduvcloud.c file, manually lower four_PI_G constant to 
    a low, but non-zero value
    d) Magnetic fields and no radiation:
    --with-problem=raduvcloud --with-gas=mhd --with-flux=hlld --with-integrator=vl --with-eos=isothermal --enable-fft --with-gravity=fft_obc --enable-starparticle --enable-mpi
    Note the change in flux from hllc to hlld
    e) Magnetic fields and radiation:
    --with-problem=raduvcloud --with-gas=mhd --enable-radiation --with-nscalars=5 --with-flux=hlld --with-integrator=vl --with-radflux=hlle_m1 --with-eos=isothermal --enable-fft --with-gravity=fft_obc --enable-starparticle --enable-mpi

3)  make all (in main directory) (May need to load openmpi and fftw first)

4)  Input files are stored in tst/new/RadUVCloud/$fileext (where $fileext is set in the bash file) 
    although this can be changed so long as the bash input file references the right input file. An example 
    input file is shown in the folder as athinput.raduvcloud

    Currently set to output 1) hst file 2) full vtk dump 3) starparticle vtk 4-6) vtk density slices 7-9) vtk surface density slices 10) restart dumps. Times are currently set to run for 4 freefall times (reset with tlim), with 10000 hst updates and 100 full vtk dumps. 
    maxout sets how many of these output blocks to use (default to just hst)

    Key new parameters different from earlier versions are:
    rho_small - factor by which the density in the surrounding ISM is lower than the mean cloud density
    rho_ffac - factor by which density floor is lower than rho_small
    rho_fdt - time intervals at which flux and PDF outputs are made
    surfd_out - This flag is used to indicated whether we want to output PDFs or not (1 or 0)
    fluxrad_out - Flag indicates whether to output flux of several quantities described below at radial shells away from either the box center or the stellar center of mass. 0 for no output or else fluxrad_out defines the number of radial bins

    Key physical parameters (described in the input file) are:
    iso_csound, M_GMC, rcloud, Psi, kappa, alpha_vir, beta

    Key numerical parameters (described in the input file) are:
    crad, theta, rsrc

5)  Bash files used as input to slurm are also stored in tst/new although again this can change. An example bash file is shown in sbrun.raduv.base

    Key sbatch submission parameters are:
    SBATCH -N #of nodes
    #SBATCH --ntasks-per-node=#cores per node
        Total number of cores must equal NGrid_x1 * NGrid_x2 * NGrid_x3 in input file
    SBATCH -t hh:mm:ss max run time
    SBATCH --mail-user=email address to send notifications (begin and/or end)

    Key parameters related to where input and output files should be stored
    basedir = Location of Hyperion main directory
    fileext = location of input file
    fileextout = Name of folder in which to store output files
    bindir = Location to which .out file is published
    outputdir = Location to which all other output files are published (usually to scratch)

6)  Can run this bash file individually for a single simulation run using sbatch sbrun sbrun.raduv.base. Alternately can use example bash file as shown in runall.bash to run a set of simulations. We can set parameters in the input file using input lists masses, radii, alphas, betas, fr, sr, psi and res. Meanwhile slurm parameters can be set using rtime.