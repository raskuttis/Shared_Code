from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_math import *
from hyp_fluxes import *
from hyp_pdf import *
from hyp_star import *
from hyp_models import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import fminbound

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
pdffile = 'sdpdfxy.dat'
starfile = 'star'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

#dflist = set_model_lookup('SDs/')
dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.1_LD_Fluxes/',
          'UV_M5.0e4_R15.0_N256_Tf4_B0.5_LD_Fluxes/',
          'UV_M5.0e4_R15.0_N256_Tf4_B1.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B2.0_LD_Fluxes/',
          'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B10.0_LD_Fluxes/',
          'UV_M5.0e4_R15.0_N256_Tf4_B20.0_LD_Fluxes/', 'UV_M5.0e4_R15.0_N256_Tf4_B50.0_LD_Fluxes/']
nfs = len(dflist)

for i in xrange(0,nfs):
    
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    
    tff = out_tff(outlines)
    time = hst_time(hstdata, tff)
        
    print i, datafolder, max(time)