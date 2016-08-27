from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_math import *
from hyp_fluxes import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import fminbound

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

# Missing 1e5r15

mlist = ['5.0e3', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '20.0', '25.0', '35.0', '15.0', '25.0', '35.0']

mlist = ['2.0e4', '5.0e4', '2.0e4', '5.0e4', '1.0e5', '2.0e4', '1.0e4', '5.0e4', '1.0e4', '1.0e5', '2.0e5', '5.0e3', '2.0e4', '5.0e4', '1.0e5', '2.0e4', '2.0e5', '1.0e4', '5.0e5', '1.0e5', '5.0e4', '2.0e5', '5.0e4', '2.0e4', '2.0e5']
rlist = ['25.0', '35.0', '20.0', '25.0', '35.0', '15.0', '10.0', '20.0', '8.0', '25.0', '35.0', '5.0', '10.0', '15.0', '20.0', '8.0', '25.0', '5.0', '35.0', '15.0', '10.0', '20.0', '8.0', '5.0', '15.0']

nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflistnf = np.core.defchararray.add(dflist, ['_N256_Tf4_NF_Fluxes/'] * nfs)
dflist = np.core.defchararray.add(dflist, ['_N256_Tf4_Fluxes/'] * nfs)

#headoutlist = ['model', 'e', 'eadj', 'ts', 't50', 't90', 'tunb', 'beta', 'tbreak', 'eff']
headoutlist = ['model', 'alpha', 'fedd']
nheads = len(headoutlist)
nconv = 20

tablefile = 'restable.dat'
tbf = open(plotdir + tablefile, 'w+')

for i in xrange(0,2):
    
    datafolder = dflist[i]
    datafoldernf = dflistnf[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    hstnfdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
    
    tff = out_tff(outlines)
    tMyr = out_tMyr(outlines)
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    mcout = mcloud / msol
    rcloud = out_rcloud(outlines)
    sigmacloud = out_sigma(outlines)
    alpha = out_alphavir(outlines)
    vrms = out_vturb(outlines)
    vesc = vrms * np.sqrt(5.0 * 2.0 / (3.0 * alpha))
    nh0 = out_rhocloud(outlines)
    mprint = float2explatex(mcout)
    mname = modelstr(mcout, rcloud)
    clight = out_clight(outlines)
    kappa = out_kappa(outlines)
    
    time = hst_time(hstdata, tff)
    timenf = hst_time(hstnfdata, tff)
    
    eff = hst_mstar(hstdata, mcloud)
    eps = hst_eff(hstdata, mcloud)
    epsnf = hst_eff(hstnfdata, mcloud)
    epsadj = eps / (1.0 - 0.12)
    ts = hst_tstar(hstdata, mcloud, tff)
    t50 = hst_teff(hstdata, mcloud, tff, 0.5)
    t90 = hst_teff(hstdata, mcloud, tff, 0.9)
    
    strout = ''
    for j in xrange(0, nheads):
        head = headoutlist[j]
        if j < nheads-1:
            fend = ' & '
        else:
            fend = '\\\\ \n'
        if head == 'model':
            strout = strout + mname + fend
        elif head == 'sigma':
            strout = strout + '{:.2f}'.format(sigmacloud) + fend
        elif head == 'r':
            strout = strout + '{:d}'.format(int(rcloud)) + fend
        elif head == 'tff':
            strout = strout + '{:.2f}'.format(tff / tMyr) + fend
        elif head == 'm':
            strout = strout + mprint + fend
        elif head == 'a':
            strout = strout + '{:.1f}'.format(alpha) + fend
        elif head == 'vesc':
            strout = strout + '{:.2f}'.format(vesc) + fend
        elif head == 'vrms':
            strout = strout + '{:.2f}'.format(vrms) + fend
        elif head == 'nh':
            strout = strout + '{:.2f}'.format(nh0) + fend
        elif head == 'e':
            strout = strout + '{:.2f}'.format(eps) + fend
        elif head == 'eadj':
            strout = strout + '{:.2f}'.format(epsadj) + fend
        elif head == 'ts':
            strout = strout + '{:.2f}'.format(ts) + fend
        elif head == 't50':
            strout = strout + '{:.2f}'.format(t50) + fend
        elif head == 't90':
            strout = strout + '{:.2f}'.format(t90) + fend
        elif head == 'alpha':
            fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
            tmini = np.abs(eff - 0.5 * eps).argmin()
            thist = time[tmini]
            alpha = flux_tableouts(fluxdata, thist, 'alpha', rcloud, tff, mcloud, 1.0)
            strout = strout + '{:.2f}'.format(alpha) + fend
        elif head == 'fedd':
            fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
            tmini = np.abs(eff - 0.9 * eps).argmin()
            thist = time[tmini]
            fout = flux_tableouts(fluxdata, thist, 'fedd', rcloud, tff, mcloud, kappa / clight)
            strout = strout + '{:.2f}'.format(fout) + fend
        else:
            strout = strout + '{:.2f}'.format(0.0) + fend

    tbf.write(strout)

    print i, datafolder, max(time), max(timenf)

tbf.close()

