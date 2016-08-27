from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
datadir = '/tigress-hsm/raskutti/tigress-hsm/Hyperion/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['UV_M2.0e4_R15.0_N256_Tf4_All/', 'UV_M2.0e5_R15.0_N256_Tf4_All/']
nds = len(dflist)

for i in xrange(0,nds):
    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    eff = hst_eff(hstdata, mcloud)
    t50 = hst_teff(hstdata, mcloud, tff, 0.5)
    m50 = hst_mteff(hstdata, mcloud, tff, 1.06) / eff
    t90 = hst_teff(hstdata, mcloud, tff, 0.9)
    m90 = hst_mteff(hstdata, mcloud, tff, 1.57) / eff
    print i, datafolder, eff, t50, m50, t90, m90


