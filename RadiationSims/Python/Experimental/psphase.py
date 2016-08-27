from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_star import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

datadir = '/Users/sudhirraskutti/Desktop/Thesis/Models/'
hstfile = 'RadParGrav.hst'
outfile = 'RadParGrav.out'
starfile = 'star'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dflist = ['Fiducial/']
nds = len(dflist)
times = []
cfracs = []
totsmasses = []

for i in xrange(0,nds):
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_local_outfile(datadir + datafolder + outfile)
    hstdata = read_local_hstfile(datadir + datafolder + hstfile)
    
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    
    stardata = read_local_allstars(datadir + datafolder + starfile,time,tff)
    
    sgra = star_grid(stardata, 256, 2.0*rcloud, 0.33*rcloud)
    sgrb = star_grid(stardata, 256, 2.0*rcloud, 0.67*rcloud)
    sgrc = star_grid(stardata, 256, 2.0*rcloud, rcloud)
    
    times.append(time)
    cfracs.append(sgra)
    cfracs.append(sgrb)
    cfracs.append(sgrc)


plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

colors = ['k','r','b','g','c']
plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
plt.plot(times[0],cfracs[0],'k',times[0],cfracs[1],'r',times[0],cfracs[2],'b')
plt.axis([0,4,0,0.1])
plt.xticks([])
plt.yticks([0,0.05,0.1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle V$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'psmasscfrac.pdf')
pp.savefig()
pp.close()


