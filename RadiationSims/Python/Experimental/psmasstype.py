from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *

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
masses = []
totsmasses = []
totnewsmasses = []
totmergesmasses = []
totaccsmasses = []

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
    
    smass = star_mass(stardata) / mcloud
    newsmass = star_new_mass(stardata) / mcloud
    mergesmass = star_merge_mass(stardata, rcloud) / mcloud
    
    times.append(time)
    masses.append(mstar)
    totsmasses.append(smass)
    totnewsmasses.append(newsmass)
    totmergesmasses.append(mergesmass)
    totaccsmasses.append(smass - newsmass)


plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

colors = ['k','r','b','g','c']
plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
plt.plot(times[0],totsmasses[0],'k',times[0],totaccsmasses[0],'r',times[0],totnewsmasses[0],'b',times[0],totmergesmasses[0],'g')
plt.axis([0,4,0,1])
plt.xticks([])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

pp = PdfPages(plotdir + 'psmasstype.pdf')
pp.savefig()
pp.close()


