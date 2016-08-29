from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tlist = [1.0,1.5,2.0,2.5]
nts = len(tlist)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
msol = out_msol(outlines)
dx = out_dx(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)

pdffile = 'vmasspdfall.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotvs = []
plotpdfs = []
for i in xrange(0,nts):

    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf_zip(pdflines,tminpdfi)
    plotv = sim_pdfx(pdflines)
    mnorm = mcloud / (dx**3)
    plotpdfs.append(plotpdf/mnorm)
    plotvs.append(plotv)

pdffile = 'vrmasspdfall.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotvrs = []
plotpdfrs = []
for i in xrange(0,nts):
    
    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(pdflines,tminpdfi)
    plotv = sim_pdfx(pdflines)
    mnorm = mcloud / (dx**3)
    plotpdfrs.append(plotpdf/mnorm)
    plotvrs.append(plotv)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.plot(plotvs[0],plotpdfs[0],'k',plotvs[1],plotpdfs[1],'r',plotvs[2],plotpdfs[2],'b',plotvs[3],plotpdfs[3],'g')
#plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([0.0,15.0,0.0,0.01])
#plt.xticks([])

plt.subplot(1,2,2)
plt.plot(plotvrs[0],plotpdfrs[0],'k',plotvrs[1],plotpdfrs[1],'r',plotvrs[2],plotpdfrs[2],'b',plotvrs[3],plotpdfrs[3],'g')
#plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([-15.0,15.0,0.0,0.01])
#plt.xticks([])

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()

