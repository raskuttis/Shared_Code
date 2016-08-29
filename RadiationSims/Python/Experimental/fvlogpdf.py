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

tlist = [0.1,0.5,0.9,0.99]
nts = len(tlist)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
msol = out_msol(outlines)
dx = out_dx(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)
nconv = 25

pdffile = 'vmasspdfall.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotvs = []
plotpdfs = []
for i in xrange(0,nts):

    tmini = np.abs(mstar - tlist[i] * eff).argmin()
    tmin = time[tmini]
    print i, tmin
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf_zip(pdflines,tminpdfi)
    plotv = sim_pdfx(pdflines)
    mnorm = mcloud / (dx**3)
    plotpdfs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotvs.append(plotv)

pdffile = 'vrmasspdfall.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotvrs = []
plotpdfrs = []
for i in xrange(0,nts):
    
    tmini = np.abs(mstar - tlist[i] * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(pdflines,tminpdfi)
    plotv = sim_pdfx(pdflines)
    mnorm = mcloud / (dx**3)
    plotpdfrs.append(np.convolve(plotpdf/mnorm, np.ones((nconv,))/nconv, mode='same'))
    plotvrs.append(plotv)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
t1, t2, t3, t4 = plt.plot(plotvs[0],plotpdfs[0],'k',plotvs[1],plotpdfs[1],'r',plotvs[2],plotpdfs[2],'b',plotvs[3],plotpdfs[3],'g')
plt.legend((t1,t2,t3,t4), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.99 \varepsilon_{\rm final}$"),prop={'size':8})
#plt.xscale('log')
plt.yscale('log')
plt.text(0.05*15,10**(0.9*3-5),r"$\displaystyle(a)$")
plt.axis([0.0,20.0,1.0e-5,1.0e-2])
plt.xticks([0.0,5.0,10.0,15.0,20.0])
plt.xlabel(r"$\displaystyle \mid v \mid / {\rm km s^{-1}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(1,2,2)
plt.plot(plotvrs[0],plotpdfrs[0],'k',plotvrs[1],plotpdfrs[1],'r',plotvrs[2],plotpdfrs[2],'b',plotvrs[3],plotpdfrs[3],'g')
#plt.xscale('log')
plt.yscale('log')
plt.text(0.05*30-15,10**(0.9*3-5),r"$\displaystyle(b)$")
plt.axis([-10.0,20.0,1.0e-5,1.0e-2])
plt.yticks([])
plt.xlabel(r"$\displaystyle v_r / {\rm km s^{-1}}$")

pp = PdfPages(plotdir + 'fvlog.pdf')
pp.savefig()
pp.close()

