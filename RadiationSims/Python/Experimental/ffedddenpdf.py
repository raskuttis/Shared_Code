from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_SDs_Alt/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tlist = [1.0,1.5,2.0,2.5]
sdlist = [0.1,1.0,10.0,100.0]
nts = len(tlist)
nsds = len(sdlist)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
msol = out_msol(outlines)
dx = out_dx(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)

pdffile = 'denpdfall.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'fedddenpdfall.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'feddsigmadenpdfall.dat'
vspdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdens = []
plotpdfs = []
nconv = 75
for i in xrange(0,nts):

    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plotpdfs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same'))
    plotdens.append(plotden)

sdtimes = []
sdvs = []
for i in xrange(0,nsds):
    
    npdfts = len(pdftime)
    sdind = np.abs(plotden - sdlist[i]).argmin()
    sdtime = []
    sdv = []
    
    for j in xrange(10,npdfts-10):
        
        plotpdf = sim_pdf_zip(vpdflines,j)
        plotdenpdf = sim_pdf(pdflines,j)
        plotpdf = plotpdf / plotdenpdf
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        sdv.append(plotpdf[sdind])
        sdtime.append(pdftime[j])
    
    sdtimes.append(sdtime)
    sdvs.append(np.convolve(sdv, np.ones((nconv,))/nconv, mode='same'))

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[1],plotpdfs[1],'r',plotdens[2],plotpdfs[2],'b',plotdens[3],plotpdfs[3],'g')
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,1.0e-1,1.0e2])
#plt.xticks([])

plt.subplot(1,2,2)
plt.plot(sdtimes[0],sdvs[0],'k',sdtimes[1],sdvs[1],'r',sdtimes[2],sdvs[2],'b',sdtimes[3],sdvs[3],'g')
#plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([0,3,0.0,100.0])
#plt.xticks([])

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()

