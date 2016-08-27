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

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_SDs/'
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
pdffile = 'vdenpdfall.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vsigmadenpdfall.dat'
vspdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdens = []
plotpdfs = []
plotps = []
plotns = []
nconv = 25
for i in xrange(0,nts):

    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    ploterr = sim_pdf_zip(vspdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    ploterr = ploterr / plotdenpdf
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.convolve(ploterr, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.sqrt(ploterr - plotpdf**2)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plotpdfs.append(plotpdf)
    plotdens.append(plotden)
    plotps.append(plotpdf + ploterr)
    plotns.append(plotpdf - ploterr)

pdffile = 'vrdenpdfall.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vsigmadenpdfall.dat'
vspdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdenrs = []
plotpdfrs = []
plotprs = []
plotnrs = []
for i in xrange(0,nts):
    
    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    ploterr = sim_pdf_zip(vspdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    ploterr = ploterr / plotdenpdf
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.convolve(ploterr, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.sqrt(ploterr - plotpdf**2)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plotpdfrs.append(plotpdf)
    plotdenrs.append(plotden)
    plotprs.append(plotpdf + ploterr)
    plotnrs.append(plotpdf - ploterr)

pdffile = 'ofdenpdfall.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vrofdenpdfall.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vrofsigmadenpdfall.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdenofs = []
plotpdfofs = []
plotpofs = []
plotnofs = []
nconv = 25
for i in xrange(0,nts):
    
    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    ploterr = sim_pdf_zip(vspdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    ploterr = ploterr / plotdenpdf
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.convolve(ploterr, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.sqrt(ploterr - plotpdf**2)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plotpdfofs.append(plotpdf)
    plotdenofs.append(plotden)
    plotpofs.append(plotpdf + ploterr)
    plotnofs.append(plotpdf - ploterr)

pdffile = 'ifdenpdfall.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vrifdenpdfall.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdffile = 'vrifsigmadenpdfall.dat'
vpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdenifs = []
plotpdfifs = []
plotpifs = []
plotnifs = []
nconv = 25
for i in xrange(0,nts):
    
    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf_zip(vpdflines,tminpdfi)
    ploterr = sim_pdf_zip(vspdflines,tminpdfi)
    plotdenpdf = sim_pdf(pdflines,tminpdfi)
    plotpdf = plotpdf / plotdenpdf
    ploterr = ploterr / plotdenpdf
    plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.convolve(ploterr, np.ones((nconv,))/nconv, mode='same')
    ploterr = np.sqrt(ploterr - plotpdf**2)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    plotpdfifs.append(plotpdf)
    plotdenifs.append(plotden)
    plotpifs.append(plotpdf + ploterr)
    plotnifs.append(plotpdf - ploterr)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[1],plotpdfs[1],'r',plotdens[2],plotpdfs[2],'b',plotdens[3],plotpdfs[3],'g')
plt.plot(plotdens[0],plotps[0],'--k',plotdens[1],plotps[1],'--r',plotdens[2],plotps[2],'--b',plotdens[3],plotps[3],'--g')
plt.plot(plotdens[0],plotns[0],'--k',plotdens[1],plotns[1],'--r',plotdens[2],plotns[2],'--b',plotdens[3],plotns[3],'--g')
plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])

plt.subplot(2,2,2)
plt.plot(plotdenrs[0],plotpdfrs[0],'k',plotdenrs[1],plotpdfrs[1],'r',plotdenrs[2],plotpdfrs[2],'b',plotdenrs[3],plotpdfrs[3],'g')
plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])

plt.subplot(2,2,3)
plt.plot(plotdenofs[0],plotpdfofs[0],'k',plotdenofs[1],plotpdfofs[1],'r',plotdenofs[2],plotpdfofs[2],'b',plotdenofs[3],plotpdfofs[3],'g')
plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,0.0,25.0])
#plt.xticks([])

plt.subplot(2,2,4)
plt.plot(plotdenifs[0],plotpdfifs[0],'k',plotdenifs[1],plotpdfifs[1],'r',plotdenifs[2],plotpdfifs[2],'b',plotdenifs[3],plotpdfifs[3],'g')
plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.axis([1.0e-2,1.0e4,-25.0,0.0])
#plt.xticks([])

pp = PdfPages(plotdir + 'ftest.pdf')
pp.savefig()
pp.close()

