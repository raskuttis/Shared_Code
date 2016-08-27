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

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_Fluxes/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tlist = [0.3, 0.6, 0.8, 1.5]
nts = len(tlist)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
msol = out_msol(outlines)
dx = out_dx(outlines)

time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
eff = hst_eff(hstdata, mcloud)

print 'Reading PDF1'
pdffile = 'denpdfall_r0_r1.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdens = []
plotpdfs = []
rawplotpdfs = []
nconv = 25
for i in xrange(0,nts):

    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    mnorm = mcloud / (dx**3)
    rawplotpdfs.append(plotpdf[1:])
    plotpdf = plotpdf * plotden
    plotpdfs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same') / np.sum(plotpdf))
    plotdens.append(plotden)

print 'Reading PDF2'
pdffile = 'denpdfof_r0_r1.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
plotdenofs = []
plotpdfofs = []
plotpdfifs = []
nconv = 25
for i in xrange(0,nts):
    
    tmini = np.abs(time - tlist[i]).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    plotpdf = sim_pdf(pdflines,tminpdfi)
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden
    mnorm = mcloud / (dx**3)
    plotpdf = plotpdf[1:]
    plotden = plotden[1:]
    plotifpdf = rawplotpdfs[i] - plotpdf
    plotpdf = plotpdf * plotden
    plotifpdf = plotifpdf * plotden
    plotpdfofs.append(np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same') / np.sum(rawplotpdfs[i] * plotden))
    plotpdfifs.append(np.convolve(plotifpdf, np.ones((nconv,))/nconv, mode='same') / np.sum(rawplotpdfs[i] * plotden))
    plotdenofs.append(plotden)

xcap = 10**(0.05*7-1)
ycap = 10**(0.9*5-6)
ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
pof, pif, pall = plt.plot(plotdenofs[0],plotpdfofs[0],'k',plotdenofs[0],plotpdfifs[0],'r',plotdens[0],plotpdfs[0],'b')
plt.xscale('log')
plt.yscale('log')
plt.legend((pof, pif, pall), (r"$\displaystyle v_r > 0$",r"$\displaystyle v_r < 0$",r"$\displaystyle All$"),prop={'size':8},loc=1)
plt.text(xcap,ycap,r"$\displaystyle(a) t = 0.3~t_{\rm ff}$")
plt.axis([1.0e-1,1.0e6,1.0e-6,1.0e-1])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6],[' ',' ',' ',' '])
plt.yticks([1.0e-6,1.0e-4,1.0e-2])
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(2,2,2)
pof, pif, pall = plt.plot(plotdenofs[1],plotpdfofs[1],'k',plotdenofs[1],plotpdfifs[1],'r',plotdens[1],plotpdfs[1],'b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b) t = 0.6~t_{\rm ff}$")
plt.axis([1.0e-1,1.0e6,1.0e-6,1.0e-1])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6],[' ',' ',' ',' '])
plt.yticks([1.0e-6,1.0e-4,1.0e-2],[' ',' ',' '])

plt.subplot(2,2,3)
pof, pif, pall = plt.plot(plotdenofs[2],plotpdfofs[2],'k',plotdenofs[2],plotpdfifs[2],'r',plotdens[2],plotpdfs[2],'b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(c) t = 1.06~t_{\rm ff}$")
plt.axis([1.0e-1,1.0e6,1.0e-6,1.0e-1])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])
plt.yticks([1.0e-6,1.0e-4,1.0e-2])
plt.ylabel(r"$\displaystyle P_M$")
plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")

plt.subplot(2,2,4)
pof, pif, pall = plt.plot(plotdenofs[3],plotpdfofs[3],'k',plotdenofs[3],plotpdfifs[3],'r',plotdens[3],plotpdfs[3],'b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(d) t = 1.57~t_{\rm ff}$")
plt.axis([1.0e-1,1.0e6,1.0e-6,1.0e-1])
plt.xticks([1.0,1.0e2,1.0e4,1.0e6])
plt.yticks([1.0e-6,1.0e-4,1.0e-2],[' ',' ',' '])
plt.xlabel(r"$\displaystyle n_H / {\rm cm^{-3}}$")

pp = PdfPages(plotdir + 'fpdfshell.pdf')
pp.savefig()
pp.close()

