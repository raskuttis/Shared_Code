from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 25

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Alt_SDs/'
pdflist = ['sdpdfxy.dat', 'sdmasspdfallcirc.dat', 'sdmasspdfmeancirc.dat']
masslist = [0,1,1,1]
nds = len(pdflist)

tlist = [0.1,0.5,0.9,0.5]
ttypelist = [1,1,1,0]
nts = len(tlist)

plotdens = []
plotpdfs = []
plotfits = []

for i in xrange(0,nds):

    pdffile = pdflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    msol = out_msol(outlines)

    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    mgas = hst_mgas(hstdata, mcloud)
    mof = mgas[0] - mgas - mstar
    eff = hst_eff(hstdata, mcloud)
    pdftime = sim_pdftime(pdflines, tff)
    
    for j in xrange(0,nts):

        if ttypelist[j] == 1:
            tmini = np.abs(mstar - tlist[j] * eff).argmin()
        else:
            tmini = np.abs(mof - tlist[j] * (1.0 - eff)).argmin()
        tmin = time[tmini]
        tminpdfi = np.abs(pdftime - tmin).argmin()
        print i, j, tmin, tminpdfi

        if (masslist[i] == 1):
            plotpdf = sim_pdf(pdflines,tminpdfi)
            plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        else:
            plotpdf = sim_pdfm(pdflines,tminpdfi)
            plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)

        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden / msol
        plotfit = sim_fitlnpdf(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)

        plotdens.append(plotden)
        plotpdfs.append(plotpdf)
        plotfits.append(plotfit)

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)


plt.subplot(2,2,1)
pfid, pall, pmean = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[4],plotpdfs[4],'r',plotdens[8],plotpdfs[8],'b')
plt.plot(plotdens[0],plotfits[0],'--k',plotdens[4],plotfits[4],'--r',plotdens[8],plotfits[8],'--b')
plt.legend((pfid, pall, pmean), (r"$\displaystyle \Sigma$",r"$\displaystyle \Sigma^{c}_{all}$",r"$\displaystyle \Sigma^{c}_{*}$"),prop={'size':8})
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(a): \varepsilon(t) = 0.1 \varepsilon_{\rm final}$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1])

#plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(2,2,2)
pfid, pall, pmean = plt.plot(plotdens[1],plotpdfs[1],'k',plotdens[5],plotpdfs[5],'r',plotdens[9],plotpdfs[9],'b')
plt.plot(plotdens[1],plotfits[1],'--k',plotdens[5],plotfits[5],'--r',plotdens[9],plotfits[9],'--b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b): \varepsilon(t) = 0.5 \varepsilon_{\rm final}$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1],[' ',' ',' ',' '])

plt.subplot(2,2,3)
pfid, pall, pmean = plt.plot(plotdens[2],plotpdfs[2],'k',plotdens[6],plotpdfs[6],'r',plotdens[10],plotpdfs[10],'b')
plt.plot(plotdens[2],plotfits[2],'--k',plotdens[6],plotfits[6],'--r',plotdens[10],plotfits[10],'--b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(c): \varepsilon(t) = 0.9 \varepsilon_{\rm final}$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(2,2,4)
pfid, pall, pmean = plt.plot(plotdens[3],plotpdfs[3],'k',plotdens[7],plotpdfs[7],'r',plotdens[11],plotpdfs[11],'b')
plt.plot(plotdens[3],plotfits[3],'--k',plotdens[7],plotfits[7],'--r',plotdens[11],plotfits[11],'--b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(d): \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of, final}$")
plt.axis([1.0e-2,1.0e4,1.0e-4,1.0e-1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1],[' ',' ',' ',' '])

plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
# plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'flnstar.pdf')
pp.savefig()
pp.close()


