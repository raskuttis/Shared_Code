from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
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

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_SDs/'
pdflist = ['sdpdfxy.dat', 'sdpdfallcirc.dat', 'sdpdfmeancirc.dat']
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
    kappa = out_kappa(outlines)

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
        tmini = np.abs(mstar - 0.1 * eff).argmin()
        tminpdfi = np.abs(pdftime - tmin).argmin()

        plotpdf = sim_pdf(pdflines,tminpdfi)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        nvals = np.sum(plotpdf)
        plotpdf = plotpdf / np.sum(plotpdf)

        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden
        plottau = plotden * kappa

        plotdens.append(plottau)
        plotpdfs.append(plotpdf)

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_NF_SDs/'
pdflist = ['sdpdfxy.dat', 'sdpdfallcirc.dat', 'sdpdfmeancirc.dat']
nds = len(pdflist)

altplotdens = []
altplotpdfs = []
altplotfits = []

for i in xrange(0,nds):
    
    pdffile = pdflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    msol = out_msol(outlines)
    kappa = out_kappa(outlines)
    
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
        tmini = np.abs(mstar - 0.1 * eff).argmin()
        tminpdfi = np.abs(pdftime - tmin).argmin()
    
        plotpdf = sim_pdf(pdflines,tminpdfi)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        nvals = np.sum(plotpdf)
        plotpdf = plotpdf / np.sum(plotpdf)
    
        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden
        plottau = plotden * kappa
    
        altplotdens.append(plottau)
        altplotpdfs.append(plotpdf)

xcap = 10**(-4.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
pall, pmean = plt.plot(plotdens[4],plotpdfs[4],'r',plotdens[8],plotpdfs[8],'b')
plt.plot(altplotdens[4],altplotpdfs[4],'--r')
plt.legend((pall, pmean), (r"$\displaystyle \Sigma^{c}_{all}$",r"$\displaystyle \Sigma^{c}_{*}$"),prop={'size':8})
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(a): \varepsilon(t) = 0.1 \varepsilon_{\rm final}$")
plt.axis([1.0e-4,1.0e2,1.0e-4,1.0e-1])
plt.xticks([1.0e-4,1.0e-2,1.0,1.0e2],[' ',' ',' ',' '])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1])

#plt.xlabel(r"$\displaystyle \tau$")
plt.ylabel(r"$\displaystyle P_{\Omega}$")

plt.subplot(2,2,2)
pall, pmean = plt.plot(plotdens[5],plotpdfs[5],'r',plotdens[9],plotpdfs[9],'b')
plt.plot(altplotdens[5],altplotpdfs[5],'--r')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b): \varepsilon(t) = 0.5 \varepsilon_{\rm final}$")
plt.axis([1.0e-4,1.0e2,1.0e-4,1.0e-1])
plt.xticks([1.0e-4,1.0e-2,1.0,1.0e2],[' ',' ',' ',' '])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1],[' ',' ',' ',' '])

#plt.xlabel(r"$\displaystyle \tau$")
#plt.ylabel(r"$\displaystyle P_{\Omega}$")

plt.subplot(2,2,3)
pall, pmean = plt.plot(plotdens[6],plotpdfs[6],'r',plotdens[10],plotpdfs[10],'b')
plt.plot(altplotdens[6],altplotpdfs[6],'--r')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(c): \varepsilon(t) = 0.9 \varepsilon_{\rm final}$")
plt.axis([1.0e-4,1.0e2,1.0e-4,1.0e-1])
plt.xticks([1.0e-4,1.0e-2,1.0,1.0e2])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1])

plt.xlabel(r"$\displaystyle \tau$")
plt.ylabel(r"$\displaystyle P_{\Omega}$")

plt.subplot(2,2,4)
pall, pmean = plt.plot(plotdens[7],plotpdfs[7],'r',plotdens[11],plotpdfs[11],'b')
plt.plot(altplotdens[7],altplotpdfs[7],'--r')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(d): \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of, final}$")
plt.axis([1.0e-4,1.0e2,1.0e-4,1.0e-1])
plt.xticks([1.0e-4,1.0e-2,1.0,1.0e2])
plt.yticks([1.0e-4,1.0e-3,1.0e-2,1.0e-1],[' ',' ',' ',' '])

plt.xlabel(r"$\displaystyle \tau$")
#plt.ylabel(r"$\displaystyle P_{\Omega}$")

pp = PdfPages(plotdir + 'ftau.pdf')
pp.savefig()
pp.close()


