from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 25

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_NF_B50.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.5_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.05_LD_SDs/']
pdflist = ['sdpdfxy.dat', 'sdpdfyz.dat', 'sdpdfallcirc.dat', 'sdpdfmeancirc.dat']
masslist = [0,0,1,1]
nps = len(pdflist)
nds = len(dflist)

plotdens = []
plotpdfs = []
plotfits = []

for i in xrange(0,nds):

    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    msol = out_msol(outlines)
    kappa = 7.228571e-03

    time = hst_time(hstdata, tff)
    mstar = hst_num(hstdata, 13) / mcloud
    
    for j in xrange(0,nps):

        pdffile = pdflist[j]
        pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
        pdftime = sim_pdftime(pdflines, tff)
        tmini = np.abs(mstar - 0.02).argmin()
        tmin = time[tmini]
        tminpdfi = np.abs(pdftime - tmin).argmin()
        print i, j, tmin, tminpdfi

        plotpdf = sim_pdf(pdflines,tminpdfi)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        nvals = np.sum(plotpdf)
        plotpdf = plotpdf / np.sum(plotpdf)

        plotlogden = sim_pdfx(pdflines)
        fnorm = plotlogden[1] - plotlogden[0]
        plotden = 10**plotlogden
        plottau = plotden * kappa
        #plotfit = sim_fitlnpdf(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)

        plotdens.append(plottau)
        plotpdfs.append(plotpdf / fnorm)
        #plotfits.append(plotfit / fnorm)

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-2.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)


plt.subplot(2,2,1)
pa, pb, pc = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[4],plotpdfs[4],'r',plotdens[8],plotpdfs[8],'b')
#plt.plot(plotdens[0],plotfits[0],'--k',plotdens[4],plotfits[4],'--r',plotdens[8],plotfits[8],'--b')
plt.legend((pa, pb, pc), (r"$\displaystyle \mu_{\Phi,0} \rightarrow \infty$",r"$\displaystyle \mu_{\Phi,0} = 9.5$",r"$\displaystyle \mu_{\Phi,0} = 3.0$"),prop={'size':8},loc=1)
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(a): \varepsilon(t) = 0.1 \varepsilon_{\rm final}$")
plt.text(xcap,ycap,r"$\displaystyle(a): \Sigma_{\rm x-y}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1])

#plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_A$")

plt.subplot(2,2,2)
pfid, pall, pmean = plt.plot(plotdens[1],plotpdfs[1],'k',plotdens[5],plotpdfs[5],'r',plotdens[9],plotpdfs[9],'b')
#plt.plot(plotdens[1],plotfits[1],'--k',plotdens[5],plotfits[5],'--r',plotdens[9],plotfits[9],'--b')
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(b): \varepsilon(t) = 0.5 \varepsilon_{\rm final}$")
plt.text(xcap,ycap,r"$\displaystyle(b): \Sigma_{\rm y-z}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1],[' ',' ',' ',' '])

plt.subplot(2,2,3)
pfid, pall, pmean = plt.plot(plotdens[2],plotpdfs[2],'k',plotdens[6],plotpdfs[6],'r',plotdens[10],plotpdfs[10],'b')
#plt.plot(plotdens[2],plotfits[2],'--k',plotdens[6],plotfits[6],'--r',plotdens[10],plotfits[10],'--b')
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(c): \varepsilon(t) = 0.9 \varepsilon_{\rm final}$")
plt.text(xcap,ycap,r"$\displaystyle(c): \Sigma_{\rm all}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1])

plt.xlabel(r"$\displaystyle \tau$")
plt.ylabel(r"$\displaystyle P_\Omega$")

plt.subplot(2,2,4)
pfid, pall, pmean = plt.plot(plotdens[3],plotpdfs[3],'k',plotdens[7],plotpdfs[7],'r',plotdens[11],plotpdfs[11],'b')
#plt.plot(plotdens[3],plotfits[3],'--k',plotdens[7],plotfits[7],'--r',plotdens[11],plotfits[11],'--b')
plt.xscale('log')
plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(d): \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of, final}$")
plt.text(xcap,ycap,r"$\displaystyle(d): \Sigma_{*}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1],[' ',' ',' ',' '])

plt.xlabel(r"$\displaystyle \tau$")
# plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'ft6-6.pdf')
pp.savefig()
pp.close()


