from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 50
nconvall = 10
nconvallh = 50
bzwrong = 1

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_B5.0_LD_Test_SDs/', 'UV_M5.0e4_R15.0_N128_Tf4_B0.05_LD_Test_SDs/']
dftlist = [1,0,0]
betas = [50.0, 5.0, 0.05]
pdffile = 'sdmasspdfallcirc.dat'
xpdffile = 'sdpdfxy.dat'
nds = len(dflist)

plottimes = []
plotsigmas = []
plotxs = []
ptimes = []
pmus = []
pbs = []
pbdbs = []
pmachs = []

for i in xrange(0,nds):

    datafolder = dflist[i]
    dftype = dftlist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    xpdflines = read_pdffile(hostname,datadir + datafolder + xpdffile)

    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    msol = out_msol(outlines)
    dx = out_dx(outlines)
    rcloud = out_rcloud(outlines)
    sigmacloud = out_sigma(outlines)
    sigmaadj = sigmacloud * (1.0 - 0.12)
    cs = out_csound(outlines)
    gauss = out_gauss(outlines)
    gravc = out_G(outlines)
    rhobar = out_rhocloud(outlines)
    ketot = hst_ke(hstdata)
    beta = betas[i]

    time = hst_time(hstdata, tff)

    bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / beta)
    if bzwrong == 1:
        bz = bz * np.sqrt(4.0 * np.pi)
    if dftype == 1:
        bztote = np.zeros(len(time))
        bztot = np.zeros(len(time))
        bz = 0.0
    else:
        bztote = hst_num(hstdata, 18)
        bztot = hst_num(hstdata, 19)

    edbe = bztote
    ebtote = 0.5 * (bz**2 / (4.0 * np.pi)) * (4.0*rcloud)**3 + edbe
    edb = bztot
    ebtot = 0.5 * (bz**2 / (4.0 * np.pi)) * 4.0 / 3.0 * np.pi * (rcloud)**3 + edb
    #ebcheck = 0.5 * (bztot / (4.0 * rcloud)**3)**2

    bmean = np.sqrt(4.0 * np.pi * 2.0 * ebtot / (4.0*rcloud)**3)
    bmeandb = np.sqrt(4.0 * np.pi * 2.0 * edb / (4.0*rcloud)**3)

    egrav = (3.0 * gravc * mcloud**2 / (5.0 * rcloud))
    eke = ketot

    mfall = 2 * np.sqrt(gravc) * mcloud / (bmean * rcloud**2)
    ptimes.append(time)
    pmus.append(mfall)
    pbs.append(ebtot / egrav)
    pbdbs.append(eke / egrav)
    
    machn = hst_mach(hstdata, cs)
    pmachs.append(machn)

    pdftime = sim_pdftime(pdflines, tff)
    nts = len(pdftime)
    
    plottime = []
    plotsigma = []
    plotx = []
    
    startflag = False
    
    print i, pdffile
    for j in xrange(75,76):
        
        tmini = np.abs(time - pdftime[j]).argmin()

        plotpdf = sim_pdf(pdflines,j)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        mtot = np.sum(plotpdf) / 1.0e9
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)
        
        xplotpdf = sim_pdfm(xpdflines,j)
        xplotpdf = np.convolve(xplotpdf, np.ones((nconv,))/nconv, mode='same')
        xmtot = np.sum(xplotpdf) / 1.0e9
        xplotpdf = xplotpdf / np.sum(xplotpdf)
        xnvals = np.sum(xplotpdf)

        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden / msol
        
        xplotlogden = sim_pdfx(xpdflines)
        xplotden = 10**xplotlogden / msol
        if startflag:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(xplotpdf,xplotlogden,xplotpdf,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
            amph, muh, sigmah, h3, h4, chisq = sim_fitlnhermpdfvals(xplotpdf,xplotlogden,xplotpdf,0.01,0.95,nvals,pin=[ampold,muold,sigmaold,0.0,0.0])
        else:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(xplotpdf,xplotlogden,xplotpdf,0.15,0.95,nvals)
            amph, muh, sigmah, h3, h4, chisqh = sim_fitlnhermpdfvals(xplotpdf,xplotlogden,xplotpdf,0.01,0.95,nvals)
        
        #print j, pdftime[j], mu, sigma, mtot

        if (sigma > 0.0):
            
            [ampold, sigmaold, muold] = [amp, sigma, mu]
            sigma = sigma / np.log10(np.exp(1))
            mu = mu / np.log10(np.exp(1))
            mean = np.exp(mu - 0.5 * sigma**2)
            
            plottime.append(pdftime[j])
            startflag = True
        # else:
            # print sigmam, h3m, h4m, chisqm
            
            ## Calculate sigma from the distribution itself
            mulogten = np.sum(plotlogden * plotpdf)
            sigmalogten = np.sqrt(np.sum(plotlogden * plotlogden * plotpdf) - mulogten**2)
            sigma = sigmalogten / np.log10(np.exp(1))
            x = np.sqrt(sigmacloud * msol * mtot / (mean))
            print j, sigma, x
            plotsigma.append(sigma)
            plotx.append(x)

    plotxs.append(plotx)
    plottimes.append(plottime)
    plotsigmas.append(np.convolve(plotsigma, np.ones((nconvall,))/nconvall, mode='same'))

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axa = plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
pa, pb, pc = plt.plot(ptimes[0],pbs[0],'k',ptimes[1],pbs[1],'r',ptimes[2],pbs[2],'b')
pa, pb, pc = plt.plot(ptimes[0],pbdbs[0],'--k',ptimes[1],pbdbs[1],'--r',ptimes[2],pbdbs[2],'--b')
plt.legend((pa, pb, pc), (r"$\displaystyle {\rm Run~I}$",r"$\displaystyle {\rm Run~MIII}$",r"$\displaystyle {\rm Run~MI}$"),prop={'size':8},loc=1)
plt.axis([0.0,3.0,1.0e-2,1.0e2])
plt.xticks([0,1,2,3])
#plt.yticks([0,20,40])
plt.yscale('log')
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle E / E_{\rm G,0}$")

#axb = plt.subplot(1,2,2)
#axb.yaxis.set_label_position("right")
#pa, pb, pc = plt.plot(ptimes[0],pmachs[0],'k',ptimes[1],pmachs[1],'r',ptimes[2],pmachs[2],'b')
#plt.axis([0.0,3.0,0,30])
#plt.xticks([0,1,2,3])
#plt.yticks([0,10,20,30])
#plt.yscale('log')
#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle \mathcal{M}$")

pp = PdfPages(plotdir + 'ft6-7b.pdf')
pp.savefig()
pp.close()


