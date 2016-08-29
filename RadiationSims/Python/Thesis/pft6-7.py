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

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_NF_B50.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.5_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.05_LD_SDs/']
pdffile = 'sdmasspdfallcirc.dat'
nds = len(dflist)

plottimes = []
plotsigmas = []

for i in xrange(0,nds):

    datafolder = dflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    msol = out_msol(outlines)
    dx = out_dx(outlines)
    rcloud = out_rcloud(outlines)
    sigmacloud = out_sigma(outlines)
    sigmaadj = sigmacloud * (1.0 - 0.12)

    time = hst_time(hstdata, tff)

    pdftime = sim_pdftime(pdflines, tff)
    nts = len(pdftime)
    
    plottime = []
    plotsigma = []
    
    startflag = False
    
    print i, pdffile
    for j in xrange(75,509):
        
        tmini = np.abs(time - pdftime[j]).argmin()

        plotpdf = sim_pdf(pdflines,j)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        mtot = np.sum(plotpdf) / 1.0e9
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)

        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden / msol
        if startflag:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
            amph, muh, sigmah, h3, h4, chisq = sim_fitlnhermpdfvals(plotpdf,plotlogden,plotpdf,0.01,0.95,nvals,pin=[ampold,muold,sigmaold,0.0,0.0])
        else:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)
            amph, muh, sigmah, h3, h4, chisqh = sim_fitlnhermpdfvals(plotpdf,plotlogden,plotpdf,0.01,0.95,nvals)
        
        #print j, pdftime[j], mu, sigma, mtot

        if (sigma > 0.0):
            
            [ampold, sigmaold, muold] = [amp, sigma, mu]
            sigma = sigma / np.log10(np.exp(1))
            
            plottime.append(pdftime[j])
            startflag = True
        # else:
            # print sigmam, h3m, h4m, chisqm
            
            ## Calculate sigma from the distribution itself
            mulogten = np.sum(plotlogden * plotpdf)
            sigmalogten = np.sqrt(np.sum(plotlogden * plotlogden * plotpdf) - mulogten**2)
            sigma = sigmalogten / np.log10(np.exp(1))
            print j, sigma
            plotsigma.append(sigma)

    plottimes.append(plottime)
    plotsigmas.append(np.convolve(plotsigma, np.ones((nconvall,))/nconvall, mode='same'))

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axsigma = plt.subplot(1,1,1)
pa, pb, pc = plt.plot(plottimes[0],plotsigmas[0],'k',plottimes[1],plotsigmas[1],'r',plottimes[2],plotsigmas[2],'b')
plt.legend((pa, pb, pc), (r"$\displaystyle \mu_{\Phi,0} \rightarrow \infty$",r"$\displaystyle \mu_{\Phi,0} = 9.5$",r"$\displaystyle \mu_{\Phi,0} = 3.0$"),prop={'size':8},loc=1)
plt.subplots_adjust(bottom=0.2)
plt.axis([0.0,2.0,0,3])
plt.xticks([0,1,2])
plt.yticks([0,1,2,3])
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")

pp = PdfPages(plotdir + 'ft6-7.pdf')
pp.savefig()
pp.close()


