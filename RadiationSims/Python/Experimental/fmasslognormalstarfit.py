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
nconvall = 50
nconvallh = 50

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_SDs/'
pdflist = ['sdpdfxy.dat', 'sdmasspdfallcirc.dat', 'sdmasspdfmeancirc.dat']
xdenomlist = [1.0, 1.0, 1.0]
masslist = [0,1,1]
nds = len(pdflist)

plottimes = []
plotsigmas = []
plotmeans = []
ploth3s = []
ploth4s = []
plotmasses = []

for i in xrange(0,nds):

    pdffile = pdflist[i]
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
    mstar = hst_mstar(hstdata, mcloud)
    eff = hst_eff(hstdata, mcloud)

    pdftime = sim_pdftime(pdflines, tff)
    nts = len(pdftime)
    
    plottime = []
    plotsigma = []
    plotmean = []
    ploth3 = []
    ploth4 = []
    plotmass = []
    xdenom = xdenomlist[i]
    
    startflag = False
    
    print i, pdffile
    for j in xrange(1250,5100):
        
        tmini = np.abs(time - pdftime[j]).argmin()
        epst = mstar[tmini]

        if (masslist[i]==1):
            plotpdf = sim_pdf(pdflines,j)
            plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
            mtot = np.sum(plotpdf) / 1.0e9
        else:
            plotpdf = sim_pdfm(pdflines,j)
            plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
            mtot = sim_pdf_mass(pdflines,j,(4.0*rcloud)**2,mcloud)
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
            mu = mu / np.log10(np.exp(1))
            mean = np.exp(mu - 0.5 * sigma**2)
            x = np.sqrt(sigmaadj * msol * (1.0 - epst) / (mean * xdenom))
            print j, x
            
            plottime.append(pdftime[j])
            plotsigma.append(sigma)
            plotmean.append(mean / msol)
            ploth3.append(x)
            ploth4.append(h4)
            plotmass.append(mtot)
            
            startflag = True
        # else:
            # print sigmam, h3m, h4m, chisqm

    plottimes.append(plottime)
    plotsigmas.append(np.convolve(plotsigma, np.ones((nconvall,))/nconvall, mode='same'))
    plotmeans.append(np.convolve(plotmean, np.ones((nconvall,))/nconvall, mode='same'))
    plottimes.append(plottime)
    ploth3s.append(np.convolve(ploth3, np.ones((nconvall,))/nconvallh, mode='same'))
    ploth4s.append(np.convolve(ploth4, np.ones((nconvall,))/nconvallh, mode='same'))
    plotmasses.append(plotmass)

plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.subplots_adjust(left=0.2)
pfid, pall, pmean = plt.plot(plottimes[0],plotmeans[0],'k',plottimes[1],plotmeans[1],'r',plottimes[2],plotmeans[2])
#pfid, pall, pmean = plt.plot(plottimes[0],plotmeans[0],'--k',plottimes[2],plotmeans[2],'--r',plottimes[4],plotmeans[4],'--b')
plt.legend((pfid, pall, pmean), (r"$\displaystyle \Sigma$",r"$\displaystyle \Sigma^{c}_{all}$",r"$\displaystyle \Sigma^{c}_{*}$"),prop={'size':8})
plt.text(0.05*3,10**(0.9*2.7)-3.0,r"$\displaystyle(a)$")
plt.yscale('log')
plt.axis([0.0,3.0,3.0,500.0])
plt.yticks([1.0e1,1.0e2])
plt.xticks([0,1,2,3],[' ',' ',' ',' '])
#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \langle \Sigma \rangle_{\rm cloud} / M_{\odot} {\rm pc^{-2}}$")

axsigma = plt.subplot(2,2,2)
pfid, pall, pmean = plt.plot(plottimes[0],plotsigmas[0],'k',plottimes[1],plotsigmas[1],'r',plottimes[2],plotsigmas[2],'b')
#pfid, pall, pmean = plt.plot(plottimes[0],plotsigmas[0],'--k',plottimes[2],plotsigmas[2],'--r',plottimes[4],plotsigmas[4],'--b')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axsigma.yaxis.set_label_position("right")
plt.text(0.05*3,0.9*3,r"$\displaystyle(b)$")
plt.axis([0.0,3.0,0,3])
plt.xticks([0,1,2,3],[' ',' ',' ',' '])
plt.yticks([0,1,2,3])
#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")

axh3 = plt.subplot(2,2,3)
pfid, pall, pmean = plt.plot(plottimes[0],ploth3s[0],'k',plottimes[1],ploth3s[1],'r',plottimes[2],ploth3s[2],'b')
plt.text(0.05*3,0.9*2,r"$\displaystyle(c)$")
plt.axis([0,3,0.0,2.0])
plt.xticks([0,1,2,3])
plt.yticks([0.0,1.0,2.0])
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle x$")

axh3 = plt.subplot(2,2,4)
pfid, pall, pmean = plt.plot(plottimes[0],plotmasses[0],'k',plottimes[1],plotmasses[1],'r',plottimes[2],plotmasses[2],'b')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axh3.yaxis.set_label_position("right")
plt.text(0.05*3,0.9*1.1,r"$\displaystyle(d)$")
plt.axis([0,3,0.0,1.1])
plt.xticks([0,1,2,3])
plt.yticks([0.0,0.5,1.0])
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl, 0}$")

pp = PdfPages(plotdir + 'flnstarfit.pdf')
pp.savefig()
pp.close()


