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

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 25
nconvall = 50

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
pdflist = ['sdpdfxy.dat', 'sdpdfallcirc.dat']
xdenomlist = [1.0, 4.0]
nds = len(pdflist)

plottimes = []
plotsigmas = []
plotmeans = []
ploth3s = []
ploth4s = []
plotxs = []

for i in xrange(0,nds):

    pdffile = pdflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    msol = out_msol(outlines)
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
    plotsigmam = []
    plotmeanm = []
    ploth3 = []
    ploth4 = []
    ploth3m = []
    ploth4m = []
    plotx = []
    plotxm = []
    startflag = False
    xdenom = xdenomlist[i]
    
    print i, pdffile
    for j in xrange(100,5100):
        
        tmini = np.abs(time - pdftime[j]).argmin()
        epst = mstar[tmini]

        plotpdf = sim_pdf(pdflines,j)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        plotpdfm = sim_pdfm(pdflines,j)
        plotpdfm = np.convolve(plotpdfm, np.ones((nconv,))/nconv, mode='same')
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)

        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden / msol
        if startflag:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdfm,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
            ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(plotpdfm,plotlogden,plotpdfm,0.15,0.95,nvals,pin=[ampmold,mumold,sigmamold])
            amph, muh, sigmah, h3, h4, chisq = sim_fitlnhermpdfvals(plotpdf,plotlogden,plotpdfm,0.01,0.95,nvals,pin=[ampold,muold,sigmaold,0.0,0.0])
            ampmh, mumh, sigmamh, h3m, h4m, chisqm = sim_fitlnhermpdfvals(plotpdfm,plotlogden,plotpdfm,0.01,0.95,nvals,pin=[ampmold,mumold,sigmamold,0.0,0.0])
        else:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdfm,0.15,0.95,nvals)
            ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(plotpdfm,plotlogden,plotpdfm,0.15,0.95,nvals)
            amph, muh, sigmah, h3, h4, chisqh = sim_fitlnhermpdfvals(plotpdf,plotlogden,plotpdfm,0.01,0.95,nvals)
            ampmh, mumh, sigmamh, h3m, h4m, chisqmh = sim_fitlnhermpdfvals(plotpdfm,plotlogden,plotpdfm,0.01,0.95,nvals)
    
        if (sigmam > 0.0 and sigma > 0.0):
            
            [ampold, sigmaold, muold] = [amp, sigma, mu]
            [ampmold, sigmamold, mumold] = [ampm, sigmam, mum]
            sigma = sigma / np.log10(np.exp(1))
            mu = mu / np.log10(np.exp(1))
            sigmam = sigmam / np.log10(np.exp(1))
            mum = mum / np.log10(np.exp(1))
            mean = np.exp(mu + 0.5 * sigma**2)
            meanm = np.exp(mum - 0.5 * sigmam**2)
            x = np.sqrt(sigmaadj * msol * (1.0 - epst) / (mean * xdenom))
            xm = np.sqrt(sigmaadj * msol * (1.0 - epst) / (meanm * xdenom))
            
            plottime.append(pdftime[j])
            plotsigma.append(sigma)
            plotmean.append(mean / msol)
            plotsigmam.append(sigmam)
            plotmeanm.append(meanm / msol)
            ploth3.append(h3)
            ploth4.append(h4)
            ploth3m.append(h3m)
            ploth4m.append(h4m)
            plotx.append(x)
            plotxm.append(xm)
            startflag = True

            print j, pdftime[j], mu, sigma, mum, sigmam, x, xm, epst

        # else:
            # print sigmam, h3m, h4m, chisqm

    plottimes.append(plottime)
    plotsigmas.append(np.convolve(plotsigma, np.ones((nconvall,))/nconvall, mode='same'))
    plotmeans.append(np.convolve(plotmean, np.ones((nconvall,))/nconvall, mode='same'))
    plottimes.append(plottime)
    plotsigmas.append(np.convolve(plotsigmam, np.ones((nconvall,))/nconvall, mode='same'))
    plotmeans.append(np.convolve(plotmeanm, np.ones((nconvall,))/nconvall, mode='same'))
    ploth3s.append(np.convolve(ploth3, np.ones((nconvall,))/nconvall, mode='same'))
    ploth4s.append(np.convolve(ploth4, np.ones((nconvall,))/nconvall, mode='same'))
    ploth3s.append(np.convolve(ploth3m, np.ones((nconvall,))/nconvall, mode='same'))
    ploth4s.append(np.convolve(ploth4m, np.ones((nconvall,))/nconvall, mode='same'))
    plotxs.append(np.convolve(plotx, np.ones((nconvall,))/nconvall, mode='same'))
    plotxs.append(np.convolve(plotxm, np.ones((nconvall,))/nconvall, mode='same'))

plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.subplots_adjust(left=0.2)
pfid, pall = plt.plot(plottimes[1],plotmeans[1],'k',plottimes[3],plotmeans[3],'r')
#pfid, pall, pmean = plt.plot(plottimes[0],plotmeans[0],'--k',plottimes[2],plotmeans[2],'--r',plottimes[4],plotmeans[4],'--b')
plt.legend((pfid, pall), (r"$\displaystyle \Sigma$",r"$\displaystyle \Sigma^{c}_{all}$"))
plt.text(0.05*2,10**(0.9*2.7)-3.0,r"$\displaystyle(a)$")
plt.yscale('log')
plt.axis([0.0,2.0,3.0,500.0])
plt.yticks([1.0e1,1.0e2])
plt.xticks([0,1,2],[' ',' ',' '])
#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \langle \Sigma \rangle_{\rm cloud} / M_{\odot} {\rm pc^{-2}}$")

axsigma = plt.subplot(2,2,2)
pfid, pall = plt.plot(plottimes[1],plotsigmas[1],'k',plottimes[3],plotsigmas[3],'r')
#pfid, pall, pmean = plt.plot(plottimes[0],plotsigmas[0],'--k',plottimes[2],plotsigmas[2],'--r',plottimes[4],plotsigmas[4],'--b')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axsigma.yaxis.set_label_position("right")
plt.text(0.05*2,0.9*3,r"$\displaystyle(c)$")
plt.axis([0.0,2.0,0,3])
plt.xticks([0,1,2],[' ',' ',' '])
plt.yticks([0,1,2,3])
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")

axh3 = plt.subplot(2,2,3)
pfid, pall = plt.plot(plottimes[1],ploth3s[1],'k',plottimes[3],ploth3s[3],'r')
plt.text(0.05*2,0.9*0.2-0.1,r"$\displaystyle(d)$")
plt.axis([0,2,-0.1,0.1])
plt.xticks([0,1,2])
plt.yticks([-0.1,0,0.1])
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle h_3$")

axh3 = plt.subplot(2,2,4)
pfid, pall = plt.plot(plottimes[1],plotxs[1],'k',plottimes[3],plotxs[3],'r')
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
axh3.yaxis.set_label_position("right")
plt.text(0.05*2,0.9*0.2-0.1,r"$\displaystyle(d)$")
plt.axis([0,2,0.5,1.5])
plt.xticks([0,1,2])
plt.yticks([0.5,1.0,1.5])
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle x$")

pp = PdfPages(plotdir + 'flnstarfitx.pdf')
pp.savefig()
pp.close()


