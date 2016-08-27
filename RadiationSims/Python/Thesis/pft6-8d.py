from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
from ..Hyperion.hyp_models import *
import matplotlib.pyplot as plt

## Plot showing evolution of best-fit lognormal width and x as a function of magnetic field strength

## Define the locations where data is located and where to plot from hyp_models
hostname, datadir, hstfile, outfile, plotdir = init_dirs('tiger')
fname = 'ft6-8d.pdf'
pdffile = 'sdmasspdfmeancirc.dat'
xpdffile = 'sdpdfxy.dat'

## List of models to plot with
dflist, betas = set_mag_model_lookup()
nds = len(dflist)

## Parameters for smoothing PDFs
nconv = 50
nconvall = 50
nconvallh = 50
bzwrong = 1

#pdffile = 'sdmpdfxyplane_8.dat'

tvs = [0.1, 0.5, 0.8]
ntvs = len(tvs)

plotmfs = np.zeros(nds)
plotmfvars = np.zeros((ntvs,nds))
plotsigmas = np.zeros((ntvs,nds))
plotxs = np.zeros((ntvs,nds))

for i in xrange(0,nds):

    datafolder = dflist[i]
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
    gravc = out_G(outlines)
    rhobar = out_rhocloud(outlines)
    beta = betas[i]
    bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / beta)
    if bzwrong == 1:
        bz = bz * np.sqrt(4.0 * np.pi)
    mf0 = 2 * np.sqrt(gravc) * mcloud / (bz * rcloud**2)
    plotmfs[i] = mf0

    time = hst_time(hstdata, tff)
    mstar = hst_num(hstdata, 13) / mcloud
    eff = np.max(mstar)
    bztote = hst_num(hstdata, 17)
    bsq = 2.0 * (bztote / (4.0 * rcloud**3)) + bz**2
    bmean = np.sqrt(bsq)
    mfall = 2 * np.sqrt(gravc) * mcloud / (bmean * rcloud**2)

    pdftime = sim_pdftime(pdflines, tff)
    nts = len(pdftime)
    
    plottime = []
    plotsigma = []
    plotx = []
    
    startflag = False
    
    print i, pdffile, mf0
    for j in xrange(125,390):
        
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
            amp, mu, sigma, chisq = sim_fitlnpdfvals(xplotpdf,xplotlogden,xplotpdf,0.25,0.95,xnvals,pin=[ampold,muold,sigmaold])
            amph, muh, sigmah, h3, h4, chisq = sim_fitlnhermpdfvals(xplotpdf,xplotlogden,xplotpdf,0.01,0.95,xnvals,pin=[ampold,muold,sigmaold,0.0,0.0])
        else:
            amp, mu, sigma, chisq = sim_fitlnpdfvals(xplotpdf,xplotlogden,xplotpdf,0.25,0.95,xnvals)
            amph, muh, sigmah, h3, h4, chisqh = sim_fitlnhermpdfvals(xplotpdf,xplotlogden,xplotpdf,0.01,0.95,xnvals)
        
        if startflag:
            samp, smu, ssigma, schisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.3,0.95,xnvals,pin=[sampold,smuold,ssigmaold])
            samph, smuh, ssigmah, sh3, sh4, schisq = sim_fitlnhermpdfvals(plotpdf,plotlogden,plotpdf,0.01,0.95,xnvals,pin=[sampold,smuold,ssigmaold,0.0,0.0])
        else:
            samp, smu, ssigma, schisq = sim_fitlnpdfvals(plotpdf,plotlogden,xplotpdf,0.3,0.95,xnvals)
            samph, smuh, ssigmah, sh3, sh4, schisq = sim_fitlnhermpdfvals(plotpdf,plotlogden,plotpdf,0.01,0.95,xnvals)
        
        #print j, pdftime[j], mu, sigma, mtot

        if (sigma > 0.0):
            
            [ampold, sigmaold, muold] = [amp, sigma, mu]
            sigma = sigma / np.log10(np.exp(1))
            mu = mu / np.log10(np.exp(1))
            mean = np.exp(mu - 0.5 * sigma**2)
            x = np.sqrt(sigmacloud * msol * mtot / (mean))
            
            [sampold, ssigmaold, smuold] = [samp, ssigma, smu]
            ssigma = ssigma / np.log10(np.exp(1))
            smu = smu / np.log10(np.exp(1))
            smean = np.exp(smu - 0.5 * ssigma**2)
            
            plottime.append(pdftime[j])
            mulogten = np.sum(plotlogden * plotpdf)
            sigmalogten = np.sqrt(np.sum(plotlogden * plotlogden * plotpdf) - mulogten**2)
            sigma = sigmalogten / np.log10(np.exp(1))
            plotsigma.append(sigma**2)
            plotx.append(x)
            #print j, sigma
            startflag = True
        # else:
            # print sigmam, h3m, h4m, chisqm

    for j in xrange(0,ntvs):

        tmini = np.abs(mstar - tvs[j] * eff).argmin()
        plotmfvars[j,i] = mfall[tmini]
        tmin = time[tmini]
        tminpdfi = np.abs(plottime - tmin).argmin()
        plotsigmas[j,i] = plotsigma[tminpdfi]
        plotxs[j,i] = plotx[tminpdfi]
        print j, tmin, plotsigma[tminpdfi], plotx[tminpdfi], np.max(time)

## Linear fits to the variance as a function of magnetic field
ma, ba = np.polyfit(np.log10(plotmfs), plotsigmas[1,:], 1)
thmfs = np.logspace(-1.0, 2.0)
thsigs = np.log10(thmfs) * ma + ba
print ma, ba

mb, bb = np.polyfit(np.log10(plotmfvars[0,:]), plotsigmas[0,:], 1)
thsigbs = np.log10(thmfs) * mb + bb
print mb, bb

## Plotting setup
ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
#plt.subplots_adjust(left=0.2)
t1, t2, t3, t4 = plt.plot(plotmfs,plotsigmas[0,:],'k.',plotmfs,plotsigmas[1,:],'r.',plotmfs,plotsigmas[2,:],'b.',
                      thmfs,thsigs,'--k')
plt.legend((t1,t2,t3), (r"$\displaystyle t_{\rm 10}$",r"$\displaystyle t_{\rm 50}$",r"$\displaystyle t_{\rm 90}$"),prop={'size':8},loc=4)
plt.text(10**(-0.3+0.05*2.3),0.9*2.5,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([0.5,1.0e2,0,2.5])
plt.yticks([0,1,2])
plt.xticks([1.0,1.0e1,1.0e2])

plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")
plt.ylabel(r"$\displaystyle \sigma^2_{{\rm ln}\Sigma}$")

axx = plt.subplot(1,2,2)
axx.yaxis.set_label_position("right")
t1, t2, t3 = plt.plot(plotmfs,plotxs[0,:],'k.',plotmfs,plotxs[1,:],'r.',plotmfs,plotxs[2,:],'b.')
plt.text(10**(-0.3+0.05*2.3),0.9*2.5,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([0.5,1.0e2,0.,2.5])
plt.yticks([0.,1,2.],[' ',' ',' '])
plt.xticks([1.0,1.0e1,1.0e2])

plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")
plt.ylabel(r"$\displaystyle x$")

pp = PdfPages(plotdir + fname)
pp.savefig()
pp.close()


