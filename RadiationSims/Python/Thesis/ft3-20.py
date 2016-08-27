from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_pdf import *
from ..Hyperion.hyp_models import *
import matplotlib.pyplot as plt

## Plot showing evolution of fits to the density PDF

## Define the locations where data is located and where to plot from hyp_models
hostname, datadir, hstfile, outfile, plotdir = init_dirs('ast')
pdffile = 'denpdfall.dat'
fname = 'ft3-20.pdf'

## Fiducial models with no feedback and no gravity
dflist = ['UV_M5.0e4_R15.0_N256_Tf4_NF/', 'UV_M5.0e4_R15.0_N256_Tf4_NG/']
nds = len(dflist)

plottimes = []
plotsigmas = []
plotmeans = []
nconv = 25

for i in xrange(0,nds):
    
    datafolder = dflist[i]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    msol = out_msol(outlines)
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, mcloud)
    pdftime = sim_pdftime(pdflines, tff)
    nlines = len(pdflines)
    nts = int(np.floor(nlines / 2))
    times = []
    sigmas = []
    means = []
    
    ## Iterate over times in the density pdf and find best fits
    for j in xrange(300, nts-1):
        
        pdf = sim_pdf(pdflines,j)
        pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
        nvals = np.sum(pdf)
        pdflogx = sim_pdfx(pdflines)
        pdfm = sim_pdfm(pdflines,j)
        pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
        #amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9)
        ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
        #sigma = sigma / np.log10(np.exp(1))
        #mu = mu / np.log10(np.exp(1))
        sigmam = sigmam / np.log10(np.exp(1))
        mum = mum / np.log10(np.exp(1))
        #mean = np.exp(mu + 0.5 * sigma**2)
        meanm = np.exp(mum - 0.5 * sigmam**2)
        
        if sigmam > 0.0:
            times.append(pdftime[j])
            sigmas.append(sigmam)
            means.append(meanm)
            print j, sigmam, meanm

    plottimes.append(times)
    plotsigmas.append(sigmas)
    plotmeans.append(means)

## Plotting setup
ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
pfid, pgrav, = plt.plot(plottimes[0],plotsigmas[0],'k',plottimes[1],plotsigmas[1],'r')
plt.legend((pfid, pgrav), (r"Fiducial No Feedback",r"No Gravity"),prop={'size':8},loc=4)
plt.text(0.05*2,0.9*3,r"$\displaystyle(a)$")
plt.axis([0,2,0,3])
plt.yticks([0,1,2,3])
#plt.xticks([0,1,2,3])
plt.xticks([0,1,2])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\rho}$")

axy = plt.subplot(1,2,2)
axy.yaxis.set_label_position("right")
#plt.subplots_adjust(left=0.2)
plt.plot(plottimes[0],plotmeans[0],'k',plottimes[1],plotmeans[1],'r',)
plt.text(0.05*2,10**(0.9*2+0.48),r"$\displaystyle(b)$")
plt.axis([0,2,3.0,300.0])
plt.yscale('log')
plt.yticks([1.0e1,1.0e2])
plt.xticks([0,1,2])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \langle \rho \rangle_{\rm cloud}$")

pp = PdfPages(plotdir + fname)
pp.savefig()
pp.close()


