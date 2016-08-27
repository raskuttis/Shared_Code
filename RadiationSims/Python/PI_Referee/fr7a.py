from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_pdf import *
from hyp_out import *
from hyp_hst import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'sdpdfxy.dat'
hostname = 'raskutti@bellona.astro.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

tff = out_tff(outlines)
msol = out_msol(outlines)
pdftime = sim_pdftime(pdflines, tff)
mcloud = out_mcloud(outlines)
time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
mgas = hst_mgas(hstdata, mcloud)

nlines = len(pdflines)
nts = int(np.floor(nlines / 2))
plottimes = []
plotsigmas = []
plotmeans = []
massmeans = []
volmeans = []
plotmassmeans = []
plotchis = []
for i in xrange(1000, nts):
    pdf = sim_pdf(pdflines,i-1)
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,i-1)
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9,nvals)
    ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    sigmam = sigmam / np.log10(np.exp(1))
    mum = mum / np.log10(np.exp(1))
    mean = np.exp(mu + 0.5 * sigma**2)
    meanm = np.exp(mum + 0.5 * sigmam**2)
    exmean, exmeanm = sim_pdf_mean(pdf,pdflogx)
    
    print i, sigmam, sigma
    
    if sigmam > 0.0:
        plottimes.append(pdftime[i-1])
        plotsigmas.append(sigma)
        plotmeans.append(mean / msol)
        plotmassmeans.append(meanm / msol)
        plotchis.append(chisq)
        volmeans.append(exmean / msol)
        massmeans.append(exmeanm / msol)

xcap = 0.05*3
ycap = 0.9*3
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
pfid, pnf = plt.plot(plottimes,plotmeans,'k',plottimes,plotmassmeans,'r')
plt.plot(plottimes,volmeans,'k--',plottimes,massmeans,'r--')
plt.legend((pfid, pnf), (r"$\displaystyle \langle \Sigma \rangle_{\rm A}$",r"$\displaystyle \langle \Sigma \rangle_{\rm M}$"),prop={'size':8})
#plt.text(xcap,10**(0.9*4),r"$\displaystyle(a)$")
plt.yscale('log')
plt.axis([0,2,1.0,10**4])
plt.yticks([1.0,1.0e1,1.0e2,1.0e3,1.0e4])
plt.xticks([0,1,2])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \langle \Sigma \rangle / M_{\odot} {\rm pc^{-2}}$")

pp = PdfPages(plotdir + 'fr7a.pdf')
pp.savefig()
pp.close()


