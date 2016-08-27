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
nconv = 25

nlines = len(pdflines)
nts = int(np.floor(nlines / 2))
plottimes = []
plotsigmas = []
plotmeans = []
massmeans = []
volmeans = []
plotmassmeans = []
plotchis = []
plotchitimes = []
for i in xrange(100, nts):
    pdforig = sim_pdf(pdflines,i-1)
    pdf = np.convolve(pdforig, np.ones((nconv,))/nconv, mode='same')
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfmorig = sim_pdfm(pdflines,i-1)
    pdfm = np.convolve(pdfmorig, np.ones((nconv,))/nconv, mode='same')
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.15,0.95,nvals)
    ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
    ampo, muo, sigmao, chisqo = sim_fitlnpdfvals(pdforig,pdflogx,pdfmorig,0.15,0.95,nvals)
    ampmo, mumo, sigmamo, chisqmo = sim_fitlnpdfvals(pdfmorig,pdflogx,pdfmorig,0.15,0.95,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    sigmam = sigmam / np.log10(np.exp(1))
    mum = mum / np.log10(np.exp(1))
    mean = np.exp(mu + 0.5 * sigma**2)
    meanm = np.exp(mum + 0.5 * sigmam**2)
    exmean, exmeanm = sim_pdf_mean(pdf,pdflogx)
    
    print i, sigma, sigmam, chisq, chisqm
    
    if sigmam > 0.0:
        plotchis.append((chisqo + chisqo))
        plotchitimes.append(pdftime[i-1])
        if pdftime[i-1] > 0.4:
            plottimes.append(pdftime[i-1])
            plotsigmas.append((sigma + sigma) / 2.0)
            plotmeans.append(mean / msol)
            plotmassmeans.append(meanm / msol)
            volmeans.append(exmean / msol)
            massmeans.append(exmeanm / msol)

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_NF/'
altpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
altpdftime = sim_pdftime(altpdflines, tff)

altnlines = len(altpdflines)
altnts = int(np.floor(altnlines / 2))
altplottimes = []
altplotsigmas = []
altplotmeans = []
altplotchis = []
altplotchitimes = []
for i in xrange(100, altnts):
    pdforig = sim_pdf(altpdflines,i-1)
    pdf = np.convolve(pdforig, np.ones((nconv,))/nconv, mode='same')
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(altpdflines)
    pdfmorig = sim_pdfm(altpdflines,i-1)
    pdfm = np.convolve(pdfmorig, np.ones((nconv,))/nconv, mode='same')
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.15,0.95,nvals)
    ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
    ampo, muo, sigmao, chisqo = sim_fitlnpdfvals(pdforig,pdflogx,pdfmorig,0.15,0.95,nvals)
    ampmo, mumo, sigmamo, chisqmo = sim_fitlnpdfvals(pdfmorig,pdflogx,pdfmorig,0.15,0.95,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    mean = np.exp(mu + 0.5 * sigma**2)
    
    print i, sigma, sigmam, chisq, chisqm
    
    if sigmam > 0.0:
        altplotchis.append((chisqo + chisqo))
        altplotchitimes.append(altpdftime[i-1])
        if altpdftime[i-1] > 0.4:
            altplottimes.append(altpdftime[i-1])
            altplotsigmas.append((sigma + sigma) / 2.0)
            altplotmeans.append(mean / msol)

xcap = 0.05*3
ycap = 0.9*3
plt.figure(figsize = [xsize,3*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(3,1,1)
plt.subplots_adjust(left=0.2)
pfid, pnf = plt.plot(plottimes,plotmeans,'k',altplottimes,altplotmeans,'r')
#plt.plot(plottimes,plotmassmeans,'b--',plottimes,volmeans,'b.',plottimes,massmeans,'b-.')
plt.legend((pfid, pnf), (r"$\displaystyle {\rm Fiducial}$",r"$\displaystyle {\rm No Feedback}$"))
plt.text(xcap,10**(0.9*2+0.48),r"$\displaystyle(a)$")
plt.yscale('log')
plt.axis([0,2,3.0,300.0])
plt.yticks([1.0e1,1.0e2])
plt.xticks([0,1,2],[' ',' ',' '])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \langle \Sigma \rangle_{\rm cloud} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(3,1,2)
pfid, pnf = plt.plot(plottimes,plotsigmas,'k',altplottimes,altplotsigmas,'r')
plt.text(xcap,ycap,r"$\displaystyle(b)$")
plt.axis([0,2,0,3])
plt.xticks([0,1,2],[' ',' ',' '])
plt.yticks([0,1,2,3])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")

plt.subplot(3,1,3)
pfid, pnf = plt.plot(plotchitimes,plotchis,'k',altplotchitimes,altplotchis,'r')
plt.text(xcap,0.9*20,r"$\displaystyle(c)$")
plt.axis([0,2,0,20])
plt.xticks([0,1,2])
plt.yticks([0,10,20])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \chi^2$")

pp = PdfPages(plotdir + 'f19.pdf')
pp.savefig()
pp.close()


