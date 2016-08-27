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
hostname = 'raskutti@bellona.astro.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
pdffile = 'sdpdfxy.dat'
print pdffile
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

tff = out_tff(outlines)
msol = out_msol(outlines)
pdftime = sim_pdftime(pdflines, tff)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
mgas = hst_mgas(hstdata, mcloud)

nlines = len(pdflines)
nts = int(np.floor(nlines / 2))
plottimes = []
altplottimes = []
denxs = []
denmxs = []
sdxs = []
sdmxs = []
for i in xrange(1000, nts-1):
    pdf = sim_pdf(pdflines,i)
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,i)
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9,nvals)
    ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    sigmam = sigmam / np.log10(np.exp(1))
    mum = mum / np.log10(np.exp(1))
    mean = np.exp(mu + 0.5 * sigma**2)
    meanm = np.exp(mum - 0.5 * sigmam**2)
    tmini = np.abs(time - pdftime[i]).argmin()
    mt = mgas[tmini] * mcloud
    meanr = mt / (np.pi * rcloud**2)
    #exmean, exmeanm = sim_pdf_mean(pdf,pdflogx)
    if sigmam > 0.0:
        plottimes.append(pdftime[i])
        sdxs.append(np.sqrt(meanr / (1.0 * mean)))
        sdmxs.append(np.sqrt(meanr / (1.0 * meanm)))

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
pdffile = 'sdpdfmeancirc.dat'
haltpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
haltpdftime = sim_pdftime(haltpdflines, tff)

haltnlines = len(haltpdflines)
haltnts = int(np.floor(haltnlines / 2))
for i in xrange(1000, haltnts-1):
    pdf = sim_pdf(haltpdflines,i)
    pdflogx = sim_pdfx(haltpdflines)
    pdfm = sim_pdfm(haltpdflines,i)
    nvals = np.sum(pdf)
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9,nvals)
    ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    sigmam = sigmam / np.log10(np.exp(1))
    mum = mum / np.log10(np.exp(1))
    mean = np.exp(mu + 0.5 * sigma**2)
    meanm = np.exp(mum - 0.5 * sigmam**2)
    tmini = np.abs(time - pdftime[i]).argmin()
    mt = mgas[tmini] * mcloud
    meanr = mt / (np.pi * rcloud**2)
    #exmean, exmeanm = sim_pdf_mean(pdf,pdflogx)
    
    if sigmam > 0.0:
        plottimes.append(pdftime[i])
        sdxs.append(np.sqrt(meanr / (4.0 * mean)))
        sdmxs.append(np.sqrt(meanr / (4.0 * meanm)))

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
pdffile = 'denpdfall.dat'
altpdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
altpdftime = sim_pdftime(altpdflines, tff)

altnlines = len(altpdflines)
altnts = int(np.floor(altnlines / 2))
for i in xrange(1000, altnts-1):
    pdf = sim_pdf(altpdflines,i)
    pdflogx = sim_pdfx(altpdflines)
    pdfm = sim_pdfm(altpdflines,i)
    nvals = np.sum(pdf)
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9,nvals)
    ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9,nvals)
    sigma = sigma / np.log10(np.exp(1))
    mu = mu / np.log10(np.exp(1))
    sigmam = sigmam / np.log10(np.exp(1))
    mum = mum / np.log10(np.exp(1))
    mean = np.exp(mu + 0.5 * sigma**2)
    meanm = np.exp(mum - 0.5 * sigmam**2)
    tmini = np.abs(time - altpdftime[i]).argmin()
    mt = mgas[tmini] * mcloud
    meanr = 3.0 * mt / (4.0 * np.pi * rcloud**3)
    #exmean, exmeanm = sim_pdf_mean(pdf,pdflogx)
    if sigmam > 0.0:
        altplottimes.append(pdftime[i])
        denxs.append((meanr / mean)**(1.0/3.0))
        denmxs.append((meanr / meanm)**(1.0/3.0))

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
#plt.plot(plottimes,sdxs,'k',altplottimes,denxs,'r')
plt.plot(plottimes,sdmxs,'k',altplottimes,denmxs,'r')
plt.axis([0,3,0,3])
plt.yticks([0,1,2,3])
plt.xticks([0,1,2,3])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle x$")

pp = PdfPages(plotdir + 'fr7.pdf')
pp.savefig()
pp.close()


