from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_math import *
from ..Hyperion.hyp_pdf import *

import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'denpdfall.dat'
starfile = 'star'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
rhocloud = out_rhocloud(outlines)
msol = out_msol(outlines)
time = hst_time(hstdata, tff)
mstar = hst_mstar(hstdata, mcloud)
tstar = time[np.argmax(mstar > 0.0)]
mgas = hst_mgas(hstdata, mcloud)

pfit, bfit = mstar_plfit(time, mstar, 0.01, 0.3)
tcut, afpl, aspl, anfpl, anspl = mstar_brokenplfit(time, mstar, 0.01, 0.3)
tcutlin, anflin, anslin = mstar_brokenlinfit(time, mstar, 0.01, 0.3)
tstar = time[np.argmax(mstar > 0.0)]
mstarfit = broken_pl(time, tstar, 4.0, pfit, 1.0, 10**bfit, 1.0)
mstarfit = 10**(pfit * np.log10(time-tstar) + bfit)
plmstarfit = broken_pl(time, tstar, tcut, afpl, aspl, anfpl, anspl)
linmstarfit = broken_pl(time, tstar, tcutlin, 1.0, 1.0, anflin, anslin)
mcut = mstar[np.argmax(time > tcut)]

# Read in PDF for calculation of mean density
pdftime = sim_pdftime(pdflines, tff)
nlines = len(pdflines)
nts = int(np.floor(nlines / 2))
times = []
besteff = 0.45
besteffbreak = 1.21
mstarkt = []
mstarktbreak = []

told = 0.0
mstarold = 0.0
mbreakold = mcut
for j in xrange(300, nts-1):
        
    pdf = sim_pdf(pdflines,j)
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,j)
    #amp, mu, sigma, chisq = sim_fitlnpdfvals(pdf,pdflogx,pdfm,0.1,0.9,nvals)
    ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.1,0.9,nvals)
    #sigma = sigma / np.log10(np.exp(1))
    #mu = mu / np.log10(np.exp(1))
    sigmam = sigmam / np.log10(np.exp(1))
    mum = mum / np.log10(np.exp(1))
    #mean = np.exp(mu + 0.5 * sigma**2)
    meanm = np.exp(mum - 0.5 * sigmam**2)
    
    if pdftime[j] < tstar:
        mstarnew = mstarold
    else:
        tmini = np.abs(time - pdftime[j]).argmin()
        mgaskt = mgas[tmini]
        mgaskt = 1.0 - mstarold
        delt = pdftime[j] - told
        dmstar = besteff * mgaskt * np.sqrt(meanm / rhocloud) * delt
        mstarnew = mstarold + dmstar
        print j, pdftime[j], mstarnew, dmstar

    if pdftime[j] < tcut:
        mbreaknew = mbreakold
        mstarktbreak.append(0.0)
    else:
        tmini = np.abs(time - pdftime[j]).argmin()
        mgaskt = mgas[tmini]
        mgaskt = 1.0 - mbreakold
        delt = pdftime[j] - told
        dmstar = besteffbreak * mgaskt * np.sqrt(meanm / rhocloud) * delt
        mbreaknew = mbreakold + dmstar
        mstarktbreak.append(mbreaknew)

    mstarkt.append(mstarnew)
    times.append(pdftime[j])
    told = pdftime[j]
    mstarold = mstarnew
    mbreakold = mbreaknew

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-4.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,2,1)
plt.plot(time, mstar, 'k',time,mstarfit,'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2],[' ',' ',' '])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(a)$")

plt.subplot(2,2,2)
plt.plot(time, mstar, 'k', time, plmstarfit, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2],[' ',' ',' '])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5],[' ',' ',' ',' ',' ',' '])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(b)$")

plt.subplot(2,2,3)
plt.plot(time, mstar, 'k', time, linmstarfit, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(c)$")

plt.subplot(2,2,4)
plt.plot(time, mstar, 'k', times, mstarkt, 'k--', times, mstarktbreak, 'r--')
plt.axis([0,2,0,0.5])
plt.xticks([0,1,2])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5],[' ',' ',' ',' ',' ',' '])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*2,0.9*0.5,r"$\displaystyle(d)$")

pp = PdfPages(plotdir + 'f31.pdf')
pp.savefig()
pp.close()


