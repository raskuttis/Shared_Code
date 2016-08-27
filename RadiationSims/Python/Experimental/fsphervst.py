from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
from hyp_hst import *
from hyp_fluxes import *
from hyp_math import *
from hyp_star import *
from hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hostname = 'raskutti@bellona.astro.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'
starfile = 'star'
nconv = 1
nconvpdf = 25

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_Fluxes/'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
clight = out_clight(outlines)
tff = out_tff(outlines)
psi = out_Psi(outlines)
kappa = out_kappa(outlines)
gravc = out_G(outlines)
rhocloud = out_rhocloud(outlines)
sigmacloud = out_sigma(outlines)
sigmaadj = sigmacloud * (1.0 - 0.12)
msol = out_msol(outlines)

time = hst_time(hstdata, tff)
eff = hst_mstar(hstdata, mcloud)
mgas = hst_mgas(hstdata, mcloud)
stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)

datafolder = 'UV_M5.0e4_R15.0_N128_Tf4_SDs_Alt/'
pdffile = 'sdmasspdfallcirc.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
nts = len(pdftime)

plottime = []
plotx = []
startflag = False
for j in xrange(100,4500):
    
    epst = eff[j]
    plotpdf = sim_pdf(pdflines,j)
    plotpdf = np.convolve(plotpdf, np.ones((nconvpdf,))/nconvpdf, mode='same')
    mtot = np.sum(plotpdf) / 1.0e9
    plotpdf = plotpdf / np.sum(plotpdf)
    nvals = np.sum(plotpdf)
        
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    if startflag:
        amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
    else:
        amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)

    if (sigma > 0.0):
            
        [ampold, sigmaold, muold] = [amp, sigma, mu]
        sigma = sigma / np.log10(np.exp(1))
        mu = mu / np.log10(np.exp(1))
        mean = np.exp(mu - 0.5 * sigma**2)
        x = np.sqrt(sigmacloud * mtot * msol / (mean * 2.0))
        startflag = True
        
        plottime.append(pdftime[j])
        plotx.append(x)

        print j, x, mtot

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
pdffile = 'sdpdfxy.dat'
pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
pdftime = sim_pdftime(pdflines, tff)
nts = len(pdftime)

aplottime = []
aplotx = []
startflag = False
for j in xrange(100,7500):
    
    epst = eff[j]
    plotpdf = sim_pdf(pdflines,j)
    plotpdf = np.convolve(plotpdf, np.ones((nconvpdf,))/nconvpdf, mode='same')
    plotpdfm = sim_pdfm(pdflines,j)
    plotpdfm = np.convolve(plotpdfm, np.ones((nconvpdf,))/nconvpdf, mode='same')
    plotpdf = plotpdf / np.sum(plotpdf)
    nvals = np.sum(plotpdf)
    
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / msol
    if startflag:
        amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdfm,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
        ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(plotpdfm,plotlogden,plotpdfm,0.15,0.95,nvals,pin=[ampmold,mumold,sigmamold])
    else:
        amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdfm,0.15,0.95,nvals)
        ampm, mum, sigmam, chisqm = sim_fitlnpdfvals(plotpdfm,plotlogden,plotpdfm,0.15,0.95,nvals)

    if (sigma > 0.0):
        
        [ampold, sigmaold, muold] = [amp, sigma, mu]
        [ampmold, sigmamold, mumold] = [ampm, sigmam, mum]
        sigma = sigma / np.log10(np.exp(1))
        mu = mu / np.log10(np.exp(1))
        sigmam = sigmam / np.log10(np.exp(1))
        mum = mum / np.log10(np.exp(1))
        mean = np.exp(mu + 0.5 * sigma**2)
        meanm = np.exp(mum - 0.5 * sigmam**2)
        x = np.sqrt(sigmaadj * msol * (1.0 - epst) / (mean * 1.0))
        xm = np.sqrt(sigmaadj * msol * (1.0 - epst) / (meanm * 1.0))
        startflag = True
        
        aplottime.append(pdftime[j])
        aplotx.append(x)

        print j, x, xm

ftime, mr = fluxintrvst(fluxdata, 5, 2.0, rcloud, tff)
ftime, rmr = fluxintrposvst(fluxdata, 5, 2.0, rcloud, tff)
rmean = rmr / (mr * rcloud)

rcom = star_mu(stardata) / rcloud
dcom = star_sigma(stardata) / rcloud
maxcom = star_max(stardata) / rcloud

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
rc,dc,mc = plt.plot(time,rcom,'k',time,dcom,'r',time,maxcom,'b')
plt.legend((rc, dc, mc), (r"$\displaystyle \mu_*$",r"$\displaystyle \sigma_*$",r"$\displaystyle r_{max,*}$"),prop={'size':8},loc=2)
plt.axis([0,3,0,2])
plt.xticks([0,1,2,3])
plt.yticks([0,1,2])
plt.text(0.9*3,0.9*2,r"$\displaystyle(a)$")
plt.ylabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")

plt.subplot(1,2,2)
x1, x2 = plt.plot(ftime,rmean,'k',plottime,plotx,'r')
plt.legend((x1, x2), (r"$\displaystyle \langle \rho r \rangle / \rho_{\rm cloud}$",r"$\displaystyle x$"),prop={'size':8},loc=4)
plt.axis([0,3,0,2])
plt.xticks([0,1,2,3])
plt.yticks([0,1,2],[' ',' ',' '])
plt.text(0.05*3,0.9*2,r"$\displaystyle(b)$")
#plt.ylabel(r"$\displaystyle r / r_{\rm cloud}$")
plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")

pp = PdfPages(plotdir + 'fsphervst.pdf')
pp.savefig()
pp.close()




