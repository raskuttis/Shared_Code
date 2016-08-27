from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import fminbound

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 50
nconvall = 50
nconvallh = 50
bzwrong = 1

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.2_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.5_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B1.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B2.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B10.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B20.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B50.0_LD_SDs/']
betas = [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
pdffile = 'sdpdfxy.dat'
nds = len(dflist)
tvs = [0.01, 0.1, 0.3]
ntvs = len(tvs)

plotmfs = np.zeros(nds)
plotmfvars = np.zeros((ntvs,nds))
plotsigmas = np.zeros((ntvs,nds))
ploteffs = np.zeros(nds)

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
    bztote = hst_num(hstdata, 17)
    bsq = 2.0 * (bztote / (4.0 * rcloud**3)) + bz**2
    bmean = np.sqrt(bsq)
    mfall = 2 * np.sqrt(gravc) * mcloud / (bmean * rcloud**2)
    ploteffs[i] = np.max(mstar) / (1.0-0.12)

    pdftime = sim_pdftime(pdflines, tff)
    nts = len(pdftime)
    
    plottime = []
    plotsigma = []
    
    startflag = False
    
    print i, pdffile, mf0
    for j in xrange(15,390):
        
        tmini = np.abs(time - pdftime[j]).argmin()

        plotpdf = sim_pdfm(pdflines,j)
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
            mu = mu / np.log10(np.exp(1))
            mean = np.exp(mu - 0.5 * sigma**2)
            
            plottime.append(pdftime[j])
            mulogten = np.sum(plotlogden * plotpdf)
            sigmalogten = np.sqrt(np.sum(plotlogden * plotlogden * plotpdf) - mulogten**2)
            x = np.sqrt(sigmacloud * msol * (1.0 - 0.12) / (mean))
            #sigma = sigmalogten / np.log10(np.exp(1))
            plotsigma.append(x)
            #print j, sigma
            startflag = True
        # else:
            # print sigmam, h3m, h4m, chisqm

    plotsigma = np.convolve(plotsigma, np.ones((nconvall,))/nconvall, mode='same')

    for j in xrange(0,ntvs):

        tmini = np.abs(mstar - tvs[j]).argmin()
        plotmfvars[j,i] = mfall[tmini]
        tmin = time[tmini]
        tminpdfi = np.abs(plottime - tmin).argmin()
        plotsigmas[j,i] = plotsigma[tminpdfi]
        if j == 0:
            print j, tmin, plotsigma[tminpdfi], np.max(time)



fitmu = np.logspace(-1.0, 3.0, num = 1000.0)
xmu = np.interp(fitmu,plotmfs,np.squeeze(plotsigmas[0,:]))
psi = 2000.0
fiteff = []
fitepseff = []
nfs = len(fitmu)
for i in xrange(0,nfs):
    x = xmu[i] * 0.84/1.12
    sigma0 = sigmacloud * (1.0 - 0.12) / x**2
    sigma = np.sqrt(0.23*np.log10(fitmu[i])+1.6)
    print i, fitmu[i], sigma, sigma0
    def min_func(eps):
        return 1.0 - eps_of(eps, x, sigma, sigma0, psi)
    eps = fminbound(min_func,1.0e-6,1.0)
    fitepseff.append(eps)
    fiteff.append(min_func(eps))

ma, ba = np.polyfit(np.log10(plotmfs), plotsigmas[0,:], 1)
thmfs = np.logspace(-1.0, 2.0)
thsigs = np.log10(thmfs) * ma + ba
print ma, ba

mb, bb = np.polyfit(np.log10(plotmfvars[0,:]), plotsigmas[0,:], 1)
thsigbs = np.log10(thmfs) * mb + bb
print mb, bb

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
t1 = plt.plot(plotmfs,ploteffs,'r.',fitmu,fiteff,'-k',fitmu,fitepseff,':k')
#plt.legend((t1,t2,t3), (r"$\displaystyle {\rm SFE} = 2~\%$",r"$\displaystyle {\rm SFE} = 10~\%$",r"$\displaystyle {\rm SFE} = 30~\%$"),prop={'size':8},loc=1)
plt.text(10**(-0.3+0.05*2.3),0.9*1,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([0.5,1.0e2,0.,1.])
plt.yticks([0.,0.5,1])
plt.xticks([1.0,1.0e1,1.0e2])

plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")
plt.ylabel(r"$\displaystyle \varepsilon$")

pp = PdfPages(plotdir + 'ft6-8c.pdf')
pp.savefig()
pp.close()


