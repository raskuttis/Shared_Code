from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
pdffile = 'sdpdfxy.dat'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

mlist = ['5.0e3', '1.0e4', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '20.0', '25.0', '35.0', '15.0', '25.0', '35.0']
nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dflist = np.core.defchararray.add(dflist, ['_N256_Tf4/'] * nfs)

sigmalow = []
sigmamid = []
sigmahigh = []
modelsigma = []
nconv = 25

for i in xrange(0,nfs):
    
    datafolder = dflist[i-1]
    print i, datafolder
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    sigmacloud = out_sigma(outlines)

    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, tff)
    eff = hst_eff(hstdata, tff)
    pdftime = sim_pdftime(pdflines, tff)

    tmini = np.abs(mstar - 0.1 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    pdf = sim_pdf(pdflines,tminpdfi)
    pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,tminpdfi)
    pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
    sigma = sigma / np.log10(np.exp(1))
    sigmalow.append(sigma)
    
    tmini = np.abs(mstar - 0.5 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    pdf = sim_pdf(pdflines,tminpdfi)
    pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,tminpdfi)
    pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
    sigma = sigma / np.log10(np.exp(1))
    sigmamid.append(sigma)
    
    tmini = np.abs(mstar - 0.9 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    pdf = sim_pdf(pdflines,tminpdfi)
    pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,tminpdfi)
    pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
    sigma = sigma / np.log10(np.exp(1))
    sigmahigh.append(sigma)
    
    modelsigma.append(sigmacloud)


alist = ['_a0.1_RT/', '_a0.2_RT/', '_a0.4_RT/', '_a0.8_RT/', '_a1.5_RT/', '/', '_a3.0_RT/', '_a6.0_RT/', '_a10.0_RT/']
nfs = len(alist)
dflist = np.core.defchararray.add(['UV_M5.0e4_R15.0_N256_Tf4'] * nfs, alist)

asigmalow = []
asigmamid = []
asigmahigh = []
modelalpha = []
nconv = 25

for i in xrange(0,nfs):
    
    datafolder = dflist[i-1]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)
    
    tff = out_tff(outlines)
    msol = out_msol(outlines)
    alphavir = out_alphavir(outlines)
    
    time = hst_time(hstdata, tff)
    mstar = hst_mstar(hstdata, tff)
    eff = hst_eff(hstdata, tff)
    pdftime = sim_pdftime(pdflines, tff)
    
    tmini = np.abs(mstar - 0.1 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    pdf = sim_pdf(pdflines,tminpdfi)
    pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
    nvals = np.sum(pdf)
    pdflogx = sim_pdfx(pdflines)
    pdfm = sim_pdfm(pdflines,tminpdfi)
    pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
    amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals)
    [ampold, muold, sigmaold] = [amp, mu, sigma]
    sigma = sigma / np.log10(np.exp(1))
    asigmalow.append(sigma)
    print i, datafolder, sigma
    
    tmini = np.abs(mstar - 0.5 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()
    
    sigma = 0
    aug = 0
    while sigma == 0 and aug <= 10:
        pdf = sim_pdf(pdflines,tminpdfi + aug)
        pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
        nvals = np.sum(pdf)
        pdflogx = sim_pdfx(pdflines)
        pdfm = sim_pdfm(pdflines,tminpdfi)
        pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
        amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
        sigma = sigma / np.log10(np.exp(1))
        aug = aug + 1
    
    asigmamid.append(sigma)
    print i, datafolder, sigma
    
    tmini = np.abs(mstar - 0.9 * eff).argmin()
    tmin = time[tmini]
    tminpdfi = np.abs(pdftime - tmin).argmin()

    sigma = 0
    aug = 0
    while sigma == 0 and aug <= 10:
        pdf = sim_pdf(pdflines,tminpdfi + aug)
        pdf = np.convolve(pdf, np.ones((nconv,))/nconv, mode='same')
        nvals = np.sum(pdf)
        pdflogx = sim_pdfx(pdflines)
        pdfm = sim_pdfm(pdflines,tminpdfi)
        pdfm = np.convolve(pdfm, np.ones((nconv,))/nconv, mode='same')
        amp, mu, sigma, chisq = sim_fitlnpdfvals(pdfm,pdflogx,pdfm,0.15,0.95,nvals,pin=[ampold,muold,sigmaold])
        sigma = sigma / np.log10(np.exp(1))
        aug = aug + 1

    asigmahigh.append(sigma)
    print i, datafolder, sigma
    
    modelalpha.append(alphavir)


xcap = 0.05*3
ycap = 0.9*4
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plow, pmid, phigh = plt.plot(modelsigma,sigmalow,'k.',modelsigma,sigmamid,'r.',modelsigma,sigmahigh,'b.')
plt.legend((plow,pmid,phigh), (r"$\displaystyle 10~\%$",r"$\displaystyle 50~\%$",r"$\displaystyle 90~\%$"))
plt.text(10**(1+0.05*2),ycap,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0e1,1.0e3,0,4])
plt.yticks([0,1,2,3,4])
plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")
plt.xlabel(r"$\displaystyle \Sigma_{\rm cl,0} / M_{\odot} {\rm pc^{-2}}$")

plt.subplot(1,2,2)
plow, pmid, phigh = plt.plot(modelalpha,asigmalow,'k.',modelalpha,asigmamid,'r.',modelalpha,asigmahigh,'b.')
plt.text(10**(-1+0.05*2),ycap,r"$\displaystyle(b)$")
plt.xscale('log')
plt.axis([1.0e-1,1.0e1,0,4])
plt.yticks([0,1,2,3,4],[' ',' ',' ',' '])
plt.xticks([1.0e-1,1.0,1.0e1])

#plt.ylabel(r"$\displaystyle \sigma_{{\rm ln}\Sigma}$")
plt.xlabel(r"$\displaystyle \alpha_{\rm vir,0}$")

pp = PdfPages(plotdir + 'f21.pdf')
pp.savefig()
pp.close()


