from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_pdf import *
from ..Hyperion.hyp_star import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_fluxes import *
import numpy as np
import scipy.integrate as spint
import scipy.interpolate as spinterp
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r0.dat'
pdffile = 'sdmasspdfallcirc.dat'
xpdffile = 'sdpdfxy.dat'
starfile = 'star'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

mlist = ['5.0e3', '1.0e4', '1.0e4', '2.0e4', '2.0e4', '2.0e4', '2.0e4', '5.0e4', '5.0e4', '5.0e4', '5.0e4', '1.0e5', '1.0e5', '1.0e5', '1.0e5', '2.0e5', '2.0e5', '2.0e5']
rlist = ['5.0', '8.0', '10.0', '5.0', '8.0', '10.0', '15.0', '8.0', '10.0', '15.0', '20.0', '15.0', '20.0', '25.0', '35.0', '15.0', '25.0', '35.0']
mlist = ['2.0e4', '5.0e4', '2.0e4', '5.0e4', '1.0e5', '2.0e4', '1.0e4', '5.0e4', '1.0e4', '1.0e5', '2.0e5', '5.0e3', '2.0e4', '5.0e4', '1.0e5', '2.0e4', '2.0e5', '1.0e4', '5.0e5', '1.0e5', '5.0e4', '2.0e5', '5.0e4', '2.0e4', '2.0e5']
rlist = ['25.0', '35.0', '20.0', '25.0', '35.0', '15.0', '10.0', '20.0', '8.0', '25.0', '35.0', '5.0', '10.0', '15.0', '20.0', '8.0', '25.0', '5.0', '35.0', '15.0', '10.0', '20.0', '8.0', '5.0', '15.0']
mlist = ['2.0e4', '5.0e4', '1.0e5', '2.0e4', '1.0e4', '5.0e4', '1.0e4', '1.0e5', '2.0e5', '5.0e3', '2.0e4', '5.0e4', '1.0e5', '2.0e4', '2.0e5', '5.0e5', '1.0e5', '5.0e4', '2.0e5', '5.0e4', '2.0e4', '2.0e5']
rlist = ['20.0', '25.0', '35.0', '15.0', '10.0', '20.0', '8.0', '25.0', '35.0', '5.0', '10.0', '15.0', '20.0', '8.0', '25.0', '35.0', '15.0', '10.0', '20.0', '8.0', '5.0', '15.0']

nfs = len(mlist)
dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
dflist = np.core.defchararray.add(dflist, rlist)
dfflist = np.core.defchararray.add(dflist, ['_N128_Tf4_Fluxes/'] * nfs)
dfpdflist = np.core.defchararray.add(dflist, ['_N128_Tf4_SDs/'] * nfs)
dfnflist = np.core.defchararray.add(dflist, ['_N128_Tf4_NF_Fluxes/'] * nfs)

modelsigma = []
modeleff = []
modeleffof = []
tplot = [0.1,0.5,0.9,0.5,3,8]
ttypelist = [1,1,1,0,2,2]
ntplot = len(tplot)
nconv = 20

allfabs = np.zeros((ntplot,nfs))
allfabscum = np.zeros((ntplot,nfs))
allalpha = np.zeros((ntplot,nfs))
allsigma = np.zeros(nfs)
allsigmas = np.zeros((ntplot,nfs))
allxs = np.zeros((ntplot,nfs))
allfedd = np.zeros((ntplot,nfs))
allpr = np.zeros((ntplot,nfs))
allprfin = np.zeros(nfs)
funbs = np.zeros(nfs)
calcpdf = 0

for i in xrange(0,nfs):
    
    datafolder = dfflist[i]
    datafoldernf = dfnflist[i]
    datafolderpdf = dfpdflist[i]
    outlines = read_outfile(hostname,datadir + datafolder + outfile)
    hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
    fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
    if calcpdf:
        pdflines = read_pdffile(hostname,datadir + datafolderpdf + pdffile)
        xpdflines = read_pdffile(hostname,datadir + datafolderpdf + xpdffile)
    #fluxnfdata = read_fluxfile(hostname,datadir + datafoldernf + fluxfile)
    
    tff = out_tff(outlines)
    tMyr = out_tMyr(outlines)
    tffMyr = tff / tMyr
    msol = out_msol(outlines)
    mcloud = out_mcloud(outlines)
    sigmacloud = out_sigma(outlines)
    psi = out_Psi(outlines)
    lsol = psi * msol
    clight = out_clight(outlines)
    gravc = out_G(outlines)
    kappa = out_kappa(outlines)
    rcloud = out_rcloud(outlines)
    vturb = out_vturb(outlines)
    sigmaadj = sigmacloud * (1.0 - 0.12)
    
    time = hst_time(hstdata, tff)
    eff = hst_mstar(hstdata, mcloud)
    effgas = hst_mgas(hstdata, mcloud)
    effof = effgas[0] - effgas - eff
    eps = hst_eff(hstdata, mcloud)
    tstar = hst_tstar(hstdata, mcloud, tff)
    ke = hst_ke(hstdata)
    pturb = np.sqrt(2.0 * effgas * mcloud * ke)
    print i, nfs, datafolder, max(time), sigmacloud

    stesc = np.max(eff) - eff[-1]
    stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
    funb = star_funboundattt(stardata, time, 4.0, gravc, mcloud)
    funbs[i] = (funb + stesc) / eps
    print i, funbs[i]

    if calcpdf:
        pdftime = sim_pdftime(pdflines, tff)
        xpdftime = sim_pdftime(xpdflines, tff)
    
    timef, frcl, routs = fluxvst(fluxdata, 6, 2.0, rcloud, tff)
    fstar = spinterp.interp1d(time, eff * mcloud)
    lstarf = psi * fstar(timef)
    lstarfint = spint.cumtrapz(lstarf, timef * tff, initial=0.0)
    fabsfcl = 1.0 - 4.0 * np.pi * (routs * rcloud)**2 * frcl / lstarf
    timef, prejout = fluxintvst(fluxdata, 12, 2.0, rcloud, tff)
    pradin = lstarfint / clight
    
    allsigma[i] = sigmacloud
    allprfin[i] = np.max(prejout) / (eps * mcloud)

    tdel = rcloud / vturb / tff
    for j in xrange(0,ntplot):
        
        if ttypelist[j] == 2:
            tmini = np.abs(time - tplot[j] / tffMyr - tstar).argmin()
            tminidel = np.abs(time - tplot[j] / tffMyr - tstar - tdel).argmin()
        elif ttypelist[j] == 1:
            tmini = np.abs(eff - tplot[j] * eps).argmin()
            tminidel = np.abs(eff - tplot[j] * eps - tdel).argmin()
        else:
            tmini = np.abs(effof - tplot[j] * (1.0 - eps)).argmin()
            tminidel = np.abs(effof - tplot[j] * (1.0 - eps) - tdel).argmin()
        thist = time[tmini]
        tminfi = np.abs(timef - thist).argmin()
        if calcpdf:
            tminpdfi = np.abs(pdftime - thist).argmin()
            tminxpdfi = np.abs(xpdftime - thist).argmin()
        
        rad, mr = fluxintvsr(fluxdata, 5, thist, rcloud, tff)
        goodinds = np.logical_and((rad > 0.5),(rad < 1.0))
        pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
        allalpha[j,i] = 3.0 - pcoeffs[0]

        rad, rhofr = fluxvsr(fluxdata, 7, thist, rcloud, tff)
        rad, rhodphidr = fluxvsr(fluxdata, 15, thist, rcloud, tff)
        fedd = rhofr * kappa / (rhodphidr * clight)
        fedd = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
        goodind = np.argmin(np.fabs(rad - 1.0))
        allfedd[j,i] = fedd[goodind]

        allfabs[j,i] = fabsfcl[tminfi]

        allpr[j,i] = (prejout[tminfi] + pturb[tmini]) / pradin[tminfi]

        if calcpdf:
            plotlogden = sim_pdfx(pdflines)
            plotden = 10**plotlogden / msol
            plotpdf = sim_pdf(pdflines,tminpdfi)
            plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
            plotpdf = plotpdf / np.sum(plotpdf)
            nvals = np.sum(plotpdf)
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)
            sigma = sigma / np.log10(np.exp(1))
            allsigmas[j,i] = sigma
            
            plotlogden = sim_pdfx(xpdflines)
            plotden = 10**plotlogden / msol
            plotpdf = sim_pdfm(xpdflines,tminxpdfi)
            plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
            plotpdf = plotpdf / np.sum(plotpdf)
            nvals = np.sum(plotpdf)
            amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)
            sigma = sigma / np.log10(np.exp(1))
            mu = mu / np.log10(np.exp(1))
            mean = np.exp(mu - 0.5 * sigma**2)
            x = np.sqrt(sigmaadj * msol * (1.0 - eff[tmini]) / (mean))
            allxs[j,i] = x


for i in xrange(0,nfs):
    print allsigma[i], allalpha[0,i], funbs[i], allfedd[2,i], allfabs[4,i], allpr[4,i], allprfin[i]

