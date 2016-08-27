from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_hst import *
from hyp_out import *
from hyp_math import *
from hyp_fluxes import *
from hyp_pdf import *
from hyp_star import *
from hyp_models import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import fminbound

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperIII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
pdffile = 'sdpdfxy.dat'
starfile = 'star'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

calclist = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
#calclist = [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

dflist = set_model_lookup('Fluxes/')
dflistsds = set_model_lookup('SDs/')
dflistnf = set_model_lookup('NF_Fluxes/')
nfs = len(dflist)

headoutlist = ['model', 'e', 'eadj', 'ts', 't50', 't90', 'tunb', 'beta', 'tbreak', 'eff']
#headoutlist = ['model', 'sigma', 'alpha', 'funb', 'fedd', 'sigln', 'x', 'fabs', 'prat', 'vout']
headcalclist = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
nheads = len(headoutlist)
nconv = 20

importflag = 0
tablefile = 'restable.dat'
tbf = open(plotdir + tablefile, 'r')
if importflag == 1:
    alllines = tbf.readlines()
tbf.close()

tablefile = 'p1table.dat'
tbf = open(plotdir + tablefile, 'w+')

for i in xrange(0,nfs):
    
    calcfile = calclist[i]
    
    if calcfile == 1:
    
        datafolder = dflist[i]
        datafoldernf = dflistnf[i]
        datafolderpdf = dflistsds[i]
        outlines = read_outfile(hostname,datadir + datafolder + outfile)
        hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
        hstnfdata = read_hstfile(hostname,datadir + datafoldernf + hstfile)
        
        tff = out_tff(outlines)
        tMyr = out_tMyr(outlines)
        tffMyr = tff / tMyr
        msol = out_msol(outlines)
        mcloud = out_mcloud(outlines)
        mcout = mcloud / msol
        rcloud = out_rcloud(outlines)
        sigmacloud = out_sigma(outlines)
        alpha = out_alphavir(outlines)
        vrms = out_vturb(outlines)
        vesc = vrms * np.sqrt(5.0 * 2.0 / (3.0 * alpha))
        nh0 = out_rhocloud(outlines)
        mprint = float2explatex(mcout)
        mname = modelstr(mcout, rcloud)
        clight = out_clight(outlines)
        kappa = out_kappa(outlines)
        psi = out_Psi(outlines)
        gravc = out_G(outlines)
        sigmaadj = sigmacloud * (1.0 - 0.12)
        
        time = hst_time(hstdata, tff)
        timenf = hst_time(hstnfdata, tff)
        
        eff = hst_mstar(hstdata, mcloud)
        effgas = hst_mgas(hstdata, mcloud)
        eps = hst_eff(hstdata, mcloud)
        epsnf = hst_eff(hstnfdata, mcloud)
        epsadj = eps / (1.0 - 0.12)
        ts = hst_tstar(hstdata, mcloud, tff)
        t50 = hst_teff(hstdata, mcloud, tff, 0.5)
        t90 = hst_teff(hstdata, mcloud, tff, 0.9)
        ke = hst_ke(hstdata)
        pturb = np.sqrt(2.0 * effgas * mcloud * ke)
        
        if importflag == 1:
            splitline = alllines[i].rstrip('\\\\ \n')
            splitline = splitline.split(' & ')
        
        strout = ''
        for j in xrange(0, nheads):
            hc = headcalclist[j]
            if hc == 1:
                head = headoutlist[j]
            else:
                head = 'dummy'
            if j < nheads-1:
                fend = ' & '
            else:
                fend = '\\\\ \n'
            if head == 'model':
                strout = strout + mname + fend
            elif head == 'sigma':
                strout = strout + '{:.2f}'.format(sigmacloud) + fend
            elif head == 'r':
                strout = strout + '{:d}'.format(int(rcloud)) + fend
            elif head == 'tff':
                strout = strout + '{:.2f}'.format(tff / tMyr) + fend
            elif head == 'm':
                strout = strout + mprint + fend
            elif head == 'a':
                strout = strout + '{:.1f}'.format(alpha) + fend
            elif head == 'vesc':
                strout = strout + '{:.2f}'.format(vesc) + fend
            elif head == 'vrms':
                strout = strout + '{:.2f}'.format(vrms) + fend
            elif head == 'nh':
                strout = strout + '{:.2f}'.format(nh0) + fend
            elif head == 'e':
                strout = strout + '{:.2f}'.format(eps) + fend
            elif head == 'eadj':
                strout = strout + '{:.2f}'.format(epsadj) + fend
            elif head == 'ts':
                strout = strout + '{:.2f}'.format(ts) + fend
            elif head == 't50':
                strout = strout + '{:.2f}'.format(t50) + fend
            elif head == 't90':
                strout = strout + '{:.2f}'.format(t90) + fend
            elif head == 'alpha':
                fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
                tmini = np.abs(eff - 0.5 * eps).argmin()
                thist = time[tmini]
                alpha = flux_tableouts(fluxdata, thist, 'alpha', rcloud, tff, mcloud, 1.0)
                strout = strout + '{:.2f}'.format(alpha) + fend
            elif head == 'fedd':
                fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
                tmini = np.abs(eff - 0.9 * eps).argmin()
                thist = time[tmini]
                fout = flux_tableouts(fluxdata, thist, 'fedd', rcloud, tff, mcloud, kappa / clight)
                strout = strout + '{:.2f}'.format(fout) + fend
            elif head == 'fabs':
                fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
                tmini = np.abs(time - 3.0 / tffMyr - ts).argmin()
                thist = time[tmini]
                fout = flux_tableouts(fluxdata, thist, 'fabs', rcloud, tff, mcloud, psi, time = time, eff = eff)
                strout = strout + '{:.2f}'.format(fout) + fend
            elif head == 'vout':
                fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
                timef, prejout = fluxintvst(fluxdata, 12, 2.0, rcloud, tff)
                fout = np.max(prejout) / (eps * mcloud)
                strout = strout + '{:.2f}'.format(fout) + fend
            elif head == 'prat':
                fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
                tmini = np.abs(time - 3.0 / tffMyr - ts).argmin()
                thist = time[tmini]
                prin = flux_tableouts(fluxdata, thist, 'pradin', rcloud, tff, mcloud, psi / clight, time = time, eff = eff)
                prej = flux_tableouts(fluxdata, thist, 'prej', rcloud, tff, mcloud, 1.0)
                fout = (prej + pturb[tmini]) / prin
                strout = strout + '{:.2f}'.format(fout) + fend
            elif head == 'funb':
                stesc = np.max(eff) - eff[-1]
                stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)
                funb = star_funboundattt(stardata, time, 4.0, gravc, mcloud)
                fout = (funb + stesc) / eps
                print funb, stesc, eps, fout
                strout = strout + '{:.2f}'.format(fout) + fend
            elif head == 'x' or head == 'sigln':
                pdflines = read_pdffile(hostname,datadir + datafolderpdf + pdffile)
                pdftime = sim_pdftime(pdflines, tff)
                ndiff = 20

                tmini = np.abs(eff - 0.1 * eps).argmin()
                thist = time[tmini]
                tminpdfi = np.abs(pdftime - thist).argmin()
                sigaout, meanout = sim_fitmeansigma(pdflines, tff, msol, thist, ndiff)
                mtot = sim_pdf_mass(pdflines,tminpdfi,(4.0*rcloud)**2,mcloud)
                xaout = np.sqrt(sigmaadj * msol * mtot / (meanout))
                
                tmini = np.abs(eff - 0.5 * eps).argmin()
                thist = time[tmini]
                tminpdfi = np.abs(pdftime - thist).argmin()
                sigbout, meanout = sim_fitmeansigma(pdflines, tff, msol, thist, ndiff)
                mtot = sim_pdf_mass(pdflines,tminpdfi,(4.0*rcloud)**2,mcloud)
                xbout = np.sqrt(sigmaadj * msol * mtot / (meanout))
                
                tmini = np.abs(eff - 0.8 * eps).argmin()
                thist = time[tmini]
                tminpdfi = np.abs(pdftime - thist).argmin()
                sigcout, meanout = sim_fitmeansigma(pdflines, tff, msol, thist, ndiff)
                mtot = sim_pdf_mass(pdflines,tminpdfi,(4.0*rcloud)**2,mcloud)
                xcout = np.sqrt(sigmaadj * msol * mtot / (meanout))
                
                if head == 'x':
                    faout = xbout - xaout
                    fbout = xbout
                    fcout = xcout - xbout
                if head == 'sigln':
                    faout = sigbout - sigaout
                    fbout = sigbout
                    fcout = sigcout - sigbout
                print fbout, faout, fcout
                strout = strout + '$' + '{:.2f}'.format(fbout) + '^{+' + '{:.2f}'.format(fcout) + '}_{-' + '{:.2f}'.format(faout) + '}$' + fend
            else:
                if importflag == 1:
                    strout = strout + splitline[j] + fend
                else:
                    strout = strout + '{:.2f}'.format(0.0) + fend

        print i, datafolder, max(time), max(timenf)

    else:

        if importflag == 1:
            strout = alllines[i]
        else:
            strout = ''
            for j in xrange(0, nheads):
                if j < nheads-1:
                    fend = ' & '
                else:
                    fend = '\\\\ \n'
                strout = strout + '{:.2f}'.format(0.0) + fend

    tbf.write(strout)

tbf.close()

