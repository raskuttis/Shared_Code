from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 25

pdflist = ['sdmasspdfmeancirc.dat', 'sdmpdfxyplane_7.dat', 'sdmpdfyzplane_7.dat']
masslist = [1,1,1]
nps = len(pdflist)

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B50000.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_SDs/']
dftlist = [0,0,0]
nds = len(dflist)

tlist = [0.2,0.5]
ttypelist = [1,1]
nts = len(tlist)

plotdens = []
plotpdfs = []
plotfits = []

for kk in xrange(0,nps):
    
    pdffile = pdflist[kk]

    for i in xrange(0,nds):

        datafolder = dflist[i]
        dftype = dftlist[i]
        outlines = read_outfile(hostname,datadir + datafolder + outfile)
        hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)
        pdflines = read_pdffile(hostname,datadir + datafolder + pdffile)

        tff = out_tff(outlines)
        mcloud = out_mcloud(outlines)
        rcloud = out_rcloud(outlines)
        msol = out_msol(outlines)

        time = hst_time(hstdata, tff)
        if dftype == 1:
            mstar = hst_mstar(hstdata, mcloud)
            mgas = hst_mgas(hstdata, mcloud)
            eff = hst_eff(hstdata, mcloud)
        else:
            mstar = hst_num(hstdata, 13) / mcloud
            mgas = hst_num(hstdata, 2) / mcloud
            eff = np.max(mstar)
        mof = mgas[0] - mgas - mstar
        pdftime = sim_pdftime(pdflines, tff)
        
        for j in xrange(0,nts):

            if ttypelist[j] == 1:
                tmini = np.abs(mstar - tlist[j] * eff).argmin()
                #tmini = np.abs(time - tlist[j]).argmin()
            else:
                tmini = np.abs(mof - tlist[j] * (1.0 - eff)).argmin()
            tmin = time[tmini]
            tminpdfi = np.abs(pdftime - tmin).argmin()

            if (masslist[kk] == 1):
                plotpdf = sim_pdf(pdflines,tminpdfi)
                plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
            else:
                plotpdf = sim_pdfm(pdflines,tminpdfi)
                plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
            print i, j, tmin, tminpdfi, pdffile, np.sum(plotpdf) / 1.0e9

            plotpdf = plotpdf / 1.0e9
            #plotpdf = plotpdf / np.sum(plotpdf)
            nvals = np.sum(plotpdf)

            plotlogden = sim_pdfx(pdflines)
            fnorm = plotlogden[1] - plotlogden[0]
            plotden = 10**plotlogden / msol

            plotdens.append(plotden)
            plotpdfs.append(plotpdf / fnorm)


xcap = 10**(-2.0+0.05*6)
ycap = 10**(-2.0+0.9*3)
plt.figure(figsize = [2*xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)


plt.subplot(2,2,1)
pfid, pall, pmean = plt.plot(plotdens[0],plotpdfs[0],'k',plotdens[6],plotpdfs[6],'r',plotdens[12],plotpdfs[12],'b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(a): {\rm Run~I}, t_{\rm 10}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1])

#plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(2,2,2)
pfid, pall, pmean = plt.plot(plotdens[1],plotpdfs[1],'k',plotdens[7],plotpdfs[7],'r',plotdens[13],plotpdfs[13],'b')
plt.legend((pfid, pall, pmean), (r"$\displaystyle \Sigma_*^c$",r"$\displaystyle \Sigma_{*,\perp}^{c}$",r"$\displaystyle \Sigma_{*,\parallel}^{c}$"),prop={'size':8},loc=1)
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(b): {\rm Run~I}, t_{\rm 50}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1],[' ',' ',' ',' '])

#plt.subplot(3,2,3)
#pfid, pall, pmean = plt.plot(plotdens[2],plotpdfs[2],'k',plotdens[8],plotpdfs[8],'r',plotdens[14],plotpdfs[14],'b')
#plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(c): {\rm Run~MIII}, t_{\rm 10}$")
#plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
#plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
#plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1])

#plt.xlabel(r"$\displaystyle \Sigma / M_{\odot} {\rm pc^{-2}}$")
#plt.ylabel(r"$\displaystyle P_M$")

#plt.subplot(3,2,4)
#pfid, pall, pmean = plt.plot(plotdens[3],plotpdfs[3],'k',plotdens[9],plotpdfs[9],'r',plotdens[15],plotpdfs[15],'b')
#plt.xscale('log')
#plt.yscale('log')
#plt.text(xcap,ycap,r"$\displaystyle(d): {\rm Run~MIII}, t_{\rm 50}$")
#plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
#plt.xticks([1.0e-2,1.0,1.0e2,1.0e4],[' ',' ',' ',' '])
#plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1],[' ',' ',' ',' '])

plt.subplot(2,2,3)
pfid, pall, pmean = plt.plot(plotdens[4],plotpdfs[4],'k',plotdens[10],plotpdfs[10],'r',plotdens[16],plotpdfs[16],'b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(e): {\rm Run~MI}, t_{\rm 10}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1])

plt.xlabel(r"$\displaystyle \Sigma / [M_{\odot} {\rm pc^{-2}}]$")
plt.ylabel(r"$\displaystyle P_M$")

plt.subplot(2,2,4)
pfid, pall, pmean = plt.plot(plotdens[5],plotpdfs[5],'k',plotdens[11],plotpdfs[11],'r',plotdens[17],plotpdfs[17],'b')
plt.xscale('log')
plt.yscale('log')
plt.text(xcap,ycap,r"$\displaystyle(f): {\rm Run~MI}, t_{\rm 50}$")
plt.axis([1.0e-2,1.0e4,1.0e-2,1.0e1])
plt.xticks([1.0e-2,1.0,1.0e2,1.0e4])
plt.yticks([1.0e-2,1.0e-1,1.0,1.0e1],[' ',' ',' ',' '])

plt.xlabel(r"$\displaystyle \Sigma / [M_{\odot} {\rm pc^{-2}}]$")
# plt.ylabel(r"$\displaystyle P_M$")

pp = PdfPages(plotdir + 'ft6-18a.pdf')
pp.savefig()
pp.close()


