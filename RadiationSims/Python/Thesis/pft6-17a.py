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

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B50000.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_SDs/']
dftlist = [0,0,0]
pdffile = 'angmasspdfmean.dat'
nds = len(dflist)

tlist = [0.02,0.1,0.5,0.8]
ttypelist = [1,1,1,1]
nts = len(tlist)

plotdens = []
plotpdfs = []
plotfits = []

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
        else:
            tmini = np.abs(mof - tlist[j] * (1.0 - eff)).argmin()
        tmin = time[tmini]
        tminpdfi = np.abs(pdftime - tmin).argmin()

        plotpdf = sim_pdf(pdflines,tminpdfi)
        print i, j, tmin, np.sum(plotpdf) / 1.0e9

        plotpdf = plotpdf / 1.0e9
        plotangle = sim_pdfx(pdflines)

        plotangle = np.sin(np.arccos(2.0 * plotangle - 1.0) - np.pi / 2.0)

        plotdens.append(plotangle)
        plotpdfs.append(plotpdf)

xcap = 10**(-2.0+0.05*6)
ycap = 10**(-2.0+0.9*3)
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
pfid, pall, pmean = plt.plot(plotdens[2],plotpdfs[2],'k',plotdens[6],plotpdfs[6],'r',plotdens[10],plotpdfs[10],'b')
plt.legend((pfid, pall, pmean), (r"$\displaystyle {\rm Run~I}$",r"$\displaystyle {\rm Run~MIII}$",r"$\displaystyle {\rm Run~MI}$"),prop={'size':8},loc=2)
#plt.text(xcap,ycap,r"$\displaystyle(a): t_{\rm 2}$")
plt.axis([-1.0,1.0,0.0,0.05])
plt.xticks([-1,0,1])
plt.yticks([0.0, 0.02, 0.04])

plt.xlabel(r"$\displaystyle {\rm sin}\theta$")
plt.ylabel(r"$\displaystyle dM / d\Omega$")

pp = PdfPages(plotdir + 'ft6-17a.pdf')
pp.savefig()
pp.close()


