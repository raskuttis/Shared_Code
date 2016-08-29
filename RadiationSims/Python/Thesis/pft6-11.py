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
nconv = 50
nconvall = 50
nconvallh = 50

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_NF_B50.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B20.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B10.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B2.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B1.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.5_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.2_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.1_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.05_LD_SDs/']
betas = [50.0, 20.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]
pdffile = 'sdpdfallcirc.dat'
nds = len(dflist)
tvs = [0.02, 0.1, 0.3]
ntvs = len(tvs)

plotmfs = np.zeros(nds)
plotmfvars = np.zeros((ntvs,nds))
plotsigmas = np.zeros((ntvs,nds))
nconv = 300

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
    kappa = 7.228571e-03
    beta = betas[i]
    bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / beta)
    mf0 = 2 * np.sqrt(gravc) * mcloud / (bz * rcloud**2)
    plotmfs[i] = mf0

    time = hst_time(hstdata, tff)
    mstar = hst_num(hstdata, 13) / mcloud
    bztote = hst_num(hstdata, 17)
    bsq = 2.0 * (bztote / (4.0 * rcloud**3)) + bz**2
    bmean = np.sqrt(bsq)
    mfall = 2 * np.sqrt(gravc) * mcloud / (bmean * rcloud**2)

    pdftime = sim_pdftime(pdflines, tff)
    nts = len(pdftime)
    
    msconv = np.convolve(mstar, np.ones((nconv,))/nconv, mode='same')
    sfr = np.diff(msconv) / np.diff(time)
    plottime = time[1:]
    plotsfr = np.convolve(sfr, np.ones((nconv,))/nconv, mode='same')
    
    startflag = False
    
    print i, pdffile, mf0
    for j in xrange(125,410):
        
        tmini = np.abs(time - pdftime[j]).argmin()

        plotpdf = sim_pdf(pdflines,j)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        mtot = np.sum(plotpdf) / 1.0e9
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)

        plotlogden = sim_pdfx(pdflines)
        plotden = 10**plotlogden
        plottau = plotden * kappa

    for j in xrange(0,ntvs):

        tmini = np.abs(mstar - tvs[j]).argmin()
        plotmfvars[j,i] = mfall[tmini]
        tmin = time[tmini]
        plotsigmas[j,i] = plotsfr[tmini-1]
        print j, tmin, plotsfr[tmini-1]

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
t1,t2,t3 = plt.plot(plotmfs,plotsigmas[0,:],'k.',plotmfs,plotsigmas[1,:],'r.',plotmfs,plotsigmas[2,:],'b.')
plt.legend((t1,t2,t3), (r"$\displaystyle {\rm SFE} = 2~\%$",r"$\displaystyle {\rm SFE} = 10~\%$",r"$\displaystyle {\rm SFE} = 30~\%$"),prop={'size':8},loc=1)
#plt.text(10**(0+0.05*2),0.9*1+0.5,r"$\displaystyle(a)$")
plt.xscale('log')
plt.axis([1.0,1.0e2,0,1])
plt.yticks([0,0.5,1])
plt.xticks([1.0,1.0e1,1.0e2])

plt.xlabel(r"$\displaystyle \mu_{\Phi,0}$")
plt.ylabel(r"$\displaystyle \varepsilon_{\rm ff}$")

pp = PdfPages(plotdir + 'ft6-11.pdf')
pp.savefig()
pp.close()


