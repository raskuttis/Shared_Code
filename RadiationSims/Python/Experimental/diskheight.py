from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_fil import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_star import *
from scipy.interpolate import RegularGridInterpolator as rgi
from astropy.io import fits
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
locdatadir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/'
locdatafolder = 'RB5.0/'

datadir = '/tigress-hsm/raskutti/tigress-hsm/Hyperion/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'
datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_R_B5.0_Out/'
starfile = 'star'

outlines = read_outfile(hostname,datadir + datafolder + outfile)
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

tff = out_tff(outlines)
mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
msol = out_msol(outlines)
dtout = out_dtout(outlines) / tff
nouts = np.floor(4.0 / dtout).astype(int)
tmax = dtout * nouts
dtouts = np.linspace(0.0,tmax,num=nouts+1)
dmouts = np.zeros(nouts)

time = hst_time(hstdata, tff)
mstar = hst_num(hstdata, 13) / mcloud
eff = np.max(mstar)

for i in xrange(0,nouts):

    ti = np.abs(time - dtouts[i]).argmin()
    dmouts[i] = mstar[ti] / eff

tmini = np.abs(dmouts - 0.1).argmin()
tmaxi = np.abs(dmouts - 0.5).argmin()
delti = tmaxi - tmini

rad = 15.
msol = 29.06
nres = 256
nconv = 10

for kk in xrange(0,delti):

    hdulist = fits.open(locdatadir + locdatafolder + projfile)
    hdu = hdulist[0]
    sd = hdu.data / msol

    conv = (4.0 * rad) / nres

    x = np.linspace(0, nres-1, nres) + 0.5
    y = np.linspace(0, nres-1, nres) + 0.5

    fn = rgi((x,y), sd)

rfil = rfil * conv
#perpr = perpr * conv

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axa = plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
cols = ['k', 'r', 'b', 'g']
for kk in xrange(0,delti):
    allrs = allrfils[kk]
    allfils = allnfils[kk]
    for i in xrange(0, nplotfils[kk]):
        plt.plot(allrs[i], allfils[i], cols[kk])
plt.axis([-15.0,15.0,0.0,300.0])
#plt.xticks([0,1,2,3])
plt.yticks([0,100,200,300])
#plt.yscale('log')
plt.xlabel(r"$\displaystyle r$")
plt.ylabel(r"$\displaystyle \Sigma$")

axa = plt.subplot(1,2,2)
plt.subplots_adjust(bottom=0.2)
cols = ['k', 'r', 'b', 'g']
for kk in xrange(0,delti):
    allrs = allrpfils[kk]
    allfils = allnpfils[kk]
    for i in xrange(0, nplotfils[kk]):
        plt.plot(allrs[i], allfils[i], cols[kk])
plt.axis([0.0,15.0,0.0,300.0])
#plt.xticks([0,1,2,3])
plt.yticks([0,100,200,300])
#plt.yscale('log')
plt.xlabel(r"$\displaystyle r$")
plt.ylabel(r"$\displaystyle \Sigma$")

pp = PdfPages(plotdir + 'ft-filprojmiii.pdf')
pp.savefig()
pp.close()


