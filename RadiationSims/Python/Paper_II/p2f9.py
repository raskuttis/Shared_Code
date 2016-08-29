from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_fluxes import *
from ..Hyperion.hyp_math import *
from ..Hyperion.hyp_star import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hostname = 'raskutti@tiger.princeton.edu'
outfile = 'RadParGrav.out'
fluxfile = 'fluxes_r1.dat'
hstfile = 'id0/RadParGrav.hst'
starfile = 'star'
nconv = 2

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Fluxes/'
print 'Reading Out'
outlines = read_outfile(hostname,datadir + datafolder + outfile)
print 'Reading Flux'
fluxdata = read_fluxfile(hostname,datadir + datafolder + fluxfile)
print 'Reading Hst'
hstdata = read_hstfile(hostname,datadir + datafolder + hstfile)

mcloud = out_mcloud(outlines)
rcloud = out_rcloud(outlines)
clight = out_clight(outlines)
tff = out_tff(outlines)
psi = out_Psi(outlines)
kappa = out_kappa(outlines)
gravc = out_G(outlines)
rhocloud = out_rhocloud(outlines)

time = hst_time(hstdata, tff)
eff = hst_mstar(hstdata, mcloud)
mgas = hst_mgas(hstdata, mcloud)

tplots = [0.6, 1.06, 1.55, 2.07]
nts = len(tplots)

fedds = []
feddcums = []
feddspecs = []
rads = []
radcums = []

for i in xrange(0,nts):

    tplot = tplots[i]

    tmini = np.abs(time - tplot).argmin()

    rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
    rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
    rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
    #rad, rhofrint = fluxintoutvsr(fluxdata, 7, tplot, rcloud, tff)
    rad, rhofrint = fluxintvsr(fluxdata, 7, tplot, rcloud, tff)
    rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
    #rad, rhodphidrint = fluxintoutvsr(fluxdata, 15, tplot, rcloud, tff)
    rad, rhodphidrint = fluxintvsr(fluxdata, 15, tplot, rcloud, tff)
    rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
    rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)

    fedd = rhofr * kappa / (rhodphidr * clight)
    feddspec = fr * kappa / (dphidr * clight)
    feddcum = rhofrint * kappa / (rhodphidrint * clight)
    
    feddinds = np.where(fedd > 0.0)
    feddcuminds = np.where(feddcum > 0.0)
    
    feddcums.append(np.convolve(feddcum[feddcuminds], np.ones((nconv,))/nconv, mode='same'))
    fedds.append(np.convolve(fedd[feddinds], np.ones((nconv,))/nconv, mode='same'))
    feddspecs.append(np.convolve(feddspec, np.ones((nconv,))/nconv, mode='same'))
    
    rads.append(rad[feddinds])
    radcums.append(rad[feddcuminds])


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
t1, t2, t3, t4 = plt.plot(rads[0],fedds[0],'k',rads[1],fedds[1],'r',rads[2],fedds[2],'b',rads[3],fedds[3],'g')
#plt.legend((t1,t2,t3,t4), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of,final}$"),prop={'size':8},loc=1)
plt.legend((t1,t2,t3,t4), (r"$\displaystyle t_{\rm 10}$",r"$\displaystyle t_{\rm 50}$",r"$\displaystyle t_{\rm 90}$",r"$\displaystyle t_{\rm of, 50}$"),prop={'size':8},loc=1)
plt.axis([0,2,1.0e-2,1.0e4])
plt.xticks([0,1,2])
#plt.yticks([0,5,10])
plt.yscale('log')
plt.ylabel(r"$\displaystyle f_{\rm Edd}$")
plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")

axrad = plt.subplot(1,2,2)
axrad.yaxis.set_label_position("right")
plt.plot(radcums[0],feddcums[0],'k',radcums[1],feddcums[1],'r',radcums[2],feddcums[2],'b',radcums[3],feddcums[3],'g')
plt.axis([0,2,1.0e-2,1.0e4])
plt.xticks([0,1,2])
#plt.yticks([0,5,10])
plt.yscale('log')
plt.ylabel(r"$\displaystyle f_{\rm Edd, cum}$")
plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")

pp = PdfPages(plotdir + 'f9.pdf')
pp.savefig()
pp.close()




