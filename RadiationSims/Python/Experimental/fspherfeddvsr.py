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
fluxfile = 'fluxes_r0.dat'
hstfile = 'id0/RadParGrav.hst'
starfile = 'star'
nconv = 20

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_Fluxes/'
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

time = hst_time(hstdata, tff)
eff = hst_mstar(hstdata, mcloud)
mgas = hst_mgas(hstdata, mcloud)
stardata = read_allstars(hostname,datadir + datafolder + starfile,time,tff)

tplot = 0.6
allstarrad, allstarmvsr = star_massvsr(stardata)
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhofrint = fluxintoutvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, rhodphidrint = fluxintoutvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
tmini = np.abs(time - tplot).argmin()
mstar = eff[tmini] * mcloud
lstar = psi * mstar
radearly = rad
rstarearly = starrad / rcloud
mstarearly = starmvsr / mstar
massratearly = mr / (mgas[tmini] * mcloud)
rhoratearly = rho / rhocloud
phiratearly = (rad*rcloud)**2 * dphidr / (gravc * (mstar + mgas[tmini] * mcloud))
phiratearlytot = phiratearly * eff[tmini]
radratearly = 4.0 * np.pi * (rad*rcloud)**2 * fr / lstar
#phiratearly = 4.0 * np.pi * (rad*rcloud)**2 * dphidr * clight / (kappa * lstar)
fedd = rhofr * kappa / (rhodphidr * clight)
feddspec = fr * kappa / (dphidr * clight)
feddcum = rhofrint * kappa / (rhodphidrint * clight)
feddcumearly = np.convolve(feddcum, np.ones((nconv,))/nconv, mode='same')
feddearly = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
feddspecearly = np.convolve(feddspec, np.ones((nconv,))/nconv, mode='same')

tplot = 1.06
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhofrint = fluxintoutvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, rhodphidrint = fluxintoutvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
tmini = np.abs(time - tplot).argmin()
mstar = eff[tmini] * mcloud
lstar = psi * mstar
radout = rad
rstarout = starrad / rcloud
mstarout = starmvsr / mstar
massrat = mr / (mgas[tmini] * mcloud)
rhorat = rho / rhocloud
phirat = (rad*rcloud)**2 * dphidr / (gravc * (mstar + mgas[tmini] * mcloud))
phirattot = phirat * eff[tmini]
radrat = 4.0 * np.pi * (rad*rcloud)**2 * fr / lstar
#phirat = 4.0 * np.pi * (rad*rcloud)**2 * dphidr * clight / (kappa * lstar)
fedd = rhofr * kappa / (rhodphidr * clight)
feddspec = fr * kappa / (dphidr * clight)
feddcum = rhofrint * kappa / (rhodphidrint * clight)
feddcumrat = np.convolve(feddcum, np.ones((nconv,))/nconv, mode='same')
feddrat = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
feddspecrat = np.convolve(feddspec, np.ones((nconv,))/nconv, mode='same')

tplot = 1.55
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhofrint = fluxintoutvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, rhodphidrint = fluxintoutvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
tmini = np.abs(time - tplot).argmin()
mstar = eff[tmini] * mcloud
lstar = psi * mstar
radmid = rad
rstarmid = starrad / rcloud
mstarmid = starmvsr / mstar
massratmid = mr / (mgas[tmini] * mcloud)
rhoratmid = rho / rhocloud
phiratmid = (rad*rcloud)**2 * dphidr / (gravc * (mstar + mgas[tmini] * mcloud))
phiratmidtot = phiratmid * eff[tmini]
radratmid = 4.0 * np.pi * (rad*rcloud)**2 * fr / lstar
#phiratmid = 4.0 * np.pi * (rad*rcloud)**2 * dphidr * clight / (kappa * lstar)
fedd = rhofr * kappa / (rhodphidr * clight)
feddspec = fr * kappa / (dphidr * clight)
feddcum = rhofrint * kappa / (rhodphidrint * clight)
feddcummid = np.convolve(feddcum, np.ones((nconv,))/nconv, mode='same')
feddmid = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
feddspecmid = np.convolve(feddspec, np.ones((nconv,))/nconv, mode='same')

tplot = 2.07
tmini = np.abs(time - tplot).argmin()
starrad = allstarrad[tmini,:]
starmvsr = allstarmvsr[tmini,:]
sortinds = np.argsort(starrad)
starrad = starrad[sortinds]
starmvsr = starmvsr[sortinds]
starmvsr = np.cumsum(starmvsr)
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhofrint = fluxintoutvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, rhodphidrint = fluxintoutvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
tmini = np.abs(time - tplot).argmin()
mstar = eff[tmini] * mcloud
lstar = psi * mstar
radlate = rad
rstarlate = starrad / rcloud
mstarlate = starmvsr / mstar
massratlate = mr / (mgas[tmini] * mcloud)
rhoratlate = rho / rhocloud
phiratlate = (rad*rcloud)**2 * dphidr / (gravc * (mstar + mgas[tmini] * mcloud))
phiratlatetot = phiratlate * eff[tmini]
radratlate = 4.0 * np.pi * (rad*rcloud)**2 * fr / lstar
#phiratlate = 4.0 * np.pi * (rad*rcloud)**2 * dphidr * clight / (kappa * lstar)
fedd = rhofr * kappa / (rhodphidr * clight)
feddspec = fr * kappa / (dphidr * clight)
feddcum = rhofrint * kappa / (rhodphidrint * clight)
feddcumlate = np.convolve(feddcum, np.ones((nconv,))/nconv, mode='same')
feddlate = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
feddspeclate = np.convolve(feddspec, np.ones((nconv,))/nconv, mode='same')

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
t1, t2, t3, t4 = plt.plot(radearly,feddearly,'k',radout,feddrat,'r',radmid,feddmid,'b',radlate,feddlate,'g')
plt.legend((t1,t2,t3,t4), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of,final}$"),prop={'size':8},loc=1)
plt.axis([0,2,1.0e-2,1.0e4])
plt.xticks([0,1,2])
#plt.yticks([0,5,10])
plt.yscale('log')
plt.ylabel(r"$\displaystyle f_{\rm edd}$")
plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")


pp = PdfPages(plotdir + 'fspherfedd.pdf')
pp.savefig()
pp.close()




