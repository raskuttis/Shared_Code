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

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_NF_Fluxes/'
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
nconv = 20

tplot = 0.57
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
aradearly = rad
arhoratearly = mr / mcloud
arhofitearly = 10**pcoeffs[1] * (rad)**pcoeffs[0]

tplot = 1.06
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
aradout = rad
arhorat = mr / mcloud
arhofit = 10**pcoeffs[1] * (rad)**pcoeffs[0]

tplot = 1.55
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
aradmid = rad
arhoratmid = mr / mcloud
arhofitmid = 10**pcoeffs[1] * (rad)**pcoeffs[0]

tplot = 2.07
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
aradlate = rad
arhoratlate = mr / mcloud
arhofitlate = 10**pcoeffs[1] * (rad)**pcoeffs[0]


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
nconv = 20

tplot = 0.57
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
radearly = rad
rhoratearly = mr / mcloud
rhofitearly = 10**pcoeffs[1] * (rad)**pcoeffs[0]

tplot = 1.06
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
radout = rad
rhorat = mr / mcloud
rhofit = 10**pcoeffs[1] * (rad)**pcoeffs[0]

tplot = 1.55
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
radmid = rad
rhoratmid = mr / mcloud
rhofitmid = 10**pcoeffs[1] * (rad)**pcoeffs[0]

tplot = 2.07
rad, mr = fluxintvsr(fluxdata, 5, tplot, rcloud, tff)
goodinds = np.logical_and((rad > 0.3),(rad < 1.0))
pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mcloud), 1)
print pcoeffs
rad, rho = fluxvsr(fluxdata, 5, tplot, rcloud, tff)
rad, rhofr = fluxvsr(fluxdata, 7, tplot, rcloud, tff)
rad, rhodphidr = fluxvsr(fluxdata, 15, tplot, rcloud, tff)
rad, fr = fluxvsr(fluxdata, 6, tplot, rcloud, tff)
rad, dphidr = fluxvsr(fluxdata, 14, tplot, rcloud, tff)
radlate = rad
rhoratlate = mr / mcloud
rhofitlate = 10**pcoeffs[1] * (rad)**pcoeffs[0]

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
plt.plot(radearly,rhofitearly,'--k',radout,rhofit,'--r',radmid,rhofitmid,'--b',radlate,rhofitlate,'--g')
t1, t2, t3, t4 = plt.plot(radearly,rhoratearly,'k',radout,rhorat,'r',radmid,rhoratmid,'b',radlate,rhoratlate,'g')
plt.legend((t1,t2,t3,t4), (r"$\displaystyle \varepsilon(t) = 0.1 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.5 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon(t) = 0.9 \varepsilon_{\rm final}$",r"$\displaystyle \varepsilon_{\rm of}(t) = 0.5 \varepsilon_{\rm of,final}$"),prop={'size':8},loc=2)
plt.text(10**(-1+0.9*1.3),10**(-4+0.9*6),r"$\displaystyle(a)$")
plt.axis([1.0e-1,2.0,1.0e-4,1.0e2])
plt.yscale('log')
plt.yticks([1.0e-4,1.0e-2,1.0,1.0e2])
plt.xscale('log')
plt.ylabel(r"$\displaystyle M(r) / M_{\rm cl,0}$")
plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")

plt.subplot(1,2,2)
plt.subplots_adjust(bottom=0.2)
plt.plot(aradearly,arhofitearly,'--k',aradout,arhofit,'--r',aradmid,arhofitmid,'--b',aradlate,arhofitlate,'--g')
t1, t2, t3, t4 = plt.plot(aradearly,arhoratearly,'k',aradout,arhorat,'r',aradmid,arhoratmid,'b',aradlate,arhoratlate,'g')
plt.text(10**(-1+0.9*1.3),10**(-4+0.9*6),r"$\displaystyle(b)$")
plt.axis([1.0e-1,2.0,1.0e-4,1.0e2])
plt.yscale('log')
plt.yticks([1.0e-4,1.0e-2,1.0,1.0e2],[' ',' ',' ',' '])
plt.xscale('log')
#plt.ylabel(r"$\displaystyle \langle \rho \rangle / \rho_{\rm cloud}$")
plt.xlabel(r"$\displaystyle r / r_{\rm cloud}$")

pp = PdfPages(plotdir + 'flogrhoprof.pdf')
pp.savefig()
pp.close()




