from matplotlib.backends.backend_pdf import PdfPages
from hyp_math import *
import numpy as np
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

dslist = [10.0, 100.0, 1000.0]
nds = len(dslist)
epsin = np.linspace(0,1,num=1000)
epsout = []
x = 1.0
sigma = 1.5
psi = 2000.0

for i in xrange(0,nds):
    sigma0 = dslist[i]
    epsof = eps_of_all(epsin, x, sigma, sigma0, psi)
    epsout.append(epsof)

dslist = [0.5, 1.5, 2.5]
nds = len(dslist)
sepsout = []
sigma0 = 100.0

for i in xrange(0,nds):
    sigma = dslist[i]
    epsof = eps_of_all(epsin, x, sigma, sigma0, psi)
    sepsout.append(epsof)

plt.figure(figsize = [xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,1,1)
plt.subplots_adjust(left=0.2)
plow, pmid, phigh, = plt.plot(epsin, epsout[0], 'k', epsin, epsout[1], 'r', epsin, epsout[2], 'b')
plt.legend((plow, pmid, phigh), (r"$\displaystyle \Sigma_{\rm cl,0} = 10~M_\odot~{\rm pc^{-2}}$",r"$\displaystyle 10^2$",r"$\displaystyle 10^3$"),prop={'size':8})
plt.axis([0,1,0,1])
plt.xticks([0,0.5,1],[' ',' ',' '])
plt.yticks([0,0.5,1])

#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \varepsilon_{\rm of}$")
plt.text(0.05*1,0.9*1,r"$\displaystyle(a)$")

plt.subplot(2,1,2)
plow, pmid, phigh, = plt.plot(epsin, sepsout[0], 'k', epsin, sepsout[1], 'r', epsin, sepsout[2], 'b')
plt.legend((plow, pmid, phigh), (r"$\displaystyle \sigma_{{\rm ln}\Sigma} = 0.5$",r"$\displaystyle 1.5$",r"$\displaystyle 2.5$"),prop={'size':8})
plt.axis([0,1,0,1])
plt.xticks([0,0.5,1])
plt.yticks([0,0.5,1])

plt.xlabel(r"$\displaystyle \varepsilon$")
plt.ylabel(r"$\displaystyle \varepsilon_{\rm of}$")
plt.text(0.05*1,10**(-1+0.9*3),r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'f27adj.pdf')
pp.savefig()
pp.close()


