from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_models import *
import matplotlib.pyplot as plt
import numpy as np

## Define the locations where data is located and where to plot from ..Hyperion.hyp_models
hostname, datadir, hstfile, outfile, plotdir = init_dirs('ast')
fname = 'ft2-1.pdf'

## Results from Heyer
mheyerlist = [3.23e5, 7.86e4, 1.41e4, 1.31e4, 1.89e4, 2.61e5, 2.77e3, 1.1e5, 1.28e5, 1.21e5,
          7.85e4, 1.89e5, 1.95e4, 2.33e5, 1.01e4, 4.92e4, 2.41e4, 1.19e5, 3.5e3, 2.02e5,
          6.19e4, 5.83e4, 4.43e4, 3.44e3, 9.2e4, 5.55e4, 3.97e4, 5.32e3, 2.08e6, 1.21e5,
          6.68e3, 4.13e5, 1.45e4, 4.18e5, 5.81e3, 2.1e4, 8.19e4, 4.55e3, 5.27e4, 1.19e5,
          2.48e4, 3.77e4, 7.69e3, 1.16e5, 4.17e4, 6.89e4, 1.13e5, 1.19e5, 9.07e4, 4.4e3,
          5.94e5, 5.52e5, 1.87e4, 3.16e4, 1.5e5, 3.55e3, 1.15e5, 2.06e5, 4.14e3, 1.02e4,
          7.04e3, 4.48e4, 1.49e5, 8.11e4, 7.88e3, 6.97e3, 3.05e4, 2.55e3, 1.42e4, 1.85e4,
          1.32e5, 6.45e3, 1.66e4, 1.88e4, 1.55e5, 2.19e4, 2.48e5, 8.32e4, 2.6e4, 2.23e5,
          2.16e5, 6.74e4, 3.15e4, 1.82e5, 9.35e3, 9.68e4, 3.74e5, 2.82e4, 5.87e4, 2.69e5,
          8.28e5, 3.45e4, 2.04e4, 1.67e4, 8.7e4, 2.69e5, 3.35e4, 3.17e4, 1.1e6, 6.9e3,
          2.61e3, 2.22e5, 3.96e4, 4.84e4, 5.06e5, 1.63e4, 8.94e4, 3.08e4, 1.57e5, 4.02e4,
          1.38e5, 2.81e4, 1.06e5, 1.4e5, 1.16e5, 3.11e4, 2.48e4, 7.76e4, 3.56e4, 1.07e5,
          7.86e3, 3.37e4, 3.6e4, 1.43e5, 7.76e3, 4.5e3, 7.36e4, 5.13e4, 5.95e4, 3.66e4,
          1.24e5, 2.95e3, 7.71e4, 3.16e4, 2.28e4, 1.19e6, 8.35e3, 1.29e4, 2.52e5, 6.72e4,
          6.33e4, 5.52e4, 5.74e4, 6.79e3, 4.59e5, 4.43e5, 5.82e4, 6.81e4, 1.95e4, 1.28e5,
          1.89e4, 7.69e3, 9.99e4, 1.9e4, 1.28e4, 2.26e4, 5.03e4, 1.12e5]
rheyerlist = [40.6, 42.2, 6.9, 9.7, 10.3, 41.3, 4.6, 32.1, 23.2, 33.8,
          22.8, 27.1, 16.7, 45.0, 15.3, 23.7, 18.1, 32.0, 7.5, 46.0,
          22.9, 18.4, 17.1, 10.6, 11.4, 16.0, 19.8, 6.4, 83.4, 28.2,
          8.8, 23.3, 11.6, 29.8, 9.7, 8.2, 22.7, 5.7, 25.4, 19.3,
          9.7, 20.2, 9.1, 26.4, 12.0, 25.8, 26.3, 23.3, 34.7, 7.3,
          52.0, 53.7, 16.3, 17.9, 42.7, 7.2, 34.0, 37.4, 9.6, 12.8,
          8.4, 16.1, 25.8, 38.5, 12.2, 12.5, 16.9, 6.1, 17.7, 14.8,
          29.0, 10.5, 9.2, 17.0, 23.3, 15.0, 37.2, 19.5, 13.8, 24.3,
          26.5, 16.6, 16.3, 31.6, 9.4, 23.8, 41.1, 14.0, 19.9, 40.8,
          34.9, 12.9, 13.2, 5.0, 26.1, 44.7, 14.4, 19.4, 41.8, 8.8,
          8.4, 36.2, 12.5, 14.7, 58.6, 13.6, 21.9, 19.0, 37.7, 15.0,
          44.4, 23.0, 34.0, 35.7, 32.9, 21.3, 19.9, 30.1, 21.5, 18.1,
          7.6, 16.1, 15.7, 36.0, 12.8, 6.4, 31.7, 17.2, 25.3, 13.9,
          31.8, 6.4, 20.8, 10.8, 10.8, 90.8, 12.5, 11.4, 42.6, 25.3,
          20.6, 15.9, 16.9, 7.1, 46.1, 47.0, 13.6, 19.8, 17.7, 32.0,
          11.9, 7.3, 31.3, 12.3, 9.0, 15.9, 12.8, 28.2]
vheyerlist = [2.6, 1.5, 1.9, 3.0, 2.2, 2.4, 1.3, 2.7, 2.8, 2.9,
          1.9, 3.1, 3.1, 4.0, 1.4, 2.0, 1.9, 1.7, 2.2, 2.8,
          2.8, 1.9, 1.1, 0.5, 3.1, 2.8, 2.3, 2.0, 2.3, 1.6,
          0.8, 4.6, 1.5, 5.9, 2.0, 2.1, 1.8, 1.6, 1.7, 2.7,
          1.7, 0.9, 1.0, 2.5, 1.9, 2.1, 3.0, 3.9, 1.8, 1.8,
          4.1, 2.5, 1.7, 1.8, 2.0, 1.3, 2.1, 2.5, 1.4, 1.5,
          1.2, 2.0, 2.5, 3.3, 1.6, 1.0, 1.7, 1.8, 1.6, 1.4,
          2.6, 3.2, 1.6, 1.0, 3.6, 1.5, 3.7, 1.3, 1.8, 3.5,
          4.1, 3.2, 0.8, 2.5, 1.9, 2.5, 2.9, 1.7, 2.4, 2.9,
          4.5, 2.1, 1.1, 1.2, 1.9, 2.9, 2.4, 0.8, 6.8, 1.8,
          0.7, 3.3, 2.6, 3.2, 3.7, 1.6, 2.6, 2.5, 2.1, 1.4,
          1.7, 1.6, 2.4, 2.2, 3.4, 1.5, 1.2, 2.2, 2.1, 2.8,
          2.3, 3.3, 1.0, 2.7, 1.2, 1.7, 3.6, 2.0, 3.1, 2.2,
          3.1, 1.6, 3.4, 2.7, 1.5, 3.8, 0.9, 1.8, 1.8, 1.4,
          3.0, 2.2, 3.8, 1.2, 4.0, 4.4, 3.4, 2.1, 1.7, 1.9,
          1.1, 2.0, 2.5, 1.6, 2.3, 1.1, 2.5, 2.9]
mheyer = np.asarray(mheyerlist)
rheyer = np.asarray(rheyerlist)
vheyer = np.asarray(vheyerlist)
alphaheyer = (vheyer * 1.0e5)**2 * rheyer * 3.09e18 * 5.0 / (3.0 * 6.67e-8 * mheyer * 2.0e33)
sigmaheyer = mheyer / (np.pi * rheyer * rheyer)
vtheyer = np.sqrt(alphaheyer * 0.6 * 6.67e-8 * mheyer * 2.0e33 / (rheyer * 3.09e18)) / 1.0e5

## Results from Duval
filename = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Data/Duval2010.txt'
with open(filename,'r') as f:
    data = np.loadtxt(f,usecols=(6,7,12,13))

rduval = data[:,0]
mduval = data[:,1]
sigmaduval = data[:,2]
alphaduval = data[:,3]
vtduval = np.sqrt(alphaduval * 0.6 * 6.67e-8 * mduval * 2.0e33 / (rduval * 3.09e18)) / 1.0e5

## List of masses and radii from our sims to plot
mlist = np.asarray([2.0e4, 5.0e4, 2.0e4, 5.0e4, 1.0e5, 2.0e4, 1.0e4, 5.0e4, 1.0e4, 1.0e5, 2.0e5, 5.0e3, 2.0e4, 5.0e4, 1.0e5, 2.0e4, 2.0e5, 1.0e4, 5.0e5, 1.0e5, 5.0e4, 2.0e5, 5.0e4, 2.0e4, 2.0e5])
rlist = np.asarray([25.0, 35.0, 20.0, 25.0, 35.0, 15.0, 10.0, 20.0, 8.0, 25.0, 35.0, 5.0, 10.0, 15.0, 20.0, 8.0, 25.0, 5.0, 35.0, 15.0, 10.0, 20.0, 8.0, 5.0, 15.0])
alpha = 2.0

cmasses = mlist
csigmas = mlist / (np.pi * rlist**2)
cvturb = np.sqrt(3.0 * 6.67e-8 * mlist * 2.0e33 * alpha / (5.0 * rlist * 3.09e18)) / 1.0e5

alist = np.asarray([0.1, 0.2, 0.4, 0.8, 1.5, 2.0, 3.0, 6.0, 10.0])

amasses = 5.0e4 * np.ones(9)
asigmas = 5.0e4 / (np.pi * 15.0**2) * np.ones(9)
avturb = np.sqrt(3.0 * 6.67e-8 * 5.0e4 * 2.0e33 * alist / (5.0 * 15.0 * 3.09e18)) / 1.0e5

## Plotting setup
ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

xcap = 10**(3+0.05*3)
ycap = 10**(0.9*3)
ycaptwo = 10**(0.9)
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
plt.subplots_adjust(bottom=0.2)
plt.plot(mheyer,sigmaheyer,'k.',mduval,sigmaduval,'r.',markersize=1.5)
plt.plot(cmasses,csigmas,'bv',amasses,asigmas,'gv')
plt.text(xcap,ycap,r"$\displaystyle(a)$")
plt.yscale('log')
plt.xscale('log')
plt.axis([1.0e3,1.0e6,1.0,1.0e3])
plt.yticks([1.0,1.0e1,1.0e2,1.0e3])
plt.xticks([1.0e3,1.0e4,1.0e5,1.0e6])

plt.ylabel(r"$\displaystyle \langle \Sigma \rangle / M_{\odot} {\rm pc^{-2}}$")
plt.xlabel(r"$\displaystyle M / M_{\odot}$")

axv = plt.subplot(1,2,2)
axv.yaxis.set_label_position("right")
plt.plot(mheyer,vtheyer,'k.',mduval,vtduval,'r.',markersize=1.5)
plt.plot(cmasses,cvturb,'bv',amasses,avturb,'gv')
plt.text(xcap,ycaptwo,r"$\displaystyle(b)$")
plt.yscale('log')
plt.xscale('log')
plt.axis([1.0e3,1.0e6,1.0,1.0e1])
plt.yticks([1.0,1.0e1])
plt.xticks([1.0e3,1.0e4,1.0e5,1.0e6])

plt.ylabel(r"$\displaystyle v_{\rm RMS} / {\rm km s^{-1}}$")
plt.xlabel(r"$\displaystyle M / M_{\odot}$")

pp = PdfPages(plotdir + fname)
pp.savefig()
pp.close()


