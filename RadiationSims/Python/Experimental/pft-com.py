from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_fil import *
from scipy.interpolate import RegularGridInterpolator as rgi
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
comfile = 'com_star.dat'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'

plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

axa = plt.subplot(1,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
pa = plt.plot(rfil,fildens,'k',perpr,perpden,'r')
plt.axis([0.0,25.0,1.0e-1,1.0e4])
#plt.xticks([0,1,2,3])
#plt.yticks([0,20,40])
plt.yscale('log')
plt.xlabel(r"$\displaystyle r$")
plt.ylabel(r"$\displaystyle n_H$")

pp = PdfPages(plotdir + 'ft-fil.pdf')
pp.savefig()
pp.close()


