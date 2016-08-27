from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
from hyp_out import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/tigress-hsm/raskutti/tigress-hsm/Hyperion/RadParGrav/'
outdatadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
outfile = 'RadParGrav.out'
outhostname = 'raskutti@bellona.astro.princeton.edu'
vtkbase = 'RadParGrav'
hostname = 'raskutti@tiger.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_AllOut/'
outdatafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'

outlines = read_outfile(outhostname,outdatadir + outdatafolder + outfile)
rhocloud = out_rhocloud(outlines)
data = AthenaDataRemote(hostname,datadir + datafolder,vtkbase,20,64)
xplot, yplot, denplot = data.get_slice_xyz('density',xcut=0.0)

lvindices = denplot < rhocloud
hvindices = denplot >= rhocloud
logdenplot = denplot
logdenplot[lvindices] = 0.0
logdenplot[hvindices] = 1.0

ysize = 0.22 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,2*ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(2,1,1)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
plt.pcolormesh(xplot, yplot, denplot, cmap='Blues', norm=LogNorm(vmin=0.01, vmax=100000))
#plt.axis([1.0e1,1.0e3,0,1])
#plt.yticks([0,0.2,0.4,0.6,0.8,1])
#plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle Y$")
plt.xlabel(r"$\displaystyle Z$")

plt.subplot(2,1,2)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
plt.pcolormesh(xplot, yplot, logdenplot, cmap='Blues')
#plt.axis([1.0e1,1.0e3,0,1])
#plt.yticks([0,0.2,0.4,0.6,0.8,1])
#plt.xticks([1.0e1,1.0e2,1.0e3])

plt.ylabel(r"$\displaystyle Y$")
plt.xlabel(r"$\displaystyle Z$")

pp = PdfPages(plotdir + 'fden.pdf')
pp.savefig()
pp.close()


