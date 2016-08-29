from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_fil import *
from scipy.interpolate import RegularGridInterpolator as rgi
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
locdatadir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/'
vtkfile = 'RadParGrav_joined.0002.vtk'
filfile = 'den02.fits_c100.up.NDskl.TRIM.ASMB.S006.BRK.a.NDskl'
locdatafolder = 'RB0.05_N128/'

nres = 128
locdata = AthenaData(locdatadir + locdatafolder + vtkfile)
dens = locdata.get_field('density')

x = np.linspace(0, nres-1, nres) + 0.5
y = np.linspace(0, nres-1, nres) + 0.5
z = np.linspace(0, nres-1, nres) + 0.5

fn = rgi((x,y,z), dens)

fillines = read_local_outfile(locdatadir + locdatafolder + filfile)
fildata = out_fil(fillines)
fillocs = find_filheaders(fildata)
filcoords = get_fil(fildata, fillocs, 0)

pts = np.transpose(filcoords)
fildens = fn(pts)

rfil, rperp, perpcoords, nlen, ndivs = extract_r(filcoords)
allperps = []
allrperps = []
for i in xrange(0,nlen-2):
    pts = np.squeeze(perpcoords[i,:,:])
    ptdens = fn(pts)
    allperps.append(ptdens)
    allrperps.append(rperp)

allperps = np.reshape(np.asarray(allperps), (nlen-2)*ndivs)
allrperps = np.reshape(np.asarray(allrperps), (nlen-2)*ndivs)

nbins = 10
perpden, perpr = np.histogram(allrperps,bins=10,weights=allperps)
perpnum, perpr = np.histogram(allrperps,bins=10)
perpr = perpr[0:nbins]
perpden = perpden / perpnum

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


