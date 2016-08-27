from matplotlib.backends.backend_pdf import PdfPages
from hyp_read import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Figures/'
datadir = '/tigress-hsm/raskutti/tigress-hsm/Hyperion/RadParGrav/'
locdatadir = '/Users/sudhirraskutti/Desktop/Thesis/PaperII/Visit/'
vtkfile = 'id0/RadParGrav.0000.vtk'
vtkbase = 'RadParGrav'
hostname = 'raskutti@tiger.princeton.edu'

datafolder = 'UV_M5.0e4_R15.0_N256_Tf4_AllOut/'
locdatafolder = 'Fiducial/'

#locdata = AthenaData(locdatadir + locdatafolder + vtkfile)
#momentum_x = locdata.get_field('momentum_x')
#print momentum_x[0,3,7]

data = AthenaDataRemote(hostname,datadir + datafolder,vtkbase,64)
density = data.get_field('density')
print np.shape(density)


