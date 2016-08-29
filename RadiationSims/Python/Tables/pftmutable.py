from matplotlib.backends.backend_pdf import PdfPages
from ..Hyperion.hyp_read import *
from ..Hyperion.hyp_hst import *
from ..Hyperion.hyp_out import *
from ..Hyperion.hyp_pdf import *
import matplotlib.pyplot as plt

plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Dissertation_2/Figures/'
datadir = '/scratch/gpfs/raskutti/RadParGrav/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@tiger.princeton.edu'

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
nconv = 50
nconvall = 50
nconvallh = 50
bzwrong = 1

dflist = ['UV_M5.0e4_R15.0_N256_Tf4_NF_B50.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B20.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B10.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B2.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B1.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.5_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.2_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.1_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_NF_B0.05_LD_SDs/']
betas = np.asarray([50.0, 20.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05])
cs = 2.0e4
mcloud = 5.0e4 * 2.0e33
rcloud = 15 * 3.09e18
rhobar = mcloud / (4.0 / 3.0 * np.pi * (rcloud**3))
bz = np.sqrt(8.0 * np.pi * cs**2 * rhobar / betas)
if bzwrong == 1:
    bz = bz * np.sqrt(4.0 * np.pi) * 0.97688287
mf = 2 * np.sqrt(6.67e-8) * mcloud / (bz * rcloud**2)
va = bz / np.sqrt(4.0 * np.pi * rhobar)

print mf
print bz * 1.0e6
print va / 1.0e5
