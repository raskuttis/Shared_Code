import numpy as np
import copy as cp

def set_model_lookup(suffix):
    
    mlist = ['2.0e4', '5.0e4', '2.0e4', '5.0e4', '1.0e5', '2.0e4', '5.0e4', '1.0e5', '2.0e5', '5.0e3', '2.0e4', '5.0e4', '1.0e5', '2.0e4', '2.0e5', '1.0e4', '1.0e5', '5.0e4', '2.0e5', '5.0e4', '2.0e4', '2.0e5']
    rlist = ['25.0', '35.0', '20.0', '25.0', '35.0', '15.0', '10.0', '25.0', '35.0', '5.0', '10.0', '15.0', '20.0', '8.0', '25.0', '5.0', '15.0', '10.0', '20.0', '8.0', '5.0', '15.0']
    reslist = ['256', '256', '256', '256', '256', '256', '256', '256', '256', '256', '128', '256', '256', '256', '256', '256', '256', '256', '256', '256', '256', '128']
    supplist = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']

    nfs = len(mlist)
    dflist = np.core.defchararray.add(['UV_M'] * nfs, mlist)
    dflist = np.core.defchararray.add(dflist, ['_R'] * nfs)
    dflist = np.core.defchararray.add(dflist, rlist)
    dflist = np.core.defchararray.add(dflist, ['_N'] * nfs)
    dflist = np.core.defchararray.add(dflist, reslist)
    dflist = np.core.defchararray.add(dflist, ['_Tf4_'] * nfs)
    dflist = np.core.defchararray.add(dflist, supplist)
    dflist = np.core.defchararray.add(dflist, [suffix] * nfs)

    return dflist
    
def init_dirs(loc = 'tiger'):

    plotdir = '/Users/sudhirraskutti/Desktop/Thesis/Figures/'
    hstfile = 'id0/RadParGrav.hst'
    outfile = 'RadParGrav.out'

    if loc == 'ast':
        hostname = 'raskutti@bellona.astro.princeton.edu'
        datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
    if loc == 'tiger':
        hostname = 'raskutti@tiger.princeton.edu'
        datadir = '/scratch/gpfs/raskutti/RadParGrav/'

    return hostname, datadir, hstfile, outfile, plotdir

def set_mag_model_lookup():

    dflist = ['UV_M5.0e4_R15.0_N256_Tf4_B200.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B100.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B50.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B20.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B10.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B5.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B2.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B1.0_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.5_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.1_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.2_LD_SDs/', 'UV_M5.0e4_R15.0_N256_Tf4_B0.05_LD_SDs/']
    betas = [200.0, 100.0, 50.0, 20.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]

    return dflist, betas