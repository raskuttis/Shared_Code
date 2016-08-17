import numpy as np
import sys
import random
import scipy.io as spio
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Code to plot the raw song features as a scatter plot
plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Figures'

# Define characteristics of feature vector to import
featdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Features'

outfxname = '_key_All.mat'
outfyname = '_tempo_All.mat'

# Import feature vector and corresponding labels
FVs = featdir + '/FV' + outfxname
LBs = featdir + '/LB' + outfxname
featx = spio.loadmat(FVs)['FV'][0]
labels = spio.loadmat(LBs)['LB'][0]
N = len(featx)
plotx = []
    
FVs = featdir + '/FV' + outfyname
LBs = featdir + '/LB' + outfyname
featy = spio.loadmat(FVs)['FV'][0]
labels = spio.loadmat(LBs)['LB'][0]
ploty = []

genres = ["Blues", "Classical", "Country", "Disco", "Hiphop", "Jazz", "Metal", "Pop", "Reggae", "Rock"]
ngenres = len(genres)

for i in range(0,N):
    plotx.append(featx[i][0][0])
    ploty.append(featy[i][0][0])

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)
    
plt.subplot(1,1,1)
plt.scatter(featx, featy,c=labels,marker='+')
plt.xlabel("Key")
plt.ylabel("Tempo")
    
pp = PdfPages(plotdir + '/KeyvsTempo.pdf')
pp.savefig()
pp.close()
