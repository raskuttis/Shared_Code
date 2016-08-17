import numpy as np
import sys
import random
import scipy.io as spio
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def plot_meanfeat(featname, prange):

    # Code to plot the raw features for on a histogram as a function of genre
    plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Figures'

    # Define characteristics of feature vector to import
    featdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Features'
    outfname = '_' + featname + '_All.mat'

    # Import feature vector and corresponding labels
    FVs = featdir + '/FV' + outfname
    LBs = featdir + '/LB' + outfname
    feat = spio.loadmat(FVs)['FV'][0]
    labels = spio.loadmat(LBs)['LB'][0]

    means = []
    frames = []
    lbs = []
    genres = ["Blues", "Classical", "Country", "Disco", "Hiphop", "Jazz", "Metal", "Pop", "Reggae", "Rock"]
    ngenres = len(genres)

    for i in range(1,ngenres+1):
        data = feat[labels == i]
        print i
        for j in range(len(data)):
            test = data[j].T
            test[np.isnan(test)] = 0.0
            N = len(np.mean(test, axis=1))
            means = means + list(np.mean(test, axis=1))
            frames = frames +  range(N)
            lbs = lbs + [i for k in range(N)]

    lbs = np.asarray(lbs)
    means = np.asarray(means)

    gmean = []
    for i in range(1,ngenres+1):
        gmean.append(np.mean(means[lbs == i]))
        print np.mean(means[lbs == i]), np.std(means[lbs == i])
    print np.std(gmean)

    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '10'
    plt.figure(figsize = [2*xsize,ngenres/2*ysize])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)

    for i in range(1,ngenres+1):
        plt.subplot(ngenres/2,2,i)
        plt.hist(means[lbs == i], bins=50, normed=1, facecolor='green')
        plt.xlim(prange[0], prange[1])
        plt.xlabel(genres[i-1])
        plt.ylabel("N")

    pp = PdfPages(plotdir + '/FeatureMean_' + featname + '.pdf')
    pp.savefig()
    pp.close()

featnamelist = ['mfc','eng','zerocross','chroma','brightness','keystrength','roughness','inharmonic','hcdf']
featx = [-0.3,-4.0,0,0,0,-8,0,0.2,0]
featy = [0.3,0,6000,1,1,8,2000,0.6,4]
nfeats = len(featnamelist)

for i in range(0,nfeats):
    featname = featnamelist[i]
    prange = [featx[i], featy[i]]
    plot_meanfeat(featname, prange)
