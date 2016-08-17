from classify import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Code to compare the error for features across a number of different classifiers
plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Figures'
featdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Features'

# Define characteristics of feature vector to import
Nclust = 3
Nex = 5
classname = 'RBFSVM'
featname = 'mfc'

outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_PCA90.mat'

FVs = featdir + '/FV' + outfname
LBs = featdir + '/LB' + outfname
feat = spio.loadmat(FVs)['FV']
labels = spio.loadmat(LBs)['LB'][0]

[fals, falserr, cfmat] = class_error(feat, labels, classname)

print 1.0-fals, falserr
print cfmat