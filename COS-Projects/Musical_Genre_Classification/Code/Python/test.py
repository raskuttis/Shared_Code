from classify import *

# Define characteristics of feature vector to import
Nclust = 3
Nex = 5
featname = 'mfc'
outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_All.mat'
#outfname = '_' + featname + '_All.mat'

ae = class_error(outfname, 'GNB')
print ae
