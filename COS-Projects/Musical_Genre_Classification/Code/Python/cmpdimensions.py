from classify import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Code to compare the error for features across a number of different classifiers
plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Figures'
featdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Features'

# Define characteristics of feature vector to import
Nclust = 3
Nex = 5
pcavals = [25,50,65,80,90,95,100]
classlist = ['GNB', 'LSVM', 'KNN', 'RF', 'RBFSVM']
collist = ['k','r','b','g','c','m']
npca = len(pcavals)
nclasses = len(classlist)
featname = 'mfc'
nfeats = [0 for x in range(0,npca)]

corrects = []
correrrs = []

for j in range(0,nclasses):

    classname = classlist[j]
    print classname

    correct = []
    correrr = []

    for i in range(0,npca):
        
        pcaval = pcavals[i]
        print pcaval
        if pcaval == 100:
            outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_All.mat'
        else:
            outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_PCA' + str(pcaval) + '.mat'

        FVs = featdir + '/FV' + outfname
        LBs = featdir + '/LB' + outfname
        feat = spio.loadmat(FVs)['FV']
        nfeats[i] = feat.shape[0]
        labels = spio.loadmat(LBs)['LB'][0]
        [fals, falserr] = class_error(feat, labels, classname)
        print 1.0 - fals
        correct.append(1.0 - fals)
        correrr.append(falserr)

    corrects.append(correct)
    correrrs.append(correrr)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [2*xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,2,1)
#plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
allps = []
for j in range(0,nclasses):
    pfid = plt.errorbar(pcavals,corrects[j],yerr=correrrs[j])
    allps.append(pfid)
plt.axis([0.0,100.0,0.0,1.0])
plt.legend(allps, classlist,loc=2,fontsize='8')
plt.xlabel(r"PCA $\%$")
plt.ylabel(r"$\%$ Accurate Classification")

plt.subplot(1,2,2)
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
allps = []
for j in range(0,nclasses):
    pfid = plt.errorbar(nfeats,corrects[j],yerr=correrrs[j])
    allps.append(pfid)
plt.axis([0.0,2000.0,0.0,1.0])
plt.xlabel(r"N Features")
#plt.ylabel(r"$\%$ Accurate Classification")

pp = PdfPages(plotdir + '/Compare_PCA_Error.pdf')
pp.savefig()
pp.close()
