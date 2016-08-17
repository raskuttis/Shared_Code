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
classlist = ['KNN']
collist = ['k','r','b','g','c','m']
nclasses = len(classlist)
featname = 'mfc'
nns = [1,3,5,7,10,15,20,50,100]

corrects = []
correrrs = []

for j in range(0,nclasses):

    classname = classlist[j]
    print classname

    correct = []
    correrr = []

    outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_PCA80.mat'

    FVs = featdir + '/FV' + outfname
    LBs = featdir + '/LB' + outfname
    feat = spio.loadmat(FVs)['FV']
    labels = spio.loadmat(LBs)['LB'][0]
    
    for i in range(0,len(nns)):
        [fals, falserr] = class_error(feat, labels, classname, pin = nns[i])
        print 1.0 - fals
        correct.append(1.0 - fals)
        correrr.append(falserr)

    corrects.append(correct)
    correrrs.append(correrr)

ysize = 0.25 * 11.69
xsize = 0.4 * 8.27
fontsize = '10'
plt.figure(figsize = [xsize,ysize])
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)

plt.subplot(1,1,1)
#plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
allps = []
for j in range(0,nclasses):
    pfid = plt.errorbar(nns,corrects[j],yerr=correrrs[j])
    allps.append(pfid)
plt.axis([0.0,100.0,0.0,1.0])
plt.legend(allps, classlist,loc=2,fontsize='8')
plt.xlabel(r"NN $\%$")
plt.ylabel(r"$\%$ Accurate Classification")

pp = PdfPages(plotdir + '/Compare_KNN_Neighbours.pdf')
pp.savefig()
pp.close()
