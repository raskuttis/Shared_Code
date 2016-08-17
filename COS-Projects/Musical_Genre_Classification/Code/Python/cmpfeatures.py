from classify import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Code to compare the error for features across a number of different classifiers
plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Figures'
featdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Features'

# Define characteristics of feature vector to import
Nclust = 3
Nex = 5
featnamelist = ['key','tempo','keystrength','roughness','hcdf','eng','zerocross','chroma','brightness','mfc']
lablist = ['Key', 'Tempo', 'KS', 'R', 'HCDF', 'E', 'ZC', 'CHR', 'SF', 'MFCC', 'ALL']
altlablist = ['Key', 'Tempo', 'KS', 'R', 'HCDF', 'E', 'ZC', 'CHR', 'SF', 'MFCC', 'MFCC + HCDF']
frmlvl = [0,0,1,1,1,1,1,1,1,1]
pcalvl = [0,0,1,1,1,1,1,1,1,1]
incallfeat = [0,0,0,0,1,0,0,0,0,1]
classlist = ['GNB', 'LSVM', 'KNN', 'RF', 'RBFSVM']
collist = ['k','r','b','g','c']
nfeats = len(featnamelist)
nclasses = len(classlist)

corrects = []
correrrs = []

for j in range(0,nclasses):

    classname = classlist[j]
    print classname

    correct = []
    correrr = []
    allfeats = []

    for i in range(0,nfeats):
        
        featname = featnamelist[i]
        print featname
        frmflag = frmlvl[i]
        if frmflag == 1:
            outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_All.mat'
        else:
            outfname = '_' + featname + '_All.mat'

        FVs = featdir + '/FV' + outfname
        LBs = featdir + '/LB' + outfname
        feat = spio.loadmat(FVs)['FV']
        labels = spio.loadmat(LBs)['LB'][0]
        if i == 0:
            allfeats = feat
        else:
            allfeats = np.concatenate((allfeats, feat))
        [fals, falserr] = class_error(feat, labels, classname)
        print 1.0 - fals
        correct.append(1.0 - fals)
        correrr.append(falserr)

    # Finally fit using all the features
    print 'all'
    print allfeats.shape
    [fals, falserr] = class_error(allfeats, labels, classname)
    print 1.0 - fals
    correct.append(1.0 - fals)
    correrr.append(falserr)

    corrects.append(correct)
    correrrs.append(correrr)

pcacorrects = []
pcacorrerrs = []

for j in range(0,nclasses):
    
    classname = classlist[j]
    print classname
    
    correct = []
    correrr = []
    allfeats = []
    incflag = 0
    
    for i in range(0,nfeats):
        
        featname = featnamelist[i]
        print featname
        frmflag = frmlvl[i]
        pcaflag = pcalvl[i]
        if frmflag == 1:
            if pcaflag == 1:
                outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_PCA80.mat'
            else:
                outfname = '_' + featname + '_GMM' + str(Nclust) + '_Nex' + str(Nex) + '_All.mat'
        else:
            outfname = '_' + featname + '_All.mat'
        
        FVs = featdir + '/FV' + outfname
        LBs = featdir + '/LB' + outfname
        feat = spio.loadmat(FVs)['FV']
        labels = spio.loadmat(LBs)['LB'][0]
        if incallfeat[i] == 1:
            if incflag == 0:
                allfeats = feat
                incflag = 1
            else:
                allfeats = np.concatenate((allfeats, feat))
        [fals, falserr] = class_error(feat, labels, classname)
        print 1.0 - fals
        correct.append(1.0 - fals)
        correrr.append(falserr)

    # Finally fit using all the features
    print 'all'
    print allfeats.shape
    [fals, falserr] = class_error(allfeats, labels, classname)
    print 1.0 - fals
    correct.append(1.0 - fals)
    correrr.append(falserr)

    pcacorrects.append(correct)
    pcacorrerrs.append(correrr)

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
    pfid = plt.errorbar(range(0,nfeats+1),corrects[j],yerr=correrrs[j])
    allps.append(pfid)
plt.axis([-0.5,nfeats+1.5,0.0,1.0])
plt.xticks(range(0,nfeats+1),lablist,rotation='vertical')
plt.legend(allps, classlist,loc=2,fontsize='8')
#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\%$ Accurate Classification")

plt.subplot(1,2,2)
plt.subplots_adjust(left=0.2)
plt.subplots_adjust(bottom=0.2)
allps = []
for j in range(0,nclasses):
    pfid = plt.errorbar(range(0,nfeats+1),pcacorrects[j],yerr=pcacorrerrs[j])
    allps.append(pfid)
plt.axis([-0.5,nfeats+1.5,0.0,1.0])
plt.xticks(range(0,nfeats+1),lablist,rotation='vertical')
#plt.legend(allps, classlist,loc=2,fontsize='8')
#plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
#plt.ylabel(r"$\%$ Accurate Classification")

pp = PdfPages(plotdir + '/Compare_Classifier_Error.pdf')
pp.savefig()
pp.close()
