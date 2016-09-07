import numpy as np
import csv as csv
from collections import OrderedDict
from collections import Counter
from datetime import datetime, timedelta, date
from sklearn.preprocessing import OneHotEncoder
from sklearn.feature_extraction import DictVectorizer
from sklearn.naive_bayes import BernoulliNB
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.cross_validation import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.decomposition import PCA
from sklearn import svm

## Extracts a set of features from the dictionary and converts them to a usable matrix format
def input_matrix(fulldata, featurelist, hkencflag = 0, pcaflag = 0):

    nfeats = len(featurelist)
    ncatalls = []
    ylist = np.asarray(fulldata['CASE_STATUS']).astype(int)
    ndps = len(ylist)
    xmat = np.zeros((ndps, nfeats))
    allcats = []
    for i in xrange(0, nfeats):
        feat = featurelist[i]
        featrep = feat.replace('_',' ')
        xdict = fulldata[feat]
        xdictfreqs = Counter(xdict)
        ncats = len(xdictfreqs)
        ncatalls.append(ncats)
        cats = xdictfreqs.keys()
        catfulls = [featrep + '-' + l for l in cats]
        allcats.append(catfulls)
        xlist = np.zeros(ndps)
        j = 0
        for (k,v) in xdictfreqs.items():
            goodinds = np.squeeze(np.where(k == xdict))
            xlist[goodinds] = j
            j = j + 1

        xmat[:,i] = xlist
    xmathk = xmat

    if (hkencflag == 1):

        enc = OneHotEncoder(n_values = ncatalls)
        enc.fit(xmat)
        xmathk = enc.transform(xmat)
    
    if (pcaflag > 0):

        xmatin = xmathk.toarray()
        pca = PCA(n_components=pcaflag)
        pca.fit(xmatin)
        xmathk = pca.transform(xmatin)
        print np.sum(pca.explained_variance_ratio_)

    return ylist, xmathk, ncatalls, allcats

def add_numfeatures(fulldata, featurelist, xmat):

    xmatin = xmat.astype(float)
    xmatin = xmatin.toarray()
    nfeats = len(featurelist)
    for i in xrange(0, nfeats):
        feat = featurelist[i]
        feats = np.asarray(fulldata[feat]).astype(float)
        ndps = np.shape(xmatin)[0]
        nexfs = np.shape(xmatin)[1]
        xmatin = np.insert(xmatin, nexfs, feats, axis=1)

    return xmatin

def un_encode(ytrain, xtrain, ncats):

    ndps = np.shape(xtrain)[0]
    ntvs = np.shape(xtrain)[1]
    ndiff = ntvs - np.sum(ncats)
    ncmin = 0
    xsmall = np.zeros((ndps,len(ncats) + ndiff))
    for i in xrange(0, len(ncats)):
        nc = ncats[i]
        xsubs = xtrain[:,ncmin:nc+ncmin]
        goodinds = np.where(ytrain == 1)
        xsubgood = xsubs[goodinds]
        xsubars = np.sum(xsubgood, axis = 0) / np.sum(xsubs, axis = 0)
        for j in xrange(0, nc):
            xsubs = xtrain[:,ncmin+j]
            goodinds = np.where(xsubs == 1)
            xsmall[goodinds,i] = xsubars[j]
    xsmall[:,len(ncats):len(ncats)+ndiff] = xtrain[:,np.sum(ncats):np.sum(ncats)+ndiff]
    
    return xsmall

def classifier(ytrain, xtrain, cname, pval = 0.0):

    pv = 0.0
    if (cname == 'nb'):
        nb = BernoulliNB()
        nb.fit(xtrain, ytrain)
        return nb, pv
    if (cname == 'gnb'):
        nb = GaussianNB()
        nb.fit(xtrain, ytrain)
        return nb, pv
    if (cname == 'lr'):
        lr = LogisticRegression()
        lr.fit(xtrain, ytrain)
        return lr, pv
    if (cname == 'rf'):
        rf = RandomForestClassifier(n_estimators=100,criterion="gini")
        rf.fit(xtrain, ytrain)
        return rf, pv
    if (cname == 'knn'):
        knn = KNeighborsClassifier(n_neighbors=50)
        knn.fit(xtrain, ytrain)
        return knn, pv
    if (cname == 'svmrbf'):
        cfsvm = svm.SVC(C=1)
        cfsvm.fit(xtrain, ytrain)
        return cfsvm, pv
    if (cname == 'ada'):
        cfada = AdaBoostClassifier(C=1)
        cfada.fit(xtrain, ytrain)
        return cfada, pv
    if (cname == 'lrc'):
        if pval == 0.0:
            cscores = {}
            for c in [0.001,0.01,0.1,1,10,100,1000]:
                lr = LogisticRegression(C=c)
                scores = cross_val_score(lr, xtrain, ytrain, scoring='roc_auc', cv=10)
                cscores[c] = scores.mean()
            cmax = max(cscores, key=cscores.get)
            pval = cmax
        lr = LogisticRegression(C=pval)
        lr.fit(xtrain, ytrain)
        return lr, pval
    if (cname == 'rfc'):
        if pval == 0.0:
            cscores = {}
            for c in [100,250,500,1000,2500]:
                print c
                rf = RandomForestClassifier(n_estimators=c,criterion="gini")
                scores = cross_val_score(rf, xtrain, ytrain, scoring='roc_auc')
                cscores[c] = scores.mean()
            print cscores, cmax
            cmax = max(cscores, key=cscores.get)
            pval = cmax
        rf = RandomForestClassifier(n_estimators=pval,criterion="gini")
        rf.fit(xtrain, ytrain)
        return rf, pval
    if (cname == 'knnc'):
        if pval == 0.0:
            cscores = {}
            for c in [10,25,50,100]:
                print c
                knn = KNeighborsClassifier(n_neighbors=c)
                scores = cross_val_score(knn, xtrain, ytrain, scoring='roc_auc')
                cscores[c] = scores.mean()
            print cscores, cmax
            cmax = max(cscores, key=cscores.get)
            pval = cmax
        knn = KNeighborsClassifier(n_neighbors=pval)
        knn.fit(xtrain, ytrain)
        return knn, pval
    if (cname == 'svmrbfc'):
        if pval == 0.0:
            cscores = {}
            for c in [0.001,0.01,0.1,1,10,100,1000]:
                print c
                cfsvm = svm.SVC(C=c)
                scores = cross_val_score(cfsvm, xtrain, ytrain, scoring='roc_auc')
                cscores[c] = scores.mean()
            cmax = max(cscores, key=cscores.get)
            pval = cmax
            print cscores
        cfsvm = svm.SVC(C=pval)
        cfsvm.fit(xtrain, ytrain)
        return cfsvm, pval
    if (cname == 'adac'):
        if pval == 0.0:
            cscores = {}
            for c in [50,100,150,200,250]:
                print c
                cfada = AdaBoostClassifier(C=c)
                scores = cross_val_score(cfada, xtrain, ytrain, scoring='roc_auc')
                cscores[c] = scores.mean()
            cmax = max(cscores, key=cscores.get)
            pval = cmax
            print cscores
        cfada = AdaBoostClassifier(C=pval)
        cfada.fit(xtrain, ytrain)
        return cfsvm, pval



