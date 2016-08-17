import numpy as np
import sys
import random
import scipy.io as spio
import scipy.stats as stats
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.metrics import zero_one_loss
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import confusion_matrix

def class_error(feat, labels, classname, nn=1, gamma=2, crbf=1):

    # Define number of songs and genres
    N = feat.shape[1]
    genres = ["Blues", "Classical", "Country", "Disco", "Hiphop", "Jazz", "Metal", "Pop", "Reggae", "Rock"]
    ngenres = len(genres)

    # Define test and training vector sets using k-fold cross validation
    k = 10
    ## Check randomness of stratification
    skf = StratifiedKFold(labels,n_folds=k,shuffle=True)
    averageCfmat = [[0 for x in range(ngenres)] for x in range(ngenres)]
    allerrors = []

    for train_index, test_index in skf:
        X_train, X_test = feat[:,train_index], feat[:,test_index]
        y_train, y_test = labels[train_index], labels[test_index]
    
        # Run classifier defined by input classname
        if classname == 'RF':
            cf = RandomForestClassifier(n_estimators=100, max_depth = 10, warm_start = False)
        elif classname == 'KNN':
            cf = KNeighborsClassifier(n_neighbors=nn)
        elif classname == 'GNB':
            cf = GaussianNB()
        elif classname == 'LSVM':
            cf = SVC(kernel="linear", C=0.025)
        elif classname == 'RBFSVM':
            cf = SVC(gamma=gamma, C=crbf)
        elif classname == 'ADAB':
            cf = AdaBoostClassifier(n_estimators=100)
        else:
            print 'Not a valid classifier'
    
        cf.fit(X_train.T,y_train)
        y_pred = cf.predict(X_test.T)
        error = zero_one_loss(y_pred,y_test)
        allerrors.append(error)
        cfmat = confusion_matrix(y_test, y_pred, labels=range(1,ngenres+1))
        averageCfmat += cfmat

    avError = np.mean(allerrors)
    stdError = np.std(allerrors)
    return [avError, stdError, averageCfmat]
