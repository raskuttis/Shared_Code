import matplotlib.pyplot as plt
import sys
import random
import scipy.io as spio
import numpy as np
from os import listdir
from matplotlib.backends.backend_pdf import PdfPages

def loadall(datadir):
    
    subdirs = listdir(datadir)
    nsubs = len(subdirs)
    
    datarr = []
    filearr = []
    labels = []
    
    for i in xrange(0,len(subdirs)):
        subdir = subdirs[i]
        if subdir[0] != '.':
            subpath = datadir + subdir + '/'
            fnames = listdir(subpath)
            for j in xrange(0,len(fnames)):
                fname = fnames[j]
                if fname[0] != '.':
                    fileN = subpath + fname
                    x = spio.loadmat(fileN)
                    datarr.append(x['DAT'])
                    filearr.append(fname)
                    lb = x['DAT']['class'][0][0][0][0]
                    labels.append(lb)

    return datarr, filearr, labels

datadir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/voxResources/data/'

d,f,lbs = loadall(datadir)

