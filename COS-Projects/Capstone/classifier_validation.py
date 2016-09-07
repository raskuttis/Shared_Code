import numpy as np
import csv as csv
from collections import OrderedDict
from collections import Counter
from datetime import datetime, timedelta
from readcsv import *
from plot_analysis import *
from data_features import *
from classify import *

## List of questions to be answered in this code run
qnums = [0]

## Read in data stored in datadir
datadir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/PERM'
plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Figures/'
inddir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/Indicators/Country/'
## Read in external data to be used as predictor variables
years, countries, indicators, inddata = read_wb_indicator(inddir, 'WB_Development')
inddata = np.squeeze(inddata[0,:,:])
alldata = read_cfs(datadir, ['PERM_FY2008'])
cleandata = clean_data(alldata, 1)
#dec_month(cleandata, monthfeature, monthformat):
change_wage(cleandata, 'WAGE_OFFER_FROM_9089', 'WAGE_OFFER_UNIT_OF_PAY_9089')
#add_feature(cleandata, 'WAIT_TIME')

featurelabel = ['PW_LEVEL_9089', 'CLASS_OF_ADMISSION', 'US_ECONOMIC_SECTOR', 'EMPLOYER_STATE', 'WAGE_OFFER_UNIT_OF_PAY_9089','COUNTRY_OF_CITZENSHIP']
numfeatures = ['WAGE_CORRECTED']
featurelist = ['Pay Level', 'Visa', 'Sector', 'State', 'Pay Type', 'Country', 'Pay Amount']

for i in xrange(0,len(numfeatures)):
    feat = cleandata[numfeatures[i]]
    normalize_indicator(cleandata, numfeatures[i])
    numfeatures[i] = numfeatures[i] + '_NORM'

clist = ['NB', 'LR', 'RF', 'KNN']
clabel = ['gnb', 'lr', 'rf', 'knn']
#plot_cmpclassifiers(cleandata, featurelabel, clist, clabel, plotdir, arflag=1, fappend = 'Acceptance')
plot_cmpclassifiers(cleandata, featurelabel, clist, clabel, plotdir, arflag=0, pcaflag=50, fappend = 'PCA50')








