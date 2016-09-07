import numpy as np
import csv as csv
from collections import OrderedDict
from collections import Counter
from datetime import datetime, timedelta
from readcsv import *
from plot_analysis import *
from data_features import *

## Read in data stored in datadir
datadir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/PERM'
plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Figures/'
inddir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/Indicators/Country/'
## Read in external data to be used as predictor variables
years, countries, indicators, inddata = read_wb_indicator(inddir, 'WB_Development')
ncountries = len(countries)
nindicators = len(indicators)
years = np.asarray(['PERM_FY' + v for v in years])
yind = np.squeeze(np.where(years == 'PERM_FY2008'))
inddata = np.squeeze(inddata[yind,:,:])
alldata = read_cfs(datadir, ['PERM_FY2008'])
cleandata = clean_data(alldata, 1)
change_wage(cleandata, 'WAGE_OFFER_FROM_9089', 'WAGE_OFFER_UNIT_OF_PAY_9089')
normalize_indicator(cleandata, 'WAGE_CORRECTED')

add_indicators(cleandata, inddata, indicators, 'COUNTRY_OF_CITZENSHIP', countries)

plot_cmpwb(cleandata, plotdir)
exit()

featurelabel = ['APPLICATION_TYPE']
featurelist = ['App Type']
featurelist = np.append(featurelist, indicators)
clabel = ['rf', 'rf']
numfeatures = indicators
plot_featimp(cleandata, featurelabel, featurelist, plotdir, arflag = 1, numfeatures = numfeatures)
#plot_cmpwb(cleandata, featurelabel, clabel, plotdir, numfeatures = indicators)
#add_feature(cleandata, 'WAIT_TIME')

#plot_corrinds(cleandata, inddata, indicators, 'COUNTRY_OF_CITZENSHIP', countries, plotdir)

