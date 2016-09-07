import numpy as np
import csv as csv
from collections import OrderedDict
from collections import Counter
from datetime import datetime, timedelta

def read_any(fullpath):

    fulldata = {}
    with open(fullpath) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            for column, value in row.items():
                fulldata.setdefault(column, []).append(value)

    return fulldata


## Function to read in csv files for a given input set of years
def read_cfs(datadir, years):

    fulldata = {}

    nyears = len(years)
    for i in xrange(0,nyears):
        year = years[i]
        print year

        datafile = '{0:s}.csv'.format(year)
        fullpath = '{0:s}/{1:s}'.format(datadir,datafile)
        keys = []
        with open(fullpath) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for column, value in row.items():
                    keys.append((column.upper()).replace(' ', '_'))
                    fulldata.setdefault((column.upper()).replace(' ', '_'), []).append(value.upper())
                fulldata.setdefault('YEAR', []).append(year)
        keys = list(set(keys))
        if (i == 0):
            intkeys = keys
        else:
            intkeys = list(set(intkeys) & set(keys))

    newdata = {}
    newdata = {k: v for (k,v) in fulldata.items() if k in intkeys}

    return newdata

def clean_data(fulldata, uflag):

    ## Keep only dictionary values
    status = np.asarray(fulldata['CASE_STATUS'])
    print np.shape(fulldata['CASE_STATUS'])
    if (uflag == 1):
        cinds = np.squeeze(np.where(status == 'CERTIFIED'))
        ceinds = np.squeeze(np.where(status == 'CERTIFIED-EXPIRED'))
        cinds = np.concatenate((cinds, ceinds))
        dinds = np.squeeze(np.where(status == 'DENIED'))
    else:
        cinds = np.squeeze(np.where(status == 'Certified'))
        dinds = np.squeeze(np.where(status == 'Denied'))
    ## Ignore those that have been withdrawn
    totinds = np.concatenate((cinds, dinds))
    ntot = len(totinds)
    accrate = float(len(cinds)) / ntot
    status[cinds] = 1
    status[dinds] = 0
    status = np.asarray(status)[totinds]
    print accrate, ntot

    alldata = {k: np.asarray(v)[totinds] for (k,v) in fulldata.items()}
    alldata['CASE_STATUS'] = status

    return alldata

def read_wb_indicator(datadir, fname):

    fulldata = {}

    datafile = '{0:s}.csv'.format(fname)
    fullpath = '{0:s}/{1:s}'.format(datadir,datafile)

    with open(fullpath) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            for column, value in row.items():
                fulldata.setdefault(column, []).append(value)

    ## Convert data dictionary to matrix of indicators
    inds = np.asarray(fulldata['Series Name'])
    goodinds = np.squeeze(np.where(inds != ''))
    fulldata = {k: np.asarray(v)[goodinds] for (k,v) in fulldata.items()}
    inds = np.asarray(fulldata['Series Name'])
    countries = np.asarray(fulldata['\xef\xbb\xbfCountry Name'])
    indfreqs = Counter(inds)
    countryfreqs = Counter(countries)
    ninds = len(indfreqs)
    ncountries = indfreqs.values()
    ncountries = ncountries[0]
    nyears = len(fulldata.keys()) - 4

    inddata = np.zeros((nyears, ncountries, ninds))
    indmean = np.zeros((nyears, ncountries, ninds))
    indstd = np.zeros((nyears, ncountries, ninds))
    years = np.sort(2014 - np.arange(nyears)).astype(str)
    indicators = inds[0:ninds]
    countries = np.unique(countries)

    ## Add indices to matrix and convert them to zero mean, standard deviation 1
    for i in xrange(0, nyears):
        key = '{0:s} [YR{0:s}]'.format(years[i])
        drow = np.asarray(fulldata[key])
        badinds = np.squeeze(np.where(drow == '..'))
        drow[badinds] = np.nan
        inddata[i,:,:] = np.reshape(drow, (ncountries,ninds))

    for i in xrange(0, ninds):
        indsub = np.squeeze(inddata[:,:,i])
        mean = np.nanmean(indsub)
        std = np.nanstd(indsub)
        inddata[:,:,i] = (indsub - mean) / std

    badinds = np.isnan(inddata)
    inddata[badinds] = 0.0

    return years, countries, indicators, inddata

