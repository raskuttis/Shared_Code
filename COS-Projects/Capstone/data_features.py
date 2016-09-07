import numpy as np
import csv as csv
from collections import OrderedDict
from collections import Counter
from datetime import datetime, timedelta, date

def weighted_pearson(x, y, w):

    meanx = np.sum(w * x) / np.sum(w)
    meany = np.sum(w * y) / np.sum(w)

    covxy = np.sum(w * (x - meanx) * (y - meany)) / np.sum(w)
    covxx = np.sum(w * (x - meanx) * (x - meanx)) / np.sum(w)
    covyy = np.sum(w * (y - meany) * (y - meany)) / np.sum(w)
    corrxy = covxy / np.sqrt(covxx * covyy)

    return corrxy

def change_wage(fulldata, wages, wagetypes):
    
    wagelist = np.asarray(fulldata[wages])
    badinds = np.where(wagelist == '')
    wagelist[badinds] = 0.0
    newwages = np.zeros(len(wagelist))
    wtlist = np.asarray(fulldata[wagetypes])
    goodinds = np.where(wtlist == 'HR')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 48.0 * 8.0
    goodinds = np.where(wtlist == 'YR')
    newwages[goodinds] = wagelist[goodinds].astype(float)
    goodinds = np.where(wtlist == 'WK')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 48.0
    goodinds = np.where(wtlist == 'MTH')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 12.0
    goodinds = np.where(wtlist == 'BI')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 12.0 * 2.0
    goodinds = np.where(wtlist == 'HOUR')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 48.0 * 8.0
    goodinds = np.where(wtlist == 'YEAR')
    newwages[goodinds] = wagelist[goodinds].astype(float)
    goodinds = np.where(wtlist == 'WEEK')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 48.0
    goodinds = np.where(wtlist == 'MONTH')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 12.0
    goodinds = np.where(wtlist == 'BI-WEEKLY')
    newwages[goodinds] = wagelist[goodinds].astype(float) * 12.0 * 2.0
    fulldata['WAGE_CORRECTED'] = newwages

def normalize_indicator(fulldata, selectfeature):

    feats = np.asarray(fulldata[selectfeature])
    badinds = np.where(feats == '')
    print len(feats), len(badinds)
    feats[badinds] = 0.0
    featurelist = feats.astype(float)
    mean = np.nanmean(featurelist)
    std = np.nanstd(featurelist)
    print mean, std
    newlist = (featurelist - mean) / std
    newlist[badinds] = 0.0
    fulldata[selectfeature + '_NORM'] = newlist

def add_indicators(fulldata, inddata, indicators, selectfeature, matchlist):

    featurelist = fulldata[selectfeature]
    featurelist = np.asarray([k.lower() for k in featurelist]).astype(str)
    ndps = len(featurelist)
    nnewvals = np.zeros((len(indicators), ndps))
    for j in xrange(0, len(matchlist)):
        match = matchlist[j].lower()
        ngood = np.shape(np.where(featurelist == match))[1]
        if ngood == 1:
            goodinds = np.where(featurelist == match)[0]
        else:
            goodinds = np.squeeze(np.where(featurelist == match))
        goodout = np.reshape(np.repeat(np.squeeze(inddata[j,:]),ngood),(len(indicators),ngood))
        nnewvals[:,goodinds] = goodout
        
    for i in xrange(0, len(indicators)):
        ind = indicators[i]
        fulldata[ind] = np.squeeze(nnewvals[i,:])

def match_indicators(fulldata, inddata, indicators, selectfeature, matchlist):
    
    featurelist = fulldata[selectfeature]
    featurelist = np.asarray([k.lower() for k in featurelist]).astype(str)
    fulldata[selectfeature + '_LOWER'] = featurelist
    matchlower = np.asarray([k.lower() for k in matchlist]).astype(str)
    featurefreqs = Counter(featurelist)
    nfeats = len(featurefreqs)
    indfeatdict = {}
    featurekeys = np.asarray(featurefreqs.keys())
    ars = acc_rate(fulldata, selectfeature + '_LOWER', 0)
    nanlist = np.empty(len(indicators))
    nanlist[:] = np.nan

    for j in xrange(0, nfeats):
        match = featurekeys[j]
        ngood = np.shape(np.where(matchlower == match))[1]
        goodind = np.where(matchlower == match)[0]
        if ngood == 1:
            indfeatdict[match] = np.squeeze(inddata[goodind,:])
        else:
            indfeatdict[match] = nanlist

    return indfeatdict, featurefreqs, ars


def acc_rate(alldata, cattitle, vthresh):
    
    ## Certified or denied status
    status = np.asarray(alldata['CASE_STATUS'])
    cinds = np.squeeze(np.where(status == '1'))
    dinds = np.squeeze(np.where(status == '0'))
    totinds = np.concatenate((cinds, dinds))
    ntot = len(totinds)

    ## Application number by categories within the column
    rawcol = np.asarray(alldata[cattitle])
    col = rawcol[totinds]
    colfreqs = Counter(col)
    colfreqs = {k: colfreqs[k] for (k,v) in colfreqs.items() if v > vthresh}
    
    ## Total acceptance rate by category within the column
    colaccs = rawcol[cinds]
    colaccfreqs = Counter(colaccs)
    ## Acceptance rate assuming from application category with more than vthresh applications
    colaccrates = {k: colaccfreqs[k] / float(v) for (k,v) in colfreqs.items() if v > vthresh}

    return colaccrates

def add_feature(fulldata, feature):

    origtime = datetime(2000, 1, 1)
    if (feature == 'START_TIME'):
        startdate = fulldata['CASE_RECEIVED_DATE']
        badinds = np.where(startdate == '')
        startdate[badinds] = '01/01/2001'
        sd = np.asarray([(datetime.strptime(a, '%m/%d/%y') - origtime).days for a in startdate])
        fulldata['START_TIME'] = sd
    if (feature == 'DEC_TIME'):
        startdate = fulldata['DECISION_DATE']
        badinds = np.where(startdate == '')
        startdate[badinds] = '01/01/2001'
        sd = np.asarray([(datetime.strptime(a, '%m/%d/%y') - origtime).days for a in startdate])
        fulldata['DEC_TIME'] = sd
    if (feature == 'DEC_TIME_1'):
        startdate = fulldata['DECISION_DATE']
        badinds = np.where(startdate == '')
        startdate[badinds] = '01-Jan-01'
        sd = np.asarray([(datetime.strptime(a, '%d-%b-%y') - origtime).days for a in startdate])
        fulldata['DEC_TIME'] = sd
    if (feature == 'WAIT_TIME'):
        add_feature(fulldata, 'START_TIME')
        add_feature(fulldata, 'DEC_TIME')
        fulldata['WAIT_TIME'] = fulldata['DEC_TIME'] - fulldata['START_TIME']

