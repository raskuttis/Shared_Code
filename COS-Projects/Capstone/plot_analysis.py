import numpy as np
import csv as csv
from collections import OrderedDict
from collections import Counter
from datetime import datetime, timedelta
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from classify import *
from data_features import *
from scipy.stats import pearsonr
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import r2_score

def classval(yout, xmat, nf, cf, reverse=False, pval=0.0, arflag=0.0, ncats=0.0):

    if reverse == True:
        yuse = -1 * yout + 1
    else:
        yuse = yout
    skf = StratifiedKFold(yuse, n_folds=nf, shuffle=True)
    
    iflag = 0
    tps = 0.0
    ps = 0.0
    fps = np.linspace(0, 1, 100)
    rs = np.linspace(0, 1, 100)
    aucs = []
    aups = []
    pv = pval
    for train_index, test_index in skf:
        y_test = yuse[test_index]
        x_test = xmat[test_index,:]
        y_train = yuse[train_index]
        x_train = xmat[train_index,:]
        if (arflag == 1):
            x_train = un_encode(y_train, x_train, ncats)
            x_test = un_encode(y_test, x_test, ncats)
        if (cf == 'rand'):
            y_predtest = np.random.rand(len(y_test))
        else:
            cfer, pv = classifier(y_train, x_train, cf, pval=pv)
            y_predtest = cfer.predict_proba(x_test)
            y_predtest = y_predtest[:,1]
        pre = average_precision_score(y_test, y_predtest)
        print pre
        p, r, th = precision_recall_curve(y_test, y_predtest)
        fp, tp, ths = roc_curve(y_test, y_predtest)
        auc = roc_auc_score(y_test, y_predtest)
        aup = average_precision_score(y_test, y_predtest)
        aucs.append(auc)
        aups.append(aup)
        ps = ps + np.interp(rs, np.flipud(r), np.flipud(p)) / float(nf)
        tps = tps + np.interp(fps, fp, tp) / float(nf)

    return ps, rs, tps, fps, aucs, aups

def plot_featimp(cleandata, featurelabel, featurelist, plotdir, reverse=False, arflag = 1, fappend = '', numfeatures = 0.0):

    yout, xmat, ncats, allcats = input_matrix(cleandata, featurelabel, 1)
    if (numfeatures != 0.0):
        xmat = add_numfeatures(cleandata, numfeatures, xmat)
    else:
        xmat = xmat.toarray()
    if reverse == True:
        yuse = -1 * yout + 1
    else:
        yuse = yout
    skf = StratifiedKFold(yuse, n_folds=10, shuffle=True)
    allcats = [v for sl in allcats for v in sl]

    fimps = []
    pimps = []
    for train_index, test_index in skf:
        y_test = yuse[test_index]
        x_test = xmat[test_index,:]
        y_train = yuse[train_index]
        x_train = xmat[train_index,:]
        if (arflag == 1):
            x_train = un_encode(y_train, x_train, ncats)
            x_test = un_encode(y_test, x_test, ncats)
        nfs = np.shape(x_train)[1]
        cfer, pv = classifier(y_train, x_train, 'rf')
        y_predtest = cfer.predict_proba(x_test)[:,1]
        acc = roc_auc_score(y_test, cfer.predict_proba(x_test)[:,1])
        pre = average_precision_score(y_test, cfer.predict_proba(x_test)[:,1])
        print acc, pre
        fimp = np.zeros(nfs)
        pimp = np.zeros(nfs)
        for i in xrange(0,nfs):
            x_t = x_test.copy()
            np.random.shuffle(x_t[:,i])
            shuff_acc = roc_auc_score(y_test, cfer.predict_proba(x_t)[:,1])
            shuff_pre = average_precision_score(y_test, cfer.predict_proba(x_t)[:,1])
            fimp[i] = (acc-shuff_acc)/acc
            pimp[i] = (pre - shuff_pre) / pre
        print fimp
        fimps.append(fimp)
        pimps.append(pimp)

    print np.shape(fimps)

    if (arflag == 1):
        cats = np.asarray(featurelist)
    else:
        cats = np.asarray(allcats)
    catnums = np.nanmean(fimps,axis=0)
    catnums = catnums / np.sum(catnums)
    caterrs = np.nanstd(fimps,axis=0)
    catpres = np.nanmean(pimps,axis=0)
    catpres = catpres / np.sum(catpres)
    catperrs = np.nanstd(pimps,axis=0)
    catsortinds = np.argsort(-1.0 * catnums)
    cats = cats[catsortinds]
    print cats[0:50]
    catnums = catnums[catsortinds]
    print catnums[0:50]
    caterrs = caterrs[catsortinds]
    catpres = catpres[catsortinds]
    catperrs = catperrs[catsortinds]
    if (arflag == 0):
        cats = cats[0:12]
        catnums = catnums[0:12]
        caterrs = caterrs[0:12]
        catpres = catnums[0:12]
        catperrs = caterrs[0:12]
    y_pos = np.arange(len(cats))

    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '10'
    plt.figure(figsize = [xsize,ysize])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)
    
    ax = plt.subplot(1,1,1)
    plt.subplots_adjust(bottom=0.2)
    #plt.subplots_adjust(left=0.25)
    
    plt.barh(y_pos, catnums, align='center', alpha=0.4, color='b', yerr=caterrs)
    plt.yticks(y_pos, y_pos+1)
    plt.axis([0.0,1.1*np.max(catnums),-1,len(cats)+1])
    #plt.xticks([0.0,0.25,0.5])
    ax.set_xlabel('AUC Importance', color='b')
    
    pp = PdfPages(plotdir + '/FeatureImportance' + fappend + '.pdf')
    pp.savefig()
    pp.close()

def plot_cmpwb(cleandata, plotdir):
    
    ps = []
    rs = []
    tps = []
    fps = []
    aucs = []
    aups = []
    pcalist = [5,15,50]
    ndps = len(cleandata['CASE_STATUS'])
    
    featurelabel = ['PW_LEVEL_9089', 'CLASS_OF_ADMISSION', 'US_ECONOMIC_SECTOR', 'EMPLOYER_STATE', 'WAGE_OFFER_UNIT_OF_PAY_9089','COUNTRY_OF_CITZENSHIP']
    numfeatures = ['WAGE_CORRECTED_NORM']
    
    arlist = [0,1]
    for i in xrange(0, len(arlist)):
        yout, xmat, ncats, allcats = input_matrix(cleandata, featurelabel, 1)
        if (numfeatures != 0.0):
            xmat = add_numfeatures(cleandata, numfeatures, xmat)
        else:
            xmat = xmat.toarray()
        p, r, tp, fp, auc, aup = classval(yout, xmat, 10, 'rf', reverse=True, pval=0.0, arflag=arlist[i], ncats = ncats)
        ps.append(p)
        rs.append(r)
        tps.append(tp)
        fps.append(fp)
        aucs.append(auc)
        aups.append(aup)
    
    featurelabel = ['PW_LEVEL_9089', 'CLASS_OF_ADMISSION', 'US_ECONOMIC_SECTOR', 'EMPLOYER_STATE', 'WAGE_OFFER_UNIT_OF_PAY_9089','COUNTRY_OF_CITZENSHIP']
    allnumfeatures = ['WAGE_CORRECTED_NORM', 'Intentional homicides (per 100,000 people)', 'Merchandise exports to low- and middle-income economies in South Asia (% of total merchandise exports)', 'Merchandise imports from low- and middle-income economies in Middle East & North Africa (% of total merchandise imports)', 'GINI index (World Bank estimate)', 'Merchandise imports from economies in the Arab World (% of total merchandise imports)', 'Import value index (2000 = 100)', 'Income share held by highest 10%', 'Fuel exports (% of merchandise exports)', 'Merchandise exports to low- and middle-income economies outside region (% of total merchandise exports)', 'Food production index (2004-2006 = 100)', 'Adjusted savings: gross savings (% of GNI)', 'Merchandise exports to low- and middle-income economies within region (% of total merchandise exports)', 'Merchandise exports to low- and middle-income economies in Sub-Saharan Africa (% of total merchandise exports)', 'Ores and metals exports (% of merchandise exports)', 'Adjusted savings: education expenditure (% of GNI)', 'Import volume index (2000 = 100)', 'Real interest rate (%)', 'Income share held by third 20%', 'Adjusted net national income (annual % growth)', 'Income share held by fourth 20%', 'Merchandise exports to high-income economies (% of total merchandise exports)', 'Youth literacy rate, population 15-24 years, both sexes (%)', 'Income share held by highest 20%', 'Merchandise exports to low- and middle-income economies in Latin America & the Caribbean (% of total merchandise exports)', 'Inflation, GDP deflator (annual %)', 'Export value index (2000 = 100)', 'Food exports (% of merchandise exports)', 'Primary completion rate, both sexes (%)', 'Manufactures exports (% of merchandise exports)', 'Rural population growth (annual %)', 'Food imports (% of merchandise imports)', 'Manufactures imports (% of merchandise imports)', 'Ores and metals imports (% of merchandise imports)', 'Merchandise imports from low- and middle-income economies within region (% of total merchandise imports)', 'Immunization, measles (% of children ages 12-23 months)', 'Real effective exchange rate index (2010 = 100)', 'Official exchange rate (LCU per US$, period average)', 'Merchandise exports to low- and middle-income economies in Middle East & North Africa (% of total merchandise exports)', 'Birth rate, crude (per 1,000 people)', 'Adjusted net savings, including particulate emission damage (% of GNI)', 'Research and development expenditure (% of GDP)', 'Consumer price index (2010 = 100)' 'Urban population growth (annual %)', 'Income share held by second 20%', 'Merchandise imports from low- and middle-income economies outside region (% of total merchandise imports)', 'Merchandise imports from low- and middle-income economies in East Asia & Pacific (% of total merchandise imports)', 'Adjusted net national income per capita (annual % growth)', 'International tourism, expenditures (% of total imports)', 'Merchandise imports from low- and middle-income economies in Sub-Saharan Africa (% of total merchandise imports)', 'Interest rate spread (lending rate minus deposit rate, %)']
    
    for i in xrange(0, len(pcalist)):
        numfeatures = allnumfeatures[0:pcalist[i]]
        if (numfeatures != 0.0):
            nfeats = len(numfeatures)
            xmat = np.zeros((ndps,nfeats))
            for j in xrange(0, nfeats):
                feat = numfeatures[j]
                feats = np.asarray(cleandata[feat]).astype(float)
                xmat[:,j] = feats
        yout = np.asarray(cleandata['CASE_STATUS']).astype(int)
        p, r, tp, fp, auc, aup = classval(yout, xmat, 10, 'rf', reverse=True, pval=0.0)
        print auc, aup
        ps.append(p)
        rs.append(r)
        tps.append(tp)
        fps.append(fp)
        aucs.append(auc)
        aups.append(aup)
    
    leglist = ['One Hot', 'Acceptance %', 'World Bank-5', 'World Bank-15', 'World Bank-50']
    plot_classareas(aucs, aups, leglist, plotdir, 'WB')
    plot_classval(ps, rs, tps, fps, leglist, plotdir, 'WB')

def plot_cmpar(cleandata, featurelabel, clabel, plotdir, numfeatures = 0.0):
    
    ps = []
    rs = []
    tps = []
    fps = []
    aucs = []
    aups = []
    arlist = [0,1]
    for i in xrange(0, len(arlist)):
        yout, xmat, ncats, allcats = input_matrix(cleandata, featurelabel, 1)
        if (numfeatures != 0.0):
            xmat = add_numfeatures(cleandata, numfeatures, xmat)
        else:
            xmat = xmat.toarray()
        p, r, tp, fp, auc, aup = classval(yout, xmat, 10, clabel[i], reverse=True, pval=0.0, arflag=arlist[i], ncats = ncats)
        ps.append(p)
        rs.append(r)
        tps.append(tp)
        fps.append(fp)
        aucs.append(auc)
        aups.append(aup)

    pcalist = [1,10,30,50,100]
    for i in xrange(0, len(pcalist)):
        yout, xmat, ncats, allcats = input_matrix(cleandata, featurelabel, 1, pcaflag = pcalist[i])
        if (numfeatures != 0.0):
            nfeats = len(numfeatures)
            for j in xrange(0, nfeats):
                feat = numfeatures[j]
                feats = np.asarray(cleandata[feat]).astype(float)
                ndps = np.shape(xmat)[0]
                nexfs = np.shape(xmat)[1]
                xmat = np.insert(xmat, nexfs, feats, axis=1)

        p, r, tp, fp, auc, aup = classval(yout, xmat, 10, clabel[0], reverse=True, pval=0.0)
        ps.append(p)
        rs.append(r)
        tps.append(tp)
        fps.append(fp)
        aucs.append(auc)
        aups.append(aup)

    leglist = ['One Hot', 'Acceptance %', 'PCA-1', 'PCA-10', 'PCA-30', 'PCA-50', 'PCA-100']
    plot_classareas(aucs, aups, leglist, plotdir, 'AR')
    plot_classval(ps, rs, tps, fps, leglist, plotdir, 'AR')

def plot_cmpdim(cleandata, pcalist, featurelabel, clabel, plotdir):
    
    ps = []
    rs = []
    tps = []
    fps = []
    aucs = []
    aups = []
    for i in xrange(0, len(pcalist)):
        yout, xmat, ncats, allcats = input_matrix(cleandata, featurelabel, 1, pcaflag = pcalist[i])
        
        p, r, tp, fp, auc, aup = classval(yout, xmat, 10, clabel, reverse=True, pval=0.0)
        ps.append(p)
        rs.append(r)
        tps.append(tp)
        fps.append(fp)
        aucs.append(auc)
        aups.append(aup)
    
    plot_classareas(aucs, aups, pcalist, plotdir, 'PCA')
    plot_classval(ps, rs, tps, fps, pcalist, plotdir, 'PCA')

def plot_cmpfeatures(cleandata, featurelist, featurelabel, clabel, plotdir, arflag=0.0, fappend = ''):

    ps = []
    rs = []
    tps = []
    fps = []
    aucs = []
    aups = []
    for i in xrange(0, len(featurelist)):
        feat = featurelabel[i]
        yout, xmat, ncats, allcats = input_matrix(cleandata, feat, 1)
        p, r, tp, fp, auc, aup = classval(yout, xmat, 10, clabel, reverse=True, arflag=arflag, ncats = ncats)
        ps.append(p)
        rs.append(r)
        tps.append(tp)
        fps.append(fp)
        aucs.append(auc)
        aups.append(aup)

    plot_classareas(aucs, aups, featurelist, plotdir, 'Feature' + fappend)
    plot_classval(ps, rs, tps, fps, featurelist, plotdir, 'Feature' + fappend)


def plot_cmpclassifiers(cleandata, featurelist, clist, clabel, plotdir, arflag=0.0, pcaflag = 0.0, fappend = '', numfeatures = 0.0):

    yout, xmat, ncats, allcats = input_matrix(cleandata, featurelist, 1, pcaflag = pcaflag)
    if (numfeatures != 0.0):
        nfeats = len(numfeatures)
        if pcaflag == 0.0:
            xmat = xmat.toarray()
        for j in xrange(0, nfeats):
            feat = numfeatures[j]
            feats = np.asarray(cleandata[feat]).astype(float)
            ndps = np.shape(xmat)[0]
            nexfs = np.shape(xmat)[1]
            xmat = np.insert(xmat, nexfs, feats, axis=1)

    ps = []
    rs = []
    tps = []
    fps = []
    aucs = []
    aups = []
    for i in xrange(0,len(clist)):
        print clabel[i]
        p, r, tp, fp, auc, aup = classval(yout, xmat, 10, clabel[i], reverse=True, arflag=arflag, ncats = ncats)
        ps.append(p)
        rs.append(r)
        tps.append(tp)
        fps.append(fp)
        aucs.append(auc)
        aups.append(aup)

    plot_classareas(aucs, aups, clist, plotdir, 'Classifier' + fappend)
    plot_classval(ps, rs, tps, fps, clist, plotdir, 'Classifier' + fappend)

def plot_classareas(aucs, aups, leglist, plotdir, fname):

    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '6'
    plt.figure(figsize = [xsize,ysize])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)
    
    ax = plt.subplot(1,1,1)
    plt.subplots_adjust(bottom=0.2)
    plt.subplots_adjust(right=0.85)
    cats = leglist
    catnums = np.nanmean(aucs,axis=1)
    catrates = np.nanmean(aups,axis=1)
    catnumerrs = np.nanstd(aucs,axis=1)
    catrateerrs = np.nanstd(aups,axis=1)
    y_pos = np.arange(len(cats))
    
    width = 0.35
    bp = plt.boxplot(aucs,positions=y_pos,widths=0.5*width,patch_artist=True,notch=True)
    plt.setp(bp['boxes'], color='b')
    plt.setp(bp['whiskers'], color='b')
    plt.setp(bp['medians'], color='b')
    plt.setp(bp['fliers'], marker='None')
    plt.xticks(y_pos, cats, rotation=90)
    #plt.axis([-1,len(cats)+1,0.0,1.0])
    #plt.yticks([0.0,0.5,1.0])
    ax.set_ylabel('AUC', color='b')
    
    axr = ax.twinx()
    bp = axr.boxplot(aups,positions=y_pos+width,widths=0.5*width,patch_artist=True,notch=True)
    plt.setp(bp['boxes'], color='r')
    plt.setp(bp['whiskers'], color='r')
    plt.setp(bp['medians'], color='r')
    plt.setp(bp['fliers'], marker='None')
    #plt.axis([-1,len(cats)+1,0.0,1.0])
    plt.xticks(y_pos, cats, rotation=90)
    #plt.yticks([0.0,0.5,1.0])
    axr.set_ylabel('Precision Score', color='r')
    
    pp = PdfPages(plotdir + '/' + fname + '_AUCs.pdf')
    pp.savefig()
    pp.close()

def plot_classval(ps, rs, tps, fps, leglist, plotdir, fname):

    nreps = len(leglist)
    cols = ['k', 'r', 'b', 'g', 'y', 'c', 'm']

    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '10'
    plt.figure(figsize = [2*xsize,ysize])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)

    ax = plt.subplot(1,2,1)
    plt.subplots_adjust(bottom=0.2)
    for i in xrange(0,nreps):
        ax.plot(rs[i],ps[i],cols[i],label=leglist[i])
    ax.legend(prop={'size':6})
    plt.text(0.05,0.9,r"$\displaystyle(a)$")
    plt.axis([0.0,1.0,0.0,1.0])
    plt.yticks([0.0,0.5,1.0])
    plt.xticks([0.0,0.5,1.0])
    plt.xlabel(r"$\displaystyle {\rm Recall}$")
    plt.ylabel(r"$\displaystyle {\rm Precision}$")

    ax2 = plt.subplot(1,2,2)
    ax2.yaxis.set_label_position("right")
    for i in xrange(0,nreps):
        plt.plot(fps[i],tps[i],cols[i])
    plt.text(0.05,0.9,r"$\displaystyle(b)$")
    plt.axis([0.0,1.0,0.0,1.0])
    plt.yticks([0.0,0.5,1.0],[' ',' ',' '])
    plt.xticks([0.0,0.5,1.0])
    plt.xlabel(r"$\displaystyle {\rm FPR}$")
    plt.ylabel(r"$\displaystyle {\rm TPR}$")

    pp = PdfPages(plotdir + '/' + fname + '_Validation.pdf')
    pp.savefig()
    pp.close()

def plot_corrinds(cleandata, inddata, indicators, colname, colvals, plotdir):

    col_inds, npercol, col_accrate = match_indicators(cleandata, inddata, indicators, colname, colvals)
    ncols = len(colvals)
    nindicators = len(indicators)

    indmat = np.zeros((ncols,nindicators))
    indacc = np.zeros(ncols)
    indfreq = np.zeros(ncols)
    i = 0
    for (k,v) in col_inds.items():
        indfreq[i] = npercol[k]
        indacc[i] = col_accrate[k]
        indmat[i,:] = v
        i = i + 1

    pcf = []
    wpcf = []
    for j in xrange(0,nindicators):
        indvals = np.squeeze(indmat[:,j])
        goodinds = np.isfinite(indvals)
        ivgood = indvals[goodinds]
        ivacc = indacc[goodinds]
        ivfreq = indfreq[goodinds]
        indcorr = pearsonr(ivgood, ivacc)
        pc = weighted_pearson(ivgood, ivacc, ivfreq)
        liminds = np.squeeze(np.where(ivfreq > 40))
        indlimcorr = pearsonr(ivgood[liminds], ivacc[liminds])
        pcf.append(indcorr[0])
        wpcf.append(pc)

    pcf = np.asarray(pcf)
    wpcf = np.asarray(wpcf)

    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '6'
    plt.figure(figsize = [xsize,ysize])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)

    nmax = 12
    ax = plt.subplot(1,1,1)
    plt.subplots_adjust(left=0.2)
    plt.subplots_adjust(right=0.85)
    plt.subplots_adjust(bottom=0.2)
    colsortinds = np.argsort(-1.0 * np.abs(wpcf))
    cats = indicators[colsortinds[0:nmax]]
    catnums = pcf[colsortinds[0:nmax]]
    catrates = wpcf[colsortinds[0:nmax]]
    colsortinds = np.argsort(catrates)
    cats = cats[colsortinds]
    catnums = catnums[colsortinds]
    catrates = catrates[colsortinds]
    print cats
    cats = np.asarray([k[0:6] for k in cats])
    y_pos = np.arange(len(cats))
    
    width = 0.35
    plt.bar(y_pos, catnums, width, align='center', alpha=0.4, color='b')
    plt.xticks(y_pos, y_pos+1)
    plt.axis([-1,len(cats)+1,-1.0,1.0])
    #plt.yticks([0.0,0.25,0.5])
    ax.set_ylabel('Pearson', color='b')
    
    axr = ax.twinx()
    axr.bar(y_pos+width, catrates, width, align='center', alpha=0.4, color='r')
    plt.axis([-1,len(cats)+1,-1.0,1.0])
    plt.yticks([0.0,0.5,1.0])
    axr.set_ylabel('Weighted Pearson', color='r')
    
    pp = PdfPages(plotdir + '/' + colname + '_IndicatorCorrelations.pdf')
    pp.savefig()
    pp.close()


def plot_certrates(alldata, cattitle, vthresh, plotdir):

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
    sortcolaccrates = sorted(colaccrates.items(), key=lambda x:x[1])

    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '6'
    plt.figure(figsize = [xsize,ysize])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)

    ax = plt.subplot(1,1,1)
    plt.subplots_adjust(bottom=0.2)
    colsortinds = np.argsort(colfreqs.values())
    cats = colfreqs.keys()
    cats = np.asarray([k[0:3] for k in cats])
    cats = cats[colsortinds]
    catnums = np.asarray(colfreqs.values()).astype(float)
    catnums = catnums[colsortinds] / float(ntot)
    catrates = np.asarray(colaccrates.values()).astype(float)
    catrates = catrates[colsortinds]
    y_pos = np.arange(len(cats))

    width = 0.35
    plt.bar(y_pos, catnums, width, align='center', alpha=0.4, color='b')
    plt.xticks(y_pos, cats, rotation=90)
    plt.axis([-1,len(cats)+1,0.0,1.1*np.max(catnums)])
    #plt.yticks([0.0,0.25,0.5])
    ax.set_ylabel('Application %', color='b')
    
    axr = ax.twinx()
    axr.bar(y_pos+width, catrates, width, align='center', alpha=0.4, color='r')
    plt.axis([-1,len(cats)+1,0.0,1.0])
    plt.yticks([0.0,0.5,1.0])
    axr.set_ylabel('Certification Rate', color='r')

    pp = PdfPages(plotdir + '/' + cattitle + '_CertRates.pdf')
    pp.savefig()
    pp.close()

def plot_violin(alldata, cattitle, plotdir):
    
    ## Certified or denied status
    status = np.asarray(alldata['CASE_STATUS'])
    cinds = np.squeeze(np.where(status == '1'))
    dinds = np.squeeze(np.where(status == '0'))
    totinds = np.concatenate((cinds, dinds))
    ntot = len(totinds)
    
    ## Application number by categories within the column
    rawcol = np.asarray(alldata[cattitle])
    col = rawcol[totinds]
    colaccs = np.asarray(rawcol[cinds]).astype(float)
    coldens = np.asarray(rawcol[dinds]).astype(float)
    
    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '6'
    plt.figure(figsize = [xsize,ysize])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)
    
    ax = plt.subplot(1,1,1)
    plt.violinplot([colaccs,coldens],showmeans=False,showmedians=True)
    #plt.axis([0,3,0.0,1000.0])
    
    pp = PdfPages(plotdir + '/' + cattitle + '_ViolinPlot.pdf')
    pp.savefig()
    pp.close()

