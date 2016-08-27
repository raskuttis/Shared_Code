import sys
import numpy as np
from hyp_math import *
from scipy.optimize import curve_fit

def sim_pdfchars(pdflines):
    line = pdflines[1]
    pdfn = int(eval(line.split()[3]))
    pdfmin = eval(line.split()[1])
    pdfdel = eval(line.split()[2])
    pdfmax = pdfmin + (pdfn - 1) * pdfdel
    return pdfmin,pdfdel,pdfmax,pdfn

def sim_pdfchars_zip(pdflines,index):
    line = pdflines[2*index+1]
    pdfn = int(eval(line.split()[3]))
    pdfmin = eval(line.split()[1])
    pdfdel = eval(line.split()[2])
    pdfmax = pdfmin + (pdfn - 1) * pdfdel
    pdfnmin = int(eval(line.split()[4]))
    delpdfn = int(eval(line.split()[5]))
    return pdfmin,pdfdel,pdfmax,pdfn,pdfnmin,delpdfn

def sim_pdftime(pdflines, tff):
    nlines = len(pdflines)
    nts = int(np.floor(nlines / 2))
    pdftime = []
    for i in xrange(1,nts):
        line = pdflines[2*i-1]
        lsplit = line.split()
        pdftime.append(eval(lsplit[0]) / tff)
    return pdftime

def sim_pdfx(pdflines):
    pdfmin, pdfdel, pdfmax, pdfn = sim_pdfchars(pdflines)
    pdfx = np.linspace(pdfmin,pdfmax,num=pdfn)
    return pdfx

def sim_pdf(pdflines,index):
    line = pdflines[2*index+2]
    pdf = np.array(map(float, line.split()))
    return pdf

def sim_pdf_zip(pdflines,index):
    line = pdflines[2*index+2]
    pdf = np.array(map(float, line.split()))
    pdfmin,pdfdel,pdfmax,pdfn,pdfnmin,delpdfn=sim_pdfchars_zip(pdflines,index)
    fullpdf = np.zeros(pdfn)
    fullpdf[pdfnmin:pdfnmin+delpdfn] = pdf
    return fullpdf

def sim_pdfm(pdflines,index):
    pdflogx = sim_pdfx(pdflines)
    pdfx = 10**pdflogx
    pdf = sim_pdf(pdflines,index)
    pdfm = pdf * pdfx
    pdfm = pdfm / np.sum(pdfm)
    return pdfm

def sim_pdf_mass(pdflines,index,atot,mconv):
    pdflogx = sim_pdfx(pdflines)
    pdfx = 10**pdflogx
    pdf = sim_pdf(pdflines,index)
    mass = np.sum(pdf * pdfx) * atot / np.sum(pdf)
    return mass / mconv

def sim_pdf_mean(pdf,pdflogx):
    pdfx = 10**pdflogx
    mean = np.sum(pdf * pdfx) / np.sum(pdf)
    meansq = np.sum(pdf * pdfx * pdfx) / np.sum(pdf)
    meanm = meansq / mean
    return mean, meanm

def sim_fitlnpdfvals(pdfin,pdflogx,pdfm,massmin,massmax,nvals,pin=None):
    pdf = pdfin / np.sum(pdfin)
    pdfcum = np.cumsum(pdfm)
    ilow = np.abs(pdfcum - massmin).argmin()
    ihigh = np.abs(pdfcum - massmax).argmin()
    muin = np.sum(pdflogx[ilow:ihigh] * pdf[ilow:ihigh]) / np.sum(pdf[ilow:ihigh])
    sigmain = np.sum(pdflogx[ilow:ihigh] * pdflogx[ilow:ihigh] * pdf[ilow:ihigh]) / np.sum(pdf[ilow:ihigh])
    sigmain = np.sqrt(sigmain - muin * muin)
    ampin = max(pdf[ilow:ihigh]) * np.sqrt(2.0 * np.pi) * sigmain
    if pin:
        [ampin,muin,sigmain] = pin
    try:
        fitvals, fitcov = curve_fit(gaussian, pdflogx[ilow:ihigh], pdf[ilow:ihigh], p0 = [ampin, muin, sigmain], sigma = np.sqrt(pdf[ilow:ihigh]/nvals))
        fita, fitmu, fitsigma = fitvals
        fitpdf = gaussian(pdflogx,fita,fitmu,fitsigma)
        fitchi = np.sum((pdf[ilow:ihigh] - fitpdf[ilow:ihigh])**2 / fitpdf[ilow:ihigh]) * nvals / (ihigh - ilow)
    except RuntimeError:
        print "Error - Curvefit Failed"
        fita = 0.0
        fitmu = 0.0
        fitsigma = 0.0
        fitchi = 0.0

    return fita, fitmu, fitsigma, fitchi

def sim_fitlnpdf(pdf,pdflogx,pdfm,massmin,massmax,nvals,pin=None):
    fitvals = sim_fitlnpdfvals(pdf,pdflogx,pdfm,massmin,massmax,nvals,pin=pin)
    fita, fitmu, fitsigma, fitchi = fitvals
    print fitsigma, fitmu
    fitpdf = gaussian(pdflogx,fita,fitmu,fitsigma)
    return fitpdf

def sim_fitlnhermpdfvals(pdfin,pdflogx,pdfm,massmin,massmax,nvals,pin=None):
    pdf = pdfin / np.sum(pdfin)
    pdfcum = np.cumsum(pdfm)
    ilow = np.abs(pdfcum - massmin).argmin()
    ihigh = np.abs(pdfcum - massmax).argmin()
    muin = np.sum(pdflogx[ilow:ihigh] * pdf[ilow:ihigh]) / np.sum(pdf[ilow:ihigh])
    sigmain = np.sum(pdflogx[ilow:ihigh] * pdflogx[ilow:ihigh] * pdf[ilow:ihigh]) / np.sum(pdf[ilow:ihigh])
    sigmain = np.sqrt(sigmain - muin * muin)
    ampin = max(pdf[ilow:ihigh]) * np.sqrt(2.0 * np.pi) * sigmain
    h3in = 0.0
    h4in = 0.0
    # ampin = pdflogx[1] - pdflogx[0]
    if pin:
        [ampin,muin,sigmain,h3in,h4in] = pin
    try:
        fitvals, fitcov = curve_fit(gausshermite, pdflogx[ilow:ihigh], pdf[ilow:ihigh], p0 = [ampin, muin, sigmain, h3in, h4in], sigma = np.sqrt(pdf[ilow:ihigh]/nvals))
        fita, fitmu, fitsigma, fith3, fith4 = fitvals
        fitpdf = gausshermite(pdflogx,fita,fitmu,fitsigma,fith3,fith4)
        fitchi = np.sum((pdf[ilow:ihigh] - fitpdf[ilow:ihigh])**2 / fitpdf[ilow:ihigh]) * nvals / (ihigh - ilow)
    except RuntimeError:
        print "Error - Curvefit Failed"
        fita = 0.0
        fitmu = 0.0
        fitsigma = 0.0
        fitchi = 0.0
        fith3 = 0.0
        fith4 = 0.0
    
    return fita, fitmu, fitsigma, fith3, fith4, fitchi

def sim_fitlnhermpdf(pdf,pdflogx,pdfm,massmin,massmax,nvals,pin=None):
    fitvals = sim_fitlnhermpdfvals(pdf,pdflogx,pdfm,massmin,massmax,nvals,pin=pin)
    fita, fitmu, fitsigma, fith3, fith4, fitchi = fitvals
    print fitsigma, fitmu, fith3, fith4
    fitpdf = gausshermite(pdflogx,fita,fitmu,fitsigma,fith3,fith4)
    return fitpdf

def sim_fitmeansigma(pdflines, tnorm, mnorm, thist, ndiff):
    
    nconv = 25
    pdftime = sim_pdftime(pdflines, tnorm)
    tminpdfi = np.abs(pdftime - thist).argmin()
    plotlogden = sim_pdfx(pdflines)
    plotden = 10**plotlogden / mnorm
    sigmas = []
    xs = []
    for var in xrange(-ndiff,ndiff):
        plotpdf = sim_pdfm(pdflines,tminpdfi+var)
        plotpdf = np.convolve(plotpdf, np.ones((nconv,))/nconv, mode='same')
        plotpdf = plotpdf / np.sum(plotpdf)
        nvals = np.sum(plotpdf)
        amp, mu, sigma, chisq = sim_fitlnpdfvals(plotpdf,plotlogden,plotpdf,0.15,0.95,nvals)
        sigma = sigma / np.log10(np.exp(1))
        mu = mu / np.log10(np.exp(1))
        mean = np.exp(mu - 0.5 * sigma**2)
        xs.append(mean)
        sigmas.append(sigma)
        
    return np.nanmean(sigmas), np.nanmean(xs)
