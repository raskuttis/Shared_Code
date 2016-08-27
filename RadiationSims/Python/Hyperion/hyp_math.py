import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import newton
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.integrate import odeint
from scipy.special import erf
from scipy.stats import lognorm
import scipy.integrate as spint
from hyp_out import *

def gaussian(x, a, mu, sigma):
    return (a/(np.sqrt(2*np.pi)*sigma)) * np.exp(-(x-mu)**2 /(2*sigma**2))

def gausshermite(x, a, mu, sigma, h3, h4):
    y = (x - mu) / sigma
    amp = (a/(np.sqrt(2*np.pi)*sigma))
    h3fac = h3 / np.sqrt(6) * (2 * np.sqrt(2) * y**3 - 3 * np.sqrt(2) * y)
    h4fac = h4 / np.sqrt(24) * (4 * y**4 - 12 * y**2 + 3)
    expont = -0.5 * y**2 * (1 + h3fac + h4fac)
    return amp * np.exp(expont)

def fabs_th(sigma,mu,kappa):
    def integrand(x):
        return np.exp(-1.0 * kappa * x)
    fabs = 1.0 - lognorm.expect(integrand, args=(sigma,), scale=np.exp(mu))
    return fabs

def eps_th(sigmae, psi, gravc, clight):
    af = psi
    bf = 2.0 * np.pi * clight * gravc
    eps = sigmae * bf / (af - sigmae * bf)
    return eps

def surf_edd(psi):
    return psi / (2.0 * np.pi * 3.0e10 * 6.67e-8) * (3.086e18)**2 / 2.0e33

def eff_simple(x, psi, sigma0):
    fedd = surf_edd(psi)
    fac = -2 * x**2 * fedd / sigma0
    eff = fac + np.sqrt(fac**2 + 1)
    return eff

def delay_pl(time, tcut, a, an, cn):
    y = an*(time-tcut)**a + cn
    return y

def broken_pl(time, tstar, tcut, a, b, an, bn):
    def f1(x):
        return delay_pl(x, tstar, a, an, 0.0)
    def f2(x):
        return delay_pl(x, tcut, b, bn, f1(tcut))
    if tcut < max(time):
        y = np.piecewise(time, [time<=tstar, ((time > tstar) & (time <= tcut)), time > tcut], [0.0, f1, f2])
    else:
        y = time * 0.0
    return y

def mstar_brokenplfit(time, mstar, mmin, mmax):
    tstar = time[np.argmax(mstar > 0.0)]
    tfit = time[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    mfit = mstar[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    def fit_func(time,tcut,a,b,an,bn):
        return broken_pl(time, tstar, tcut, a, b, an, bn)
    tcutguess = time[np.argmax(mstar > 0.1)]
    popt, pcov = curve_fit(fit_func, tfit, mfit, p0 = [tcutguess,10.0,10.0,0.3,0.3])
    return popt

def mstar_brokenlinplfit(time, mstar, mmin, mmax):
    tstar = time[np.argmax(mstar > 0.0)]
    tfit = time[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    mfit = mstar[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    def fit_func(time,tcut,an,bn):
        return broken_pl(time, tstar, tcut, 1.0, 1.0, an, bn)
    tcutguess = time[np.argmax(mstar > 0.1)]
    popt, pcov = curve_fit(fit_func, tfit, mfit, p0 = [tcutguess, 0.3, 0.3])
    return popt

def mstar_brokenlinfit(time, mstar, mmin, mmax):
    tstar = time[np.argmax(mstar > 0.0)]
    tfit = time[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    mfit = mstar[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    nts = len(tfit)
    chisqlow = 1.0e20
    for i in xrange(1,nts-1):
        tlow = tfit[0:i]
        mlow = mfit[0:i]
        thigh = tfit[i:nts-1]
        mhigh = mfit[i:nts-1]
        p = np.polyfit(tlow,mlow,1)
        anflin = p[0]
        p = np.polyfit(thigh,mhigh,1)
        anslin = p[0]
        tcutguess = tfit[i]
        popt = [tcutguess,anflin,anslin]
        chisq = np.sum((mfit - broken_pl(tfit, tstar, tfit[nts/2], 1.0, 1.0, anflin, anslin))**2)
        if chisq <= chisqlow:
            poptout = popt
            chisqlow = chisq
    return poptout

def mstar_plfit(time, mstar, mmin, mmax):
    tstar = time[np.argmax(mstar > 0.0)]
    tfit = time[np.argmax(mstar > mmin):np.argmax(mstar > mmax)] - tstar
    mfit = mstar[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    m,b = np.polyfit(np.log10(tfit), np.log10(mfit), 1)
    return m,b

def lsqerror(mstar, mfit, mmin, mmax):
    mbound = mstar[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    mfitbound = mfit[np.argmax(mstar > mmin):np.argmax(mstar > mmax)]
    return np.sum((mbound - mfitbound)**2)

def ye(epsin, x, sigmain, sigma0, psi):
    sigmae = surf_edd(psi)
    eps = epsin - 1.0e-10
    sigma = sigmain + 1.0e-10
    y = np.log(4.0 * sigmae * x**2 * eps / (sigma0 * (1.0 - eps**2))) / (np.sqrt(2.0) * sigma) - sigma / np.sqrt(8.0)
    return y

def eps_of(eps, x, sigma, sigma0, psi):
    y = ye(eps, x, sigma, sigma0, psi)
    return 0.5 * (1.0 - eps) * (1.0 + erf(y))

def eps_of_all(eps, x, sigma, sigma0, psi):
    y = ye(eps, x, sigma, sigma0, psi)
    return 0.5 * (1.0 + erf(y))

def dmwdmskt(mw,ms,psi,eff,sig0,sigln):
    mg = 1.0 - ms - mw
    fg = mg / (mg + ms)
    sigcgs = sig0 * (mg + ms) * 0.000209
    nf = psi / (4.0 * np.pi * 6.67e-8 * 3.0e10)
    xcrit = np.log(((1.0 - fg) / fg) * nf / sigcgs)
    xsi = 0.5 * (1.0 - erf((-2.0 * xcrit + sigln**2) / (2.0 * np.sqrt(2.0) * sigln)))
    return xsi / eff

def mwkt(mstar,psi,eff,sig0,sigln):
    def integrand(mw, ms):
        return dmwdmskt(mw, ms, psi, eff, sig0, sigln)
    mwind = odeint(integrand,0.0,mstar)
    return mwind

def mtotkt(mstar,psi,eff,sig0,sigln):
    def integrand(mw, ms):
        return dmwdmskt(mw, ms, psi, eff, sig0, sigln)
    mwind = odeint(integrand,0.0,mstar)
    mstar = np.asarray(mstar)
    mwind = np.squeeze(np.asarray(mwind))
    return mstar + mwind

def kteff(psi,eff,sig0,sigln):
    mstar = np.linspace(0.0,1.0,num=10000)
    mtot = mtotkt(mstar,psi,eff,sig0,sigln)
    xind = np.argmin(abs(mtot - 1.0))
    eff = mstar[xind]
    return eff

def ktvals(akt, machkt):
    rkt = 0.5*((3.0-akt) / (2.0-akt)) * ((1.0 - machkt**(2.0 * (2.0-akt))) / (1.0 - machkt**(2.0 * (3.0-akt))))
    sigktsq = np.log(1.0 + rkt * machkt**2 / 4.0)
    sigkt = np.sqrt(sigktsq)
    return sigkt

def vsigmat_out(sigmain, tin, psicgs, epsff, tbreak, epsbreak, tbrmax, outlines, rmfac = 2.0, ctmax=0):
    
    tff = out_tff(outlines)
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    gravc = out_G(outlines)
    kappa = out_kappa(outlines)
    
    sigmaedd = surf_edd(psicgs)
    vnorm = gravc * mcloud / (4.0 * rcloud**2) * tff * epsff

    sigmap = sigmain / sigmaedd
    breakfrac = epsbreak / epsff
    epsstart = sigmap / (1.0 - sigmap)
    tstart = tbreak + (epsstart - epsbreak) / epsff
    rmax = rmfac * rcloud
    nsds = np.size(sigmain)

    def vtsout(t, sigma):
        sigmap = sigma / sigmaedd
        sigfrac = sigmap / (epsff * (1.0 - sigmap))
        highinds = np.where(t > tbrmax)
        lowinds = np.where(t <= tbrmax)
        y = vnorm * (1.0 - np.exp(-159 * sigmap)) / sigmap * ((t-tbreak)**2 + (t-tbreak)*epsbreak - sigfrac**2 + breakfrac**2)
        if (np.size(highinds) > 1):
            y[highinds] = vnorm * (1.0 - np.exp(-159 * sigmap)) / sigmap * ((tbrmax-tbreak)**2 + 2.0*(tbrmax-tbreak)*(t[highinds]-tbrmax) + (t[highinds]-tbreak)*epsbreak - sigfrac**2 + breakfrac**2)
        return y
    def rtsout(t, sigma):
        sigmap = sigma / sigmaedd
        sigfrac = sigmap / (epsff * (1.0 - sigmap))
        highinds = np.where(t > tbrmax)
        y = vnorm * tff * (1.0 - np.exp(-159 * sigmap)) / sigmap * ((t-tbreak)**3 / 3.0 + (t-tbreak)**2*epsbreak / 2.0 - (sigfrac**2 - breakfrac**2) * (t - tbreak)) + rcloud - rmax
        if (np.size(highinds) > 1):
            y[highinds] = vnorm * tff * (1.0 - np.exp(-159 * sigmap)) / sigmap * ((tbrmax-tbreak)*(t[highinds]-tbreak)**2 - 2.0/3.0*(tbrmax-tbreak)**3 - (tbrmax-tbreak)**2*(t[highinds]-tbrmax) + (t[highinds]-tbreak)**2*epsbreak / 2.0 - (sigfrac**2 - breakfrac**2) * (t[highinds] - tbreak)) + rcloud - rmax
        return y

    vstart = vtsout(tstart, sigmain)
    rstart = rtsout(tstart, sigmain) + rmax

    if (nsds > 1):
        tmax = np.zeros(nsds)
        for i in xrange(0, nsds):
            sigma = sigmain[i]
            def rtout(t):
                return rtsout(t, sigma)
            rsol = brentq(rtout, tstart[i] + 0.5)
            #print i, sigma, rsol.x
            tmax[i] = rsol.x
    else:
        def rtout(t):
            return rtsout(t, sigmain)
        rsol = root(rtout, tstart + 0.5)
        tmax = rsol.x

    vmax = vtsout(tmax, sigmain)
    
    vout = vtsout(tin, sigmain)
    rout = rtsout(tin, sigmain) + rmax

    vout = np.where(rout < rcloud, 0.0, vout)
    vout = np.where(rout > rmax, vmax, vout)
    rout = np.where(rout < rcloud, rcloud, rout)
    rout = np.where(rout > rmax, rmax, rout)

    if ctmax == 1:
        return rmax, vmax, tmax
    else:
        return rout, vout, tmax

def vacc_out(mstar, time, t0, v0, vnorm):
    goodinds = np.where(time > t0)
    minteg = spint.cumtrapz(mstar[goodinds], time[goodinds], initial=0)
    vout = vnorm * minteg + v0
    tout = time[goodinds]
    return tout, vout

def racc_out(mstar, time, t0, v0, r0, vnorm):
    tout, vout = vacc_out(mstar, time, t0, v0, vnorm)
    rout = spint.cumtrapz(vout, tout, initial = 0)
    rout = rout + r0
    return tout, rout, vout

def acc_out_lims(mstar, time, t0, v0, r0, rmax, vnorm):
    tout, rout, vout = racc_out(mstar, time, t0, v0, r0, vnorm)
    goodinds = np.where(rout < rmax)
    return tout[goodinds], rout[goodinds], vout[goodinds]
