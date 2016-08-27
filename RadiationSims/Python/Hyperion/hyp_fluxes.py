import numpy as np
import copy as cp
import scipy.integrate as spint
import scipy.interpolate as spinterp

def flux_time(fluxdata, tnorm):
    data = fluxdata[:,:,0] / tnorm
    return data

def flux_radius(fluxdata, rnorm):
    data = fluxdata[:,:,1] / rnorm
    return data

def flux_full(fluxdata, num):
    data = fluxdata[:,:,num]
    return data

def fluxvst(fluxdata, num, r, rnorm, tnorm):
    allrad = flux_radius(fluxdata, rnorm)
    rinds = np.abs(allrad - r).argmin(axis=1)
    routs = allrad[range(0,np.size(rinds)),rinds]
    alltime = flux_time(fluxdata, tnorm)
    time = alltime[range(0,np.size(rinds)),rinds]
    numdata = fluxdata[:,:,num]
    data = numdata[range(0,np.size(rinds)),rinds]
    return time, data, routs

def fluxvsr(fluxdata, num, t, rnorm, tnorm):
    alltime = flux_time(fluxdata, tnorm)
    time = alltime[:,0]
    tind = np.abs(time - t).argmin()
    allrad = flux_radius(fluxdata, rnorm)
    rad = allrad[tind,:]
    data = fluxdata[tind,:,num]
    return rad, data

def fluxintvsr(fluxdata, num, t, rnorm, tnorm):
    alltime = flux_time(fluxdata, tnorm)
    time = alltime[:,0]
    tind = np.abs(time - t).argmin()
    allrad = flux_radius(fluxdata, rnorm)
    rad = allrad[tind,:]
    data = fluxdata[tind,:,num]
    data = np.cumsum(data * (rad * rnorm)**2) * 4.0 * np.pi * (rad[1] - rad[0]) * rnorm
    return rad, data

def fluxintoutvsr(fluxdata, num, t, rnorm, tnorm):
    alltime = flux_time(fluxdata, tnorm)
    time = alltime[:,0]
    tind = np.abs(time - t).argmin()
    allrad = flux_radius(fluxdata, rnorm)
    rad = allrad[tind,:]
    radrev = np.flipud(rad)
    data = fluxdata[tind,:,num]
    datarev = np.flipud(data)
    datarev = np.cumsum(datarev * (radrev * rnorm)**2) * 4.0 * np.pi * (rad[1] - rad[0]) * rnorm
    data = np.flipud(datarev)
    return rad, data

def fluxintvst(fluxdata, num, r, rnorm, tnorm):
    allrad = flux_radius(fluxdata, rnorm)
    rad = allrad[0,:]
    rind = np.abs(rad - r).argmin()
    alltime = flux_time(fluxdata, tnorm)
    time = alltime[:,rind]
    data = cp.deepcopy(fluxdata[:,rind,num])
    delt = np.ediff1d(time * tnorm)
    data[1:] = np.cumsum(data[1:] * delt)
    data[0] = 0.0
    return time, data

def fluxintsphervst(fluxdata, num, r, rnorm, tnorm):
    time, data, routs = fluxvst(fluxdata, num, r, rnorm, tnorm)
    datasph = 4.0 * np.pi * (routs * rnorm)**2 * data
    delt = np.ediff1d(time * tnorm)
    datasph[1:] = np.cumsum(datasph[1:] * delt)
    datasph[0] = 0.0
    return time, datasph, routs

def fluxintrvst(fluxdata, num, r, rnorm, tnorm):
    allrad = flux_radius(fluxdata, rnorm)
    rad = allrad[0,:]
    rind = np.abs(rad - r).argmin()
    alltime = flux_time(fluxdata, tnorm)
    time = alltime[:,rind]
    alldata = fluxdata[:,:,num]
    data = np.sum(alldata[:,0:rind] * (allrad[:,0:rind] * rnorm)**2,axis=1) * 4.0 * np.pi * (rad[1] - rad[0]) * rnorm
    return time, data

def fluxintrposvst(fluxdata, num, r, rnorm, tnorm):
    allrad = flux_radius(fluxdata, rnorm)
    rad = allrad[0,:]
    rind = np.abs(rad - r).argmin()
    alltime = flux_time(fluxdata, tnorm)
    time = alltime[:,rind]
    alldata = fluxdata[:,:,num]
    data = np.sum(alldata[:,0:rind] * (allrad[:,0:rind] * rnorm)**3,axis=1) * 4.0 * np.pi * (rad[1] - rad[0]) * rnorm
    return time, data

def flux_tableouts(fluxdata, thist, fname, rnorm, tnorm, mnorm, fnorm, time = 0.0, eff = 0.0):
    nconv = 20
    if fname == 'alpha':
        rad, mr = fluxintvsr(fluxdata, 5, thist, rnorm, tnorm)
        mrnorm = mr / np.max(mr)
        goodinds = np.logical_and((mrnorm > 0.05),(mrnorm < 0.95))
        pcoeffs = np.polyfit(np.log10(rad[goodinds]), np.log10(mr[goodinds] / mnorm), 1)
        fout = 3.0 - pcoeffs[0]
    elif fname == 'fedd':
        # fnorm = kappa / clight
        rad, rhofr = fluxvsr(fluxdata, 7, thist, rnorm, tnorm)
        rad, rhodphidr = fluxvsr(fluxdata, 15, thist, rnorm, tnorm)
        fedd = rhofr * fnorm / rhodphidr
        fedd = np.convolve(fedd, np.ones((nconv,))/nconv, mode='same')
        goodind = np.argmin(np.fabs(rad - 1.0))
        fout = fedd[goodind]
    elif fname == 'fabs':
        # fnorm = psi
        timef, frcl, routs = fluxvst(fluxdata, 6, 2.0, rnorm, tnorm)
        fstar = spinterp.interp1d(time, eff * mnorm, bounds_error=False, fill_value='extrapolate')
        lstarf = fnorm * fstar(timef)
        fabsfcl = 1.0 - 4.0 * np.pi * (routs * rnorm)**2 * frcl / lstarf
        tminfi = np.abs(timef - thist).argmin()
        fout = fabsfcl[tminfi]
    elif fname == 'facum':
        # fnorm = psi
        timef, frclint, routs = fluxintsphervst(fluxdata, 6, 2.0, rnorm, tnorm)
        fstar = spinterp.interp1d(time, eff * mnorm, bounds_error=False, fill_value='extrapolate')
        lstarf = fnorm * fstar(timef)
        lstarfint = spint.cumtrapz(lstarf, timef * tnorm, initial=0.0)
        fabsfclint = 1.0 - frclint / lstarfint
        tminfi = np.abs(timef - thist).argmin()
        fout = fabsfclint[tminfi]
    elif fname == 'pradin':
        # fnorm = psi / clight
        timef, frcl, routs = fluxvst(fluxdata, 6, 2.0, rnorm, tnorm)
        fstar = spinterp.interp1d(time, eff * mnorm, bounds_error=False, fill_value='extrapolate')
        lstarf = fnorm * fstar(timef)
        lstarfint = spint.cumtrapz(lstarf, timef * tnorm, initial=0.0)
        tminfi = np.abs(timef - thist).argmin()
        fout = lstarfint[tminfi]
    elif fname == 'prej':
        timef, prejout = fluxintvst(fluxdata, 12, 2.0, rnorm, tnorm)
        tminfi = np.abs(timef - thist).argmin()
        fout = prejout[tminfi]
    else:
        fout = 0.0

    return fout

