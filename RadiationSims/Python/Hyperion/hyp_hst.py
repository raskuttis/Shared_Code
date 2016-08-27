import numpy as np

def hst_num(hstdata, num):
    data = hstdata[:,num]
    return data

def hst_time(hstdata, tff):
    time = hstdata[:,0] / tff
    return time

def hst_deltime(hstdata, tff):
    deltime = hstdata[:,1] / tff
    return deltime

def hst_mstar(hstdata, mcloud):
    mstar = hstdata[:,15] / mcloud
    return mstar

def hst_tstar(hstdata, mcloud, tff):
    mstar = hst_mstar(hstdata, mcloud)
    time = hst_time(hstdata, tff)
    tmini = np.min(np.where(mstar > 0.0))
    tstar = time[tmini]
    return tstar

def hst_teff(hstdata, mcloud, tff, frac):
    mstar = hst_mstar(hstdata, mcloud)
    time = hst_time(hstdata, tff)
    eff = hst_eff(hstdata, mcloud)
    tmini = np.min(np.where(mstar > frac * eff))
    tstar = time[tmini]
    return tstar

def hst_mteff(hstdata, mcloud, tff, thist):
    mstar = hst_mstar(hstdata, mcloud)
    time = hst_time(hstdata, tff)
    eff = hst_eff(hstdata, mcloud)
    tmini = np.min(np.where(time > thist))
    mout = mstar[tmini]
    return mout

def hst_mflux(hstdata, mcloud):
    mflux = hstdata[:,20] / mcloud
    return mflux

def hst_eff(hstdata, mcloud):
    mstar = hstdata[:,15]
    eff = max(mstar) / mcloud
    return eff

def hst_mgas(hstdata, mcloud):
    mgas = hstdata[:,2] / mcloud
    return mgas

def hst_mofsimple(hstdata, mcloud):
    mgas = hstdata[:,2] / mcloud
    mstar = hst_mstar(hstdata, mcloud)
    mof = mgas[0] - mgas - mstar
    return mof

def hst_mof(hstdata, mcloud, tff):
    time = hst_time(hstdata, tff)
    deltatime = np.diff(time)
    deltatime = np.append(0, deltatime)
    deltime = hst_deltime(hstdata, tff)
    mflux = hst_mflux(hstdata, mcloud)
    mof = np.cumsum(mflux * deltatime / deltime)
    return mof

def hst_ke(hstdata):
    kex = hstdata[:,6]
    key = hstdata[:,7]
    kez = hstdata[:,8]
    ke = kex + key + kez
    return ke

def hst_pegas(hstdata):
    pegas = hstdata[:,17] / 2.0
    return pegas

def hst_alphavir(hstdata):
    ke = hst_ke(hstdata)
    pe = hst_pegas(hstdata)
    alphavir = -2.0 * ke / pe
    return alphavir

def hst_mach(hstdata, csound):
    ke = hst_ke(hstdata)
    mgas = hstdata[:,2]
    vt = np.sqrt(2.0*ke/mgas)
    return vt/csound

def hst_reffs(hstdata):
    mgas = hstdata[:,2]
    xeff = hstdata[:,34] / mgas
    yeff = hstdata[:,35] / mgas
    zeff = hstdata[:,36] / mgas
    return xeff, yeff, zeff

def hst_rsqeffs(hstdata):
    xeff, yeff, zeff = hst_reffs(hstdata)
    neff = len(xeff)
    mgas = hstdata[:,2]
    xxeff = hstdata[:,37] / mgas - xeff**2
    yyeff = hstdata[:,38] / mgas - yeff**2
    zzeff = hstdata[:,39] / mgas - zeff**2
    xyeff = hstdata[:,40] / mgas - xeff*yeff
    xzeff = hstdata[:,41] / mgas - xeff*zeff
    yzeff = hstdata[:,42] / mgas - yeff*zeff
    cov = np.zeros((neff,3,3))
    cov[:,0,0] = xxeff
    cov[:,1,1] = yyeff
    cov[:,2,2] = zzeff
    cov[:,0,1] = xyeff
    cov[:,1,0] = xyeff
    cov[:,0,2] = xzeff
    cov[:,2,0] = xzeff
    cov[:,1,2] = yzeff
    cov[:,2,1] = yzeff
    return cov

def hst_rgauss(hstdata):
    cov = hst_rsqeffs(hstdata)
    w, v = np.linalg.eig(cov)
    return np.sort(np.sqrt(5.0 * w),1)