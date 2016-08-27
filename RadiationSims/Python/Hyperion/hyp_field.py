import numpy as np
from hyp_star import *

def rvhist(weight, x, y, z, rnorm, wnorm):
    xv, yv, zv = np.meshgrid(x, y, z)
    rv = np.sqrt(xv**2 + yv**2 + zv**2)
    wr, rad = np.histogram(rv, bins=100, weights = weight)
    wr = wr / wnorm
    dr = rad[1] - rad[0]
    rad = (rad[0:-1] + dr) / rnorm
    return rad, wr

def rrhohist(weight, x, y, z, rnorm, wnorm):
    
    xv, yv, zv = np.meshgrid(x, y, z)
    rv = np.sqrt(xv**2 + yv**2 + zv**2)
    wr, rad = np.histogram(rv, bins=100, weights = weight)
    nr, rad = np.histogram(rv, bins=100)
    rvflat = np.ndarray.flatten(rv)
    weightflat = np.ndarray.flatten(weight)
    inds = np.digitize(rvflat, rad)
    
    nrcount = np.zeros(100,dtype=int)
    wrmeans = np.zeros(100)
    wrstds = np.zeros(100)
    wrlow = np.zeros(100)
    wrhigh = np.zeros(100)

    for i in xrange(0,100):
        goodinds = np.where(inds == i+1)
        nrcount[i] = np.size(goodinds)
        goodweights = weightflat[goodinds]
        wrmeans[i] = np.mean(goodweights)
        wrstds[i] = np.std(goodweights)
        lowinds = np.where(goodweights < wrmeans[i])
        highinds = np.where(goodweights > wrmeans[i])
        wrlow[i] = np.size(lowinds) * 1.0 / (np.size(goodinds) * 1.0)
        wrhigh[i] = np.size(highinds) * 1.0 / (np.size(goodinds) * 1.0)
    
    dr = rad[1] - rad[0]
    rad = (rad[0:-1] + 0.5 * dr) / rnorm

    return rad, wrmeans / wnorm, wrstds / wnorm, wrlow, wrhigh

def rfrhist(weight, x, y, z, rnorm, wnorm):
    xv, yv, zv = np.meshgrid(x, y, z)
    rv = np.sqrt(xv**2 + yv**2 + zv**2)
    ## Calculate radial component of field
    wx = weight[:,:,:,0]
    wy = weight[:,:,:,1]
    wz = weight[:,:,:,2]
    wrad = (wx * xv + wy * yv + wz * zv) / rv
    wr, rad = np.histogram(rv, bins=100, weights = wrad)
    wr = wr / wnorm
    dr = rad[1] - rad[0]
    rad = (rad[0:-1] + dr) / rnorm
    return rad, wr

def field_centre(weight, x, y, z):
    xv, yv, zv = np.meshgrid(x, y, z)
    mtot = np.sum(weight)
    xcom = np.sum(weight * xv) / mtot
    ycom = np.sum(weight * yv) / mtot
    zcom = np.sum(weight * zv) / mtot
    return xcom, ycom, zcom

def field_genrecentre(x, y, z, xcom, ycom, zcom):
    xout = x - xcom
    yout = y - ycom
    zout = z - zcom
    return xout, yout, zout

def field_starrecentre(x, y, z, stardata, tout, time):
    xeff, yeff, zeff = star_pos(stardata)
    tmini = np.abs(time - tout).argmin()
    xcom = xeff[tmini]
    ycom = yeff[tmini]
    zcom = zeff[tmini]
    print xcom, ycom, zcom
    xout, yout, zout = field_genrecentre(x, y, z, xcom, ycom, zcom)
    return xout, yout, zout

def field_recentre(weight, x, y, z):
    xcom, ycom, zcom = field_centre(weight, x, y, z)
    print xcom, ycom, zcom
    xout, yout, zout = field_genrecentre(x, y, z, xcom, ycom, zcom)
    return xout, yout, zout