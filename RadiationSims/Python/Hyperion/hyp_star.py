import re
import sys
import copy as cp
import numpy as np
import scipy as sp

def star_mass(stardata):
    mass = np.sum(stardata[:,:,1],axis=1)
    return mass

def star_new_mass(stardata):
    smass = stardata[:,:,1]
    nts = np.shape(smass)[0]
    nstars = np.shape(smass)[1]
    newsmass = np.zeros((nts,nstars))
    for i in xrange(0,nstars):
        sts = np.squeeze(smass[:,i])
        fs = np.squeeze(np.where(sts > 0.0))
        newsmass[fs[0],i] = sts[fs[0]]
    mass = np.sum(newsmass,axis=1)
    mass = np.cumsum(mass)
    return mass

def star_grid(stardata, nz, rbox, rin):
    z = np.linspace(-1.0*rbox, 1.0*rbox, num=nz,endpoint=False)
    delz = z[1] - z[0]
    xg, yg, zg = np.meshgrid(z,z,z)
    rg = np.sqrt(xg**2 + yg**2 + zg**2)
    rglin = np.reshape(rg,nz**3)
    smass = stardata[:,:,1]
    xst = stardata[:,:,2]
    yst = stardata[:,:,3]
    zst = stardata[:,:,4]
    nts = np.shape(smass)[0]
    nstars = np.shape(smass)[1]
    ist = np.floor((xst + rbox) / delz)
    jst = np.floor((yst + rbox) / delz)
    kst = np.floor((zst + rbox) / delz)
    indst = (kst * nz * nz + jst * nz + ist).astype(int)
    for i in xrange(-1,2):
        for j in xrange(-1,2):
            for k in xrange(-1,2):
                lindst = ((kst + k) * nz * nz + (jst + j) * nz + (ist + i)).astype(int)
                indst = np.hstack((indst,lindst))

    stargrid = np.zeros(nts)
    goodinds = np.squeeze(np.where(rglin < rin))
    ngood = len(goodinds)
    lstargrid = np.zeros(nz**3)
    for i in xrange(0,nts):
        starinds = indst[i,:]
        lstargrid[starinds] = 1.0
        ff = np.sum(lstargrid[goodinds]) / (1.0 * ngood)
        stargrid[i] = ff

    return stargrid

def star_merge_mass(stardata, rcloud):
    smass = stardata[:,:,1]
    sx = stardata[:,:,2]
    sy = stardata[:,:,3]
    sz = stardata[:,:,4]
    nts = np.shape(smass)[0]
    nstars = np.shape(smass)[1]
    newsmass = np.zeros((nts,nstars))
    for i in xrange(0,nstars):
        sts = np.squeeze(smass[:,i])
        fs = np.squeeze(np.where(sts > 0.0))
        if (fs[-1] < nts-1):
            dmax = np.max([np.abs(sx[fs[-1],i]), np.abs(sy[fs[-1],i]), np.abs(sz[fs[-1],i])])
            if (dmax < 0.99 * rcloud):
                newsmass[fs[-1],i] = sts[fs[-1]]
    mass = np.sum(newsmass,axis=1)
    mass = np.cumsum(mass)
    return mass

def star_pos(stardata):
    mass = stardata[:,:,1]
    tmass = star_mass(stardata)
    xeff = np.sum(stardata[:,:,2] * mass,axis=1) / tmass
    yeff = np.sum(stardata[:,:,3] * mass,axis=1) / tmass
    zeff = np.sum(stardata[:,:,4] * mass,axis=1) / tmass
    naninds = np.isnan(xeff)
    xeff[naninds] = 0.0
    yeff[naninds] = 0.0
    zeff[naninds] = 0.0
    return xeff, yeff, zeff

def star_vel(stardata):
    mass = stardata[:,:,1]
    tmass = star_mass(stardata)
    vxeff = np.sum(stardata[:,:,5] * mass,axis=1) / tmass
    vyeff = np.sum(stardata[:,:,6] * mass,axis=1) / tmass
    vzeff = np.sum(stardata[:,:,7] * mass,axis=1) / tmass
    naninds = np.isnan(vxeff)
    vxeff[naninds] = 0.0
    vyeff[naninds] = 0.0
    vzeff[naninds] = 0.0
    return vxeff, vyeff, vzeff

def star_mu(stardata):
    xeff, yeff, zeff = star_pos(stardata)
    mu = np.sqrt(xeff**2 + yeff**2 + zeff**2)
    return mu

def star_vcom(stardata):
    vxeff, vyeff, vzeff = star_vel(stardata)
    vcom = np.sqrt(vxeff**2 + vyeff**2 + vzeff**2)
    return vcom

def star_sigma(stardata):
    xeff, yeff, zeff = star_pos(stardata)
    tmass = star_mass(stardata)
    mass = stardata[:,:,1]
    xstar = (stardata[:,:,2].transpose() - xeff).transpose()
    ystar = (stardata[:,:,3].transpose() - yeff).transpose()
    zstar = (stardata[:,:,4].transpose() - zeff).transpose()
    stardist = np.sqrt(xstar**2 + ystar**2 + zstar**2)
    dist = np.sum(stardist * mass, axis=1) / tmass
    return dist

def star_max(stardata):
    xeff, yeff, zeff = star_pos(stardata)
    tmass = star_mass(stardata)
    mass = stardata[:,:,1]
    xstar = (stardata[:,:,2].transpose() - xeff).transpose()
    ystar = (stardata[:,:,3].transpose() - yeff).transpose()
    zstar = (stardata[:,:,4].transpose() - zeff).transpose()
    stardist = np.sqrt(xstar**2 + ystar**2 + zstar**2)
    dist = np.max(stardist, axis=1)
    return dist

def star_dists(stardata):
    xeff, yeff, zeff = star_pos(stardata)
    tmass = star_mass(stardata)
    mass = stardata[:,:,1]
    xstar = (stardata[:,:,2].transpose() - xeff).transpose()
    ystar = (stardata[:,:,3].transpose() - yeff).transpose()
    zstar = (stardata[:,:,4].transpose() - zeff).transpose()
    stardist = np.asarray([xstar.transpose(), ystar.transpose(), zstar.transpose()])
    return stardist, mass

def star_locs(stardata):
    mass = stardata[:,:,1]
    xstar = stardata[:,:,2]
    ystar = stardata[:,:,3]
    zstar = stardata[:,:,4]
    stardist = np.asarray([xstar.transpose(), ystar.transpose(), zstar.transpose()])
    return stardist, mass

def star_locsatt(stardata, time, tint):
    starlocs, mass = star_locs(stardata)
    tmini = np.abs(time - tint).argmin()
    tpos = np.squeeze(starlocs[:,:,tmini])
    tmass = np.squeeze(mass[tmini,:])
    return tpos, tmass

def star_velcs(stardata):
    vxeff, vyeff, vzeff = star_vel(stardata)
    mass = stardata[:,:,1]
    vxstar = (stardata[:,:,5].transpose() - vxeff).transpose()
    vystar = (stardata[:,:,6].transpose() - vyeff).transpose()
    vzstar = (stardata[:,:,7].transpose() - vzeff).transpose()
    return vxstar, vystar, vzstar, mass

def star_distsatt(stardata, time, tint):
    stardist, mass = star_dists(stardata)
    tmini = np.abs(time - tint).argmin()
    dmatrix = np.squeeze(stardist[:,:,tmini]).transpose()
    ms = np.squeeze(mass[tmini,:])
    allds = sp.spatial.distance.squareform(sp.spatial.distance.pdist(dmatrix, 'euclidean')) + 1.0e-10
    return allds, ms

def star_peatt(stardata, time, tint, gravc):
    stardist, mass = star_distsatt(stardata, time, tint)
    nstars = np.shape(stardist)[0]
    pes = np.zeros(nstars)
    for i in xrange(0,nstars):
        tmass = cp.deepcopy(mass)
        tmass[i] = 0.0
        pes[i] = -1.0 * mass[i] * np.sum(gravc * tmass / stardist[i,:])
    return pes, mass

def star_keatt(stardata, time, tint):
    vxstar, vystar, vzstar, mass = star_velcs(stardata)
    tmini = np.abs(time - tint).argmin()
    kes = 0.5 * np.squeeze(mass[tmini,:]) * np.squeeze(vxstar[tmini,:]**2 + vystar[tmini,:]**2 + vzstar[tmini,:]**2)
    return kes

def star_funboundattt(stardata, time, tint, gravc, mcloud):
    kes = star_keatt(stardata, time, tint)
    pes, mass = star_peatt(stardata, time, tint, gravc)
    unbinds = np.where(kes + pes > 0.0)
    funb = np.sum(mass[unbinds]) / mcloud
    return funb

def star_massvsr(stardata):
    xeff, yeff, zeff = star_pos(stardata)
    tmass = star_mass(stardata)
    mass = stardata[:,:,1]
    xstar = (stardata[:,:,2].transpose() - xeff).transpose()
    ystar = (stardata[:,:,3].transpose() - yeff).transpose()
    zstar = (stardata[:,:,4].transpose() - zeff).transpose()
    stardist = np.sqrt(xstar**2 + ystar**2 + zstar**2)
    return stardist, mass

def star_massvsratt(stardata, time, tint):
    allstarrad, allstarmvsr = star_massvsr(stardata)
    tmini = np.abs(time - tint).argmin()
    starrad = allstarrad[tmini,:]
    starmvsr = allstarmvsr[tmini,:]
    sortinds = np.argsort(starrad)
    starrad = starrad[sortinds]
    starmvsr = starmvsr[sortinds]
    starmvsr = np.cumsum(starmvsr)
    starrad, starradinds = np.unique(starrad,return_index=True)
    starmvsr = starmvsr[starradinds]
    return starrad, starmvsr

def star_quartileatt(stardata, qval, time, tint):
    starrad, starmvsr = star_massvsratt(stardata, time, tint)
    starmvsr = starmvsr / np.max(starmvsr)
    goodind = np.argmin(np.fabs(starmvsr - qval))
    return starrad[goodind]

def star_massvstatr(stardata, rint):
    allstarrad, allstarmvsr = star_massvsr(stardata)
    goodstarms = np.where(allstarrad < rint, allstarmvsr, 0.0)
    starmvst = np.sum(goodstarms,axis=1)
    return starmvst

def star_masssigma(stardata,time,tplot):
    xeff, yeff, zeff = star_pos(stardata)
    tmass = star_mass(stardata)
    mass = stardata[:,:,1]
    xstar = (stardata[:,:,2].transpose() - xeff).transpose()
    ystar = (stardata[:,:,3].transpose() - yeff).transpose()
    zstar = (stardata[:,:,4].transpose() - zeff).transpose()
    stardist = np.sqrt(xstar**2 + ystar**2 + zstar**2)
    tmini = np.abs(time - tplot).argmin()
    starrad = stardist[tmini,:]
    starmvsr = mass[tmini,:] / tmass[tmini]
    sortinds = np.argsort(starrad)
    starrad = starrad[sortinds]
    starmvsr = starmvsr[sortinds]
    starmvsr = np.cumsum(starmvsr)
    tmini = np.abs(starmvsr - 0.68).argmin()
    stardist = starrad[tmini]
    return stardist