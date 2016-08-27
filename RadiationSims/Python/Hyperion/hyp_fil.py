import re
import sys
import numpy as np

def out_fil(outlines):
    i = 0
    for line in outlines:
        if re.search('FILAMENTS', line, re.I):
            filstart = i + 1
            break
        i = i + 1
        
    i = 0
    for line in outlines:
        if re.search('CRITICAL POINTS DATA', line, re.I):
            filend = i
            break
        i = i + 1
    return outlines[filstart:filend]

def find_filheaders(fildata):

    nfils = fildata[0]
    nfils = int(nfils)
    nlen = len(fildata)
    
    filstarts = np.zeros((2,nfils)).astype(int)
    ncount = 0
    for i in xrange(1,nlen):
        line = fildata[i]
        if line[0] != ' ':
            filstarts[0,ncount] = i+1
            if ncount != 0:
                filstarts[1,ncount-1] = i-1
            ncount = ncount + 1
    filstarts[1,nfils-1] = nlen-1

    return filstarts, nfils

def get_fil(fildata, fillocs, nfil, ndims = 3):

    filstart = fillocs[0,nfil]
    filend = fillocs[1,nfil]
    nlen = filend - filstart + 1

    filvals = fildata[filstart:filend+1]
    filcoords = np.zeros((ndims,nlen))

    for i in xrange(0,nlen):

        strcoords = filvals[i]
        coords = strcoords.split(' ')
        filcoords[0,i] = float(coords[1])
        filcoords[1,i] = float(coords[2])
        if (ndims > 2):
            filcoords[2,i] = float(coords[3])

    return filcoords

def extract_r(filcoords):

    dims = np.shape(filcoords)
    nlen = dims[1]
    rfil = np.zeros(nlen)
    
    nrad = 10
    ndiv = 35
    part = np.linspace(-1.0*nrad,nrad,ndiv)
    pars = np.linspace(-1.0*nrad,nrad,ndiv)
    allpart, allpars = np.meshgrid(part, pars)
    nallpars = ndiv * ndiv
    linpart = np.reshape(allpart, nallpars)
    linpars = np.reshape(allpart, nallpars)
    linrs = np.sqrt(linpart**2 + linpars**2)
    
    perpcoords = np.zeros((nlen-1, nallpars, 3))

    for i in xrange(0,nlen-1):
        
        xm = np.squeeze(filcoords[:,i])
        xp = np.squeeze(filcoords[:,i+1])
        xpos = (xm + xp) / 2.0
        xgrad = xp - xm
        delr = np.sqrt(xgrad[0]**2 + xgrad[1]**2 + xgrad[2]**2)
        rfil[i+1] = rfil[i] + delr
        
        xveca = [-1 * xgrad[1] / xgrad[0],1,0]
        xveca = xveca / np.sqrt(xveca[0]**2 + xveca[1]**2 + xveca[2]**2)
        xvecb = [-1 * xgrad[2] / xgrad[0],0,1]
        xvecb = xvecb / np.sqrt(xvecb[0]**2 + xvecb[1]**2 + xvecb[2]**2)
        
        perpcoords[i,:,0] = xveca[0] * linpart + xvecb[0] * linpars + xpos[0]
        perpcoords[i,:,1] = xveca[1] * linpart + xvecb[1] * linpars + xpos[1]
        perpcoords[i,:,2] = xveca[2] * linpart + xvecb[2] * linpars + xpos[2]
        
    return rfil, linrs, perpcoords, nlen, nallpars

def extract_rproj(filcoords):
    
    dims = np.shape(filcoords)
    nlen = dims[1]
    rfil = np.zeros(nlen)
    
    nrad = 50
    ndiv = 1000
    linpart = np.linspace(-1.0*nrad,nrad,ndiv)
    linrs = np.abs(linpart)
    
    perpcoords = np.zeros((nlen-1, ndiv, 2))
    
    for i in xrange(0,nlen-1):
        
        xm = np.squeeze(filcoords[:,i])
        xp = np.squeeze(filcoords[:,i+1])
        xpos = (xm + xp) / 2.0
        xgrad = xp - xm
        #print xm, xp, xpos, xgrad
        delr = np.sqrt(xgrad[0]**2 + xgrad[1]**2)
        rfil[i+1] = rfil[i] + delr
        
        if xgrad[0] > 0.0:
            xveca = [-1 * xgrad[1] / xgrad[0],1]
        else:
            xveca = [1.0, 0.0]
        xveca = xveca / np.sqrt(xveca[0]**2 + xveca[1]**2)
        
        perpcoords[i,:,0] = xveca[0] * linpart + xpos[0]
        perpcoords[i,:,1] = xveca[1] * linpart + xpos[1]
    
        #print xveca[0] * linpart + xpos[0], xveca[1] * linpart + xpos[1]
    
    return rfil, linrs, perpcoords, nlen, ndiv


