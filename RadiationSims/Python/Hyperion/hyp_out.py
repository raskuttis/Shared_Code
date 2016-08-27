import re
import sys
import numpy as np

def out_tff(outlines):
    for line in outlines:
        if re.search('t_ff', line, re.I):
            tff = eval(line.split()[5])
            return tff

def out_dtout(outlines):
    i = 0
    for line in outlines:
        i = i + 1
        if re.search('output2', line, re.I):
            line = outlines[i+1]
            dtout = eval(line.split()[2])
            return dtout

def out_Psi(outlines):
    for line in outlines:
        if re.search('Psi', line, re.I):
            Psi = eval(line.split()[5])
            return Psi

def out_tsec(outlines):
    for line in outlines:
        if re.search('units.s', line, re.I):
            tsec = eval(line.split()[2])
            return tsec

def out_gram(outlines):
    for line in outlines:
        if re.search('units.g', line, re.I):
            gram = eval(line.split()[2])
            return gram

def out_cm(outlines):
    for line in outlines:
        if re.search('units.cm', line, re.I):
            cm = eval(line.split()[2])
            return cm

def out_gauss(outlines):
    gram = out_gram(outlines)
    cm = out_cm(outlines)
    tsec = out_tsec(outlines)
    gauss = np.sqrt(gram / cm) / tsec
    return gauss

def out_tMyr(outlines):
    for line in outlines:
        if re.search('units.s', line, re.I):
            tsec = eval(line.split()[2])
            return tsec * 3.16e13

def out_mcloud(outlines):
    for line in outlines:
        if re.search('M_GMC in code', line, re.I):
            mcloud = eval(line.split()[5])
            return mcloud

def out_Ekin(outlines):
    for line in outlines:
        if re.search('Ekin in code', line, re.I):
            Ekin = eval(line.split()[5])
            return Ekin

def out_Egrav(outlines):
    for line in outlines:
        if re.search('Egrav in code', line, re.I):
            Egrav = eval(line.split()[5])
            return Egrav

def out_alphavir(outlines):
    Ekin = out_Ekin(outlines)
    Egrav = out_Egrav(outlines)
    alpha = 2.0 * Ekin / Egrav
    return alpha

def out_vturb(outlines):
    for line in outlines:
        if re.search('v_turb in code', line, re.I):
            vturb = eval(line.split()[5])
            return vturb

def out_kappa(outlines):
    for line in outlines:
        if re.search('kappa', line, re.I):
            kappa = eval(line.split()[5])
            return kappa

def out_rcloud(outlines):
    for line in outlines:
        if re.search('rcloud', line, re.I):
            rcloud = eval(line.split()[5])
            return rcloud

def out_clight(outlines):
    for line in outlines:
        if re.match('c in code', line, re.I):
            clight = eval(line.split()[5])
            return clight

def out_rhocloud(outlines):
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    return 3.0 * mcloud / (4.0 * np.pi * rcloud**3)

def out_ljeans(outlines):
    for line in outlines:
        if re.search('L_Jeans', line, re.I):
            ljeans = eval(line.split()[5].replace(',',''))
            return ljeans

def out_G(outlines):
    cm = out_cm(outlines)
    sec = out_tsec(outlines)
    gram = out_gram(outlines)
    return 6.67e-8 * cm**3 / (gram * sec**2)


def out_csound(outlines):
    rhocloud = out_rhocloud(outlines)
    grav = out_G(outlines)
    ljeans = out_ljeans(outlines)
    return np.sqrt(4.0 * np.pi * grav * rhocloud) * ljeans / (2.0 * np.pi)

def out_dx(outlines):
    for line in outlines:
        if re.search('L_Jeans', line, re.I):
            ljeans = eval(line.split()[5].replace(',',''))
            ljeansdx = eval(line.split()[8])
            dx = float(ljeans) / float(ljeansdx)
            return dx

def out_msol(outlines):
    for line in outlines:
        if re.search('units.g', line, re.I):
            msol = eval(line.split()[2]) * 2.0e33
            return msol

def out_sigma(outlines):
    mcloud = out_mcloud(outlines)
    rcloud = out_rcloud(outlines)
    msol = out_msol(outlines)
    sigma = mcloud / (np.pi * rcloud**2 * msol)
    return sigma