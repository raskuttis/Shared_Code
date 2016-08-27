import numpy as np
import matplotlib.pylab as plt
<<<<<<< HEAD
import time

from ath_units import *
from scipy.integrate import odeint

CLOSURE = 'M1'
#CLOSURE = 'linear'
#CLOSURE = 'P1'

if   CLOSURE=='M1':
    flim = 2.0*np.sqrt(3.0)/5.0
elif CLOSURE=='linear':
    flim = 1.0
elif CLOSURE=='P1':
    flim = 0.5
Tlim = flim**-0.25
print flim,Tlim

# Initialize physical units:
u = Units_LMT()
=======

from ath_units import Units
from scipy.integrate import odeint

#CLOSURE = 'M1'
CLOSURE = 'linear'
#CLOSURE = 'P1'

# Initialize physical units:
u = Units()
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575

# Dimensionless input parameters:
#f_Eddstar = 0.3
#tau_star  = 1.0
f_Eddstar = 0.03
tau_star  = 10.0
beta_s    = 1.8e-6
k_0       = 10.0**0.5

# Dimensional parameters:
mu        = 2.33 * u.mH
kappa_R0  = 10.0**-1.5 * u.cm**2/u.g
T_0       = 10.0 * u.K

# Derived quantities:
c_star = u.c*beta_s
T_star = mu*c_star**2/u.kB
F_0 = u.c*u.aR*T_star**4
kappa_Rstar = kappa_R0*(T_star/T_0)**2
kappa_Pstar = k_0*kappa_Rstar
Sigma_tot = tau_star/kappa_Rstar
g = kappa_Rstar*F_0/(u.c*f_Eddstar)
h_star = c_star**2/g
rho_star = Sigma_tot/h_star

print 'c_*       = {} [km s^-1]'.format(c_star/u.kms)
print 'T_*       = {} [K]'.format(T_star)
print 'F_0       = {:e} [Lsun kpc^-2]'.format(F_0/(u.Lsun/u.kpc**2))
print 'kappa_R*  = {} [cm^2 g^-1]'.format(kappa_Rstar)
print 'kappa_P*  = {} [cm^2 g^-1]'.format(kappa_Pstar)
print 'Sigma_tot = {} [g cm^-2]'.format(Sigma_tot)
print 'g         = {:e} [dyne g^-1]'.format(g)
print 'h_star    = {:e} [pc]'.format(h_star/u.pc)
print 'rho_star  = {:e} [g cm^-3]'.format(rho_star)
#input()

# Use scipy's ODE solver:
<<<<<<< HEAD
#def chi(f):
#    if   CLOSURE=='M1':
#        return (5.0 - 2.0*np.sqrt(4.0-3.0*f**2))/3.0
#    elif CLOSURE=='linear':
#        return (1.0+2.0*f)/3.0
#    elif CLOSURE=='P1':
#        return 1.0/3.0
#
#def chiprime(f):
#    if   CLOSURE=='M1':
#        return 2.0*f/np.sqrt(4.0-3.0*f**2)
#    elif CLOSURE=='linear':
#        return 2.0/3.0
#    elif CLOSURE=='P1':
#        return 0.0
#
#    
##def lam(T):
##    f = T**-4.0
###    print T,f,chi(f),chiprime(f),chi(f) - f*chiprime(f)   
###    print T,f,chiprime(f),(chi(f)-chi((1.0-1.0e-6)*f))/(1.0e-6*f)
##    return chi(f) - f*chiprime(f)    
#    
#def func(y,t):
#    P = y[0]
#    T = y[1]
#    f = 1/T**4
#    lam = chi(f) - f*chiprime(f)
#    if (T < tol_T):
#        critical = True
##    if (abs(lam) < 1.0e-7):
##        print 'LAM TOO SMALL',lam,f,y,t,critical
###        input()
##    lam = max(1.0e-7,abs(lam))*np.sign(lam)
#    
#    return [-P/T*(1.0-f_Eddstar*T**2), 
#            -tau_star*P/(4.0*lam*T**2)]


def f(T):
    return T**-4.0
    
=======
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
def chi(f):
    if   CLOSURE=='M1':
        return (5.0 - 2.0*np.sqrt(4.0-3.0*f**2))/3.0
    elif CLOSURE=='linear':
        return (1.0+2.0*f)/3.0
    elif CLOSURE=='P1':
        return 1.0/3.0

def chiprime(f):
    if   CLOSURE=='M1':
        return 2.0*f/np.sqrt(4.0-3.0*f**2)
    elif CLOSURE=='linear':
        return 2.0/3.0
    elif CLOSURE=='P1':
        return 0.0

<<<<<<< HEAD
#def chidoubleprime(f):
#    if   CLOSURE=='M1':
#        return 8.0*(4.0-3.0*f**2)**-1.5
#    elif CLOSURE=='linear':
#        return 0.0
#    elif CLOSURE=='P1':
#        return 0.0

def lam(T):
    ff = f(T)
    return chi(ff) - ff*chiprime(ff)
    
#def lamprime(f):
#    return -f*chidoubleprime(f)
    
def myodeint(y0,t):
    dz = 1.0e-3
    z = t[0]
    zmax = t[1]
    d0   = y0[0]
    T0   = y0[1]
    lam0 = lam(T0)
    P0   = d0*T0
    dT0  = -tau_star*P0/(4.0*lam0*T0**2)
    d_soln = [d0]
    T_soln = [T0]
    critical = False
    while z < zmax:    
        z += dz
        P1 = P0*(1.0 - dz/T0*(1.0-f_Eddstar*T0**2))
        T1 = T0 - dz*tau_star*P0/(4.0*lam0*T0**2)
        lam1 = lam(T1)
        dT1 = (lam0/lam1)*(T0/T1)**2*(P1/P0)*dT0
        d_soln.append(P1/T1)
        T_soln.append(T1)
        P0 = P1
        T0 = T1
        lam0 = lam1
        dT0 = dT1
        P1_max = (4.0*lam1/tau_star)*T1**2*(1.0-f_Eddstar*T1**2)
        print z,P1_max,P1,P0,T1,T0,lam1,lam0
        time.sleep(1.0)
        if P1 > P1_max:
            critical = True
            break
        if abs(dT1) < 1.0e-5:
            break
    y1 = np.array([d_soln,T_soln])
    return critical,y1.T
=======
    
#def lam(T):
#    f = T**-4.0
##    print T,f,chi(f),chiprime(f),chi(f) - f*chiprime(f)   
##    print T,f,chiprime(f),(chi(f)-chi((1.0-1.0e-6)*f))/(1.0e-6*f)
#    return chi(f) - f*chiprime(f)    
    
def func(y,t):
    P = y[0]
    T = y[1]
    f = 1/T**4
    lam = chi(f) - f*chiprime(f)
    if (T < tol_T):
        critical = True
    if (abs(lam) < 1.0e-5):
        print 'LAM TOO SMALL',lam,f,y,t,critical
#        input()
    lam = max(1.0e-5,abs(lam))*np.sign(lam)
    
    return [-P/T*(1.0-f_Eddstar*T**2), 
            -tau_star*P/(4.0*lam*T**2)]

>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575

# Solve the coupled ODE system via double-iteration
zmax = 20.0
nz = 200
tol_T = 1.0e-3
<<<<<<< HEAD
tol_intd = 1.0e-2
=======
tol_intd = 2.0e-2
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
maxiters = 1000
critical = False

z = np.linspace(0,zmax,nz+1)
dz = z[1]-z[0]
P_soln = np.zeros(z.shape)
T_soln = np.zeros(z.shape)

T_soln[0] = 4.5  # guess for T(0)
P_soln[0] = 0.2*T_soln[0]  # guess for P(0)

outer_iters = 0
converged_intd = False
<<<<<<< HEAD
converged_T = False
x = []
y = []
=======
converged_T = False    
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
while (not (converged_intd and converged_T)):
    outer_iters += 1

    inner_iters = 0
    converged_T = False    
    T_soln[0] = max(T_soln[0],1.0)
    T_soln[1:] = 0.0
    P_soln[1:] = 0.0
    T_end = 0.0
    while (not converged_T):
        inner_iters += 1
        T_soln[0] = max(T_soln[0],1.0)
        T_soln[1:] = 0.0
        P_soln[1:] = 0.0
        y0 = [P_soln[0], T_soln[0]]
                
        intd = 0.0
        for i in range(nz):
<<<<<<< HEAD
#            dP,dT = func(y0,z[i])
#            if (abs(dT) < tol_T):
##                print 'dT~0'
#                break
#            if (y0[1] < tol_T):
#                print 'T->0', y0[1]
#                break
#            y1 = odeint(func,y0,[z[i],z[i+1]])
            critical,y1 = myodeint(y0,[z[i],z[i+1]])
=======
            dP,dT = func(y0,z[i])
            if (abs(dT) < tol_T):
#                print 'dT~0'
                break
            if (y0[1] < tol_T):
                print 'T->0', y0[1]
                break
            y1 = odeint(func,y0,[z[i],z[i+1]])
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
            if (critical):
                print 'CRITICAL POINT DETECTED'
                critical = False
                break
#            y1 = odeint(func,y0,[z[i],z[i+1]],h0=1.0e-3)
#            dP,dT = func(y1[-1,:],z[i+1])
#            if (np.isnan(dP) or np.isnan(dT)):
#                dP,dT = func(y0,z[i])
#                f = 1.0/y0[1]**4                
#                print 'NAN',dP,dT,f,chi(f),y0
#                input('blah')

<<<<<<< HEAD
#            print critical,y1.shape,type(y1),y1
            input('')
=======
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
            y0 = y1[-1,:]
            P_soln[i+1],T_soln[i+1] = y0
            intd += dz*P_soln[i+1]/T_soln[i+1]
#            T_soln[i+1] = max(y1[-1,0],1.0)
#            P_soln[i+1] = max(y1[-1,0],0.0)

#            dP,dT = func(y0,z[i+1])
#            if (np.isnan(dP) or np.isnan(dT)):
#                f = 1.0/y0[1]**4                
#                print 'NAN',dP,dT,f,chi(f)
#                input('blah')
#            if (abs(dT) < tol or ):
#                break

        T_end_last = T_end
        T_end = T_soln[i-1]
        d_0 = P_soln[0]/T_soln[0]
<<<<<<< HEAD
        x.append(d_0)
        y.append(T_soln[0])
=======
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
#        d_soln = P_soln[:i]/T_soln[:i]
#        intd = dz*np.sum(d_soln[:i])
#        intd = dz*np.sum(d_soln)
#        print d_soln[0],intd,T_soln[0],T_soln[-1],T_end
#        print i,d_soln[0],intd,T_soln[0],T_soln[i],T_end
        
<<<<<<< HEAD
        relerr = T_soln[i] - Tlim
        if (T_end < T_end_last and T_end < Tlim):
=======
        relerr = T_soln[i] - 1.0
        if (T_end < T_end_last and T_end < 1.0):
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
            print 'T TOO SMALL AND NOT GROWING.  EXITING',converged_T
            break
        elif (T_end < 1.0):
            print 'T TOO SMALL.  GROWING...'
#            T_soln[0] *= 1.0 + 0.01
            T_soln[0] *= 1.0 + 0.1
        elif (abs(relerr) < tol_T):
            print '*** T CONVERGED!',T_end,intd
            converged_T = True
            break
        elif (inner_iters > maxiters):
            print 'INNER CONVERGENCE EXCEEDED MAX ITERS!'
            break
        else:            
            T_soln[0] -= relerr/10.0
#        P_soln[0] = d_soln[0]*T_soln[0]        
        P_soln[0] = d_0*T_soln[0]        

    
    
    relerr = intd - 1.0
    if (abs(relerr) < tol_intd):
        converged_intd = True
        break
    elif (outer_iters > maxiters):
        print 'OUTER CONVERGENCE EXCEEDED MAX ITERS!'
        break
    if (intd > 1.0):
        print 'INTD TOO BIG.  SHRINKING...',intd
#        d_soln[0] *= 1.0 - 0.08
        d_0 *= 1.0 - 0.1
    else:
        print 'INTD TOO SMALL.  GROWING...',intd
#        d_soln[0] -= relerr/10.0
#        d_0 -= relerr/10.0
#        d_soln[0] *= 1.0 + 0.08
        d_0 *= 1.0 + 0.01
#    P_soln[0] = d_soln[0]*T_soln[0]     
    P_soln[0] = d_0*T_soln[0]     
#    print 'new d=',d_soln[0],'new P=',P_soln[0],'T=',T_soln[0]
    print 'new d=',d_0,'new P=',P_soln[0],'T=',T_soln[0]
<<<<<<< HEAD
    if (outer_iters%100==0):
        plt.plot(x,y)
        plt.show()
=======
>>>>>>> 1c095849e7aeecf3cee1f352d8ae2ba6336ff575
    

print 'Converged = ',converged_intd and converged_T
#print y1

d_soln = P_soln[:i]/T_soln[:i]
print 'd0 =',d_soln[0],'T0=',T_soln[0]
plt.semilogy(z[:i],d_soln,'b',z[:i],T_soln[:i],'r')
plt.xlim([0,14])
plt.ylim([0.01,5])
#plt.semilogy(z,y1)

print 'done'
plt.show()



