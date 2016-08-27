from matplotlib.backends.backend_pdf import PdfPages
import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt

def read_outfile(hostname,outfile):
    sshout = subprocess.Popen(['ssh', hostname, 'cat', outfile],
                              stdout=subprocess.PIPE)
    for line in sshout.stdout.readlines():
        if re.search('t_ff', line, re.I):
            global tff
            tff = eval(line.split()[5])
        if re.search('M_GMC in code', line, re.I):
            global mcloud
            mcloud = eval(line.split()[5])


plotdir = '/Users/sudhirraskutti/Desktop/Thesis/PaperI/Figures/'
datadir = '/u/raskutti/PhD/Hyperion/Tests/RadParGrav/'
datafolder = 'UV_M5.0e4_R15.0_N256_Tf4/'
hstfile = 'id0/RadParGrav.hst'
outfile = 'RadParGrav.out'
hostname = 'raskutti@bellona.astro.princeton.edu'

read_outfile(hostname,datadir + datafolder + outfile)

sshhst = subprocess.Popen(['ssh', hostname, 'cat', datadir + datafolder + hstfile],
                       stdout=subprocess.PIPE)
data = np.loadtxt(sshhst.stdout)

time = data[:,0] / tff
mstar = data[:,15] / mcloud
mgas = data[:,31] / mcloud
mof = mgas[0] - mgas - mstar
ke = data[:,19]
pegas = data[:,17]
alphavir = -2.0 * ke / pegas


plt.figure()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='12')

plt.subplot(2,2,1)
pmstar, pmgas, pmof, = plt.plot(time, mstar, time, mgas, time, mof, linewidth = 2.0)
plt.legend((pmstar, pmgas, pmof), (r"$\displaystyle M_*$",r"$\displaystyle M_{\rm gas}$",r"$\displaystyle M_{\rm of}$"))
plt.axis([0,3,0,1.5])
plt.xticks([0,1,2,3])
plt.yticks([0,0.5,1,1.5])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle M / M_{\rm cl,0}$")
plt.text(0.05*3,0.9*1.5,r"$\displaystyle(a)$")

plt.subplot(2,2,2)
palphavir = plt.plot(time, alphavir, linewidth = 2.0)
plt.yscale('log')
plt.axis([0,3,0.1,100])
plt.xticks([0,1,2,3])
plt.yticks([0.1,1,10,100])

plt.xlabel(r"$\displaystyle t / t_{\rm ff}$")
plt.ylabel(r"$\displaystyle \alpha_{\rm vir}$")
plt.text(0.05*3,50.0,r"$\displaystyle(b)$")

pp = PdfPages(plotdir + 'f1.pdf')
pp.savefig()
pp.close()


