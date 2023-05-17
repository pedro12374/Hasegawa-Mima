import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import use
from matplotlib import cm
from scipy.stats import gaussian_kde
from scipy import signal
from numpy import linspace
from numpy.fft import rfft, rfftfreq, ifft
import itertools
use('GTK3Cairo')

plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif') # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

fig, ax = plt.subplots(2,1)
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(20*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """


MSD_MEAN = []

for j in range(100):
    t,x,y = np.loadtxt("dat/traj_"+str(j)+".dat",delimiter="\t",unpack=True)
    aux = 0
    MSD = [0]
    for i in range(1,len(t)):
        aux+= ((y[i]-y[0])*(y[i]-y[0])+(x[i]-x[0])*(x[i]-x[0]))
        MSD.append(aux/i)
    MSD_MEAN.append(MSD)
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")
    ax[0].plot(t,MSD,t,t,'k--')
    #ax.legend(['','t'])
SUME = np.sum(MSD_MEAN,axis=0)
SUME = SUME/len(t)
ax[1].plot(t,SUME,'r',t,t,'k--')
ax[1].set_yscale("log")
ax[1].set_xscale("log")
#ax.set_ylim(0,5)
plt.savefig("MSD.png") 
