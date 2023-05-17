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
import os
#use('GTK3Cairo')
import sys 
plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif') # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

fig, ax = plt.subplots(2,1)
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(20*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """
pars = sys.argv[1]

MSD_MEAN = []

a = os.getcwd()

path1 = os.path.join(a,str(pars)+"/trajs")
Nf = len(os.listdir(path1))

H = []
for j in range(Nf):
    t,x,y = np.loadtxt(path1+"/traj_"+str(j)+".dat",delimiter="\t",unpack=True)
    aux = 0
    MSD = [0]
    #H = [0]
    for i in range(1,len(t)):
        aux+= ((y[i]-y[0])*(y[i]-y[0])+(x[i]-x[0])*(x[i]-x[0]))
        MSD.append(aux/i)
        #H.append((MSD[i]-MSD[i-1])/(t[1]-t[0]))  
    MSD_MEAN.append(MSD)
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")
    #ax[0].set_ylim(0,40)
    ax[0].plot(t,MSD,t,t,'k--')
    #ax.legend(['','t'])
SUME = np.sum(MSD_MEAN,axis=0)
SUME = SUME/len(t)
ax[1].plot(t,SUME,'r',t,t,'k--')
ax[1].set_title(r"Média")
ax[1].set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_ylim(0,2)
fig.suptitle(r"$\gamma_1=\gamma_3="+pars+"$")
plt.savefig("grafs/MSD_"+pars+".png") 


