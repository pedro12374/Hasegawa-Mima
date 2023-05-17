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

fig, ax = plt.subplots()
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(20*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """


kx_range = np.arange(1, 7, 1)
ky_range = np.arange(1, 7, 1)
par=[]
kx = [0,0,0]
ky = [0,0,0]
# Gerar todas as combinações de kx e ky
combinations = itertools.product(kx_range, ky_range, repeat=3)
# Percorrer as combinações e verificar as condições
# for kx1, ky1, kx2, ky2, kx3, ky3 in combinations:
#     if kx2 + kx3 == kx1 and ky2 + ky3 == ky1 and not kx1-ky1==0 and not kx2-ky2==0 and not kx3-ky3==0 and (not(kx2*ky3==kx3*ky2)  and not(kx3*ky1==kx1*ky3) and not(kx1*ky2==kx2*ky1)    ):

#         par.append([kx1,kx2,kx3,ky1,ky2,ky3])

#file = "dat/kx_{:.1f}_{:.1f}_{:.1f}_ky_{:.1f}_{:.1f}_{:.1f}_.dat".format(par[i][0],par[i][1],par[i][2],par[i][3],par[i][4],par[i][5])
t_t,psi1_t,psi2_t,psi3_t,x_t,y_t= np.loadtxt("teta.dat",delimiter="\t",unpack=True,usecols=range(6))
x = np.split(x_t,36)
t = np.split(t_t,36)
y = np.split(y_t,36)
#a,b = np.loadtxt("testa.dat",delimiter="\t",unpack=True)
#print(a,b)

#c = np.split(a,2)
#d = np.split(b,2)
#print(c,d)


for j in range(36):
    aux = 0
    MSD = [0]
    for i in range(1,len(x[0])):
        aux+= ((x[j][i]-x[j][0])*(x[j][i]-x[j][0])+(y[j][i]-y[j][0])*(y[j][i]-y[j][0]))
        MSD.append(aux/i)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(t[0],MSD,t[0],t[0],'k--')
    #ax.legend(['','t'])
plt.savefig("MSDx.png") 
