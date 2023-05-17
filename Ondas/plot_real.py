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

fig, ax = plt.subplots(2,2)
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
for kx1, ky1, kx2, ky2, kx3, ky3 in combinations:
    if kx2 + kx3 == kx1 and ky2 + ky3 == ky1 and not kx1-ky1==0 and not kx2-ky2==0 and not kx3-ky3==0 and (not(kx2*ky3==kx3*ky2)  and not(kx3*ky1==kx1*ky3) and not(kx1*ky2==kx2*ky1)    ):

        par.append([kx1,kx2,kx3,ky1,ky2,ky3])

for i in range(len(par)):
    file = "dat/kx_{:.1f}_{:.1f}_{:.1f}_ky_{:.1f}_{:.1f}_{:.1f}_.dat".format(par[i][0],par[i][1],par[i][2],par[i][3],par[i][4],par[i][5])
    t_t,psi1_t,psi2_t,psi3_t,x_t,y_t= np.loadtxt(file,delimiter="\t",unpack=True,usecols=range(6))
    t = np.split(t_t,36)
    psi1 = np.split(psi1_t,36)
    psi2 = np.split(psi2_t,36)
    psi3 = np.split(psi3_t,36)
    x = np.split(x_t,36)
    y = np.split(y_t,36)

    psi = np.add(psi1,psi2)
    psi = np.add(psi,psi3)
    fig, ax = plt.subplots(2,2)
    plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
    fig.set_size_inches(20*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """
    for j in range(36):
        ax[0,0].plot(t[j],psi[j])
        ax[1,0].plot(t[j],x[j])
        ax[1,1].plot(t[j],y[j])
        ax[0,1].scatter(x[j],y[j],s=0.01)
    
    ax[0,0].set_xlabel('t')
    ax[0,0].set_ylabel('phi(t)')
    ax[1,0].set_xlabel('t')
    ax[1,0].set_ylabel('x(t)')
    ax[1,1].set_xlabel('t')
    ax[1,1].set_ylabel('y(t)')
    ax[0,1].set_xlabel('y')
    ax[0,1].set_ylabel('x')

    ax[0,1].set_xlim(0,np.pi)
    ax[0,1].set_ylim(0,2*np.pi)
        #ax[1,0].set_xlim(10000,15000)
        #ax[1,1].set_xlim(10000,15000)
        #ax[0,0].set_xlim(10000,15000)
    plt.savefig("graf/kx_{:.1f}-{:.1f}-{:.1f}ky_{:.1f}-{:.1f}-{:.1f}.png".format(par[i][0],par[i][1],par[i][2],par[i][3],par[i][4],par[i][5])) 
    plt.close("all")

#plt.savefig("TSAA.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img


# ax[1].plot(t_t,t_ps1_r,'-',t_t,t_ps2_r,'-',t_t,t_ps3_r,'-',linewidth=0.4)
# ax[1].legend(['psi 1','psi 2','psi 3'])
# plt.savefig("psi.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img

# plt.close()


# fig, ax = plt.subplots(2,2)
# plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
# fig.set_size_inches(30*0.393, 30*0.393) # esse fatir 0.393 é p converter polegadas p cm """






# ax[0,0].scatter(dx, x,  marker='o', s=(10./fig.dpi)**2)
# #plt.savefig("X_R.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img




# ax[0,1].scatter(x, y,  marker='o', s=(10./fig.dpi)**2)
# ax[0,1].set_xlim(0,2*np.pi)
# ax[0,1].set_ylim(0,np.pi)
# #plt.savefig("Y_R.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img



# ax[1,0].plot(t_t,x,"-",linewidth=0.4)
# a = x+y
# ax[1,1].plot(t_t,y,"-",linewidth=0.4)


# plt.savefig("PS.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img
