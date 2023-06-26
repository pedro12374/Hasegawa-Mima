import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import use
from matplotlib import cm
from scipy.stats import gaussian_kde
from scipy import signal
from numpy import linspace
from numpy.fft import rfft, rfftfreq, ifft
use('GTK3Cairo')

plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif',size=17) # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

fig, ax = plt.subplots(1,3)
plt.tight_layout(pad=0.2) # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(50*0.393, 15*0.393) # esse fatir 0.393 é p converter polegadas p cm """



t_t,x,y,dx,dy= np.loadtxt("an_k.dat",delimiter="\t",unpack=True)#vermelho

ax[0].set_xlabel(r"$x$")
ax[0].set_ylabel(r"$y$")

ax[0].scatter(x, y,  marker='o', s=(5./fig.dpi)**2)
#plt.savefig("PS_Ak.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img

ax[1].set_xlabel(r"$x$")
ax[1].set_ylabel(r"$\frac{dx}{dt}$")


ax[1].scatter(x, dx,  marker='o', s=(5./fig.dpi)**2)
#plt.savefig("PS_dxA.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img


t_t,x,y,dx,dy= np.loadtxt("An_k_t.dat",delimiter="\t",unpack=True)#vermelho
ax[2].set_xlabel(r"$t$")
ax[2].set_ylabel(r"$q$")

ax[2].plot(t_t,x,'-',t_t,y,'-',linewidth=0.4)
ax[2].legend(['x','y'])
#plt.savefig("X_R.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img

plt.savefig("PS_AK.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img

