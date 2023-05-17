import matplotlib.pyplot as plt
import numpy as np
from matplotlib import use
import os
import matplotlib.gridspec as gridspec
import sys


plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif',size=20) # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

fig, ax = plt.subplots()
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(15*0.393, 15*0.393) # esse fatir 0.393 é p converter polegadas p cm


#gs = gridspec.GridSpec(2,2)

#plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco

file = "kx_17.5_-7.5_-10.0_ky_5.0_-2.5_-2.5_.dat"
t_t,psi1_t,psi2_t,psi3_t,x,y= np.loadtxt(file,delimiter="\t",unpack=True,usecols=range(6))

y = y%(2*np.pi)
ax.plot(y,x,"k,")
ax.set_xlabel(r"$y$")
ax.set_ylabel(r"$x$")

#ax.set_title(r"$\gamma_1=\gamma_3="+args+"$")
#ax.set_xlim(0,2*np.pi)

plt.savefig("PS.png",bbox_inches='tight',dpi=300) 

