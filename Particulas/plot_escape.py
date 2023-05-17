import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import colors
from matplotlib import cm
import sys

plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif',size=20) # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png



par = sys.argv[1]
atr = str(par)+'/escape.dat'


fig, ax = plt.subplots()
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(15*0.393, 15*0.393) # esse fatir 0.393 é p converter polegadas p cm

x,y,t_e,y_e = np.loadtxt(atr,delimiter="\t",unpack=True)#vermelho

#ax.spines['top'].set_color('#FC0FC0') 
#ax.spines['top'].set_linewidth(5)
cmap = colors.ListedColormap(['#E3C567','#607744'])
bounds=[0.0,np.pi,2*np.pi]
norm = colors.BoundaryNorm(bounds, cmap.N)
ax.set_xlim(0,2*np.pi)
ax.set_ylim(0,1)
#ax.set_title(r"$\gamma_1=\gamma_3="+par+"$")
b = ax.scatter(y, x, c=y_e,  marker='o', s=(100./fig.dpi)**2, edgecolors='none', cmap=cmap)
ax.set_xlabel(r"$y$")
ax.set_ylabel(r"$x$")
plt.savefig('grafs/atr_'+str(par)+'.png',bbox_inches='tight',dpi=300)
b.remove()



fig.set_size_inches(17*0.393, 15*0.393)
im = ax.scatter(y, x, c=t_e,  marker='o', s=(100./fig.dpi)**2, edgecolors='none', cmap="jet",norm=colors.LogNorm())
cb = fig.colorbar(im,label = r"$\log(t)$")
ax.set_xlim(0,2*np.pi)
ax.set_ylim(0,1)
#ax.spines['top'].set_color('#FC0FC0') 
ax.set_xlabel(r"$y$")
ax.set_ylabel(r"$x$")
#ax.set_title(r"$\gamma_1=\gamma_3="+par+"$")
plt.savefig('grafs/te_'+str(par)+'.png',bbox_inches='tight',dpi=300)

