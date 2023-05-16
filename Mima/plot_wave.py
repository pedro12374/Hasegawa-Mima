import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import colors
from matplotlib import cm
from matplotlib import use
use('GTK3Cairo')
plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif',size=20) # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png


fig, ax = plt.subplots(2,2)
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(30*0.393, 30*0.393) # esse fatir 0.393 é p converter polegadas p cm

#ax.set_ylabel(r"$y$") # Legenda, p renderizar direito precisa do r"$blablabla$"
#ax.set_xlabel(r"$x$")



t,x_1,y_1 = np.loadtxt("dat/wave_0.1.dat",delimiter="\t",unpack=True)#vermelho
t,x_2,y_2 = np.loadtxt("dat/wave_0.5.dat",delimiter="\t",unpack=True)#vermelho
t,x_3,y_3 = np.loadtxt("dat/wave_1.5.dat",delimiter="\t",unpack=True)#vermelho
t,x_4,y_4 = np.loadtxt("dat/wave_2.0.dat",delimiter="\t",unpack=True)#vermelho

for i in range(0,2):
    for j in range(0,2):
        ax[i,j].set_xlim(0,np.pi)
        ax[i,j].set_ylim(0,2*np.pi)
#ax[0,0].scatter(x_1, y_1,  marker='o', s=(0.5/fig.dpi)**2)
#ax[0,1].scatter(x_2, y_2,  marker='o', s=(0.5/fig.dpi)**2)
#ax[1,0].scatter(x_3, y_3,  marker='o', s=(0.5/fig.dpi)**2)
#ax[1,1].scatter(x_4, y_4,  marker='o', s=(0.5/fig.dpi)**2)

ax[0,0].plot(x_1, y_1,  'k,')
ax[0,1].plot(x_2, y_2,  'k,')
ax[1,0].plot(x_3, y_3,  'k,')
ax[1,1].plot(x_4, y_4,  'k,')

ax[0,0].title.set_text(r"$A_2=0.1$")
ax[0,1].title.set_text(r"$A_2=0.5$")
ax[1,0].title.set_text(r"$A_2=1.5$")
ax[1,1].title.set_text(r"$A_2=2.0$")
#ax.plot(x,y,"k,")
plt.savefig('2waves.png',bbox_inches='tight',dpi=300)
