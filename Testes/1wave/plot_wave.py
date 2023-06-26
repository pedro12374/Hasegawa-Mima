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


fig, ax = plt.subplots(1,1)
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(15*0.393, 15*0.393) # esse fatir 0.393 é p converter polegadas p cm

#ax.set_ylabel(r"$y$") # Legenda, p renderizar direito precisa do r"$blablabla$"
#ax.set_xlabel(r"$x$")



t,x,y = np.loadtxt("ana.dat",delimiter="\t",unpack=True)#vermelho
#t,x_a,y_a = np.loadtxt("onda_an.dat",delimiter="\t",unpack=True)#vermelho
ax.set_xlim(0,np.pi)
ax.set_ylim(0,2*np.pi)
#ax[1].set_xlim(0,np.pi)
#ax[1].set_ylim(0,2*np.pi)
ax.scatter(x, y,  marker='o', s=(1./fig.dpi)**2)
#ax.plot(x,y,"k,")
#ax[1].scatter(x_a, y_a,  marker='o', s=(1./fig.dpi)**2)

#ax[0].title.set_text("Numerico")
#ax[1].title.set_text("Semi Analtico")
#ax.plot(x,y,"k,")
plt.savefig('1wave.png',bbox_inches='tight',dpi=300)
