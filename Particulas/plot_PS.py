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
args = sys.argv[-1] 
a = os.getcwd()

path = os.path.join(a,str(args))

#plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco

t,x,y= np.loadtxt(path+"/PS.dat",delimiter="\t",unpack=True)
x_i,y_i= np.loadtxt(path+"/CI.dat",unpack=True)


y = y%(2*np.pi)
ax.plot(y,x,"k,")
ax.plot(y_i,x_i,"ro",markersize=3)
ax.set_xlabel(r"$y$")
ax.set_ylabel(r"$x$")

#ax.set_title(r"$\gamma_1=\gamma_3="+args+"$")
ax.set_xlim(0,2*np.pi)
ax.set_ylim(0,1)
plt.savefig("grafs/PS_"+args+".png",bbox_inches='tight',dpi=300) 

