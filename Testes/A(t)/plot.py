import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import use
from matplotlib import cm
from scipy.stats import gaussian_kde
from numpy import linspace
use('GTK3Cairo')

plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif') # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

t,x,y = np.loadtxt("teste2.dat",delimiter="\t",unpack=True)#vermelho



fig, ax = plt.subplots()
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(10*0.393, 10*0.393) # esse fatir 0.393 é p converter polegadas p cm """


ax.plot(y,x,"k,")

ax.set_xlim(0,2*np.pi)
ax.set_ylim(0,np.pi)
x_i=[]
y_i=[]
for i in np.arange(0,np.pi,np.pi/10):
    for j in np.arange(0,2*np.pi,2*np.pi/10):
        x_i.append(i)
        y_i.append(j)


ax.plot(y_i,x_i,'r.')

plt.savefig("2wave.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img







