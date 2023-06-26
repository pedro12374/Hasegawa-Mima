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

t_t,t_ps1_r,t_ps1_i,t_ps2_r,t_ps2_i,t_ps3_r,t_ps3_i = np.loadtxt("test.dat",delimiter="\t",unpack=True)#vermelho



fig, ax = plt.subplots()
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(10*0.393, 10*0.393) # esse fatir 0.393 é p converter polegadas p cm """


t = t_t[t_t>7500]
ps1_r = t_ps1_r[t_t>7500]
ps2_r = t_ps2_r[t_t>7500]
ps3_r = t_ps3_r[t_t>7500]
ps1_i = t_ps1_i[t_t>7500]
ps2_i = t_ps2_i[t_t>7500]
ps3_i = t_ps3_i[t_t>7500]


ps1_m = np.sqrt(t_ps1_r*t_ps1_r + t_ps1_i*t_ps1_i)
ps2_m = np.sqrt(t_ps2_r*t_ps2_r + t_ps2_i*t_ps2_i)
ps3_m = np.sqrt(t_ps3_r*t_ps3_r + t_ps3_i*t_ps3_i)

ps1 = np.sqrt(ps1_r*ps1_r + ps1_i*ps1_i)
ps2 = np.sqrt(ps2_r*ps2_r + ps2_i*ps2_i)
ps3 = np.sqrt(ps3_r*ps3_r + ps3_i*ps3_i)


ax.plot(t,ps1,"r-",t,ps2,"b--",t,ps3,"g:",linewidth=1)

ax.set_xlim(8000,9500)

ax.set_xticks(np.arange(8000,9600,500))



#plt.savefig("graph.pdf",bbox_inches='tight') ## salva como pdf

plt.savefig("TS.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img

plt.cla()

ax.plot(ps1_m,ps2_m,"k-",linewidth=0.4)




plt.savefig("PS.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img

