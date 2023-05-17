import matplotlib.pyplot as plt
import numpy as np
from matplotlib import use
import matplotlib.gridspec as gridspec
use('GTK3Cairo')

plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif') # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

fig = plt.figure()
fig.set_size_inches(20*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """


#gs = gridspec.GridSpec(2,2)

#plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco

#t,x,y= np.loadtxt("traj.dat",delimiter="\t",unpack=True)
ax = plt.subplot()
#ax.plot(y,x,"k,")
#y = y%(2*np.pi)
#ax.set_xlim(0,2*np.pi)
#ax.set_ylim(0,1)
"""
ax = plt.subplot(gs[1,0])
ax.plot(t[2],x[50])

ax = plt.subplot(gs[1,1])
ax.plot(t[2],y[50])
"""

#fig.savefig("ps.png")
for j in range(100):
    t,x,y = np.loadtxt("dat/traj_"+str(j)+".dat",delimiter="\t",unpack=True)
    y = y%(2*np.pi)
    ax.plot(y,x,'k,')
    #ax.legend(['','t'])
ax.set_xlim(0,2*np.pi)
ax.set_ylim(0,1)
plt.savefig("trajs.png") 
