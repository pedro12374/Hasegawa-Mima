import matplotlib.pyplot as plt
import numpy as np
from matplotlib import use
import os
import matplotlib.gridspec as gridspec
from matplotlib import colors
import sys


plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif') # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

fig = plt.figure()
fig.set_size_inches(20*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """


#gs = gridspec.GridSpec(2,2)
args = sys.argv[-1] 
a = os.getcwd()

path = os.path.join(a,str(args))

#plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco

t,x,y= np.loadtxt(path+"/PS.dat",delimiter="\t",unpack=True)
x_i,y_i= np.loadtxt(path+"/CI.dat",unpack=True)
y = y%(2*np.pi)
ax = plt.subplot()
"""

y = y%(2*np.pi)
ax.plot(y,x,"k,")
ax.set_xlim(0,2*np.pi)
ax.set_ylim(0,1)
plt.savefig("grafs/PS_"+args+".png") 
"""
im = ax.scatter(y, x, c=t,  marker='o', s=(80./fig.dpi)**2, edgecolors='none', cmap="jet")
cb = fig.colorbar(im)
ax.set_xlim(0,2*np.pi)
ax.set_title(r"$\gamma_1=\gamma_3="+args+"$")
ax.set_ylim(0,1)
plt.savefig("grafs/PST_"+args+".pdf") 



"""
ax = plt.subplot(gs[1,0])
ax.plot(t[2],x[50])

ax = plt.subplot(gs[1,1])
ax.plot(t[2],y[50])
"""
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
"""
