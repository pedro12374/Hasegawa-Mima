import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import use
from matplotlib import cm
from scipy.stats import gaussian_kde
from scipy import signal
from numpy import linspace
from numpy.fft import rfft, rfftfreq, ifft
import sys
plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif',size=24) # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

#fig, ax = plt.subplots()
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco

fig = plt.figure()
fig.set_size_inches(45*0.393, 25*0.393) # esse fatir 0.393 é p converter polegadas p cm """

gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True)
#fig.suptitle('Sharing both axes')


#par = sys.argv[1]
par = ["-0.1","-0.211","-0.32"]
for i in range(len(par)):
    t,phi1,phi2,phi3 = np.loadtxt("phi/phi_"+par[i]+".dat",delimiter="\t",unpack=True)#vermelho
    #ax.set_xlim(0,2000)

    axs[i].plot(t,phi1,label=r"$\phi_1$")
    axs[i].plot(t,phi2,label=r"$\phi_2$")
    axs[i].plot(t,phi3,label=r"$\phi_3$")

    axs[i].text(0.08,0.7,r"$\gamma="+par[i]+"$", horizontalalignment='center',
     verticalalignment='center', transform=axs[i].transAxes)
    #axs[0].plot(x, y ** 2)
    #axs[1].plot(x, 0.3 * y, 'o')
    #axs[2].plot(x, y, '+')
    
    # Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
    ax.set_xlim(0,2000)


fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center', bbox_to_anchor=(0.5, 0.95),
          fancybox=True, shadow=True, ncol=3)
axs[1].set_ylabel(r"$\phi$")
axs[2].set_xlabel(r"$t$")
plt.savefig("grafs/Waves.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img


