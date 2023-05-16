import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import use
from matplotlib import cm
from scipy.stats import gaussian_kde
from scipy import signal
from numpy import linspace
from numpy.fft import rfft, rfftfreq, ifft
use('GTK3Cairo')

plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
plt.rc('font', family='serif') # fonte tipo serif, p fica paredico com latex msm
plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

fig, ax = plt.subplots()
plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
fig.set_size_inches(10*0.393, 10*0.393) # esse fatir 0.393 é p converter polegadas p cm """



t_t,t_ps1_r,t_ps1_i,t_ps2_r,t_ps2_i,t_ps3_r,t_ps3_i = np.loadtxt("test.dat",delimiter="\t",unpack=True)#vermelho

psi = np.add(t_ps1_r,t_ps2_r)
psi = np.add(psi,t_ps3_r)
psi = psi*0.5*10
t_o = t_t/3.82e7

ax.set_xlim(0.00028,0.0004)

f_t = t_o[t_o>0.00028]
f_t = f_t[f_t<=0.0004]

fil_psi = psi[np.isin(t_o,f_t)]




ax.plot(t_o,psi,'k-',linewidth=0.4)
plt.savefig("TSR.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img

#nbins = 100
#n,bins =np.histogram(psi,nbins)
#pdfx = np.zeros(n.size)
#pdfy = np.zeros(n.size)
#for k in range(n.size):
#    pdfx[k] = 0.5*(bins[k]+bins[k+1])
#    pdfy[k] = n[k]

#ax.plot(pdfx,pdfy)
ax.cla()
kde = gaussian_kde(psi)
dist_space = linspace(min(psi),max(psi),100)
ax.plot(dist_space,kde(dist_space))



plt.savefig("dist.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img


""" mod_psi = np.abs(psi)*np.abs(psi)

"""
N = len(f_t)
T = t_o[1]-t_o[0]
x = linspace(f_t[0],N*T,N)

xf = rfftfreq(len(f_t),d=T)

f_signal = rfft(fil_psi)
ft_psi = np.abs(f_signal)


ax.cla()
ax.set_yscale("log")
ax.set_xlim(0,1000000)
ax.plot(xf,ft_psi)
plt.savefig("fft.png",bbox_inches='tight',dpi = 300) ## salva como png, dpi eh a resolução da img
