import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
#use('GTK3Cairo')
import gc



def gridstart(N):
    X = np.linspace(0,10,10)
    Y = np.linspace(0,10,N)
    A = []
    for x in X:
        for y in Y:
            A.append([x,y])
    A = np.array(A)
    return A 

A = gridstart(35)

N = 11**2


pars = [[5.6,5.0],[6.7,5.0],[7.8,5.0],[8.9,5.0],[10.0,5.0]]

for i in pars:
	plt.rcParams["mathtext.fontset"] = "cm" # Fonte matemática pro latex
	plt.rc('font', family='serif') # fonte tipo serif, p fica paredico com latex msm
	plt.rc('text', usetex=False) # esse vc deixa True e for salvar em pdf e False se for p salvar png

	fig = plt.figure()

	gs = GridSpec(2, 3, figure=fig)
	fig.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
	fig.set_size_inches(40*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """

	config = i
	file = "dat/Ea_"+str(round(config[0],1))+"-Eb_"+str(round(config[1],1))+".csv"
    
    #file = "kx_17.5_-7.5_-10.0_ky_5.0_-2.5_-2.5_.csv"
	df = pd.read_csv(file)
	df.columns
	tt = df['time'].to_numpy()
	tft = tt/3.82e7
	x_r = df['x_r'].to_numpy()
	y_r = df['y_r'].to_numpy()

	dx_r = df['dx_r'].to_numpy()

	y_r = y_r%(2*np.pi)
	t = np.split(tt,N)
	tf = np.split(tft,N)
	
	x = np.split(x_r,N)
	dx = np.split(dx_r,N)
	y = np.split(y_r,N)

	del df
	
	
	ax1 = fig.add_subplot(gs[0,0])
	ax2 = fig.add_subplot(gs[1,0])


	ax3 = fig.add_subplot(gs[:,1])
	ax4 = fig.add_subplot(gs[:,2])

	
	E_x = np.linspace(0,np.max(x),100)
	E_y = 0
	E_y = 1.0*config[0]*E_x**2+1*config[1]*E_x   
	for j in range(N):
		ax1.plot(t[j],x[j])
		ax2.plot(t[j],dx[j])
		ax3.plot(y[j],x[j],"k,")
		ax4.plot(E_y,E_x)
	ax4.invert_xaxis()
	ax4.set_ylim(0,np.max(x))
	ax3.set_ylim(0,np.max(x))

	ax1.set_xlabel(r'$t$')
	ax1.set_ylabel(r'$x$')

	ax2.set_xlabel(r'$t$')
	ax2.set_ylabel(r'$dx/dt$')

	ax3.set_xlabel(r'$y$')
	ax3.set_ylabel(r'$x$')

	ax4.set_xlabel(r'$\phi(x)$')
	ax4.set_ylabel(r'$x$')

	"""
    ax[0,0].set_xlabel('t')
    ax[0,0].set_ylabel('phi(t)')
    ax[1,0].set_xlabel('t')
    ax[1,0].set_ylabel('x(t)')
    ax[1,1].set_xlabel('t')
    ax[1,1].set_ylabel('y(t)')
    ax[0,1].set_xlabel('y')
    ax[0,1].set_ylabel('x')
   
	"""
	plt.savefig("graf/Ea_"+str(round(config[0],1))+"-Eb_"+str(round(config[1],1))+".png",dpi=300, bbox_inches='tight')
	print("Ea_"+str(round(config[0],1))+"-Eb_"+str(round(config[1],1)))
	ax1.clear()
	ax2.clear()
	ax3.clear()
	ax4.clear()
	plt.close()
	gc.collect()



"""
file = "Ea_0.0-Eb_4.1"
read_f = "dat/"+file+".csv"
df = pd.read_csv(read_f)








for j in range(N):
    #ax[0,0].plot(tf[j],phi[j])
    ax[1,0].plot(t[j],x[j])
    ax[1,1].plot(t[j],dx[j])
    ax[0,1].scatter(y[j],x[j],s=0.01,c="black")
ax[0,0].set_xlabel('t')
ax[0,0].set_ylabel('phi(t)')
ax[1,0].set_xlabel('t')
ax[1,0].set_ylabel('x(t)')
ax[1,1].set_xlabel('t')
ax[1,1].set_ylabel('y(t)')
ax[0,1].set_xlabel('y')
ax[0,1].set_ylabel('x')

plt.savefig("graf/"+file+".png")
"""
