import numpy as np                                                                                                        
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import use
from matplotlib import cm
use('GTK3Cairo')


def func(t, x,omega,lbd,gamma):
    phi1,phi2,phi3 = x
    omega1,omega2,omega3 = omega
    lambda1,lambda2,lambda3 = lbd
    gamma1,gamma2,gamma3 = gamma
    dphi1 = -1.0j*omega1+lambda1*np.conj(phi2)*np.conj(phi3)+gamma1*phi1
    dphi2 = -1.0j*omega2+lambda2*np.conj(phi1)*np.conj(phi3)+gamma2*phi2
    dphi3 = -1.0j*omega3+lambda3*np.conj(phi1)*np.conj(phi2)+gamma3*phi3

    return [dphi1,dphi2,dphi3]
          
# Definir as condições iniciais e o intervalo de tempo
x0 = [0.1+0j,0.1+0j,0.1+0j]
lbd = [0.04,-0.5,0.4]
omega = [0.00131,0.00131,0.00131]
gamma = [-0.0211,0.01,-0.0211]
t_span = [0, 15000]

# Resolver a equação diferencial
sol = solve_ivp(func, t_span, x0,args=(omega,lbd,gamma),dense_output=True)
          
t = np.linspace(0,15000,1000000)
z = sol.sol(t) 
          
# Plotar a solução
#phi = 0.5*(z[0].T+z[1].T+z[2].T+np.conj(z[0].T)+np.conj(z[3].T)+np.conj(z[2].T)  )
fig,ax = plt.subplots()
ax.plot(np.abs(z[0].T), np.abs(z[1].T))
#ax.set_xlim(8000,9500)
#plt.plot(t, np.imag(z[0].T), label='Imaginário')
ax.set_xlabel('Tempo')
ax.set_ylabel('x(t)')
#plt.legend()
plt.savefig("testes.png") 
