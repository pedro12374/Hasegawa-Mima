import os
import sys
from numpy import arange
from numpy import loadtxt
import subprocess
import io
import numpy as np
import time

#como rodar: p3 run.py escape.cpp nucleos pars.txt


# Grade de pontos iniciais
def gridstart(N):
    X = np.linspace(0,1,N)
    Y = np.linspace(0,2*3.1415,N)
    A = []
    for x in X:
        for y in Y:
            A.append([x,y])
    A = np.array(A)
    return A 

def my_makedirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)


comp = "g++ " + sys.argv[1] + " -o "+ sys.argv[1].removesuffix(".cpp")+".out"
print(comp)
os.system(comp)
#time.sleep(5) # tempo p cancelar caso de probelma na compilaca

A = gridstart(8)
#vars = [0]
args = sys.argv[-1] # contem os valores de phi
pars = loadtxt(args)
print(pars)

a = os.getcwd()

path1 = os.path.join(a,str(pars[1]))
my_makedirs(path1)

if(len(os.listdir(path1))!=0):
    os.system("rm "+path1+"/*.dat")

path2 = os.path.join(a,str(pars[1])+"/grafs")
my_makedirs(path2)


Nrun = int(sys.argv[2])
tmax = pars[0]# Tempo total de integração

Nsim = len(A[:,0])  # numero de simulaçoes
Nfull = int(Nsim/Nrun) # Numero de rodadas cheias
Nfinal = Nsim-Nfull*Nrun # Quantidade de programas paralelos caso Nsim n seja multiplo de Nrun

# printa umas groselhas sobre o numero de simulações
print("Nrun,Nsim,Nfull,Nfinal")
print(Nrun,Nsim,Nfull,Nfinal)

n_f = 1
if Nfinal == 0:
    n_f = 0
Npar = Nrun
#
t = []
t0 = time.time()
for i in range(0,Nfull+n_f): 
    if i == Nfull:
        Npar = Nfinal
    #
    run_string = ""
    for j in range(0,Npar): # executa em rodadas de 5
        index = Nrun*i+j
        run_string +=  "./"+sys.argv[1].removesuffix(".cpp")+ ".out  "+str(pars[0]) +" " + str(pars[1])+ " " + str(A[index,0]) + " " + str(A[index,1]) + " " + ">> "+path1+"/escape.dat"+" & "
    run_string += "wait "
    # de fato roda o os prgramas
    os.system(run_string)
    print(run_string+"\n")

trun = time.time()-t0
print(trun)
