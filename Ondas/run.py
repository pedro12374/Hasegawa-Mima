import os
import sys
from numpy import arange
from numpy import savetxt
import subprocess
import io
import numpy as np
import itertools


#como rodar: p3 run.py prog.cpp nucleos extra.txt


def my_makedirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)


Ea_range = np.arange(1, 7, 1)
Eb_range = np.arange(1, 7, 1)
par=[]
kx = [0,0,0]
ky = [0,0,0]
# Gerar todas as combinações de kx e ky
combinations = itertools.product(kx_range, ky_range, repeat=3)
# Percorrer as combinações e verificar as condições
for kx1, ky1, kx2, ky2, kx3, ky3 in combinations:
    if kx2 + kx3 == kx1 and ky2 + ky3 == ky1 and not kx1-ky1==0 and not kx2-ky2==0 and not kx3-ky3==0 and (not(kx2*ky3==kx3*ky2)  and not(kx3*ky1==kx1*ky3) and not(kx1*ky2==kx2*ky1)    ):

        par.append([kx1,kx2,kx3,ky1,ky2,ky3])
#print(par)


#for i in par[:]:
#    print(i[0],i[1],i[2],i[3],i[4],i[5])
print(len(par))


Nr = int(sys.argv[2])    
rad = len(par)//Nr
rest = len(par)%Nr
print(rest)
print(rad)
k = 0
for i in range(0,rad):
    run = ""
    for j in range(0,Nr):
        run = run + "./"+sys.argv[1].removesuffix(".cpp")+ " " + str(par[k][0]) + " " + str(par[k][1]) + " " + str(par[k][2]) + " " + str(par[k][3]) + " " + str(par[k][4]) + " " + str(par[k][5]) + " "  
        k = k+1
        if(k%Nr!=0):
            run = run + " & "
    #print(run)
    #os.system(run)
    

print("AAAAAAA")
if(rest !=0):
    run = ""
    for i in range(Nr*(rad),len(par)):
        
        run = run + "./"+sys.argv[1].removesuffix(".cpp")+ " " + str(par[k][0]) + " " + str(par[k][1]) + " " + str(par[k][2]) + " " + str(par[k][3]) + " " + str(par[k][4]) + " " + str(par[k][5]) + " "  
        if(i%(len(par)-1) !=0):
            run = run + " & "
    #print(run)
    #os.system(run)


