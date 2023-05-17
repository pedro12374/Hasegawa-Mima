import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import os
import sys

folder = sys.argv[1] # pasta com os dados
data = np.loadtxt(folder )

#print(data)
index1 = np.argsort(data[:,0]) #indice da primiera coluna
data = data[index1,:]
#print("------")
#print("sorted first column")
#print(data)

L = 0
val0 = data[0,0]
for i in range(0,len(data[:,0])):
    if data[i,0] == val0:
        L +=1
    else:
        break

for i in range(0,L):
    #print(i)
    s = i*L
    temp = data[s:s+L,:]
    #print(temp)
    index2 = np.argsort(temp[:,1])
    data[s:s+L,:] = temp[index2,:]
#print("------")
#print(data)

#os.system("rm escape.dat")
np.savetxt("escape_f.dat",data)



