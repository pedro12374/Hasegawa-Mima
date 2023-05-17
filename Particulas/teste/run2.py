import numpy as np
import os
import time
import sys

# Grade de pontos iniciais
def gridstart(N):
    #X = np.linspace(0.2,0.2,N)
    Y = np.linspace(0,2*3.1415,N)
    A = []
    #for x in X:
        #for y in Y:
            #A.append([x,y])
    A = np.array(A)
    return Y 

A = gridstart(100)
vars = [0]

# compilação
t_all = time.time()
fie = sys.argv[1]
program =  fie+".out" # nome do programa
os.system("g++ "+fie+".cpp -o " + program)
time.sleep(5) # tempo p cancelar caso de probelma na compilaca
#os.system("g++ arrumado.cpp -lm -lgsl -o " + program)

 # carrega as cond. inicias num array
Nrun = 12 # numero máximo de programas simultanios
tmax = 500# Número de pontos no arquivo final
#vars = np.hstack((np.linspace(-1,-0.25,25),np.linspace(-0.25,0.25,301),np.linspace(0.25,1,25))) # array do parametro a ser variavel
rootname = "dat" # Nome principal da rodada de experimentos
############################

# this flag indicates if we are doing a large batch of simulations and the results should be
# transfered to another folder. 1 = True, 0 = False
batch_bool = 0  # Basicamente separar os resultados
############################

    ## Algumas maracutais pro processamento paralelo funcionar direito
Nsim = len(A[:])  # numero de simulaçoes
Nfull = int(Nsim/Nrun) # Numero de rodadas cheias
Nfinal = Nsim-Nfull*Nrun # Quantidade de programas paralelos caso Nsim n seja multiplo de Nrun

# printa umas groselhas sobre o numero de simulações
print("Nrun,Nsim,Nfull,Nfinal")
print(Nrun,Nsim,Nfull,Nfinal)

# Renomeia as coisas positivas e negativas
# Oraganiza as pastas do experimento, e cria a pasta com as trajetorias

# Arquivo com os registros do tempo de simulação (Talvez tirar isso aq)
timefile = open( "timelog.dat","w")
timefile.write("# #sim_paralela \t parametro \t tempo(s) \n")
# Coisas pra dar certo o paralelismp
# n_f é p garantir que caso seja um inteiro, o loop n rode a parte final 2x
n_f = 1
if Nfinal == 0:
    n_f = 0
Npar = Nrun
#
t = []
#os.system("rm "+fie+".dat")
## Numero de vezes q vai ter q rodar o role completo, N rodadas cheias + 1 finall
for i in range(0,Nfull+n_f): 
    if i == Nfull:
        Npar = Nfinal
    t0 = time.time()
    run_string = ""
    for j in range(0,Npar): # executa em rodadas de 5
        index = Nrun*i+j
        #print(index,start[index,0],start[index,1],rn)
        run_string += "./"+ program \
        + " 0.2 " + str(A[index]) \
        + " " + str(tmax) \
        + " >> dat/"+fie+"_"+str(index)+".dat" \
        + " & "
    run_string += "wait "
    # de fato roda o os prgramas
    os.system(run_string)
    print(run_string+"\n")
    # Informações a respeito do tempo de computação
    trun = time.time()-t0
    timefile.write(str(Npar) +  "\t" + str(trun) + "\n")
    t.append(trun)
   # print(i,"/",Nfull+n_f,round(trun,3),"T_sim: ",round(avgt/60,2),"m T_batch: ",round(exct/60,2),"m T_all: ",round(exct*len(vars)/(60*60),3),"h")
timefile.close()

time.sleep(1)

#os.system("python3 cat.py escape.dat")


#os.system("shutdown"
