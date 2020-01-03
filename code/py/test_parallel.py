import multiprocessing as mp
import numpy as np
import time



####FONCTIONS
def prod_scal_process(A,B):
	len=A.size
	nbrproc=mp.cpu_count() #Nombre de processus dispo
	pool = mp.Pool(nbrproc) # On utilise nbrproc processus
	granul=len//nbrproc #Granularité

	#Aschedul=np.copy(A)
	#Bschedul=np.copy(B)
	#Aschedul=Aschedul.reshape(nbrproc,-1) #len doit être divisible par nbrproc
	#Bschedul=Bschedul.reshape(nbrproc,-1)
	#print("Aschedul="+str(Aschedul[0]))

	#Partie parallèle
	T0=time.time()
	results=[pool.apply_async(prod_scal, args=(A[x*granul:(x+1)*granul],B[x*granul:(x+1)*granul])) for x in range(nbrproc)] # Il faut utiliser apply en mode asynchrone car cela permet de gagner beaucoup de temps
	T1=time.time()


	output=[p.get() for p in results]
	print("Temps réelement parallèle: "+str(T1-T0)+"s")
	return np.array(output).sum()
		


def prod_scal(A,B):
	len=A.size
	S=0
	for i in range(len):
		S=S+A[i]*B[i]
	return S


#####CODE

A=np.ones(10000000) #A partir de 1 000 000 il devient avantageux d'utiliser le multi proc
B=np.ones(10000000)

T0=time.time()
print(prod_scal(A,B))
T1=time.time()
print("Temps séquentiel avec prod_scal(A,B) : "+str(T1-T0)+"s.")


T0=time.time()
print(prod_scal_process(A,B))
T1=time.time()
print("Temps multi process avec prod_scal_process(A,B): "+str(T1-T0)+"s.")
#Ce qui prend énormément de temps en multi-process c'est le lancement des processus (la partie parallèle en elle même est très rapide)