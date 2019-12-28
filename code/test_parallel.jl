using Printf
using PyPlot
using SpecialFunctions
using LinearAlgebra 
using Base.Threads
using Distributed

include("./polar.jl")
include("./diffraction.jl")

@everywhere using SharedArrays

#Pour fixer le nombre de threads, il faut faire "export JULIA_NUM_THREADS=4" dans le terminal
#Pour fixer le nombre de processus locaux il faut faire: "julia nom_program -p n" où n est le nombre de process
#Pour fixer les processus sur des machines distantes 
 
addprocs(3) # Ajoute 3 processus en plus du processus principal qui lance julia. On a donc en tout 4 processus
println("indices des processus:\n")
println(procs()) #Affiche un vecteur conenant les indices des processus
println("Indice des ouvriers (procesus autres que le processus principal): \n")
println(workers()) #Affiche un vecteur conenant les indices des ouvriers (les processus autres que le processus principal)


function Prod_Scalaire_Distributed(X,Y)
	
	Taille_vect=size(X,1)
	output = SharedArray{Int, 1}(Taille_vect); #Array auquel tous les process ont accès
	
	result=@distributed (+) for i = 1:Taille_vect
  		output[i] = X[i]*Y[i];
	end

	return result

end

function Prod_Scalaire_Thread(X,Y)
	nbrThread=nthreads()
	
	Taille_vect=size(X,1)
	S=zeros(nbrThread)
	#println(S)
	#println(nbrThread)
	@threads for i in 1:Taille_vect
		IDThread=threadid()
		#println(IDThread)
		#println(S)
		S[IDThread]=S[IDThread] + X[i]*Y[i]
	end

	sum=Atomic{Float64}(0)    #Créer un objet atomic de type Float64 et initialisé à 0
	
	@threads for j in 1:nbrThread
		atomic_add!(sum,S[j])
		println(sum)
		#IDThread=threadid()
		#println("IDThread=",IDThread,", S[j]=",S[j])
	end

	return sum[] # Pour acceder à la valeur d'un type atomic, il faut lui ajouter des crochets [] à la fin 


end




println("Produit scalaire avec multithreads:\n")
@time res=Prod_Scalaire_Thread(ones(50),ones(50)) #Le macro @time permet d'afficher le temps (ainsi que d'autres info) d'execution de l'instruct suivante

println(res)

println("Produit scalaire avec multiprocess:\n")
@time res=Prod_Scalaire_Distributed(ones(10000000),ones(10000000))
println(res)



#println(Sys.CPU_CORES)