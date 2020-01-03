using Printf
using PyPlot
using SpecialFunctions
using LinearAlgebra 
using Base.Threads
using Distributed
using InteractiveUtils
#using ArrayFire

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
	
	#@sync begin
	#@async begin
	@threads for i in 1:Taille_vect
		#IDThread=threadid()
		#println("S[$(threadid())]=$(S[threadid()])")
	    @inbounds S[threadid()] += X[i]*Y[i]
	end
	#end
	#end

	
	#@threads for i in 1:floor(Int, Taille_vect/4)
		#IDThread=threadid()
		#println("S[$(threadid())]=$(S[threadid()])")
	#	S[threadid()] += X[i]*Y[i]
	#end
	
	#@threads for i in floor(Int, Taille_vect/4):floor(Int, Taille_vect/2)
	#	#IDThread=threadid()
	#	#println("S[$(threadid())]=$(S[threadid()])")
	#	S[threadid()] += X[i]*Y[i]
	#end
	
	#@threads for i in floor(Int, Taille_vect/2):floor(Int, 3*Taille_vect/4)
		#println("S[$(threadid())]=$(S[threadid()])")
	#	S[threadid()] += X[i]*Y[i]
	#end
	
	#@threads for i in floor(Int, 3*Taille_vect/4):Taille_vect
		#println("S[$(threadid())]=$(S[threadid()])")
	#	S[threadid()] += X[i]*Y[i]
	#end  

	sum=Atomic{Float64}(0)    #Créer un objet atomic de type Float64 et initialisé à 0
	
	@threads for j in 1:nbrThread
		atomic_add!(sum,S[j])
		#println(sum)
		#IDThread=threadid()
		#println("IDThread=",IDThread,", S[j]=",S[j])
	end

	return sum[] # Pour acceder à la valeur d'un type atomic, il faut lui ajouter des crochets [] à la fin 

end


#erreure:
#function Prod_Scalaire_Thread_simple(X,Y)
#	Taille_vect=size(X,1)
#	Res=0
#	@threads for i in 1:Taille_vect
#		Res= Res + X[i]*Y[i]
#	end
#	return Res
#end


function Prod_Scalaire_Sequentiel(X,Y)
	Taille_vect=size(X,1)
	Res=0
	for i in 1:Taille_vect
		@inbounds Res= Res + X[i]*Y[i]
	end
	return Res
end



const X=ones(Float64, 50000000)
const Y=ones(Float64, 50000000)
println("Produit scalaire avec multithreads:\n")
@time res=Prod_Scalaire_Thread(X,Y) #Le macro @time permet d'afficher le temps (ainsi que d'autres info) d'execution de l'instruct suivante
@time res=Prod_Scalaire_Thread(X,Y) # Pour avoir un temps correct on lance 2 fois la fonct car la première fois, @time prend également en compte le temps de compilation. On ne regarde donc que le second temps donné

println(res)

println("\nProduit scalaire avec multiprocess:\n")
#@time res=Prod_Scalaire_Distributed(X,Y)
#@time res=Prod_Scalaire_Distributed(X,Y)
#println(res)

println("\nProduit Scalaire en sequentiel: \n")
@time res=Prod_Scalaire_Sequentiel(X,Y)
@time res=Prod_Scalaire_Sequentiel(X,Y)
println(res)

#println(Sys.CPU_CORES)