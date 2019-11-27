using Printf
using PyPlot
using SpecialFunctions
using LinearAlgebra 
using Base.Threads

include("./polar.jl")
include("./diffraction.jl")





function Prod_Scalaire(X,Y)
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

	sum=Atomic{Float64}(0)
	
	@threads for j in 1:nbrThread
		atomic_add!(sum,S[j])
		println(sum)
		#IDThread=threadid()
		#println("IDThread=",IDThread,", S[j]=",S[j])
	end

	return sum


end


res=Prod_Scalaire(ones(50),ones(50))
println(res)