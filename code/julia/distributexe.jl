using Printf
using PyPlot
using Base.Threads
using Distributed
@everywhere using SpecialFunctions
@everywhere using LinearAlgebra 



@everywhere include("./polar.jl")
@everywhere include("./diffraction.jl")
include("./initboules.jl")

# Le macro @everywhere sert à dire que l'élément qui suit sera chargé dans tous les processus. Sans cela, seul le processus appelant les aurra chargé
# Pour lancer le code: julia distributexe.jl -p "nombre de processus"    Dans notre Cas on a prit 4 prcessus

#using Polar
#import Polar


addprocs(3) #On ajoute 3 ouvriers. Avec le processus courrant, cela fait 4 processus en tout. 
Granular=7 # Granularité

k              = 2*pi
lambda         = 2*pi/k # longueur d'onde
h              = lambda/30  #pas de la grille
taille_espace  = 10 # taille total de la grille
taille_matrice = convert(Int64, taille_espace*1/h)

beta = pi #Angle de l'onde incidente
e    = 10^(-12)
NbrBoulles = 100
	
#####Function########
@everywhere function Image_Calcul_Ligne(Idx_Ligne, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)
	#println("test")
	#println(Idx_Ligne)
	res= zeros(Float64, taille_matrice)
	for j=1:taille_matrice
		x,y = coordonnees(Idx_Ligne, j, h, taille_espace)
		r, lambda = conversion_polaire(x, y)

		if !Is_inDisk(x,y,Obstacle, M)
			res[j] = Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, M, Beta)
			#println("Image [",i,",",j,"] = ", Image[i,j],"\n")
		end

	end

	return res

end

@everywhere function Image_Calcul_Portion(Idx, nbr_ligne, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)

	res= zeros(Float64, nbr_ligne, taille_matrice)
	for i=Idx:(Idx+nbr_ligne-1)
		for j=1:taille_matrice
			x,y = coordonnees(i, j, h, taille_espace)
			r, lambda = conversion_polaire(x, y)

			if !Is_inDisk(x,y,Obstacle, M)
				res[i-Idx+1,j] = Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, M, Beta)
				#println("Image [",i,",",j,"] = ", Image[i,j],"\n")
			end

		end
	end

	return res

end


function Image_Mulit(obstacle, taille_matrice, taille_espace,Cm,Dm, k, M,Beta, h, Granular)

	np = nprocs() #Nombre de processus
	#println(np)

# declaration de la matrice
	Image = zeros(Float64, taille_matrice, taille_matrice)

	i = 1
    #nextidx() = (idx = i; i += 1; idx) #Granularité: A chaque fois qu'un ouvrier demande du travaille on lui passe 1 ligne 
	nextidx() = (idx = i; i+=Granular; idx)
	@sync begin # See below the discussion about all this part.
        for p = 1:np
        	if p != myid() || np == 1 #Seul les 3 processus ouvrier vont vraiment travailler
        		@async begin
        			#println(p)
                    while true
                        idx = nextidx()
                        if idx > taille_matrice
                            break
                        end
                        #temp = remotecall(Image_Calcul_Ligne, p, idx , Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h )

                        longueur=Granular

                        if taille_matrice - idx < Granular
                        	longueur=taille_matrice -idx +1
                        end
                        #println(p)
                        #println(idx)
                        temp = remotecall(Image_Calcul_Portion, p, idx, longueur, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)
                        #println(size(fetch(temp)))
                        Image[idx:(idx+longueur-1),:]=fetch(temp)
                        


                    end
                end

        	end
        end
    end

	# Affichage graphique
	
	imshow(transpose(Image), vmin=-2.5, vmax=2.5, extent=(-5, 5, -5, 5))
	colorbar()
	savefig("resMultProc.svg")
end




##########CODE############

Np = floor(Int64, k*1 + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*1)^(1/3) +1)


#Obstacle=[[0,0,1,Np], [4,4,0.01,Np], [2,3,1,Np], [1,2,1,Np], [2,2,1,Np], [1,1,1,Np], [3,3,1,Np],[1,4,1,Np], [3,1,1,Np], [3,4,1,Np]]
#println(Np,"\n")
(Obstacle, )= initBoulesGrid(NbrBoulles, taille_espace, Np )
#println(Obstacle)	
Dm = Extraire_Dm(NbrBoulles, Obstacle, beta, k)
println("Temps calcule B: ")
B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)
B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)
println("\nTemps calcule A: ")
A  = @time Calcule_A(NbrBoulles, Obstacle, k)
A  = @time Calcule_A(NbrBoulles, Obstacle, k)
C  = Calcule_C(A,B)


# println(Dm)
#println(A)
# println(B)
# println(C)

println("\nTemps calcule C: ")
Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)
Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)
# println(Cm)

@everywhere include("./polar.jl")
@everywhere include("./diffraction.jl") #On réecrit les @everywhere avant d'appeler la fonction Image_Multi sans quoi les autres ...
@everywhere using SpecialFunctions      # process ne reconnaissent pas les fonctions (coordonnees(), Is_inDisk(), besselh...) 
@everywhere using LinearAlgebra 

println("\nTemps calcule Image: ")
@time Image_Mulit(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular)
#@timev Image_Mulit(Obstacle,Cm,Dm,NbrBoulles,beta)
@time Image_Mulit(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular)
