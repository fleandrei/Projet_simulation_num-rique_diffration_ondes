using Printf
using PyPlot
using Base.Threads
using Distributed
@everywhere using SpecialFunctions
@everywhere using LinearAlgebra 



@everywhere include("./polar.jl")
@everywhere include("./diffraction_para.jl")
include("./initboules.jl")

#######INFO GENERALES
# Le macro @everywhere sert à dire que l'élément qui suit sera chargé dans tous les processus. Sans cela, seul le processus appelant les aurra chargé

# Le macro @time permet d'afficher le temps (ainsi que d'autres info) d'execution de l'instruct devant laquelle il est placé.
# Pour avoir un temps correct on lance 2 fois la fonct car la première fois, @time prend également en compte le temps de compilation. On ne regarde donc que le second temps donné

# Pour lancer le code: julia distributexe.jl -p "nombre de processus"    Dans notre Cas on a prit 4 prcessus
# Pour lancer le code sur des machines distantes: julia distributexe.jl --machinefile hostfile.jl
#######

#using Polar
#import Polar


addprocs(3) #On ajoute 3 ouvriers. Avec le processus courrant, cela fait 4 processus en tout. 
Granular_Image=5 # Granularité pour le calcul de l'image i.e. nombre de lignes qui vont être attribué aux processus ui va demander du travail
Granular_A=2 # Granularité pour le calcul de la matrice *a 

k              = 2*pi
lambda         = 2*pi/k # longueur d'onde
h              = lambda/30  #pas de la grille
taille_espace  = 20 # taille total de la grille
taille_matrice = floor(Int64, taille_espace*1/h) #convert(Int64, taille_espace*1/h)

beta = pi #Angle de l'onde incidente
e    = 10^(-12)


#######Création des Boulles########""

println("Si vous voulez générer des boules aléatoirement tapez 'A'. \n Si vous voulez générer des boulles alignées par grille tapez: 'G' \n")
TypInitBoule=readline();

if TypInitBoule=="A"
	boules_min=5
	boules_max=5
	xmin=-taille_espace/2
	xmax=taille_espace/2
	ymin=xmin
	ymax=xmax
	rmin=0.1
	rmax=1.5
	Obstacle,NbrBoulles = initBoulesAlea(boules_min,boules_max,xmin,xmax,ymin,ymax,rmin,rmax)
	println("Nombre de boulles: $NbrBoulles")
end


if TypInitBoule=="G"
	NbrBoulles = 50
	(Obstacle, )= initBoulesGrid(NbrBoulles, taille_espace)
end

#####Function########
# Calcule la ligne num 'Idx'_Ligne de l'image.
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


# Calcule une portion de l'image: On calcule 'nbr_ligne' lignes à partir de la ligne d'indice 'Idx'
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


function Image_Multi_Proc(obstacle, taille_matrice, taille_espace,Cm,Dm, k, M,Beta, h, Granular)

	np = nprocs() #Nombre de processus
	#println(np)

# declaration de la matrice
	Image = zeros(Float64, taille_matrice, taille_matrice)

	i = 1
    #nextidx() = (idx = i; i += 1; idx) #Granularité: A chaque fois qu'un ouvrier demande du travaille on lui passe 1 ligne 
	
	nextidx() = (idx = i; i+=Granular; idx) #Fonction commune à toutes les tasks lancées dans Image_Multi_Proc. Renvoie l'indice de la prochaine ligne à calculer. 

	@sync begin # On entre dans une Section synchrone: Permet de s'assurer que toutes les tasks sont terminées avant de sortir de cette zone

        for p = 1:np # On parcourt les processus 

        	if p != myid() || np == 1 #Seul les 3 processus ouvrier vont vraiment travailler (sauf si il n'y a qu'un seul processus): Modèle Patron/Ouvrier

        		@async begin # On rentre dans une Section Asynchrone. Le code contenu dans cette zone correspond à une Task. On lance ainsi 'np' task différentes de manière asynchrones MAIS pas dans des 'np' processus différent. Toutes les task s'exécutent sur le processus principal i.e. le processus p=1 

        			#println(p)
                    while true
                        idx = nextidx() # Indice de la prochaine ligne à calculer
                        if idx > taille_matrice # Si il n'y a plus rien à calculer, on sort de la boucle (et la task est terminé)
                            break
                        end

                        #temp = remotecall(Image_Calcul_Ligne, p, idx , Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h )

                        longueur=Granular # Granularité: nombre de lignes que l'on va confier au processus lorsqu'il viendra demander du travail

                        if taille_matrice - idx < Granular  # Si il reste moins de 'Granular' lignes à calculer
                        	longueur=taille_matrice -idx +1
                        end

                        temp = remotecall(Image_Calcul_Portion, p, idx, longueur, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)  # On appelle le processus 'p' et on lui demande d'executer la fonction Image_Calcul_Portion avec les arguments (idx, longueur, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)

                        #println(size(fetch(temp)))
                        Image[idx:(idx+longueur-1),:]=fetch(temp) # le résultat renvoyé par remotecall n'est pas directement exploitable c'est pourquoi on utilise fetch()
                        


                    end
                end # Sortie de @async

        	end
        end
    end # Sortie de @sync: Attend que toutes les taches soient finit avant de sortir

	# Affichage graphique
	scale_min = 0
	scale_max = taille_espace
	
	#imshow(transpose(Image), vmin=-2.5, vmax=2.5, extent=(-5, 5, -5, 5))
	imshow(transpose(Image), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("resMultProc.svg")
end




##########CODE############

Np = floor(Int64, k*1 + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*1)^(1/3) +1)


#Obstacle=[[0,0,1,Np], [4,4,0.01,Np], [2,3,1,Np], [1,2,1,Np], [2,2,1,Np], [1,1,1,Np], [3,3,1,Np],[1,4,1,Np], [3,1,1,Np], [3,4,1,Np]]
#println(Np,"\n")

 
#println(Obstacle)	
Dm = Extraire_Dm(NbrBoulles, Obstacle, beta, k)
println("Temps calcule B: ")
B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)
B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)
println("\nTemps calcule A: ")
@everywhere include("./diffraction_para.jl")
@everywhere using SpecialFunctions
@everywhere using LinearAlgebra
A  = @time Calcule_Parallel_A(NbrBoulles, Obstacle, k)
A  = @time Calcule_Parallel_A(NbrBoulles, Obstacle, k)
println("Taille A=$(size(A)),  B=$(length(B)) \n")
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
@everywhere include("./diffraction_para.jl") #On réecrit les @everywhere avant d'appeler la fonction Image_Multi_Proc ... 
@everywhere using SpecialFunctions # ... sans quoi les autre process ne reconnaissent pas les fonctions ...					@everywhere using LinearAlgebra    #... (coordonnees(),Is_inDisk(), besselh...)

println("\nTemps calcule Image: ")
@time Image_Multi_Proc(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular_Image)
#@timev Image_Multi_Proc(Obstacle,Cm,Dm,NbrBoulles,beta)
@time Image_Multi_Proc(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular_Image)
