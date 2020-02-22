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


addprocs(length(Sys.cpu_info()) - 1) #On ajoute 3 ouvriers. Avec le processus courrant, cela fait 4 processus en tout. 
Granular_Image=9 # Granularité pour le calcul de l'image i.e. nombre de lignes qui vont être attribué aux processus qui va demander du travail
Granular_A=2 # Granularité pour le calcul de la matrice *a 
println("On a $(nprocs()) processus\n")




########## Global variables ############


k              = 2*pi
lambda         = 2*pi/k # longueur d'onde
h              = lambda/30  #pas de la grille
taille_espace  = 20 # taille total de la grille
taille_matrice = floor(Int64, taille_espace*1/h) #convert(Int64, taille_espace*1/h)

beta = pi #Angle de l'onde incidente
e    = 10^(-12)

##################################
########## Run Script ############
##################################



####### Création des Boulles ########

println("--------------------------------------------------\n")
println("---------- Initialisation du Probleme ------------\n")
println("--------------------------------------------------\n")

println("Si vous voulez générer des boules aléatoirement tapez 'A'.\nSi vous voulez générer des boulles alignées par grille tapez: 'G' \n")
TypInitBoule=readline();


#-------------------------
# -- Random disposition --
#-------------------------
if (TypInitBoule == "A" || TypInitBoule == "a")
	boules_min = 5
	boules_max = 10

	xmin = -taille_espace/2
	xmax =  taille_espace/2
	ymin = xmin
	ymax = xmax
	rmin = 0.1
	rmax = 1.5

	Obstacle,NbrBoulles = initBoulesAlea(boules_min,boules_max,xmin,xmax,ymin,ymax,rmin,rmax)
	println("Nombre de disques: $NbrBoulles")
end

#-------------------------
# -- grid structure --
#-------------------------
if (TypInitBoule == "G" || TypInitBoule == "g")
	println("Entrez le nombre de disques: \n")

	NbrBoulles = parse(Int, readline()) 
	if sqrt(NbrBoulles) > taille_espace #Si le nombre de boules est trop grand par rapport à la taille de l'espace, on augmente ce dernier

		taille_espace = floor(Int64, sqrt(NbrBoulles))+1

		print("Changement de la taille du côté du carré à :")
		print(taille_espace)
		print("\n \n")

		taille_matrice = floor(Int64, taille_espace*1/h)
	end

	(Obstacle, ) = initBoulesGrid(NbrBoulles, taille_espace) #Génération de la grille de disques
end



#--------------------------
# -- custom distribution --  #under work
#--------------------------
# Np = floor(Int64, k*1 + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*1)^(1/3) +1)


#Obstacle=[[0,0,1,Np], [4,4,0.01,Np], [2,3,1,Np], [1,2,1,Np], [2,2,1,Np], [1,1,1,Np], [3,3,1,Np],[1,4,1,Np], [3,1,1,Np], [3,4,1,Np]]
#println(Np,"\n")



####### Calculs de la solution ########
println("Si vous voulez utiliser le parallelisme entrez 'P',\nSinon entrez 'S'\n")
P=readline()
if(P == "P" || P == "p")
	Parallel = true

elseif (P == "S" || P == "s")
	Parallel = false

else
	println("Erreur veuillez entrer 'S' ou 'P' \n")
end

if(Parallel)
	println("--------------------------------------------------\n")
	println("----------- Construction du Probleme -------------\n")
	println("--------------------------------------------------\n")

	Dm = Extraire_Dm(NbrBoulles, Obstacle, beta, k)
	println("Temps calcule B: ")
	B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)

	println("\nTemps calcule A: ")
	@everywhere include("./diffraction_para.jl")
	@everywhere using SpecialFunctions
	@everywhere using LinearAlgebra
	
	A  = @time Calcule_Parallel_A(NbrBoulles, Obstacle, k) #Ici on n'a pas donné la Granularité (param facultatif) : Elle sera définie par la fonction

	println("\n------- Solving the system -------")
	C  = Calcule_C(A,B)


	println("\nTemps calcule C: ")
	Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)


	@everywhere include("./polar.jl")
	@everywhere include("./diffraction_para.jl") #On réecrit les @everywhere avant d'appeler la fonction Image_Multi_Proc ... 
	@everywhere using SpecialFunctions # ... sans quoi les autre process ne reconnaissent pas les fonctions ...					@everywhere using LinearAlgebra    #... (coordonnees(),Is_inDisk(), besselh...)


	println("\n------- Creating Images -------")
	println("\nTemps calcule Image: ")

	Image = @time Image_Multi_Proc(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular_Image)

	save_image_para(Image)

else

	println("--------------------------------------------------\n")
	println("----------- Construction du Probleme -------------\n")
	println("--------------------------------------------------\n")

	vv
	Dm = Extraire_Dm(NbrBoulles, Obstacle, beta, k)
	println("Temps calcule B: ")
	B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)

	println("\nTemps calcule A: ")
	@everywhere include("./diffraction_para.jl")
	@everywhere using SpecialFunctions
	@everywhere using LinearAlgebra

	A  = @time Calcule_A(NbrBoulles, Obstacle, k) #Ici on n'a pas donné la Granularité (param facultatif) : Elle sera définie par la fonction

	println("\n------- Solving the system -------")
	C  = Calcule_C(A,B)




	println("\nTemps calcule C: ")
	Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)


	@everywhere include("./polar.jl")
	@everywhere include("./diffraction_para.jl") #On réecrit les @everywhere avant d'appeler la fonction Image_Multi_Proc ... 
	@everywhere using SpecialFunctions # ... sans quoi les autre process ne reconnaissent pas les fonctions ...					@everywhere using LinearAlgebra    #... (coordonnees(),Is_inDisk(), besselh...)

	println("\n------- Creating Images -------")
	println("\nTemps calcule Image: ")

	@time Image_Mulit(Obstacle, Cm,Dm, NbrBoulles,beta)
end

plot_RCS(Cm, Obstacle, NbrBoulles)






