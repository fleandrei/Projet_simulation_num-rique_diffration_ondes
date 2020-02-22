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





k              = 2*pi
lambda         = 2*pi/k # longueur d'onde
h              = lambda/30  #pas de la grille
taille_espace  = 20 # taille total de la grille
taille_matrice = floor(Int64, taille_espace*1/h) #convert(Int64, taille_espace*1/h)

beta = pi #Angle de l'onde incidente
e    = 10^(-12)


#######Création des Boulles########""

println("--------------------------------------------------\n")
println("---------- Initialisation du Probleme ------------\n")
println("--------------------------------------------------\n")

println("Si vous voulez générer des boules aléatoirement tapez 'A'.\nSi vous voulez générer des boulles alignées par grille tapez: 'G' \n")
TypInitBoule=readline();

if TypInitBoule=="A"
	boules_min=5
	boules_max=10
	xmin=-taille_espace/2
	xmax=taille_espace/2
	ymin=xmin
	ymax=xmax
	rmin=0.1
	rmax=1.5
	Obstacle,NbrBoulles = initBoulesAlea(boules_min,boules_max,xmin,xmax,ymin,ymax,rmin,rmax)
	println("Nombre de disques: $NbrBoulles")
end


if TypInitBoule=="G"
	println("Entrez le nombre de disques: \n")

	NbrBoulles = parse(Int, readline()) 
	if sqrt(NbrBoulles) > taille_espace #Si le nombre de boules est trop grand par rapport à la taille de l'espace, on augmente ce dernier

		taille_espace=floor(Int64, sqrt(NbrBoulles))+1
		print("Changement de la taille du côté du carré à :")
		print(taille_espace)
		print("\n \n")
		taille_matrice = floor(Int64, taille_espace*1/h)
	end

	(Obstacle, )= initBoulesGrid(NbrBoulles, taille_espace) #Génération de la grille de disques
end

#####Function########
# Calcule la ligne num 'Idx'_Ligne de l'image.
@everywhere function Image_Calcul_Ligne(Idx_Ligne, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)

	#println(Idx_Ligne)
	line_inc  = zeros(Float64, taille_matrice)
	line_diff = zeros(Float64, taille_matrice)
	line_tot  = zeros(Float64, taille_matrice)

	for j=1:taille_matrice
		x,y = coordonnees(Idx_Ligne, j, h, taille_espace)
		r, lambda = conversion_polaire(x, y)

		if !Is_inDisk(x,y,Obstacle, M)
			line_inc[j], line_diff[j], line_tot[j]= Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, M, Beta)
			#println("Image [",i,",",j,"] = ", Image[i,j],"\n")
		end

	end

	return line_inc, line_diff, line_tot

end


# Calcule une portion de l'image: On calcule 'nbr_ligne' lignes à partir de la ligne d'indice 'Idx'
@everywhere function Image_Calcul_Portion(Idx, nbr_ligne, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)


	bloc_inc  = zeros(Float64, nbr_ligne, taille_matrice)
	bloc_diff = zeros(Float64, nbr_ligne, taille_matrice)
	bloc_tot  = zeros(Float64, nbr_ligne, taille_matrice)


	for i=Idx:(Idx+nbr_ligne-1)
		for j=1:taille_matrice
			x,y = coordonnees(i, j, h, taille_espace)
			r, lambda = conversion_polaire(x, y)

			if !Is_inDisk(x,y,Obstacle, M)
				bloc_inc[i-Idx+1,j], bloc_diff[i-Idx+1,j], bloc_tot[i-Idx+1,j] = Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, M, Beta)
				#println("Image [",i,",",j,"] = ", Image[i,j],"\n")
			end

		end
	end

	return bloc_inc, bloc_diff, bloc_tot

end


function Image_Multi_Proc(obstacle, taille_matrice, taille_espace,Cm,Dm, k, M,Beta, h, Granular)

	np = nprocs() #Nombre de processus
	#println(np)

# declaration de la matrice
	# Image = zeros(Float64, taille_matrice, taille_matrice)
	Image_inc  = zeros(Float64, taille_matrice, taille_matrice)
	Image_diff = zeros(Float64, taille_matrice, taille_matrice)
	Image_tot  = zeros(Float64, taille_matrice, taille_matrice)
	println("wsh")

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

                        longueur = Granular # Granularité: nombre de lignes que l'on va confier au processus lorsqu'il viendra demander du travail

                        if taille_matrice - idx < Granular  # Si il reste moins de 'Granular' lignes à calculer
                        	longueur = taille_matrice - idx + 1
                        end

                        temp_inc, temp_diff, temp_tot = remotecall(Image_Calcul_Portion, p, idx, longueur, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)  # On appelle le processus 'p' et on lui demande d'executer la fonction Image_Calcul_Portion avec les arguments (idx, longueur, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)
                        println("pute")
                        #println(size(fetch(temp)))
                        Image_inc[idx:(idx+longueur-1),:]  = fetch(temp_inc) # le résultat renvoyé par remotecall n'est pas directement exploitable c'est pourquoi on utilise fetch()
                        Image_diff[idx:(idx+longueur-1),:] = fetch(temp_diff)
                        Image_tot[idx:(idx+longueur-1),:]  = fetch(temp_tot)


                    end
                end # Sortie de @async

        	end
        end
    end # Sortie de @sync: Attend que toutes les taches soient finit avant de sortir

	return Image_inc, Image_diff, Image_tot
end


function Image_save_para(Image_inc, Image_diff, Image_tot)

# declaration de la matrice
	Image_inc2  = zeros(Float64, taille_matrice, taille_matrice)
	Image_diff2 = zeros(Float64, taille_matrice, taille_matrice)
	Image_tot2  = zeros(Float64, taille_matrice, taille_matrice)
	

	# Image_inc, Image_diff, Image_tot = Image_Mulit(obstacle,Cm,Dm, M,Beta)
	# Parcoure et remplissage de la matrice
	for i = 1:taille_matrice

		for j = 1:taille_matrice
			Image_inc2[i,j]  = real(Image_inc[i,j])
			Image_diff2[i,j] = abs(Image_diff[i,j])
			Image_tot2[i,j]  = abs(Image_tot[i,j])
		end
	end
	
	# Affichage graphique

	scale_min = - taille_espace/2
	scale_max =   taille_espace/2
	# println(Image_tot2)
	# imshow(transpose(Image), vmin=0.0, vmax=2.0, extent=(scale_min, scale_max, scale_min, scale_max))
	imshow(transpose(Image_inc2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("../results/resMult_inc.svg")
	println("----- incindent saved in '../results/resMult_inc.svg' -----:\n")
	PyPlot.clf()
	imshow(transpose(Image_diff2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("../results/resMult_diff.svg")
	println("----- diffraction saved in '../results/resMult_diff.svg' -----:\n")
	PyPlot.clf()
	imshow(transpose(Image_tot2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("../results/resMult_tot.svg")
	println("----- total saved in '../results/resMult_tot.svg' -----:\n")

	return Image_inc2, Image_diff2, Image_tot2
end

function Gener_Sequence(Matrice_Image_Spatial,Pas_Temps, Temps_Final)
	Nbr_Image=0
	Image=zeros(Float64, taille_matrice, taille_matrice)
	scale_min = 0
	scale_max = taille_espace

	for t=0:Pas_Temps:Temps_Final
		for i=1:taille_matrice
			for j=1:taille_matrice
				Image[i,j]=abs(exp(- im*2*pi/lambda*t) * Matrice_Image_Spatial[i,j])
			end
		end
		imshow(transpose(Image), extent=(scale_min, scale_max, scale_min, scale_max))
		#colorbar()
		savefig("resMultProc_t$(t).png")
		Nbr_Image=Nbr_Image+1
	end

end


function RCS(Obstacle, numdisque, Cm, R, Beta, k)
	if numdisque>length(Obstacle)
		println("Erreur le numéro de disque que vous avez entré est trop grand: max=$(length(Obstacle))")
	else

		x=Obstacle[numdisque][1]
		y=Obstacle[numdisque][2]
		origin=beta - pi
		Ysource, Xsource= conversion_cart(R, origin)
		Ei=calculUinc_exact(Xsource,Ysource,Beta,k)
		Deg=zeros(361)
		RCS=zeros(361)

		for i=0:360
			Deg[i+1]=i
			Es=calculUp(R, (origin+Deg2Rad(i))%(2*pi), Cm[numdisque][1], Obstacle[numdisque][4], k)
			println(abs(Es))
			println(real(Ei))
			RCS[i+1]=4*pi*R*abs(Es)*abs(Es)/(real(Ei)*real(Ei))
		end

		plot(Deg,RCS)
		show()
	end

end


########## Run Script ############

Np = floor(Int64, k*1 + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*1)^(1/3) +1)


#Obstacle=[[0,0,1,Np], [4,4,0.01,Np], [2,3,1,Np], [1,2,1,Np], [2,2,1,Np], [1,1,1,Np], [3,3,1,Np],[1,4,1,Np], [3,1,1,Np], [3,4,1,Np]]
#println(Np,"\n")

println("Si vous voulez utiliser le parallelisme entrez 'P',\nSinon entrez 'S'\n")
P=readline()
if(P=="P")
	Parallel=true
elseif P=="S"
	Parallel=false
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
	#B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)
	println("\nTemps calcule A: ")
	@everywhere include("./diffraction_para.jl")
	@everywhere using SpecialFunctions
	@everywhere using LinearAlgebra
	A  = @time Calcule_Parallel_A(NbrBoulles, Obstacle, k) #Ici on n'a pas donné la Granularité (param facultatif) : Elle sera définie par la fonction
	#A  = @time Calcule_Parallel_A(NbrBoulles, Obstacle, k)
	#println("Taille A=$(size(A)), Taille B=$(length(B)) \n")
	println("\n------- Solving the system -------\n")
	C  = Calcule_C(A,B)


	# println(Dm)
	#println(A)
	# println(B)
	# println(C)

	println("\nTemps calcule C: ")
	Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)
	#Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)
	# println(Cm)

	@everywhere include("./polar.jl")
	@everywhere include("./diffraction_para.jl") #On réecrit les @everywhere avant d'appeler la fonction Image_Multi_Proc ... 
	@everywhere using SpecialFunctions # ... sans quoi les autre process ne reconnaissent pas les fonctions ...					@everywhere using LinearAlgebra    #... (coordonnees(),Is_inDisk(), besselh...)

	println("\nTemps calcule Image: ")
	#@time Image_Multi_Proc(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular_Image)
	#@timev Image_Multi_Proc(Obstacle,Cm,Dm,NbrBoulles,beta)
	Image_inc, Image_diff, Image_tot = @time Image_Multi_Proc(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular_Image)
	Image_save_para(Image_inc, Image_diff, Image_tot)

else
	println("--------------------------------------------------\n")
	println("----------- Construction du Probleme -------------\n")
	println("--------------------------------------------------\n")

	Dm = Extraire_Dm(NbrBoulles, Obstacle, beta, k)
	println("Temps calcule B: ")
	B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)
	#B  = @time Calcule_B(NbrBoulles ,Obstacle, beta,k,Dm)
	println("\nTemps calcule A: ")
	@everywhere include("./diffraction_para.jl")
	@everywhere using SpecialFunctions
	@everywhere using LinearAlgebra
	A  = @time Calcule_A(NbrBoulles, Obstacle, k) #Ici on n'a pas donné la Granularité (param facultatif) : Elle sera définie par la fonction
	#A  = @time Calcule_Parallel_A(NbrBoulles, Obstacle, k)
	#println("Taille A=$(size(A)),  B=$(length(B)) \n")
	println("\n------- Solving the system -------\n")
	C  = Calcule_C(A,B)


	# println(Dm)
	#println(A)
	# println(B)
	# println(C)

	println("\nTemps calcule C: ")
	Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)
	#Cm = @time Extraire_Cm(C,NbrBoulles,Obstacle)
	# println(Cm)

	@everywhere include("./polar.jl")
	@everywhere include("./diffraction_para.jl") #On réecrit les @everywhere avant d'appeler la fonction Image_Multi_Proc ... 
	@everywhere using SpecialFunctions # ... sans quoi les autre process ne reconnaissent pas les fonctions ...					@everywhere using LinearAlgebra    #... (coordonnees(),Is_inDisk(), besselh...)

	println("\nTemps calcule Image: ")
	#@time Image_Multi_Proc(Obstacle, taille_matrice, taille_espace, Cm,Dm, k, NbrBoulles,beta, h, Granular_Image)
	#@timev Image_Multi_Proc(Obstacle,Cm,Dm,NbrBoulles,beta)
	@time Image_Mulit(Obstacle, Cm,Dm, NbrBoulles,beta)
	
end
	
RCS(Obstacle, 1, Cm, 30, beta, k)

#Gener_Sequence(Image_Spatiale,1,5)