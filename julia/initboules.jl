using Printf
using PyPlot

# a partir de l'indice d'une matrice donne les coordonnées cartésiennes
function coordonnees(i,j,h) 
	x=i*h
	y=-j*h
	return(x,y)
end


function distance(x1,y1,x2,y2)
	d=sqrt((x2-x1)^(2)+(y2-y1)^(2))
	return(d)
end


function initBoulesAlea(Np) # Crée des boulles dont le nombre, la position et le rayon sont aléatoires
	nb_boules_min = 2
	nb_boules_max = 20
	nb_boules = rand(nb_boules_min:nb_boules_max)
	boules = zeros(Float64, nb_boules, 4)

	breaking = 0

	xmin = 0
	xmax = 10
	ymin = 0
	ymax = 10
	rmin = 0.5
	rmax = 3

	i = 0
	while i < nb_boules
		breaking = 0
		r = rand()*(rmax-rmin)+rmin
		xxmin = xmin + r
		xxmax = xmax - r
		x = rand()*(xxmax-xxmin)+xxmin
		yymin = ymin + r
		yymax = ymax - r
		y = rand()*(yymax-yymin)+yymin
		if (i != 0)
			for k = 1:i
				if (sqrt((boules[k,1] - x)^2 + (boules[k,2] - y)^2) < r + boules[k,3] && breaking == 0)
					breaking = 1
				end
			end
		end
		if (breaking == 0)
			boules[i+1, 1] = x
			boules[i+1, 2] = y
			boules[i+1, 3] = r
			boules[i+1, 4] = Np
			i = i+1
		end
	end

	#println(boules)
	return (boules, nb_boules)
end


function initBoulesGrid(NbrBoules, TailleMat, Np) # Crée NbrBoules boulles allignées et espacé à interval constant
	boules = [[] for i=1:NbrBoules]
	nbrBoule_cote=floor(Int64, sqrt(NbrBoules))  #Nombre de boules le long de chaque dimenssion de la matrice
	#nbrBoule_cote=nbrBoule_cote + floor(Int64, (NbrBoules - nbrBoule_cote*nbrBoule_cote + 1)/2) 

	if (NbrBoules - nbrBoule_cote*nbrBoule_cote !=0) 
		nbrBoule_cote=nbrBoule_cote + 1
	end
	
	espace_Boule=TailleMat/nbrBoule_cote    # Longueur nécéssaire pour représenter une boule (la boule en elle même plus son espace vital) sur la matrice.
	espace_Vital=espace_Boule*0.1   # L'espace vital de chaque boulle correspond à 10% de l'espace_Boule
	Rayon= (espace_Boule-espace_Vital)/2

	count=0
	for i in 1:nbrBoule_cote #Parcourt les lignes de la matrice
		j=0
		while (count<NbrBoules && j<nbrBoule_cote) #Parcourt les colonnes de la matrice
			push!(boules[count+1], (i-1)*espace_Boule + espace_Boule/2 - TailleMat/2, j*espace_Boule + espace_Boule/2 - TailleMat/2, Rayon, Np)
			j=j+1
			count =count+1
		end

	end 

	return (boules, count)

end




#######Code exemple

Np=14
#(boules, nb_boules) = initBoulesAlea(Np)
#taille_espace = 10
#(boules, nb_boules)= initBoulesGrid(100, taille_espace, Np )
#h = 10
#breaking = 0
#println(boules)
#taille_matrice = taille_espace * h
#M = zeros(Float64, taille_matrice, taille_matrice)

#for i = 1:taille_matrice
#	for j = 1:taille_matrice
#		breaking = 0
#		for k = 1:nb_boules
#			if (distance(i/h, j/h, boules[k,1], boules[k,2]) < boules[k,3] && breaking == 0)
#				breaking = 1
#			end
#		end
#		if (breaking == 0)
#			M[i,j] = 1
#		end
#	end
#end

# Affichage graphique
	
#imshow(M, extent=(0, 10, 0, 10))
#savefig("multiboules.svg")
