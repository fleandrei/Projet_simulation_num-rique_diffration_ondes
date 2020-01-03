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


function initBoules()
	nb_boules_min = 2
	nb_boules_max = 20
	nb_boules = rand(nb_boules_min:nb_boules_max)
	boules = zeros(Float64, nb_boules, 3)

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
			i = i+1
		end
	end

	#println(boules)
	return (boules, nb_boules)
end

(boules, nb_boules) = initBoules()

taille_espace = 10
h = 10
breaking = 0

taille_matrice = taille_espace * h
M = zeros(Float64, taille_matrice, taille_matrice)

for i = 1:taille_matrice
	for j = 1:taille_matrice
		breaking = 0
		for k = 1:nb_boules
			if (distance(i/h, j/h, boules[k,1], boules[k,2]) < boules[k,3] && breaking == 0)
				breaking = 1
			end
		end
		if (breaking == 0)
			M[i,j] = 1
		end
	end
end

# Affichage graphique
	
imshow(M, extent=(0, 10, 0, 10))
savefig("multiboules.svg")

