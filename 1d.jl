using Printf

########### VARIABLE #############
k = 2*pi
lambda = 2*pi/k # longueur d'onde
h = lambda/10  #pas de la grille
alpha = pi # angle de l'onde incidente
a = 1 # rayon de l'obstacle
taille_espace=6 # taille total de la grille
taille_matrice=convert(Int64, taille_espace*1/h)

########## FONCTION ##############

# a partir de l'indice d'une matrice donne les coordonnées cartésiennes
function coordonnees(i,j,h,taille_espace) 
	x=i*h-taille_espace/2
	y=j*h+taille_espace/2
	return(x,y)
end


# convertit les cartésiennes en polaires
function conversion_polaire(x,y)
	r=sqrt(x*x+y*y)
	lambda=atan(y/x)
	return(r,lambda)
end


########### CODE #############


# declaration de la matrice
M = ones(Float64,taille_matrice, taille_matrice)

for i = 1:taille_matrice
	for j=1:taille_matrice
		x,y=coordonnees(i,j,h,taille_espace)
		r,lambda=conversion_polaire(x,y)
		if r<1
			M[i,j]=0
		else
			M[i,j]=1
		end
	end
end

#println(M[30, 34])
#println(conversion_polaire(1, 1))












