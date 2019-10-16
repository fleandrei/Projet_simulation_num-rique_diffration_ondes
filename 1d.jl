using Printf

########### VARIABLE #############
k = 2*pi
lambda = 2*pi/k # longueur d'onde
h = lambda/10
alpha = pi # angle de l'onde incidente
a = 1 # rayon de l'obstacle
taille_espace=6 # taille total de la grille

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



M = ones(Float64, 60, 60)
println(M[30, 34])

println(conversion_polaire(1, 1))












