using Printf
using PyPlot
using SpecialFunctions
using LinearAlgebra 

include("./polar.jl")
include("./diffraction.jl")

#using Polar
#import Polar

k              = 2*pi
lambda         = 2*pi/k # longueur d'onde
h              = lambda/60  #pas de la grille
taille_espace  = 6 # taille total de la grille
taille_matrice = convert(Int64, taille_espace*1/h)

beta = pi #Angle de l'onde incidente
e    = 10^(-12)
NbrObstacle = 2
	
#####Function########	
function Image_Mulit(obstacle,Cm,Dm, NbrObstacle)

# declaration de la matrice
	M = zeros(Float64, taille_matrice, taille_matrice)
	
	# Parcoure et remplissage de la matrice
	for i = 1:taille_matrice

		for j = 1:taille_matrice

			x,y      = coordonnees(i, j, h, taille_espace)
			r,lambda = conversion_polaire(x, y)

			if !Is_inDisk(x,y,Obstacle, NbrObstacle)

				M[i,j] = Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, NbrObstacle)
				println("M[",i,",",j,"] = ", M[i,j],"\n")
			end
		end
	end
	
	# Affichage graphique
	
	imshow(M, extent=(-3, 3, -3, 3))
	savefig("resMult.svg")
end




##########CODE############

Np = floor(Int64, k*1 + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*1)^(1/3) +1)



Obstacle=[[-2,0,1,Np], [2,0,1,Np]]

Dm = Extraire_Dm(NbrObstacle, Obstacle, beta, k)
B  = Calcule_B(NbrObstacle ,Obstacle, beta,k,Dm)
A  = Calcule_A(NbrObstacle, Obstacle, k)
C  = Calcule_C(A,B)


# println(Dm)
#println(A)
# println(B)
# println(C)


Cm = Extraire_Cm(C,NbrObstacle,Obstacle)
# println(Cm)



Image_Mulit(Obstacle,Cm,Dm,NbrObstacle)
