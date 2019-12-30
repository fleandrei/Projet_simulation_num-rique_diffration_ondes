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
taille_espace  = 10 # taille total de la grille
taille_matrice = convert(Int64, taille_espace*1/h)

beta = pi #Angle de l'onde incidente
e    = 10^(-12)
M = 2
	
#####Function########	
function Image_Mulit(obstacle,Cm,Dm, M,Beta)

# declaration de la matrice
	Image = zeros(Float64, taille_matrice, taille_matrice)
	
	# Parcoure et remplissage de la matrice
	for i = 1:taille_matrice

		for j = 1:taille_matrice

			x,y      = coordonnees(i, j, h, taille_espace)
			r,lambda = conversion_polaire(x, y)

			if !Is_inDisk(x,y,Obstacle, M)

				Image[i,j] = Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, M,Beta)
				#println("Image [",i,",",j,"] = ", Image[i,j],"\n")
			end
		end
		println("Image [",i,"]","\n")
	end
	
	# Affichage graphique

	scale_min = - taille_espace/2
	scale_max =   taille_espace/2
	
	# imshow(transpose(Image), vmin=-2.5, vmax=2.5, extent=(scale_min, scale_max, scale_min, scale_max))
	imshow(transpose(Image), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("resMult.svg")
end




##########CODE############

Np = floor(Int64, k*1 + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*1)^(1/3) +1)



Obstacle=[[0,0,1,Np], [4,4,0.01,1]]
println(Np,"\n")

Dm = Extraire_Dm(M, Obstacle, beta, k)
B  = Calcule_B(M ,Obstacle, beta,k,Dm)
A  = Calcule_A(M, Obstacle, k)
C  = Calcule_C(A,B)

# println("----- C -----:\n")
# println(C)

println("----- A -----:\n")
println(A)

println("----- B -----:\n")
println(B)


Cm = Extraire_Cm(C,M,Obstacle)
println("----- Cm -----:\n")
println(Cm)



Image_Mulit(Obstacle,Cm,Dm,M,beta)
