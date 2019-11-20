using Printf
using PyPlot
using SpecialFunctions
using LinearAlgebra 

include("./polar.jl")
include("./diffraction.jl")
#using Polar
#import Polar

k = 2*pi
lambda = 2*pi/k # longueur d'onde
h = lambda/60  #pas de la grille
taille_espace=6 # taille total de la grille
taille_matrice=convert(Int64, taille_espace*1/h)

beta=pi #Angle de l'onde incidente
e=10^(-12)
M=2
	
Np=floor(Int64, k*1 + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*1)^(1/3) +1)



Obstacle=[[2,0,1,Np], [0,2,1,Np]]


B=Calcule_B(M ,Obstacle, beta,k)
A=Calcule_A(M, Obstacle, k)
C=Calcule_C(A,B)

println(C)

Cm=Extraire_Cm(C,M,Obstacle)
println(Cm[2])
