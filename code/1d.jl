using Printf
using PyPlot
using SpecialFunctions
using LinearAlgebra 

include("./polar.jl")
include("./diffraction.jl")
#using Polar
#import Polar

########### VARIABLE #############
k = 2*pi
lambda = 2*pi/k # longueur d'onde
h = lambda/60  #pas de la grille
a = 1 # rayon de l'obstacle
taille_espace=6 # taille total de la grille
taille_matrice=convert(Int64, taille_espace*1/h)


alpha=pi #Angle de l'onde incidente
alphap=0 #Angle de l'obstacle
bp=0
rp=1
e=10^(-12)
Np= floor(Int64,k*a + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*a)^(1/3) +1)


########## FONCTION ##############

# a partir de l'indice d'une matrice donne les coordonnées cartésiennes
# function coordonnees(i,j,h,taille_espace) 
# 	x=i*h-taille_espace/2
# 	y=-j*h+taille_espace/2
# 	return(x,y)
# end

# # convertit les cartésiennes en polaires
# function conversion_polaire(x,y)
# 	r=sqrt(x*x+y*y)
# 	if (x > 0 && y >= 0)
# 		lambda = atan(y/x)

# 	elseif(x > 0 && y < 0)
# 		lambda = atan(y/x) + 2*pi

# 	elseif(x < 0)
# 		lambda = atan(y/x) + pi

# 	elseif(x == 0 && y > 0)
# 		lambda = pi/2

# 	else
# 		lambda = 3*pi /2
	    
# 	end

# 	return(r,lambda)
# end

# # convertit les polaires en cartésiennes
# function conversion_cart(r,lambda)
# 	x=r*cos(lambda)
# 	y=r*sin(lambda)
# 	return(x,y)
# end

# # convertit les cartésiens en coordonnées matricielles
# function coordonnees_from_cart(x,y,h,taille_espace)
# 	i=convert(Int64,(x+taille_espace/2)/h)
# 	j=convert(Int64,-(y-taille-espace/2)/h)
# 	return(i,j)
# end


# # Calcule l'onde diffractée complexe U
# function calculUp(r,teta, Cm, Np)
# 	#N=length(Cm)
# 	U=0.0
# 	for m=1:2*Np + 1
# 		U= U + Cm[m]*besselh(m-Np,k*r)*exp(im*(m-Np)*teta)
# 	end	

# 	return U
# end


# # Calcule l'onde incidente complexe U
# function calculUinc(r,teta, Dm, Np)
# 	#N=length(Cm)
# 	U=0.0
# 	for m=1:2*Np + 1
# 		U= U + Dm[m]*besselj(m-Np,k*r)*exp(im*(m-Np)*teta)
# 	end	

# 	return U
# end


# function calculRCS(U)
# 	return 10*log10(2*pi*abs(U)^2)
# end



# function CDH(bp,alphap,alpha,Np, k) #Calcule Cm: coef Fourrier onde réfléchie, Dm: Coef Fourrier onde incidente, Hm: vecteur Besselh
# Cm=ones(Complex{Float32},2*Np+1,1)
# Dm=ones(Complex{Float32},2*Np+1,1)
# Hm=ones(Complex{Float32},2*Np+1,1)
# temp=exp(im*k*cos(alpha-alphap)*bp)

# for m=1:2*Np+1
# 	Dm[m]=temp*exp(im*(pi*0.5 - alpha)*m -Np)
# 	Hm[m]=besselh(m-Np,a*k)
# 	Cm[m]=besselj(m-Np,a*k) * Hm[m]^(-1) * Dm[m]

# end

# return (Cm,Dm,Hm)
# end


# function CoeFourrier_OndeInc(bp,alphap,alpha,Np, k)
# Dm=ones(Complex{Float32},2*Np+1,1)
# temp=exp(im*k*cos(alpha-alphap)*bp)
# for m=1:2*Np+1
# 	Dm[m]=temp*exp(im*(pi*0.5 - alpha)*m -Np)
# end

# return Dm
# end


# function BesselHVector(Np,a,k)
# Hm=ones(Complex{Float32},2*Np+1,1)
# temp=a*k
# for m=1:2*Np+1
# 	Hm[m]=besselh(m-Np,temp)
# end

# return Hm
# end


# function Calcule_b(bp,alphap,alpha,Np, k, a) 
# Dm=CoeFourrier_OndeInc(bp,alphap,alpha,Np, k)
# b=ones(Complex{Float32},2*Np+1,1)
# temp=a*k
# for m=1:2*Np+1
# 	b[m]=besselj(m-Np, temp)*Dm[m]
# end
# return b
# end


# function CoeFourrier_OndeRefracte(bp,alphap,alpha,Np, k, a)
# Cm=ones(Complex{Float32},2*Np+1,1)
# b=Calcule_b(bp,alphap,alpha,Np, k, a)
# Hm=BesselHVector(Np,a,k)
# #A=zeros(Complex{Float32},2*Np+1, 2*Np+1)
# A=Matrix(I,2*Np+1, 2*Np+1)
# A=A.*Hm
# Cm=A\b

# return Cm
# end
########### CODE #############


#Cm,Dm,Hm=CDH(bp,alphap,alpha,Np,k)
Cm=CoeFourrier_OndeRefracte(bp,alphap,alpha,Np, k, a)
Dm=CoeFourrier_OndeInc(bp,alphap,alpha,Np, k)

# declaration de la matrice
M = ones(Float64,taille_matrice, taille_matrice)


# Parcoure et remplissage de la matrice
for i = 1:taille_matrice
	for j=1:taille_matrice
		x,y=coordonnees(i,j,h,taille_espace)
		r,lambda=conversion_polaire(x,y)
		#println("x=",x,"  y=",y,"  r=",r)
		if r < a
			#println(" r=",r)
			M[i,j]=0
		else
			M[i,j]=abs(calculUp(r,lambda, Cm, Np)+calculUinc(r,lambda, Dm, Np))
			#println(" autre=",r)
		end
	end
end

# Affichage graphique

#figimage(M)
imshow(M, extent=(-3, 3, -3, 3))
savefig("res.svg")












