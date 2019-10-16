using Printf
using Gadfly
using SpecialFunctions

########### VARIABLE #############
k = 2*pi
lambda = 2*pi/k # longueur d'onde
h = lambda/50  #pas de la grille
alpha = pi # angle de l'onde incidente
a = 1 # rayon de l'obstacle
taille_espace=6 # taille total de la grille
taille_matrice=convert(Int64, taille_espace*1/h)


alpha=pi #Angle de l'onde incidente
alphap=0 #Angle de l'obstacle
bp=0
a=1
rp=1
e=10^(-12)
Np= floor(Int64,k*a + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*a)^(1/3) +1)


########## FONCTION ##############

# a partir de l'indice d'une matrice donne les coordonnées cartésiennes
function coordonnees(i,j,h,taille_espace) 
	x=i*h-taille_espace/2
	y=-j*h+taille_espace/2
	return(x,y)
end

# convertit les cartésiennes en polaires
function conversion_polaire(x,y)
	r=sqrt(x*x+y*y)
	lambda=atan(abs(y/x))
	return(r,lambda)
end

# convertit les polaires en cartésiennes
function conversion_cart(r,lambda)
	x=r*cos(lambda)
	y=r*sin(lambda)
	return(x,y)
end

# convertit les cartésiens en coordonnées matricielles
function coordonnees_from_cart(x,y,h,taille_espace)
	i=convert(Int64,(x+taille_espace/2)/h)
	j=convert(Int64,-(y-taille-espace/2)/h)
	return(i,j)
end


function calculU(r,teta, Cm)
	N=length(Cm)
	U=0.0
	for m=1:N
		U= U + Cm[m]*besselh(m-N/2,k*r)*exp(im*(m-N/2)*teta)
	end	

	return U
end




########### CODE #############


Cm=ones(Complex{Float32},2*Np+1,1)
Dm=ones(Complex{Float32},2*Np+1,1)
Hm=ones(Complex{Float32},2*Np+1,1)
temp=exp(im*k*cos(alpha-alphap)*bp)
u=0.0

for m=1:2*Np+1
	Dm[m]=temp*exp(im*(pi*0.5 - alpha)*m)
	Hm[m]=besselh(m-Np,a*k)
	Cm[m]=besselj(m-Np,a*k) * Hm[m]^(-1) * Dm[m]

end



# declaration de la matrice
M = ones(Float64,taille_matrice, taille_matrice)


# Parcoure et remplissage de la matrice
for i = 1:taille_matrice
	for j=1:taille_matrice
		x,y=coordonnees(i,j,h,taille_espace)
		r,lambda=conversion_polaire(x,y)
		#println("x=",x,"  y=",y,"  r=",r)
		if r<1
			M[i,j]=0
		else
			M[i,j]=real(calculU(r,lambda, Cm))
		end
	end
end

# Affichage graphique
p = spy(M)
img = SVG("output.svg", 8inch, 8inch)
draw(img, p)










