#module Polar
import Base.angle

# a partir de l'indice d'une matrice donne les coordonnées cartésiennes
function coordonnees(i,j,h,taille_espace) 
	x=i*h-taille_espace/2
	y=-j*h+taille_espace/2
	return(x,y)
end

# convertit les cartésiennes en polaires
function conversion_polaire(x,y)
	r = sqrt(x*x + y*y)
	
	# lambda = angle(x+y)
	lambda = atan(y,x) 

	return(r,lambda)
end

# convertit les polaires en cartésiennes
function conversion_cart(r,lambda)
	x = r*cos(lambda)
	y = r*sin(lambda)
	return(y,x)
end

# convertit les cartésiens en coordonnées matricielles
function coordonnees_from_cart(x,y,h,taille_espace)
	i = convert(Int64, (x+taille_espace/2)/h)
	j = convert(Int64,-(y-taille-espace/2)/h)
	return(i,j)
end

function distance(x1,y1,x2,y2)
	d = sqrt((x2-x1)^(2) + (y2-y1)^(2))
	return(d)
end

function angle_2p(x1,y1,x2,y2)
	
	x = x2 - x1
	y = y2 - y1

	
	alpha = atan(y,x)
	# alpha = angle(x+y) 
	return alpha
end

	#export coordonnees

#end

function Rad2Deg(Rad)
	return Rad*360/(2*pi)
end

function  Deg2Rad(Deg)
	return Deg*2*pi/360
end

push!(LOAD_PATH, pwd())