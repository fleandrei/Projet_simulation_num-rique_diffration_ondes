#module Polar

# a partir de l'indice d'une matrice donne les coordonnées cartésiennes
function coordonnees(i,j,h,taille_espace) 
	x=i*h-taille_espace/2
	y=-j*h+taille_espace/2
	return(x,y)
end

# convertit les cartésiennes en polaires
function conversion_polaire(x,y)
	r=sqrt(x*x+y*y)
	if (x > 0 && y >= 0)
		lambda = atan(y/x)

	elseif(x > 0 && y < 0)
		lambda = atan(y/x) + 2*pi

	elseif(x < 0)
		lambda = atan(y/x) + pi

	elseif(x == 0 && y > 0)
		lambda = pi/2

	else
		lambda = 3*pi /2
	    
	end

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

function distance(x1,y1,x2,y2)
	d=sqrt((x2-x1)^(2)+(y2-y1)^(2))
	return(d)
end

function angle(x1,y1,x2,y2)
	#alpha=atan(abs(y2-y1)/abs(x2-x1))
	alpha=atan(y2-y1/x2-x1)
	return(alpha)
end

	#export coordonnees

#end

push!(LOAD_PATH, pwd())