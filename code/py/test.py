from math import *
import cmath as cm
import numpy as np
#import matplotlib.pyplot as plt
from scipy import special

########### VARIABLE #############
k = 2*pi
lbda = 2*pi/k # longueur d'onde
h = lbda/60  #pas de la grille
alpha = pi # angle de l'onde incidente
a = 1 # rayon de l'obstacle
taille_espace=6 # taille total de la grille
taille_grille=int(taille_espace*1/h)


alpha=pi #Angle de l'onde incidente
alphap=0 #Angle de l'obstacle
bp=0
a=1
rp=1
eps=10**(-12)
Np= floor(k*a + (np.cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*eps))) ** 2 * np.cbrt(k*a) +1)


########## FONCTION ##############

# a partir de l'indice d'une matrice donne les coordonnées cartésiennes
def coordonnees(i,j,h,taille_espace):
	x=i*h-taille_espace/2
	y=-j*h+taille_espace/2
	return(x,y)

# convertit les cartésiennes en polaires
def conversion_polaire(x,y):
	r=sqrt(x*x+y*y)
	if (x > 0 and y >= 0):
		lbda = atan(y/x)

	elif(x > 0 and y < 0):
		lbda = atan(y/x) + 2*pi

	elif(x < 0):
		lbda = atan(y/x) + pi

	elif(x == 0 and y > 0):
		lbda = pi/2

	else:
		lbda = 3*pi /2


	return(r,lbda)

# convertit les polaires en cartésiennes
def conversion_cart(r,lbda):
	x=r*cos(lbda)
	y=r*sin(lbda)
	return(x,y)

# convertit les cartésiens en coordonnées matricielles
def coordonnees_from_cart(x,y,h,taille_espace):
	i=convert(Int64,(x+taille_espace/2)/h)
	j=convert(Int64,-(y-taille-espace/2)/h)
	return(i,j)


# Calcule l'onde diffractée complexe U
def calculUp(r,teta, Cm, Np):
	#N=length(Cm)
	U=0.0
	for m in range(-Np, Np+1):
		U += Cm[m+Np] * special.hankel1(m, k*r) * cm.exp(1j*m*teta)

	return (U)


# Calcule l'onde incidente complexe U
def calculUinc(r,teta, Dm, Np):
	#N=length(Cm)
	U=0.0
	for m in range(-Np, Np+1):
		U += Dm[m+Np] * special.jv(m, k*r) * cm.exp(1j*m*teta)

	return (U)


def calculRCS(theta):
	return (10*log10(2*pi*abs(theta)^2))





########### CODE #############


Cm=np.zeros((2*Np+1,1), dtype=np.complex_)
Dm=np.zeros((2*Np+1,1), dtype=np.complex_)
Hm=np.zeros((2*Np+1,1), dtype=np.complex_)
temp=cm.exp(1j*k*cm.cos(alpha-alphap)*bp)
#print(temp)
u=0.0

for m in range(-Np, Np+1):
	Dm[m+Np]=temp*cm.exp(1j*(pi*0.5 - alpha)*m)
	Hm[m+Np]=special.hankel1(m,a*k)
	Cm[m+Np]=special.jv(m,a*k) * Dm[m+Np] / Hm[m+Np]




# declaration de la grille
M = np.zeros((taille_grille, taille_grille))
taille_mrcs = floor(2*pi*10)
MRCS = np.zeros(taille_mrcs)


# Parcours et remplissage de la grille
for i in range(taille_grille):
	for j in range(taille_grille):
		(x,y)=coordonnees(i,j,h,taille_espace)
		(r,lbda)=conversion_polaire(x,y)
		#print("x=",x,"  y=",y,"  r=",r)
		if (r >= a):
			M[i,j]=abs(calculUp(r, lbda, Cm, Np) + calculUinc(r, lbda, Dm, Np))
			#print(" autre=",r)

# Parcours et remplissage de MRCS
for i in range(taille_mrcs):
	MRCS[i] = calculRCS(i/10)


# Affichage graphique

#plt.imshow(M, extent=(-3, 3, -3, 3))
#plt.savefig("res.svg")

# Affichage RCS

plt.imshow(MRCS)
plt.savefig("rcs.svg")