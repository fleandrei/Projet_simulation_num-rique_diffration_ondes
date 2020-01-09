from math import *
import cmath as cm
import numpy as np
import matplotlib.pyplot as plt
from scipy import special

########### VARIABLE #############
k = 2*pi
lbda = 2*pi/k # longueur d'onde
h = lbda/30  #pas de la grille
alpha = pi # angle de l'onde incidente
a = 1 # rayon de l'obstacle
taille_espace=10 # taille total de la grille
taille_grille=int(taille_espace*1/h)


beta=pi #Angle de l'onde incidente
alphap=0 #Angle de l'obstacle
bp=0
rp=1
eps=10**(-12)
M = 2





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


def distance(x1,y1,x2,y2):
	d=sqrt((x2-x1)**2+(y2-y1)**2)
	return(d)

def angle_2p(x1,y1,x2,y2):
	
	x = x2 - x1
	y = y2 - y1

	
	alpha = atan(y,x)
	# alpha = angle(x+y) 
	return(alpha)


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

# Calcule l'onde incidente complexe U
def calculUinc_exact(x,y,Beta,k):

	prod_scal = x*cos(Beta) + y*sin(Beta)
	U = cm.exp(1j*k*prod_scal)

	return (U)


def calculRCS(theta):
	return (10*log10(2*pi*abs(theta)^2))



def CDH(bp,alphap,alpha,Np, k): #Calcule Cm: coef Fourrier onde réfléchie, Dm: Coef Fourrier onde incidente, Hm: vecteur Besselh
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

	return (Cm,Dm,Hm)


def CoeFourrier_OndeInc(bp,alphap,alpha,Np, k):

	Dm=np.zeros((2*Np+1,1), dtype=np.complex_)
	temp=cm.exp(1j*k*cm.cos(alpha-alphap)*bp)

	for m in range(-Np, Np+1):
		Dm[m+Np]=temp*cm.exp(1j*(pi*0.5 - alpha)*m)
	

	return Dm



def BesselHVector(Np,a,k):

	Hm=np.zeros((2*Np+1,1), dtype=np.complex_)

	for m in range(-Np, Np+1):
		Hm[m+Np] = besselh(m, a*k)
	
	return Hm



def Calcule_b(bp,alphap,alpha,Np, k, a):

	Dm   = CoeFourrier_OndeInc(bp, alphap, alpha, Np, k)
	b    = np.zeros((2*Np+1,1), dtype=np.complex_)
	temp = a * k

	for m in range(-Np, Np+1):

		b[m + Np] = - special.jv(m, temp) * Dm[m + Np]
	

	return b




def CoeFourrier_OndeDiffracte(bp,alphap,alpha,Np, k, a):

	Cm = np.zeros((2*Np+1,1), dtype=np.complex_)
	b  = Calcule_b(bp, alphap, alpha, Np, k, a)
	Hm = BesselHVector(Np, a, k)
	A  = np.eye(2*Np+1, dtype=np.complex_)
	for i in range(2*Np+1):
		for j in range(2*Np+1):
			A[i, j] *= Hm[i]

	Cm = np.linalg.solve(A, b)

	print(Cm)

	return Cm



def Phi(m,rp, teta,k):
	return special.hankel1(m, k*rp)*cm.exp(1j*m*teta)



def Phi_Chap(m,rp,teta,k):
	return special.jv(m, rp*k)*cm.exp(1j*m*teta)



def Smn(m,n,b,teta, k):
	return Phi(m-n,b, teta,k)




def dmp(p, Obstacle, Beta,k): # Calcule les coeff Fourrier de l'onde incidente pour la boule p
	
	Np      = Obstacle[p][3]
	Np      = floor(Np) #convert to int
	print(p)
	print(Np)
	Dm      = np.zeros((2*Np+1,1), dtype=np.complex_)
	bp,teta = conversion_polaire(Obstacle[p][0], Obstacle[p][1])
	temp    = cm.exp(1j*k * cos(Beta-teta) * bp)

	for m in range(-Np, Np+1):

		Dm[m + Np] = temp * cm.exp(1j * (pi*0.5 - Beta)* m)
	

	return Dm



def Extraire_Dm(M,Obstacle,Beta,k): # Renvoie Tableau de tableau Dm tel que Dm[P]="coeff Fourrier onde incidente pour la boule P"
	Dm = [[] for i in range(M)]

	for i in range(M):
		Dm[i].append(dmp(i,Obstacle,Beta,k))
	

	return Dm



def Calcule_B(M, Obstacle, Beta,k,Dm): #Calcule le vecteur b du système "Ac=b"

	N = 0 
	for i in range(M):
		N += 2*Obstacle[i][3] + 1
	

	N = floor(N) #convert to int
	print(N)


	B = np.zeros((N,1), dtype=np.complex_)
	for i in range(M):

		Np = Obstacle[i][3]
		dm = Dm[i][0]
		ap = Obstacle[i][2]

		#print(size(dm))   ### "size(dm)" ? ###

		Np = floor(Np)
		#dm = floor(dm)
		#ap = floor(ap)

		
		for m in range(-Np, Np+1):
			if(i > 1):
				Np_prec = Obstacle[i-1][3]
			else:
				Np_prec = 0 #arbritraire
			
			Np_prec = floor(Np_prec) #convert to int

			temp   = m + Np
			idx    = temp + i*(2*Np_prec +1)
			print(str(temp) + "," + str(idx))
			B[idx] = -(special.jv(m,ap*k)/special.hankel1(m,k*ap)) * dm[temp]
		
	

	return B



def Matrix_S(Np,Nq,b,teta,Obstacle,k):
	
	S = np.zeros((2*Nq+1,2*Np+1), dtype=np.complex_)

	for n in range(2*Nq+1):
		for m in range(2*Np+1):
				S[n,m] = Smn(n,m,b,teta,k)
		
	
	return S


def Matrix_D(Np,k,ap):
	D = np.zeros((2*Np+1,2*Np+1), dtype=np.complex_)

	#N = min(Np,Nq)
	for m in range(2*Np+1):
		D[m,m] = special.jv(m,ap*k)/special.hankel1(m, ap*k)
	

	return D





def Apq(p,q,k, Obstacle): #Calcule la sous-matrice d'indices p,q de la matrice A 
	
	Np = Obstacle[p][3]
	Nq = Obstacle[q][3]

	Np = floor(Np)
	Nq = floor(Nq)

	if (p == q):
		A = np.eye(2*Np+1, dtype=np.complex_)
		#print(A)
		return A
		
	else:

		xp = Obstacle[p][0]
		yp = Obstacle[p][1]
		ap = Obstacle[p][2]

		xq = Obstacle[q][0]
		yq = Obstacle[q][1]

		xp = floor(xp)
		yp = floor(yp)
		xq = floor(xq)
		yq = floor(yq)

		b    = distance(xp, yp, xq, yq)
		teta = angle_2p(xp, yp, xq, yq)
		A    = np.zeros((2*Np+1,2*Nq+1), dtype=np.complex_)

		D = Matrix_D(Np,k,ap)
		S = Matrix_S(Np,Nq,b,teta,Obstacle,k)

		print("----------")
		print(D.shape)
		print(S.transpose().shape)
		print("----------")

		A = D * S.transpose()

		

		# for m =1:2*Np+1
		# 	D = besselj(m,ap*k)/besselh(m, ap*k)
		# 	for n= 1:2*Nq+1
		# 			A[m,n] = D * Smn(m,n,b,teta,k)
		# 	end
		# end
		# println(A,"\n b = ",b," teta = ",teta,"\n")

		print("----------")
		print(A.shape)
		print("----------")

		return A
		



def Calcule_A(M, Obstacle, k): #Calcule la matrice A du système "Ac=b"
	
	
	A=Apq(1,1,k, Obstacle) #Correspond à la première boucle de for j=1
	

	for j in range(1, M):

		A=np.concatenate((A, Apq(1,j,k, Obstacle)), axis=1)
	
	
	
	#println("\n")

	for i in range(1, M):

		Al=Apq(i,1,k, Obstacle) #Correspond à la première boucle de for j=1
		for j in range(1, M):

			Al=np.concatenate((Al, Apq(i,j,k, Obstacle)), axis=1)
			
		
		#println(Al)
		#println("\n")
		A=np.concatenate((A,Al), axis=0)
	
	

	#println(A,"\n")

	return A



def Calcule_C(A,b): #Calcule le vecteur c du système: "Ac=b".  Il s'agit du vecteur des coeff de Fourrier des ondes diffractées
	return np.linalg.solve(A, b)



def Extraire_Cm(C,M,Obstacle): #A partir du vecteur C du système 
	Cm=[[] for i in range(M)]

	curseur=0
	for i in range(M):
		
		Np=Obstacle[i][3]
		Np = floor(Np)

		Cm[i].append(C[curseur:curseur+2*Np])
		curseur=curseur+2*Np+1
		
	
	print(Cm)

	return Cm



def Boule_Proche(Obstacle,x, y,M):
	
	MinDist=1000000
	p=0
	for i in range(M):
		Dist=distance(x,y,Obstacle[i][0], Obstacle[i][1])
		if (Dist<MinDist):
			MinDist=Dist
			p=i
		
	
	return p



# function CalculeUq(Obstacle, r, teta, k, p,Cm,M)
	
# 	Np=Obstacle[p][4]

# 	somme_m=0

# 	for m=-Np:Np
# 		somme_q=0
# 		for q=1:M
# 			if q != p
# 				Nq=Obstacle[q][4]
# 				somme_n=0
# 				x1=Obstacle[p][1]
# 				y1=Obstacle[p][2]
# 				x2=Obstacle[q][1]
# 				y2=Obstacle[q][2]
# 				bpq=distance(x1,y1,x2,y2)
# 				angle_2ppq=angle_2p(x1,y1,x2,y2)
# 				Cq=Cm[q][1]
# 				for n=-Nq:Nq
# 					somme_n=somme_n+ Smn(n,m,bpq,angle_2ppq, k)*Cq[n+Nq+1]
# 				end
# 				somme_q=somme_q + somme_n
# 			end
# 		end
# 		somme_m= somme_m +	somme_q * Phi_Chap(m,r,teta,k)

# 	end

# 	return somme_m

# end


def CalculeUq(Obstacle, x,y, k,Cm,M):

	somme_q = 0

	for q in range(M):
		#somme_n=0

		Nq = Obstacle[q][3]
		Nq = floor(Nq)
		
		x1 = Obstacle[q][0]
		y1 = Obstacle[q][1]
		x2 = x
		y2 = y

		x1 = floor(x1)
		y1 = floor(y1)

		b  = distance(x1,y1,x2,y2)
		angle_2ppq = angle_2p(x1,y1,x2,y2)
		Cq = Cm[q, 1]
		#println(Cq,"\n")

		somme_q = somme_q + calculUp(b, angle_2ppq, Cq, Nq)
	

	return somme_q




# function Calcule_UDiff_MultiDisk(Obstacle,x,y,Cm,k,M)

# 	#           m m Up = calculUp(r,teta, Cm[p][1], Np)
# 	Uq = CalculeUq(Obstacle, x, y, k, Cm, M)
# 	U  = Uq 
# 	# println("Up = ",Up, " , Uq = ",Uq, ", U = ",U,"\n")
# 	return U
# end

def Calcule_Utot_MultiDisk(Obstacle,x,y,Cm,Dm,k,M,Beta):

	# p     = Boule_Proche(Obstacle,x,y,M)
	# Np    = Obstacle[p][4]
	# xp    = Obstacle[p][1]
	# yp    = Obstacle[p][2]

	# r     = distance(xp, yp, x, y)
	# teta  = angle_2p(xp, yp, x, y)

	# Uinc  = calculUinc(r,teta, Dm[p][1], Np)
	Uinc  = calculUinc_exact(x,y,Beta,k)

	Udiff = CalculeUq(Obstacle, x, y, k, Cm, M)
	

	Utot = abs( Uinc + Udiff )
	# Utot = real(Uinc)
	# Utot = real(Udiff)

	return Utot



def Is_inDisk(x,y,Obstacle,M):
	for i in range(M):

		if (distance(x,y,Obstacle[i][0], Obstacle[i][1]) <= Obstacle[i][2]):
			return true

	return false




def Image_Mulit(obstacle,Cm,Dm, M,Beta):

# declaration de la matrice
	Image = np.zeros((taille_matrice, taille_matrice))
	
	# Parcoure et remplissage de la matrice
	for i in range(taille_matrice):

		for j in range(taille_matrice):

			x,y    = coordonnees(i, j, h, taille_espace)
			r,lbda = conversion_polaire(x, y)

			if (Is_inDisk(x,y,Obstacle, M) == False):

				Image[i, j] = Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, M,Beta)
				#println("Image [",i,",",j,"] = ", Image[i,j],"\n")


		print("Image [" + str(i) + "]")

	
	# Affichage graphique
	
	plt.imshow(transpose(Image), vmin=-2.5, vmax=2.5, extent=(-5, 5, -5, 5))
	plt.colorbar()
	plt.savefig("resMult.svg")



########### CODE #############


Np= floor(k*a + (np.cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*eps))) ** 2 * np.cbrt(k*a) +1)

Obstacle=[[0,0,1,Np], [4,4,0.01,1]]
print(Np)

Dm = Extraire_Dm(M, Obstacle, beta, k)
B  = Calcule_B(M ,Obstacle, beta,k,Dm)
A  = Calcule_A(M, Obstacle, k)
C  = Calcule_C(A,B)


# println(Dm)
#println(A)
# println(B)
# println(C)


Cm = Extraire_Cm(C,M,Obstacle)
# println(Cm)



Image_Mulit(Obstacle,Cm,Dm,M,beta)

