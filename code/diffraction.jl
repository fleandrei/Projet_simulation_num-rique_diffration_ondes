include("./polar.jl")
using LinearAlgebra

# Calcule l'onde diffractée complexe U
function calculUp(r,teta, Cm, Np)

	U = 0.0

	for m = -Np:Np

		idx = m + Np + 1
		U   = U + (Cm[idx] * Phi(m,r, teta,k))
	end	



	# ordre = -14
 #    U   =  Cm[ordre+Np+1] * Phi(ordre,r, teta,k)
 #    if(real(U) > 2)
 #    	println("U = ",U,"Ordre ",ordre," : c=", Cm[ordre+Np+1],"; phi = ",Phi(ordre,r, teta,k), "\n")
 #    end
	return U
end


# Calcule l'onde incidente complexe U
function calculUinc(r,teta, Dm, Np)

	U = 0.0
	for m = -Np:Np

		idx = m + Np + 1
		U   = U + Dm[idx] * Phi_Chap(m,r,teta,k)
	end	

	return U
end

# Calcule l'onde incidente complexe U
function calculUinc_exact(x,y,Beta,k)

	prod_scal = x*cos(Beta) + y*sin(Beta)
	U = exp(im*k*prod_scal)

	return U
end


function calculRCS(U)
	return 10 * log10(2*pi * abs(U)^2)
end


function CDH(bp,alphap,alpha,Np, k) #Calcule Cm: coef Fourrier onde réfléchie, Dm: Coef Fourrier onde incidente, Hm: vecteur hankelh1
	Cm   = zeros(Complex{Float64}, 2*Np + 1, 1)
	Dm   = zeros(Complex{Float64}, 2*Np + 1, 1)
	Hm   = zeros(Complex{Float64}, 2*Np + 1, 1)
	temp = exp(im*k * cos(alpha-alphap) * bp)

	for m = -Np:Np

		idx     = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - alpha) * m)
		Hm[idx] = hankelh1( m, a*k )
		Cm[idx] = besselj( m, a*k ) * Hm[idx]^(-1) * Dm[idx]

	end

	return (Cm,Dm,Hm)
end


function CoeFourrier_OndeInc(bp,alphap,alpha,Np, k) #cas 1d

	Dm   = zeros(Complex{Float64}, 2*Np +1, 1)
	temp = exp(im*k * cos(alpha-alphap) * bp)

	for m = -Np:Np

		idx     = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - alpha)* m)
	end

	return Dm
end


# function hankelh1Vector(Np,a,k) #cas 1d #inutile

# 	Hm   = zeros(Complex{Float64},2*Np+1,1)
# 	temp = a * k

# 	for m = -Np:Np

# 		idx     = m + Np + 1
# 		Hm[idx] = hankelh1(m, temp)
# 	end

# 	return Hm
# end


function Calcule_b(bp,alphap,alpha,Np, k, a) #cas 1d

	Dm   = CoeFourrier_OndeInc(bp, alphap, alpha, Np, k)
	b    = zeros(Complex{Float64}, 2*Np + 1, 1)
	temp = a * k

	for m = -Np:Np

		idx    = m + Np + 1
		b[idx] = - besselj(m, temp)/hankelh1(m,temp) * Dm[idx]
	end

	return b
end



function CoeFourrier_OndeDiffracte(bp,alphap,alpha,Np, k, a) # cas 1d

	Cm = zeros(Complex{Float64},2*Np+1,1)
	b  = Calcule_b(bp, alphap, alpha, Np, k, a)
	# Hm = hankelh1Vector(Np, a, k)
	A  = Matrix(I, 2*Np+1, 2*Np+1)
	# A  = A.*Hm
	Cm = A\b

	println("----- A -----:\n")
	println(A)

	println("----- B -----:\n")
	println(b)

	println("----- Cm -----:\n")
	println(Cm)

	return Cm
end


function Phi(m,rp, teta,k)
	# println("Ordre ",m," : Hankel=", hankelh1(m, k*rp),"; exp = ",exp(im*m*teta), "\n")
	return hankelh1(m, k*rp)*exp(im*m*teta)
end


function Phi_Chap(m,rp,teta,k)
	return besselj(m, rp*k)*exp(im*m*teta)
end


function Smn(m,n,b,teta, k)
	return Phi(m-n,b, teta,k)
end



function dmp(p, Obstacle, Beta,k) # Calcule les coeff Fourrier de l'onde incidente pour la boule p
	
	Np      = Obstacle[p][4]
	Np      = floor(Int, Np) #convert to int
	# println(p)
	# println(Np)
	Dm      = zeros(Complex{Float64}, 2*Np +1, 1)
	bp,teta = conversion_polaire(Obstacle[p][1], Obstacle[p][2])
	temp    = exp(im*k * cos(Beta-teta) * bp)

	for m = -Np:Np

		idx     = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - Beta)* m)
	end

	return Dm
end


function Extraire_Dm(M,Obstacle,Beta,k) # Renvoie Tableau de tableau Dm tel que Dm[P]="coeff Fourrier onde incidente pour la boule P"
	Dm = [[] for i=1:M]

	for i=1:M
		push!(Dm[i], dmp(i,Obstacle,Beta,k))
	end

	return Dm
end


function construct_Bp(p,Obstacle,k,dm)

	Np = Obstacle[p][4]
	Np = floor(Int, Np)
	ap = Obstacle[p][3]

	Bp = zeros(Complex{Float64}, 2*Np +1, 1)

	for m = -Np:Np
		idx = m + Np +1
		Bp[idx]  = -(besselj(m,ap*k)/hankelh1(m,k*ap))* dm[idx]
	end

	return Bp
end

function Calcule_B(M, Obstacle, Beta,k,Dm)
	#iteration 1
		B = construct_Bp(1,Obstacle,k,Dm[1][1])


	for i = 2:M
		dm = Dm[i][1]
		B_p = construct_Bp(i,Obstacle,k,dm)

		B = vcat(B,B_p)
	end

	return B
end


function Matrix_S(Np,Nq,b,teta,k)
	
	S = zeros(Complex{Float64}, 2*Nq +1, 2*Np +1)

	for n = 1:2*Nq+1
		for m = 1:2*Np+1
				S[n,m] = Smn(m-(Np+1),n-(Nq+1),b,teta, k)
		end
	end
	return S
end

function Matrix_D(Np,k,ap)
	D = zeros(Complex{Float64}, 2*Np +1, 2*Np +1)

	#N = min(Np,Nq)
	for m = 1:2*Np+1
		D[m,m] = besselj(m-(Np+1),ap*k)/hankelh1(m-(Np+1), ap*k)
	end

	return D
end




function Apq(p,q,k, Obstacle) #Calcule la sous-matrice d'indices p,q de la matrice A 
	
	Np = Obstacle[p][4]
	Nq = Obstacle[q][4]

	Np = floor(Int, Np)
	Nq = floor(Int, Nq)

	if p == q
		A = 1*Matrix(I,2*Np+1,2*Np+1)
		#println(A,"\n")
		return A
		
	else

		xp = Obstacle[p][1]
		yp = Obstacle[p][2]
		ap = Obstacle[p][3]

		xq = Obstacle[q][1]
		yq = Obstacle[q][2]

		xp = floor(Int, xp)
		yp = floor(Int, yp)
		xq = floor(Int, xq)
		yq = floor(Int, yq)

		b    = distance(xp, yp, xq, yq)
		teta = angle_2p(xp, yp, xq, yq)
		A    = zeros(Complex{Float64}, 2*Np+1, 2*Nq+1)

		println("----- angle (",p,",",q,") ----------")
		println(teta)
		println("-------------------------------")

		D = Matrix_D(Np,k,ap)
		S = Matrix_S(Np,Nq,b,teta,k)

		println("----------\n")
		#  println(D)
		println(size(transpose(S)),"\n")
		println("----------\n")

		A = D * transpose(S)

		

		println("----------\n")
		println(size(A),"\n")
		println("----------\n")

		return A
	end	
end


function Calcule_A(M, Obstacle, k) #Calcule la matrice A du système "Ac=b"
	
	
	A=Apq(1,1,k, Obstacle) #Correspond à la première boucle de for j=1
	

	for j = 2:M

		A=hcat(A, Apq(1,j,k, Obstacle))
	end
	
	#println("\n")

	for i = 2:M

		Al=Apq(i,1,k, Obstacle) #Correspond à la première boucle de for j=1
		for j = 2:M

			Al=hcat(Al, Apq(i,j,k, Obstacle))
			
		end
		#println(Al)
		#println("\n")
		A=vcat(A,Al)
		println("----------\n")
		println(size(A))
		println("----------\n") 
	end
	#println(A,"\n")

	return A
end


function Calcule_C(A,b) #Calcule le vecteur c du système: "Ac=b".  Il s'agit du vecteur des coeff de Fourrier des ondes diffractées
	
	return A\b
end


function Extraire_Cm(C,M,Obstacle) #A partir du vecteur C du système 
	Cm=[[] for i=1:M]
	# Crev = reverse(C,dims=2)
	curseur=0
	for i=1:M
		
		Np=Obstacle[i][4]
		Np = floor(Int, Np)

		push!(Cm[i],C[curseur+1:curseur+2*Np+1])

		curseur=curseur+2*Np+1
		
	end
	# println(Cm)

	return Cm

end


function CalculeUq(Obstacle, x,y, k,Cm,M)

	somme_q = 0

	for q = 1:M
		#somme_n=0

		Nq = Obstacle[q][4]
		Nq = floor(Int, Nq)
		
		x1 = Obstacle[q][1]
		y1 = Obstacle[q][2]
		x2 = x
		y2 = y


		b  = distance(x1,y1,x2,y2)
		angle_2ppq = angle_2p(x1,y1,x2,y2)
		Cq = Cm[q][1]
		#println(Cq,"\n")

		somme_q = somme_q + calculUp(b, angle_2ppq, Cq, Nq)
	end

	return somme_q

end



function Calcule_Utot_MultiDisk(Obstacle,x,y,Cm,Dm,k,M,Beta)

	Uinc  = calculUinc_exact(x,y,Beta,k)

	Udiff = CalculeUq(Obstacle, x, y, k, Cm, M)
	

	Utot = Uinc + Udiff 
	# Utot = real(Uinc)
	# Utot = abs(Udiff)

	return Uinc, Udiff, Utot

end


function Is_inDisk(x,y,Obstacle,M) 
	for i = 1:M

		if distance(x,y,Obstacle[i][1], Obstacle[i][2]) <= Obstacle[i][3]
			return true
		end
	end

	return false
end









#------------------------------------------------------------------------------------------------

#fonctions liées aux test et visualisation et non à la resolution mathematiques du probleme

#------------------------------------------------------------------------------------------------



# --------- genere boule aleatoirement----------
function initBoulesAlea(nb_boules_min,nb_boules_max,xmin,xmax,ymin,ymax,rmin,rmax) # Crée des boulles dont le nombre, la position et le rayon sont aléatoires
	# nb_boules_min = 2
	# nb_boules_max = 20
	nb_boules = rand(nb_boules_min:nb_boules_max)
	boules = [[] for p=1:nb_boules]

	# breaking = 0

	# xmin = 0
	# xmax = 10
	# ymin = 0
	# ymax = 10
	# rmin = 0.5
	# rmax = 3

	i = 0
	iter = 0
	while (i < nb_boules && iter < 10000)
		breaking = 0
		r = rand()*(rmax-rmin)+rmin
		xxmin = xmin + r
		xxmax = xmax - r
		x = rand()*(xxmax-xxmin)+xxmin
		yymin = ymin + r
		yymax = ymax - r
		y = rand()*(yymax-yymin)+yymin
		if (i != 0)
			for k = 1:i
				if (sqrt((boules[k][1] - x)^2 + (boules[k][2] - y)^2) < r + boules[k][3] && breaking == 0)
					breaking = 1
				end
			end
		end
		if (breaking == 0)
			Np = floor(Int64, k*r + cbrt( (1/(2*sqrt(2))) *log(2*sqrt(2)*pi*r*k*e) )^(2) * cbrt(k*r) + 1)
			push!(boules[i+1],x,y,r,Np)

			i = i+1
		end
		iter += 1 
	end

	#println(boules)
	return (boules, nb_boules)
end




#------------ Calcule et sauve image ------------

function Image_Mulit(obstacle,Cm,Dm, M,Beta)

# declaration de la matrice
	Image_inc  = zeros(Complex{Float64}, taille_matrice, taille_matrice)
	Image_diff = zeros(Complex{Float64}, taille_matrice, taille_matrice)
	Image_tot  = zeros(Complex{Float64}, taille_matrice, taille_matrice)
	
	# Parcoure et remplissage de la matrice
	for i = 1:taille_matrice

		for j = 1:taille_matrice

			x,y      = coordonnees(i, j, h, taille_espace)
			r,lambda = conversion_polaire(x, y)

			if !Is_inDisk(x,y,Obstacle, M)

				Image_inc[i,j], Image_diff[i,j], Image_tot[i,j] = Calcule_Utot_MultiDisk(Obstacle, x, y, Cm, Dm, k, M,Beta)
				#println("Image [",i,",",j,"] = ", Image[i,j],"\n")
			end
		end
		println("Images [",i,"]","\n")
	end
	

	return Image_inc, Image_diff, Image_tot
end








function Image_save(obstacle,Cm,Dm, M,Beta)

# declaration de la matrice
	Image_inc2  = zeros(Float64, taille_matrice, taille_matrice)
	Image_diff2 = zeros(Float64, taille_matrice, taille_matrice)
	Image_tot2  = zeros(Float64, taille_matrice, taille_matrice)
	

	Image_inc, Image_diff, Image_tot = Image_Mulit(obstacle,Cm,Dm, M,Beta)
	# Parcoure et remplissage de la matrice
	for i = 1:taille_matrice

		for j = 1:taille_matrice
			Image_inc2[i,j]  = real(Image_inc[i,j])
			Image_diff2[i,j] = abs(Image_diff[i,j])
			Image_tot2[i,j]  = abs(Image_tot[i,j])
		end
	end
	
	# Affichage graphique

	scale_min = - taille_espace/2
	scale_max =   taille_espace/2
	# println(Image_tot2)
	# imshow(transpose(Image), vmin=0.0, vmax=2.0, extent=(scale_min, scale_max, scale_min, scale_max))
	imshow(transpose(Image_inc2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("../results/resMult_inc.svg")
	PyPlot.clf()
	imshow(transpose(Image_diff2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("../results/resMult_diff.svg")
	PyPlot.clf()
	imshow(transpose(Image_tot2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("../results/resMult_tot.svg")

	return Image_inc, Image_diff, Image_tot
end





     

# function Video_Mult(obstacle,Cm,Dm, M,Beta,T) # ne pas utiliser
# 	Sequence = []
# 	println("tarace")
# 	Image_inc, Image_diff, Image_tot = Image_Mulit(obstacle,Cm,Dm, M,Beta)
# 	for t in 1:T
# 		println(". - - -- -- - - - Sequence [",t,"]. -- - - - - - --  -","\n")
		
# 		for i = 1:taille_matrice

# 			for j = 1:taille_matrice
# 				Image_tot[i,j]  = real(exp(im*w*(t-1)*0.1)*Image_tot[i,j])
# 			end
# 		end
# 		push!(Sequence,Image_tot)
# 	end


# 	return Sequence
# end

		





