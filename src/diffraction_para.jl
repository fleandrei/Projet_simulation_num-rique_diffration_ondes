include("./polar.jl")
using LinearAlgebra

# Calcule l'onde diffractée complexe U
function calculUp(r,teta, Cm, Np, k)

	U = 0.0

	for m = -Np:Np

		idx = m + Np + 1
		U   = U + (Cm[idx] * Phi(m,r, teta,k))
	end	

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
	A  = Matrix(I, 2*Np+1, 2*Np+1)

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

		D = Matrix_D(Np,k,ap)
		S = Matrix_S(Np,Nq,b,teta,k)

		A = D * transpose(S)

		return A
	end	
end



function Calcule_A(M, Obstacle, k) #Calcule la matrice A du système "Ac=b"
	
	
	A=Apq(1,1,k, Obstacle) #Correspond à la première boucle de for j=1

	for j = 2:M

		A=hcat(A, Apq(1,j,k, Obstacle))
	end
	

	for i = 2:M

		Al=Apq(i,1,k, Obstacle) #Correspond à la première boucle de for j=1

		for j = 2:M

			Al=hcat(Al, Apq(i,j,k, Obstacle))
			
		end
		A=vcat(A,Al)
	end

	return A
end




function Calcule_PortionLigne_SousMatrice(LigneDebut,NbrLigne,k, Obstacle, M)
	
	A=Apq(LigneDebut,1,k, Obstacle)

	for j = 2:M

		A=hcat(A, Apq(LigneDebut,j,k, Obstacle))
	end

	if NbrLigne>1
		for i = LigneDebut+1:LigneDebut+NbrLigne - 1

			Al=Apq(i,1,k, Obstacle) #Correspond à la première boucle de for j=1
			for j = 2:M

				Al=hcat(Al, Apq(i,j,k, Obstacle))
			
			end
			A=vcat(A,Al)
		end
	end

	return A
end



function Calcule_Parallel_A(M, Obstacle, k, Granular::Int=0) #Calcule la matrice A du système "Ac=b" en multi process
	np = nprocs()
	i=1
	nextidx() = (idx = i; i+=Granular; idx)
	Dico_Mat=Dict()

	if Granular==0
		coef=0.3
		Granular=floor(Int, coef*M/(np-1))
		if Granular==0
			Granular=1
		end
	end

	@sync begin # On entre dans une Section synchrone: Permet de s'assurer que toutes les tasks sont terminées avant de sortir de cette zone

		for p = 1:np # On parcourt les processus 

        	if p != myid() || np == 1 #Seul les 3 processus ouvrier vont vraiment travailler (sauf si il n'y a qu'un seul processus): Modèle Patron/Ouvrier
				
				@async begin # On rentre dans une Section Asynchrone. Le code contenu dans cette zone correspond à une Task. On lance ainsi 'np' task différentes de manière asynchrones MAIS pas dans des 'np' processus différent. Toutes les task s'exécutent sur le processus principal i.e. le processus p=1 

                    while true
                        idx = nextidx() # Indice de la prochaine ligne de sous-matrices à calculer

                        if idx > M # Si il n'y a plus rien à calculer, on sort de la boucle (et la task est terminé)
                            break
                        end

                        longueur=Granular # Granularité: nombre de lignes de sous-matrices que l'on va confier au processus lorsqu'il viendra demander du travail

                        if M - idx < Granular  # Si il reste moins de 'Granular' lignes à calculer
                        	longueur = M -idx +1
                        end

                        temp = remotecall(Calcule_PortionLigne_SousMatrice, p, idx, longueur, k, Obstacle, M)  # On appelle le processus 'p' et on lui demande d'exécuter la fonction Calcule_PortionLigne_SousMatrice avec les arguments (idx, longueur, k, Obstacle, M)

                        temp2=fetch(temp)
                        Dico_Mat[idx] = temp2# le résultat renvoyé par remotecall n'est pas directement exploitable c'est pourquoi on utilise fetch()
                                                
                    end
                end # Sortie de @async
            end
        end
    end # Sortie de @sync: Attend que toutes les taches soient finit avant de sortir

    Cle=sort(collect(keys(Dico_Mat)))


    A=Dico_Mat[Cle[1]]
    for k = 2:length(Cle)
    	A=vcat(A,Dico_Mat[Cle[k]])

    end

    return A
end


function VCAT(Dico_Mat, Cle, Indmin, IndMax)
	Res=Dico_Mat[Cle[Indmin]]
	for i=(Indmin+1):IndMax
		Res=vcat(Res, Dico_Mat[Cle[i]])
	end
	return Res

end


function Calcule_C(A,b) #Calcule le vecteur c du système: "Ac=b".  Il s'agit du vecteur des coeff de Fourrier des ondes diffractées
	
	return A\b
end


function Extraire_Cm(C,M,Obstacle) #A partir du vecteur C du système 

	Cm=[[] for i=1:M]

	curseur=0
	for i=1:M
		
		Np=Obstacle[i][4]
		Np = floor(Int, Np)

		push!(Cm[i],C[curseur+1:curseur+2*Np+1])

		curseur=curseur+2*Np+1
		
	end

	return Cm

end


function CalculeUq(Obstacle, x,y, k,Cm,M)

	somme_q = 0

	for q = 1:M

		Nq = Obstacle[q][4]
		Nq = floor(Int, Nq)
		
		x1 = Obstacle[q][1]
		y1 = Obstacle[q][2]
		x2 = x
		y2 = y

		b  = distance(x1,y1,x2,y2)
		angle_2ppq = angle_2p(x1,y1,x2,y2)
		Cq = Cm[q][1]

		somme_q = somme_q + calculUp(b, angle_2ppq, Cq, Nq, k)
	end

	return somme_q

end



function Calcule_Utot_MultiDisk(Obstacle,x,y,Cm,Dm,k,M,Beta)

	Uinc  = calculUinc_exact(x,y,Beta,k)

	Udiff = CalculeUq(Obstacle, x, y, k, Cm, M)
	

	Utot = Uinc + Udiff 

	return Uinc, Udiff, Utot

end


function Calcule_Utot_MultiDisk_para(Obstacle,x,y,Cm,Dm,k,M,Beta)

	Uinc  = calculUinc_exact(x,y,Beta,k)

	Udiff = CalculeUq(Obstacle, x, y, k, Cm, M)
	

	Utot = Uinc + Udiff 

	return abs(Utot)

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

	nb_boules = rand(nb_boules_min:nb_boules_max)
	boules = [[] for p=1:nb_boules]

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
			
			end
		end
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

	imshow(transpose(Image_inc2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("results/resMult_inc.png")
	println("----- incindent saved in 'results/resMult_inc.png' -----:\n")
	PyPlot.clf()

	imshow(transpose(Image_diff2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("results/resMult_diff.png")
	println("----- diffraction saved in 'results/resMult_diff.png' -----:\n")
	PyPlot.clf()

	imshow(transpose(Image_tot2), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("results/resMult_tot.png")
	println("----- total saved in 'results/resMult_tot.png' -----:\n")

	return Image_inc, Image_diff, Image_tot
end


#------------------------------------------------------------------------------------------------

#fonctions liées aux test et visualisation et non à la resolution mathematiques du probleme
#									PARALELLE

#------------------------------------------------------------------------------------------------

# Calcule la ligne num 'Idx'_Ligne de l'image.
@everywhere function Image_Calcul_Ligne(Idx_Ligne, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)

	res= zeros(Float64, taille_matrice)

	for j=1:taille_matrice
		x,y = coordonnees(Idx_Ligne, j, h, taille_espace)
		r, lambda = conversion_polaire(x, y)

		if !Is_inDisk(x,y,Obstacle, M)
			res[j] = Calcule_Utot_MultiDisk_para(Obstacle, x, y, Cm, Dm, k, M, Beta)
		end

	end

	return res

end


# Calcule une portion de l'image: On calcule 'nbr_ligne' lignes à partir de la ligne d'indice 'Idx'
@everywhere function Image_Calcul_Portion(Idx, nbr_ligne, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)

	res= zeros(Float64, nbr_ligne, taille_matrice)
	for i=Idx:(Idx+nbr_ligne-1)
		for j=1:taille_matrice
			x,y = coordonnees(i, j, h, taille_espace)
			r, lambda = conversion_polaire(x, y)

			if !Is_inDisk(x,y,Obstacle, M)
				res[i-Idx+1,j] = Calcule_Utot_MultiDisk_para(Obstacle, x, y, Cm, Dm, k, M, Beta)
			end

		end
	end

	return res

end


function Image_Multi_Proc(obstacle, taille_matrice, taille_espace,Cm,Dm, k, M,Beta, h, Granular)

	np = nprocs() #Nombre de processus

	# declaration de la matrice
	Image = zeros(Float64, taille_matrice, taille_matrice)

	i = 1
	
	nextidx() = (idx = i; i+=Granular; idx) #Fonction commune à toutes les tasks lancées dans Image_Multi_Proc. Renvoie l'indice de la prochaine ligne à calculer. 

	@sync begin # On entre dans une Section synchrone: Permet de s'assurer que toutes les tasks sont terminées avant de sortir de cette zone

        for p = 1:np # On parcourt les processus 

        	if p != myid() || np == 1 #Seul les 3 processus ouvrier vont vraiment travailler (sauf si il n'y a qu'un seul processus): Modèle Patron/Ouvrier

        		@async begin # On rentre dans une Section Asynchrone. Le code contenu dans cette zone correspond à une Task. On lance ainsi 'np' task différentes de manière asynchrones MAIS pas dans des 'np' processus différent. Toutes les task s'exécutent sur le processus principal i.e. le processus p=1 

        			#println(p)
                    while true
                        idx = nextidx() # Indice de la prochaine ligne à calculer

                        if idx > taille_matrice # Si il n'y a plus rien à calculer, on sort de la boucle (et la task est terminé)
                            break
                        end

                        longueur=Granular # Granularité: nombre de lignes que l'on va confier au processus lorsqu'il viendra demander du travail

                        if taille_matrice - idx < Granular  # Si il reste moins de 'Granular' lignes à calculer
                        	longueur=taille_matrice -idx +1
                        end

                        temp = remotecall(Image_Calcul_Portion, p, idx, longueur, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)  # On appelle le processus 'p' et on lui demande d'executer la fonction Image_Calcul_Portion avec les arguments (idx, longueur, Obstacle, taille_matrice, taille_espace, Cm, Dm, k, M, Beta, h)

                        Image[idx:(idx+longueur-1),:]=fetch(temp) # le résultat renvoyé par remotecall n'est pas directement exploitable c'est pourquoi on utilise fetch()
                       
                    end
                end # Sortie de @async

        	end
        end
    end # Sortie de @sync: Attend que toutes les taches soient finit avant de sortir

	return Image
end


function save_image_para(Image)
	# Affichage graphique
	scale_min = 0
	scale_max = taille_espace
	
	imshow(transpose(Image), extent=(scale_min, scale_max, scale_min, scale_max))
	colorbar()
	savefig("results/para_resMult_tot.png")
	println("\n----- total saved in 'results/para_resMult_tot.png' -----:\n")
	PyPlot.clf()

end


function champ_lointain(theta, Cm, Obstacle,M)


	sum_disk = 0
	for p in 1:M

		Np = Obstacle[p][4]
		x  = Obstacle[p][1]
		y  = Obstacle[p][2]

		Cp = Cm[p][1]

		b, alpha = conversion_polaire(x, y)


		sum_Cn = 0
		for n in -Np:Np
			index = n + Np+1

			sum_Cn += exp( im*n*(theta - pi/2) ) * Cp[index]
		end

		sum_disk += exp(-im*b*k* cos(theta - alpha)) * sum_Cn
	end


	res = exp(-im * pi/4) * sqrt(2/(k*pi)) * sum_disk

	return res
end


function calculRCS(theta, Cm, Obstacle,M)
	F = champ_lointain(theta, Cm, Obstacle,M)
	return 10 * log10(2*pi * abs(F)^2)
end


function plot_RCS(Cm, Obstacle,M)
	Deg=zeros(361)
	RCS=zeros(361)

	for i=0:360
		theta = Deg2Rad(i)

		Deg[i+1] = i
		RCS[i+1] = calculRCS(theta, Cm, Obstacle,M)
	end

	plot(Deg,RCS)
	savefig("results/RCS.png")
	println("\n----- total saved in 'results/RCS.png' -----:\n")
	PyPlot.clf()
end


