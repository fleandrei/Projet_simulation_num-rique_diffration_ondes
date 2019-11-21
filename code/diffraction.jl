include("./polar.jl")
using LinearAlgebra

# Calcule l'onde diffractée complexe U
function calculUp(r,teta, Cm, Np)
	#N=length(Cm)
	U = 0.0
	for m = -Np:Np
		idx = m + Np + 1
		#println("m = ",m)
		#println("\nk=",k)
		#println("\nr=",r)
		#println("\nk*r",k*r)
		#println(besselh(m, k*r))
		U   = U + Cm[idx] * besselh(m, k*r) * exp(im * m * teta)
		#println("U vaut ", U, "  Cm[i] = ",Cm[idx],"\n")
	end	

	return U
end


# Calcule l'onde incidente complexe U
function calculUinc(r,teta, Dm, Np)
	#N=length(Cm)
	U = 0.0
	for m = -Np:Np
		idx = m + Np + 1
		
		U   = U + Dm[idx] * besselj(m, k*r) * exp(im * m * teta)
	end	

	return U
end


function calculRCS(U)
	return 10 * log10(2*pi * abs(U)^2)
end


function CDH(bp,alphap,alpha,Np, k) #Calcule Cm: coef Fourrier onde réfléchie, Dm: Coef Fourrier onde incidente, Hm: vecteur Besselh
	Cm   = zeros(Complex{Float32}, 2*Np + 1, 1)
	Dm   = zeros(Complex{Float32}, 2*Np + 1, 1)
	Hm   = zeros(Complex{Float32}, 2*Np + 1, 1)
	temp = exp(im*k * cos(alpha-alphap) * bp)

	for m = -Np:Np
		idx     = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - alpha) * m)
		Hm[idx] = besselh( m, a*k )
		Cm[idx] = besselj( m, a*k ) * Hm[idx]^(-1) * Dm[idx]
	end

	return (Cm,Dm,Hm)
end


function CoeFourrier_OndeInc(bp,alphap,alpha,Np, k)
	Dm   = zeros(Complex{Float32}, 2*Np +1, 1)
	temp = exp(im*k * cos(alpha-alphap) * bp)

	for m = -Np:Np
		idx   = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - alpha)* m)
	end

	return Dm
end


function BesselHVector(Np,a,k)

	Hm   = zeros(Complex{Float32},2*Np+1,1)
	temp = a * k

	for m = -Np:Np
		idx     = m + Np + 1
		Hm[idx] = besselh(m, temp)
	end

	return Hm
end


function Calcule_b(bp,alphap,alpha,Np, k, a) 

	Dm   = CoeFourrier_OndeInc(bp, alphap, alpha, Np, k)
	b    = zeros(Complex{Float32}, 2*Np + 1, 1)
	temp = a * k

	for m = -Np:Np
		idx    = m + Np + 1
		b[idx] = besselj(m, temp) * Dm[idx]
	end

	return b
end



function CoeFourrier_OndeDiffracte(bp,alphap,alpha,Np, k, a)

	Cm = zeros(Complex{Float32},2*Np+1,1)
	b  = Calcule_b(bp, alphap, alpha, Np, k, a)
	Hm = BesselHVector(Np, a, k)
	#A=zeros(Complex{Float32},2*Np+1, 2*Np+1)
	A  = Matrix(I, 2*Np+1, 2*Np+1)
	A  = A.*Hm
	Cm = A\b

	return Cm
end

function Phi(m,rp, teta,k)
	return besselh(m, k*rp)*exp(im*m*teta)
end

function Phi_Chap(m,rp,teta,k)
	return besselj(m, rp*k)*exp(im*m*teta)
end

function Smn(m,n,b,teta, k)
	return Phi(m-n,b, teta,k)
end



function dmp(p, Obstacle, Beta,k) # Calcule les coeff Fourrier de l'onde incidente pour la boule p
	Np=Obstacle[p][4]
	Dm = ones(Complex{Float32}, 2*Np +1, 1)
	println(p)
	bp,teta=conversion_polaire(Obstacle[p][1], Obstacle[p][2])
	temp = exp(im*k * cos(Beta-teta) * bp)

	for m = -Np:Np
		idx   = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - Beta)* m)
	end

	return Dm
end


function Extraire_Dm(M,Obstacle,Beta,k) # Renvoie Tableau de tableau Dm tel que Dm[P]="coeff Fourrier onde incidente pour la boule P"
	Dm=[[] for i=1:M]

	for i=1:M
		push!(Dm[i], dmp(i,Obstacle,Beta,k))
	end

	return Dm
end


function Calcule_B(M, Obstacle, Beta,k,Dm) #Calcule le vecteur b du système "Ac=b"
	B=zeros(Complex{Float32}, (2*Np +1)*M, 1)
	for i=1:M
		Np=Obstacle[i][4]
		dm=Dm[i][1]
		ap=Obstacle[i][3]
		
		for m=-Np:Np
			temp=m + Np + 1
			idx=  temp+ (i-1)*(2*Np +1)
			B[idx]=-besselj(m,ap*k)/besselh(m,k*ap)* dm[temp]
		end
	end
	return B
end

function Apq(p,q,k, Obstacle) #Calcule la sous-matrice d'indices p,q de la matrice A 
	Np=Obstacle[p][4]
	Nq=Obstacle[q][4]
	if p==q
		return 1*Matrix(I,2*Np+1,2*Np+1)
	else
		ap=Obstacle[p][3]
		b=distance(Obstacle[p][1], Obstacle[p][2], Obstacle[q][1], Obstacle[q][2])
		teta=angle(Obstacle[p][1], Obstacle[p][2], Obstacle[q][1], Obstacle[q][2])
		A=zeros(Complex{Float32}, 2*Np+1, 2*Nq+1)

		for m =1:2*Np+1
			for n= 1:2*Nq+1
				A[m,n]=besselj(m,ap*k)/besselh(m, ap*k) * Smn(m,n,b,teta,k)
			end
		end
		println(A,"\n b = ",b," teta = ",teta,"\n")

		return A
	end	
end


function Calcule_A(M, Obstacle, k) #Calcule la matrice A du système "Ac=b"
	
	
	A=Apq(1,1,k, Obstacle) #Correspond à la première boucle de for j=1
	for j = 2:M

		A=vcat(A, Apq(1,j,k, Obstacle))
	end
	#println(A)
	#println("\n")

	for i = 2:M

		Al=Apq(i,1,k, Obstacle) #Correspond à la première boucle de for j=1
		for j = 2:M

			Al=vcat(Al, Apq(i,j,k, Obstacle))
			
		end
		#println(Al)
		#println("\n")
		A=hcat(A,Al)
	
	end
	return A
end


function Calcule_C(A,b) #Calcule le vecteur c du système: "Ac=b".  Il s'agit du vecteur des coeff de Fourrier des ondes diffractées
	
	return A\b
end


function Extraire_Cm(C,M,Obstacle) #A partir du vecteur C du système 
	Cm=[[] for i=1:M]

	curseur=0
	for i=1:M
		
		Np=Obstacle[i][4]
		push!(Cm[i],C[curseur+1:curseur+2*Np+1])
		curseur=curseur+2*Np+1
		
	end

	return Cm

end

function Boule_Proche(Obstacle,x, y,M)
	
	MinDist=1000000
	p=0
	for i=1:M
		Dist=distance(x,y,Obstacle[i][1], Obstacle[i][2])
		if Dist<MinDist
			MinDist=Dist
			p=i
		end
	end
	return p
end


function CalculeUq(Obstacle, r, teta, k, p,Cm,M)
	
	Np=Obstacle[p][4]

	somme_m=0

	for m=-Np:Np
		somme_q=0
		for q=1:M
			if q != p
				Nq=Obstacle[q][4]
				somme_n=0
				x1=Obstacle[p][1]
				y1=Obstacle[p][2]
				x2=Obstacle[q][1]
				y2=Obstacle[q][2]
				bpq=distance(x1,y1,x2,y2)
				anglepq=angle(x1,y1,x2,y2)
				Cq=Cm[q][1]
				for n=-Nq:Nq
					somme_n=somme_n+ Smn(n,m,bpq,anglepq, k)*Cq[n+Nq+1]
				end
				somme_q=somme_q + somme_n
			end
		end
		somme_m= somme_m +	somme_q * Phi_Chap(m,r,teta,k)

	end

	return somme_m

end



function Calcule_UDiff_MultiDisk(Obstacle,r,teta,Cm,k,M,p)

	Up = calculUp(r,teta, Cm[p][1], Np)
	Uq = CalculeUq(Obstacle, r, teta, k, p,Cm,M)
	U= Up + Uq 
	# println("Up = ",Up, " , Uq = ",Uq, ", U = ",U,"\n")
	return U
end

function Calcule_Utot_MultiDisk(Obstacle,x,y,Cm,Dm,k,M)
	p=Boule_Proche(Obstacle,x,y,M)
	Np=Obstacle[p][4]
	r=distance(Obstacle[p][1], Obstacle[p][2], x,y)
	teta=angle(Obstacle[p][1], Obstacle[p][2], x,y)
	Udiff = Calcule_UDiff_MultiDisk(Obstacle,r,teta,Cm,k,M,p)
	Uint  = calculUinc(r,teta, Dm[p][1], Np)
	Utot= abs(Udiff + Uint)
	#println("Udiff = ",Udiff, " , Uint = ",Uint, ", Utot = ",Utot,"\n")

	return Utot

end


function Is_inDisk(x,y,Obstacle,M)
	for i=1:M
		if distance(x,y,Obstacle[i][1], Obstacle[i][2]) <= Obstacle[i][3]
			return true
		end
	end

	return false
end