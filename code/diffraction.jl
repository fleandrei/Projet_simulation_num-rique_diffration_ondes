# Calcule l'onde diffractée complexe U
function calculUp(r,teta, Cm, Np)
	#N=length(Cm)
	U = 0.0
	for m = -Np:Np
		idx = m + Np + 1
		U   = U + Cm[idx] * besselh(m, k*r) * exp(im * m * teta)
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
	Cm   = ones(Complex{Float32}, 2*Np + 1, 1)
	Dm   = ones(Complex{Float32}, 2*Np + 1, 1)
	Hm   = ones(Complex{Float32}, 2*Np + 1, 1)
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
	Dm   = ones(Complex{Float32}, 2*Np +1, 1)
	temp = exp(im*k * cos(alpha-alphap) * bp)

	for m = -Np:Np
		idx   = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - alpha)* m)
	end

	return Dm
end


function BesselHVector(Np,a,k)

	Hm   = ones(Complex{Float32},2*Np+1,1)
	temp = a * k

	for m = -Np:Np
		idx     = m + Np + 1
		Hm[idx] = besselh(m, temp)
	end

	return Hm
end


function Calcule_b(bp,alphap,alpha,Np, k, a) 

	Dm   = CoeFourrier_OndeInc(bp, alphap, alpha, Np, k)
	b    = ones(Complex{Float32}, 2*Np + 1, 1)
	temp = a * k

	for m = -Np:Np
		idx    = m + Np + 1
		b[idx] = besselj(m, temp) * Dm[idx]
	end

	return b
end



function CoeFourrier_OndeDiffracte(bp,alphap,alpha,Np, k, a)

	Cm = ones(Complex{Float32},2*Np+1,1)
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



function dmp(p, Obstacle, Beta, Np,k)
	Dm = ones(Complex{Float32}, 2*Np +1, 1)
	bp=Obstacle[p,1]
	teta=Obstacle[p,2]
	temp = exp(im*k * cos(Beta-teta) * bp)

	for m = -Np:Np
		idx   = m + Np + 1
		Dm[idx] = temp * exp(im * (pi*0.5 - Beta)* m)
	end

	return Dm
end



function Calcule_B(M,Np Obstacle, Beta,k)
	B=zeros(Complex{Float32}, (2*Np +1)*M, 1)
	for i=1:M
		dm=dmp(i,Obstacle,Beta,k)
		ap=Obstacle[p,3]
		
		for m=-Np:Np
			temp=m + Np + 1
			idx=  temp+ (i-1)*(2*Np +1)
			B[idx]=-besselj(ap*k)/besselh(k*ap)* dm[temp]
		end
	end
	return B
end

function Apq(p,q,Np,Nq,k, Obstacle)
	if p==q
		return eye(2*Np+1)
	else
		ap=Obstacle[p,3]
		b=distance(Obstacle[p,1], Obstacle[p,2], Obstacle[q,1], Obstacle[q,2])
		teta=angle(Obstacle[p,1], Obstacle[p,2], Obstacle[q,1], Obstacle[q,2])
		A=zeros(2*Np+1, 2*Nq+1)

		for m =1:2*Np+1
			for n= 1:2*Nq+1
				A[m,n]=besselj(ap*k)/besselh(ap*k) * Smn(m,n,b,teta,k)
			end
		end

		return A
	end	
end


function Calcule_A(M,Np, Obstacle, k)
	
	for i = 1:M

		for j = 1:M
			if j==1
				Al=Apq(i,j,Obstacle[i,4],Obstacle[j,4],k, Obstacle)
			else

				Al=vcat(Al, Apq(i,j,Obstacle[i,4],Obstacle[j,4],k, Obstacle))
			end
		end

		if i==1
			A=Al
		else 
			A=hcat(A,Al)
		end
	end
	return A
end