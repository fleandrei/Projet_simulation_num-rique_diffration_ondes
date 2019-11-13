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


function CoeFourrier_OndeRefracte(bp,alphap,alpha,Np, k, a)

	Cm = zeros(Complex{Float32},2*Np+1,1)
	b  = Calcule_b(bp, alphap, alpha, Np, k, a)
	Hm = BesselHVector(Np, a, k)
	#A=zeros(Complex{Float32},2*Np+1, 2*Np+1)
	A  = Matrix(I, 2*Np+1, 2*Np+1)
	A  = A.*Hm
	Cm = A\b

	return Cm
end