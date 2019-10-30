# Calcule l'onde diffractée complexe U
function calculUp(r,teta, Cm, Np)
	#N=length(Cm)
	U=0.0
	for m=1:2*Np + 1
		U= U + Cm[m]*besselh(m-Np,k*r)*exp(im*(m-Np)*teta)
	end	

	return U
end


# Calcule l'onde incidente complexe U
function calculUinc(r,teta, Dm, Np)
	#N=length(Cm)
	U=0.0
	for m=1:2*Np + 1
		U= U + Dm[m]*besselj(m-Np,k*r)*exp(im*(m-Np)*teta)
	end	

	return U
end


function calculRCS(U)
	return 10*log10(2*pi*abs(U)^2)
end