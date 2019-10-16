using SpecialFunctions

alpha=pi #Angle de l'onde incidente
alphap=0 #Angle de l'obstacle
bp=0
a=1
rp=1
e=10^(-12)
k=2*pi
Np= floor(Int64,k*a + cbrt(1/(2*sqrt(2))*log(2*sqrt(2)*pi*k*e))^(2) * (k*a)^(1/3) +1)

function calculU(r,teta, Cm)
	N=length(Cm)
	U=0.0
	for m=1:N
		U= U + Cm[m]*besselh(m-N/2,k*r)*exp(im*(m-N/2)*teta)
	end	

	return U
end


Cm=ones(Complex{Float32},2*Np+1,1)
Dm=ones(Complex{Float32},2*Np+1,1)
Hm=ones(Complex{Float32},2*Np+1,1)
temp=exp(im*k*cos(alpha-alphap)*bp)
u=0.0

for m=1:2*Np+1
	Dm[m]=temp*exp(im*(pi*0.5 - alpha)*m)
	Hm[m]=besselh(m-Np,a*k)
	Cm[m]=besselj(m-Np,a*k) * Hm[m]^(-1) * Dm[m]

end
 println(Cm)

 println("\nU=",calculU(40,3*pi/2, Cm))

 