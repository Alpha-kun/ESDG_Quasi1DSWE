using LinearAlgebra
using Plots

#solving Q1D SWE on [0,1]

#mesh parameter
L = 1
N = 100
dx = L/N

#physics constants
g = 1.0

#channel width profile
function a(x)
    return 1 .+ 0.5*sin.(2π*x)
end

#initial height profile
function ψ(x)
    return 0.1 .+ 0.1*exp.(-64*(x.-0.25).*(x.-0.25))
end

function total_entropy(V)
    AH = V[1,:]'
    AHU = V[2,:]'
    H = AH./A
    U = AHU./AH
    η = 0.5A.*(g.*H.^2+H.*U.^2)
    return sum(η)
end

#midpoint of intervals
X=[(i-0.5)*dx for i in 1:N]'

#left and right boundary of intervals
X₋=X.-0.5dx;
X₊=X.+0.5dx;
Xₛ=[X₋; X₊]

#the width vector
A = a(X)
#vector of first conserved quantity (volume)
AH = ψ(X).*A
#vector of second conserved quantity (momentum)
AHU = 0*X

#solution vectors
Vec=[AH; AHU]

Vcf=[AH; AHU]

#compute the central difference derivative of a row vector
function Q(x)
    return (circshift(x,(1,-1))-circshift(x,(1,1))) ./ (2dx)
end

#time derivatives of conserved quantities
function dVdt(V)
    AH = V[1,:]'
    AHU = V[2,:]'
    H = AH./A
    U = AHU./AH
    AH′ = -Q(AHU)
    AHU′ = -0.5*(Q(AHU.*U)+AHU.*Q(U)+U.* Q(AHU)) - g*(AH.*Q(H))
    V′ = [AH′; AHU′];
    return V′
end

function dVcfdt(V)
    AH = V[1,:]'
    AHU = V[2,:]'
    H = AH./A
    U = AHU./AH
    AH′ = -Q(AHU)
    AHU′ = -Q(AHU.*U) -0.5g*Q(AH.*H) + 0.5g*(H.^2 .* Q(A))
    V′ = [AH′; AHU′];
    return V′
end

#time parameter
t=0
dt = 0.0025

ηec=[]
ηcf=[]
anime = @animate for i in 1:800
    #Energy Conservative scheme, rk4
    k1 = dVdt(Vec)
    k2 = dVdt(Vec + 0.5dt*k1)
    k3 = dVdt(Vec + 0.5dt*k2)
    k4 = dVdt(Vec + dt*k3)
    Vec = Vec + dt*(k1+2k2+2k3+k4)/6
    push!(ηec, total_entropy(Vec))

    k1 = dVcfdt(Vcf)
    k2 = dVcfdt(Vcf + 0.5dt*k1)
    k3 = dVcfdt(Vcf + 0.5dt*k2)
    k4 = dVcfdt(Vcf + dt*k3)
    Vcf = Vcf + dt*(k1+2k2+2k3+k4)/6
    push!(ηcf, total_entropy(Vcf))
    #plot
    AH = Vec[1,:]
    H = AH./(A')
    plot(X', H, color="blue",ylims=(0,0.25),lab="Energy Conservative")

    AH = Vcf[1,:]
    H = AH./(A')
    p=plot!(X', H, color="red",ylims=(0,0.25),lab="Central Flux")
    display(p)
    t+=dt
end



plot(1:800, ηec,label="EC")
plot!(1:800, ηcf,label="CF")

gif(anime, "D:\\Rice\\spring2022\\Q1DSWE_project\\EC_Q1DSWE_FDM_1psin.gif", fps = 30)
