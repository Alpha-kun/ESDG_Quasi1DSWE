using LinearAlgebra
using Plots
using LaTeXStrings
#solving Q1D SWE on [0,1]

#mesh parameter
L = 1
N = 100
dx = L/N

#physics constants
g = 1.0

#channel width profile
function a(x)
    return 0.25 .+ 0.125*sin.(2π*x)
end

#bottom topography
function b(x)
    return 0.05*exp.(-64*(x.-0.5).*(x.-0.5))
end

#initial height profile
function ψ(x)
    return 0.1 .+ 0.1*exp.(-64*(x.-0.25).*(x.-0.25))
end


#midpoint of intervals
X=[(i-0.5)*dx for i in 1:N]

#left and right boundary of intervals
X₋=X.-0.5dx;
X₊=X.+0.5dx;
Xₛ=[X₋ X₊]

#the width vector
A = a(X)
#the Bottom vector
B = b(X)
#vector of first conserved quantity (volume)
AH = (ψ(X)-B).*A
#vector of second conserved quantity (momentum)
AHU = 0*X


#solution vectors
Vec=[AH AHU]

Vcf=[AH AHU]

function total_entropy(V)
    AH = V[:,1]
    AHU = V[:,2]
    H = AH./A
    U = AHU./AH
    η = A.*(0.5g.*H.^2 + 0.5H.*U.^2 + g*H.*B)
    return sum(η)
end

#compute the central difference derivative of a row vector
function Q(x)
    return (circshift(x,-1)-circshift(x,1)) ./ (2dx)
end

#time derivatives of conserved quantities

#EC scheme
function dVdt(V)
    AH = V[:,1]
    AHU = V[:,2]
    H = AH./A
    U = AHU./AH
    AH′ = -Q(AHU)
    AHU′ = -0.5*(Q(AHU.*U)+AHU.*Q(U)+U.* Q(AHU)) - g*(AH.*(Q(H)+Q(B)))
    V′ = [AH′ AHU′]
    return V′
end

#central difference scheme
function dVcfdt(V)
    AH = V[:,1]
    AHU = V[:,2]
    H = AH./A
    U = AHU./AH
    AH′ = -Q(AHU)
    AHU′ = -Q(AHU.*U) -0.5g*Q(AH.*H) + 0.5g*(H.^2 .* Q(A)) - g*AH .* Q(B)
    V′ = [AH′ AHU′]
    return V′
end

#time parameter
t=0
dt = 0.001

ηec=[]
ηcf=[]

total_steps = 1000
anime = @animate for i in 1:total_steps
    global Vec, Vcf, t
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

    #plot bottom
    plot(X, B, color="black",ylims=(0,0.25),lab="bottom")

    #plot EC height
    AH = Vec[:,1]
    H = AH./(A)
    p = plot!(X, B+H, color="blue",ylims=(0,0.25),lab="EC")

    #plot central flux height
    AH = Vcf[:,1]
    H = AH./(A)
    p=plot!(X, B+H, color="red",ylims=(0,0.25),lab="Central Flux")

    display(p)
    t+=dt
end


plot(1:total_steps, ηec,label="ours")
plot!(1:total_steps, ηcf,label="central flux",xlabel="time step", ylabel=L"\eta",legend=:bottomleft)

#savefig("D:\\Rice\\spring2022\\RURS\\eta_comparison.png")
#gif(anime, "D:\\Rice\\spring2022\\Q1DSWE_project\\EC_Q1DSWE_FDM_1psin.gif", fps = 30)
