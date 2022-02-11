using FastGaussQuadrature
using LinearAlgebra
using Plots

#Domain: [-1,1]

order=4
N=16 #number of element
h=1/N #half of element size

#create quadrature points
nodes, weights = gausslobatto(order)

#compute interpolation vector
Lₗ=zeros(order)
Lᵣ=zeros(order)
for i in 1:order
    denom = prod(nodes[i] .- nodes[filter(x -> x != i, 1:order)])
    Lₗ[i]=prod(-1 .- nodes[filter(x -> x != i, 1:order)]) / denom
    Lᵣ[i]=prod(1 .- nodes[filter(x -> x != i, 1:order)]) / denom
end

#compute mass matrix
M=h*Diagonal(weights)

#compute stiffness matrix
K=zeros(order,order)
for i in 1:order
    denom = prod(nodes[i] .- nodes[filter(x -> x != i, 1:order)])
    for j in 1:order
        if i==j
            L̇ᵢ=sum(1 ./ (nodes[i] .- nodes[filter(x -> x != i, 1:order)]))
            K[i,j]=weights[j]*L̇ᵢ
        else
            L̇ᵢ=prod(nodes[j] .- nodes[filter(x -> (x != i) & (x != j), 1:order)]) / denom
            K[i,j]=weights[j]*L̇ᵢ
        end
    end
end

#global mass matrix
Mglb = Diagonal(reduce(vcat, [weights for i in 1:N]))

#global stiffness matrix
Qglb = zeros(order*N, order*N)
for i in 1:N
    Qglb[(order*(i-1)+1):(order*i),(order*(i-1)+1):(order*i)]+=K
    if i < N
        Qglb[(order*i):(order*i+1),(order*i):(order*i+1)]+=[-1/2 -1/2; 1/2 1/2]
    end
end
Qglb[1,1]+=1/2
Qglb[1,order*N]+=1/2
Qglb[order*N,1]+=-1/2
Qglb[order*N,order*N]+=-1/2

#########################################################
################ problem specific set up ################
#########################################################
#physics constants
g = 1.0

#channel width profile
function a(x)
    return 1 #.+ 0.5*sin.(2π*x)
end

#initial height profile
function ψ(x)
    return 0.1 .+ 0.1*exp.(-64*(x.-0.25).*(x.-0.25))
end


function total_entropy(V)
    AH = V[:,1]
    AHU = V[:,2]
    H = AH./A
    U = AHU./AH
    η = 0.5A.*(g.*H.^2+H.*U.^2)
    return sum(Mglb*η)
end

#midpoint of intervals
X=reduce(vcat, [(-1+(2*i-1)*h) .+ h*nodes for i in 1:N])
#the width vector
A = a(X)
#vector of first conserved quantity (volume)
AH = ψ(X).*A
#vector of second conserved quantity (momentum)
AHU = 0*X

#solution vectors
Vec=[AH AHU]

#compute the central difference derivative of a row vector
function Q(x)
    return -Mglb\(Qglb*x)
end

#time derivatives of conserved quantities
function dVdt(V)
    AH = V[:,1]
    AHU = V[:,2]
    H = AH./A
    U = AHU./AH
    AH′ = -Q(AHU)
    AHU′ = -0.5*(Q(AHU.*U)+AHU.*Q(U)+U.* Q(AHU)) - g*(AH.*Q(H))
    V′ = [AH′ AHU′];
    return V′
end

#time parameter
t=0
dt = 0.002

ηec=[]
anime = @animate for i in 1:800
    #Energy Conservative scheme, rk4
    global Vec, t#WTF is going on (didn't need this line previously)
    k1 = dVdt(Vec)
    k2 = dVdt(Vec + 0.5dt*k1)
    k3 = dVdt(Vec + 0.5dt*k2)
    k4 = dVdt(Vec + dt*k3)
    Vec = Vec + dt*(k1+2k2+2k3+k4)/6
    push!(ηec, total_entropy(Vec))
    #plot
    AH = Vec[:,1]
    H = AH./(A')
    p=plot(X, H, color="blue",ylims=(0,0.25),lab="Energy Conservative")
    display(p)
    t+=dt
end

total_entropy(Vec)


plot(1:800, ηec,label="Entropy")

gif(anime, "D:\\Rice\\spring2022\\Q1DSWE_project\\EC_Q1DSWE_DG_1.gif", fps = 30)