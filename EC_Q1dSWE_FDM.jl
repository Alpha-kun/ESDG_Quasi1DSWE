using LinearAlgebra
using Plots

#solving SWE on [0,1]

#mesh parameter
L = 1
N = 400
dx = L/N

#physics constants
g = 1.0

#initial height profile
function ψ(x)
    return 0.1.+0.1*exp.(-64*(x.-0.25).*(x.-0.25))
end

#midpoint of intervals
X=[(i-0.5)*dx for i in 1:N]'
#left and right boundary of intervals
X₋=X.-0.5dx;
X₊=X.+0.5dx;
Xₛ=[X₋; X₊]
#height vector
H = ψ(X)
#momentum vector
HU = 0*X

Vlax=[H; HU]
Vec=[H; HU]

function f(H,U,i)
     return [H[i]*U[i]; H[i]*U[i]^2+(g/2)*H[i]^2]
end


function Flax(V,l,r)
    H=V[1,:]
    U=V[2,:]./H
    if l==0
        return [0; f(H,U,r)[2]]
    end
    if r==N+1
        return [0;f(H,U,l)[2]]
    end
    fₗ = f(H,U,l)
    fᵣ = f(H,U,r)
    d = 0.2
    return (fₗ+fᵣ)/2 + 0.5d*(V[:,l]-V[:,r])
end

function Fec(V,l,r)
    H=V[1,:]
    U=V[2,:]./H
    if l==0
        return [0; 0.5g*(H[r]^2)]
    end
    if r==N+1
        return [0; 0.5g*(H[l]^2)]
    end
    f₁ = 0.25*(H[l]+H[r])*(U[l]+U[r])
    f₂ = 0.125*(H[l]+H[r])*(U[l]+U[r])^2 + 0.25g*(H[l]^2+H[r]^2)
    return [f₁;f₂]
end

function dUdt(U, flx)
    U′=0*U;
    for i in 1:N
        if flx==1
            U′[:,i]=(Flax(U, i-1, i) - Flax(U, i, i+1))/dx
        else
            U′[:,i]=(Fec(U, i-1, i) - Fec(U, i, i+1))/dx
        end
    end
    return U′
end

#time parameter
t=0
dt = 0.0025

Elax=[]
Eec=[]
anime = @animate for i in 1:500
    #Lax-Friedrich Flux
    k1=dUdt(Vlax,1)
    k2=dUdt(Vlax+0.5dt*k1,1)
    k3=dUdt(Vlax+0.5dt*k2,1)
    k4=dUdt(Vlax+dt*k3,1)
    Vlax=Vlax+dt*(k1+2k2+2k3+k4)/6
    push!(Elax,0.5*sum((Vlax[2,:].^2)./Vlax[1,:]+g*Vlax[1,:].^2))
    #Energy Conservative Flux
    k1=dUdt(Vec,2)
    k2=dUdt(Vec+0.5dt*k1,2)
    k3=dUdt(Vec+0.5dt*k2,2)
    k4=dUdt(Vec+dt*k3,2)
    Vec=Vec+dt*(k1+2k2+2k3+k4)/6
    push!(Eec,0.5*sum((Vec[2,:].^2)./Vec[1,:]+g*Vec[1,:].^2))
    #plot
    plot(Xₛ,[Vlax[1,:]';Vlax[1,:]'],color="red",ylims=(0,0.25),lab="Lax-Friedrich",primary=onlyone(true))
    p=plot!(Xₛ,[Vec[1,:]';Vec[1,:]'],color="blue",ylims=(0,0.25),lab="Energy Conservative",primary=onlyone(true))
    display(p)
    #sleep(0.1)
    t+=dt
end

plot(1:500, Elax,label="Lax-Friedrich")
plot!(1:500, Eec,label="Entropy Conservative")

gif(anime, "D:\\side_projects\\EC_DGFV\\FV_SWE_Lax_ec_finer.gif", fps = 30)
