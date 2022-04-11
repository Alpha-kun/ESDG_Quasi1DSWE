using PlotlyJS

function a(x)
    return 1 .+ 0.5*sin.(π*x)
end

function b(x)
    return 0.5*exp.(-36*(x.+0.25).*(x.+0.25))
end

function h(x)
    return 0.6 +0.2*exp.(-36*(x.-0.25).*(x.-0.25))
end

X(u,v) = u
Y(u,v) = 0.5a(u)
Z(u,v) = 1 .- v .* (1 .- b(u))

U = X
V = range( 0, 1, length=51)

xs = [X(u, v) for u in U, v in V]
ys = [Y(u, v) for u in U, v in V]
zs = [Z(u, v) for u in U, v in V]

trace1=surface(x=xs, y=ys, z=zs,
    lighting_specular=0.15,
    opacity = 0.7,
    lightposition=attr(x=20, y=10, z=100),
    colorscale=[[0, "rgb(100,100,100)"], [1, "rgb(100,100,100)"]],
    showscale=false)

xs = [X(u, v) for u in U, v in V]
ys = [-Y(u, v) for u in U, v in V]
zs = [Z(u, v) for u in U, v in V]

trace2=surface(x=xs, y=ys, z=zs,
    lighting_specular=0.15,
    opacity = 0.7,
    lightposition=attr(x=20, y=10, z=100),
    colorscale=[[0, "rgb(100,100,100)"], [1, "rgb(100,100,100)"]],
    showscale=false)

U = range(-1, 1,   length=51)
V = range(-1, 1, length=51)

xs = [u for u in U, v in V]
ys = [0.5v*a(u) for u in U, v in V]
zs = [b(u) for u in U, v in V]

trace3 = surface(x=xs, y=ys, z=zs,
    lighting_specular=0.15,
    opacity = 0.7,
    lightposition=attr(x=20, y=10, z=100),
    colorscale=[[0, "rgb(100,100,100)"], [1, "rgb(100,100,100)"]],
    showscale=false)

U = range(-1, 1,   length=51)
V = range(-1, 1, length=51)

xs = [u for u in U, v in V]
ys = [0.5v*a(u) for u in U, v in V]
zs = [h(u) for u in U, v in V]

trace4 = surface(x=xs, y=ys, z=zs,
    lighting_specular=0.05,
    lightposition=attr(x=20, y=10, z=1e3),
    colorscale=[[0, "rgb(10,70,150)"], [1, "rgb(210,230,245)"]])

plot([trace1, trace2, trace3, trace4])
