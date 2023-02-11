using PlotlyJS

####################################################
Y(u,v) = 0.5a(u)
Z(u,v) = 0.25 .- v .* (0.25 .- b(u))

U = X
V = range(0, 1, length=51)

#trace1: the left wall
xs = [u for u in U, v in V]
ys = [Y(u, v) for u in U, v in V]
zs = [Z(u, v) for u in U, v in V]

trace1=PlotlyJS.surface(x=xs, y=ys, z=zs,
    lighting_specular=0.15,
    opacity = 0.7,
    lightposition=attr(x=20, y=10, z=100),
    colorscale=[[0, "rgb(100,100,100)"], [1, "rgb(100,100,100)"]],
    showscale=false)

#trace2: the right wall
xs = [u for u in U, v in V]
ys = [-Y(u, v) for u in U, v in V]
zs = [Z(u, v) for u in U, v in V]

trace2=PlotlyJS.surface(x=xs, y=ys, z=zs,
    lighting_specular=0.15,
    opacity = 0.7,
    lightposition=attr(x=20, y=10, z=100),
    colorscale=[[0, "rgb(100,100,100)"], [1, "rgb(100,100,100)"]],
    showscale=false)

#trace3: bottom
U = X
V = range(-1, 1, length=51)

xs = [u for u in U, v in V]
ys = [0.5v*a(u) for u in U, v in V]
zs = [b(u) for u in U, v in V]

trace3 = PlotlyJS.surface(x=xs, y=ys, z=zs,
    lighting_specular=0.15,
    opacity = 0.7,
    lightposition=attr(x=20, y=10, z=100),
    colorscale=[[0, "rgb(100,100,100)"], [1, "rgb(100,100,100)"]],
    showscale=false)

#trace4: water surface
U = X
V = range(-1, 1, length=51)

AH = Vec[:,1]
H = AH./(A)

xs = [u for u in U, v in V]
ys = [0.5v*a(u) for u in U, v in V]

zs = B+H
for i in 1:50
    zs=hcat(zs,B+H)
end

trace4 = PlotlyJS.surface(x=xs, y=ys, z=zs,
    lighting_specular=0.05,
    lightposition=attr(x=20, y=10, z=1e3),
    colorscale=[[0, "rgb(10,70,150)"], [1, "rgb(210,230,245)"]])


# set aspect ratio, title, and other things as needed
layout = PlotlyJS.Layout(scene=attr(aspectmode="data"),
                title=attr(text="t=0.8"))

# plot everything
PlotlyJS.plot([trace1, trace2, trace3, trace4],layout)
