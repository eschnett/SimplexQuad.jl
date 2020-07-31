using Makie
using SimplexQuad

vertices = [
    0.0 0.0
    1.0 0.0
    0.0 1.0
]
# Need at least 4 vertices for plotting?
vertices1 = [vertices; vertices[1:1, :]]

for npoints in [1, 2, 3, 4, 10, 100]
    X, W = simplexquad(npoints, vertices)

    scene = Scene()
    lines!(scene, vertices1; color = :black, linewidth = 2)
    scatter!(scene, X[:, 1], X[:, 2]; color = :red, markersize = 0.3 * sqrt.(W))
    text!(scene, "N=$npoints"; textsize = 0.15, position = (0.5, 0.9))
    scale!(scene, 1, 1)

    Makie.save("figures/gauß-points-$npoints.png", scene;
               resolution = (300, 300))
end
