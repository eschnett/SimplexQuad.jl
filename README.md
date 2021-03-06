# Gaussian Quadrature for an n-dimensional simplex

* [GitHub](https://github.com/eschnett/SimplexQuad.jl): Source code
  repository
* [![GitHub
  CI](https://github.com/eschnett/SimplexQuad.jl/workflows/CI/badge.svg)](https://github.com/eschnett/SimplexQuad.jl/actions)

## Provenance of this package

This code was originally published by Greg von Winckel (Contact:
gregvw(at)math(dot)unm(dot)edu, <http://math.unm.edu/~gregvw>) on the
MathWorks File Exchange
<https://www.mathworks.com/matlabcentral/fileexchange/9435-n-dimensional-simplex-quadrature>.
(The given email and web addresses seems now defunct; however, a web
search for Greg von Winckel easily finds up-to-date pointers.) The
code in this package is a fairly literal translation from Matlab to
Julia.

## Description

Construct Gauss points and weights for an `n`-dimensional simplex
domain with vertices specified by the `n*(n-1)` matrix `vert`, where
each row contains the coordinates `(x1,...,xn)` for a vertex. The
order of the quadrature scheme in each dimension must be the same in
this implementation.

## Sample usage

```Julia
X, W = simplexquad(n, vert)     # Specify the vertices
X, W = simplexquad(T, n, dim)   # Specify the dimension and use unit simplex
```

`X` will be a matrix for which the `j`th column are the grid points in
each coordinate `xj`.

To integrate a function `f`, use e.g.

```Julia
sum(W[i] * f(X[i,:]) for i in 1:length(W))
```

I tested the package for up to `D=5` dimensions and order `N=10`, and
found the integration error for polynomials of order `P≤N` (which
should have only floating-point round-off error) to be less than
`10eps`. This is tested by the test suite.

## Results for a 2d triangle

Integration points and weights for various numbers of points:

![N=1](figures/gauß-points-1.png)
![N=2](figures/gauß-points-2.png)
![N=3](figures/gauß-points-3.png)
![N=4](figures/gauß-points-4.png)
![N=10](figures/gauß-points-10.png)
![N=100](figures/gauß-points-100.png)

The area of the red disks corresponds to their integration weights. It
is evident that the location of the integration points is not a tensor
product and does not respect the symmetry of the domain.

## Trivia

The first four simplexes are

n | Domain
--|------------
1 | Interval
2 | Triangle
3 | Tetrahedron
4 | Pentatope

## TODO

Examine `JuAFEM.jl` and its `generate_quadrature.jl` file
<https://kristofferc.github.io/JuAFEM.jl/dev/reference/quadrature/>.
